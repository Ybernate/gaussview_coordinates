"""MolGauss — Streamlit front-end.

Generate, edit and export 3D molecular structures ready for GaussView / Gaussian.

Run with:
    streamlit run app.py
"""

from __future__ import annotations

import io
import traceback

import pandas as pd
import streamlit as st

from molgauss import (
    InputError,
    GeometryError,
    from_smiles,
    from_inchi,
    from_name,
    from_fasta,
    from_pdb_block,
    from_xyz_block,
    xyz_perceive_bonds,
    from_cif,
    build_polymer,
    embed_and_optimize,
    write_gjf,
    write_xyz,
    RouteConfig,
    build_genecp_tail,
    suggest_resources,
    mode_route,
    ROUTE_PRESETS,
    EXPORT_MODES,
    HEAVY_ELEMENTS_FOR_ECP,
    editor,
)
from molgauss.molmodel import Atom, Molecule3D, symbol_to_atomic_number

# --------------------------------------------------------------------------- #
# Page config
# --------------------------------------------------------------------------- #
st.set_page_config(
    page_title="MolGauss — 3D structures for Gaussian",
    page_icon="🧪",
    layout="wide",
)

CPK_COLORS = {
    "H": "#FFFFFF", "C": "#909090", "N": "#3050F8", "O": "#FF0D0D",
    "F": "#90E050", "Cl": "#1FF01F", "Br": "#A62929", "I": "#940094",
    "S": "#FFFF30", "P": "#FF8000", "Si": "#F0C8A0", "Fe": "#E06633",
}


# --------------------------------------------------------------------------- #
# Session state
# --------------------------------------------------------------------------- #
def _init_state() -> None:
    st.session_state.setdefault("mol", None)          # Molecule3D or None
    st.session_state.setdefault("method_info", "")    # how it was generated
    st.session_state.setdefault("source_note", "")


_init_state()


# --------------------------------------------------------------------------- #
# 3D viewer
# --------------------------------------------------------------------------- #
def render_3d(
    mol: Molecule3D,
    style: str = "stick",
    show_indices: bool = False,
    mark_frozen: bool = True,
) -> None:
    """Render an interactive py3Dmol view from the molecule's XYZ block.

    ``show_indices`` draws each atom's table index on top of its sphere so the
    user can match the 3D view to the CRUD table row by row. ``mark_frozen``
    tints frozen atoms so a partial-optimization selection is visible in 3D.
    """
    try:
        import py3Dmol
        from stmol import showmol
    except Exception as exc:  # pragma: no cover
        st.warning(f"3D viewer unavailable ({exc}). Showing coordinates instead.")
        st.code(mol.to_xyz_block(), language="text")
        return

    view = py3Dmol.view(width=520, height=420)
    view.addModel(mol.to_xyz_block(), "xyz")
    if style == "stick":
        view.setStyle({"stick": {"radius": 0.12}, "sphere": {"scale": 0.22}})
    elif style == "ball-and-stick":
        view.setStyle({"stick": {"radius": 0.15}, "sphere": {"scale": 0.30}})
    elif style == "spacefill":
        view.setStyle({"sphere": {"scale": 1.0}})
    else:  # line
        view.setStyle({"line": {}})

    # Tint frozen atoms (partial optimization) so they stand out.
    if mark_frozen and any(a.frozen for a in mol.atoms):
        for i, a in enumerate(mol.atoms):
            if a.frozen:
                view.setStyle(
                    {"serial": i},
                    {"stick": {"radius": 0.12},
                     "sphere": {"scale": 0.34, "color": "0x1f77b4"}},
                )

    # Draw the table index above each atom. The label text includes the symbol
    # so "#14 (C)" is unambiguous when matching the editor table.
    if show_indices:
        for i, a in enumerate(mol.atoms):
            frozen_tag = " ❄" if a.frozen else ""
            view.addLabel(
                f"{i}:{a.symbol}{frozen_tag}",
                {
                    "position": {"x": a.x, "y": a.y, "z": a.z},
                    "backgroundColor": "0x111111",
                    "backgroundOpacity": 0.65,
                    "fontColor": "white",
                    "fontSize": 11,
                    "inFront": True,
                },
            )

    view.setBackgroundColor("0xf7f7f9")
    view.zoomTo()
    showmol(view, height=420, width=520)


# --------------------------------------------------------------------------- #
# Molecule info panel
# --------------------------------------------------------------------------- #
def show_summary(mol: Molecule3D) -> None:
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Atoms", len(mol))
    c2.metric("Formula", mol.formula)
    c3.metric("Charge", mol.charge)
    c4.metric("Multiplicity", mol.multiplicity)
    if st.session_state.method_info:
        st.caption(f"Geometry method: {st.session_state.method_info}")


# --------------------------------------------------------------------------- #
# Sidebar — input
# --------------------------------------------------------------------------- #
st.sidebar.title("🧪 MolGauss")
st.sidebar.caption("3D molecular structures ready for GaussView / Gaussian")

input_mode = st.sidebar.radio(
    "Input type",
    [
        "SMILES",
        "InChI",
        "Name / IUPAC (PubChem)",
        "Polymer (monomer + n)",
        "Peptide (FASTA)",
        "Upload file (PDB / CIF / XYZ / MOL)",
    ],
)

ff_choice = st.sidebar.selectbox(
    "Force field (pre-optimization)",
    ["AUTO (MMFF94, UFF for metals)", "MMFF94", "UFF"],
    help="AUTO uses MMFF94 for organics and switches to UFF when metals or "
    "heavy elements are present.",
)
FF = {"AUTO (MMFF94, UFF for metals)": "AUTO", "MMFF94": "MMFF94", "UFF": "UFF"}[ff_choice]

align_axes = st.sidebar.checkbox(
    "Align to principal axes (PCA)", value=True,
    help="Centre the molecule at the origin and put its longest axis on X. "
    "Gives GaussView a tidy orientation and helps Gaussian detect point-group "
    "symmetry, which speeds up the calculation.",
)

charge_override = st.sidebar.number_input(
    "Charge (override)", value=0, step=1,
    help="Set the total molecular charge. Leave at the auto value for neutral "
    "closed-shell molecules; set explicitly for ions / complexes.",
)
use_charge_override = st.sidebar.checkbox("Apply charge override", value=False)
mult_override = st.sidebar.number_input("Multiplicity (override)", value=1, min_value=1, step=1)
use_mult_override = st.sidebar.checkbox("Apply multiplicity override", value=False)


def _finalize(mol3d: Molecule3D, info: str, note: str) -> None:
    if use_charge_override:
        mol3d.charge = int(charge_override)
    if use_mult_override:
        mol3d.multiplicity = int(mult_override)
    st.session_state.mol = mol3d
    st.session_state.method_info = info
    st.session_state.source_note = note


def _handle_error(exc: Exception) -> None:
    if isinstance(exc, (InputError, GeometryError, ValueError, IndexError)):
        st.sidebar.error(str(exc))
    else:
        st.sidebar.error(f"Unexpected error: {exc}")
        st.sidebar.code(traceback.format_exc(), language="text")


# ------------------------------- input forms ------------------------------- #
if input_mode == "SMILES":
    smiles = st.sidebar.text_input("SMILES", value="CCO", help="e.g. CCO for ethanol")
    name = st.sidebar.text_input("Molecule name (for filenames)", value="molecule")
    if st.sidebar.button("Generate 3D", type="primary"):
        try:
            mol = from_smiles(smiles)
            m3d, info = embed_and_optimize(
                mol, force_field=FF, name=name, title=name, align=align_axes
            )
            _finalize(m3d, info, f"SMILES: {smiles}")
        except Exception as exc:
            _handle_error(exc)

elif input_mode == "InChI":
    inchi = st.sidebar.text_area("InChI", value="InChI=1S/CH4/h1H4", height=100)
    name = st.sidebar.text_input("Molecule name", value="molecule")
    if st.sidebar.button("Generate 3D", type="primary"):
        try:
            mol = from_inchi(inchi)
            m3d, info = embed_and_optimize(
                mol, force_field=FF, name=name, title=name, align=align_axes
            )
            _finalize(m3d, info, "InChI input")
        except Exception as exc:
            _handle_error(exc)

elif input_mode == "Name / IUPAC (PubChem)":
    cname = st.sidebar.text_input("Common / IUPAC name", value="aspirin")
    st.sidebar.caption("Resolved online via PubChem.")
    if st.sidebar.button("Look up & generate 3D", type="primary"):
        try:
            with st.spinner("Querying PubChem…"):
                mol, resolved = from_name(cname)
            m3d, info = embed_and_optimize(
                mol, force_field=FF, name=cname, title=cname, align=align_axes
            )
            _finalize(m3d, info, f"PubChem name → SMILES: {resolved}")
        except Exception as exc:
            _handle_error(exc)

elif input_mode == "Polymer (monomer + n)":
    monomer = st.sidebar.text_input(
        "Monomer SMILES (mark two [*] attachment points)",
        value="[*]CC[*]",
        help="e.g. [*]CC[*] = polyethylene, [*]CC(c1ccccc1)[*] = polystyrene",
    )
    n_units = st.sidebar.number_input("Degree of polymerization (n)", 1, 200, value=5, step=1)
    cap = st.sidebar.text_input("End-cap element", value="H")
    poly_conf = st.sidebar.selectbox(
        "Backbone conformation",
        ["all-trans (extended)", "tgt (trans-gauche helix, PEG-like)", "random coil"],
        help="Extended / regular backbones are far better DFT starting points "
        "than random coils and converge faster in Gaussian.",
    )
    conf_key = {"all-trans (extended)": "all-trans",
                "tgt (trans-gauche helix, PEG-like)": "tgt",
                "random coil": None}[poly_conf]
    name = st.sidebar.text_input("Polymer name", value="polymer")
    if st.sidebar.button("Build polymer & generate 3D", type="primary"):
        try:
            with st.spinner(f"Building chain of {n_units} units…"):
                mol = build_polymer(monomer, int(n_units), cap=cap)
                m3d, info = embed_and_optimize(
                    mol, force_field=FF, name=name, title=name,
                    conformation=conf_key, align=align_axes,
                )
            _finalize(m3d, info, f"Polymer {monomer} × {n_units} ({poly_conf})")
        except Exception as exc:
            _handle_error(exc)

elif input_mode == "Peptide (FASTA)":
    seq = st.sidebar.text_area(
        "Amino-acid sequence (one-letter)", value="ACDEFG", height=100,
        help="Standard 20 residues (ACDEFGHIKLMNPQRSTVWY). Header lines with > are ignored.",
    )
    pep_conf = st.sidebar.selectbox(
        "Initial secondary structure",
        ["extended (β strand)", "alpha helix", "polyproline-II", "random coil"],
        help="β-strand (φ≈-135, ψ≈135) or α-helix (φ≈-57, ψ≈-47) give a regular, "
        "physically sensible backbone for the optimizer to start from.",
    )
    pconf_key = {"extended (β strand)": "extended", "alpha helix": "alpha",
                 "polyproline-II": "polypro", "random coil": None}[pep_conf]
    name = st.sidebar.text_input("Peptide name", value="peptide")
    if st.sidebar.button("Build peptide & generate 3D", type="primary"):
        try:
            with st.spinner("Building peptide backbone…"):
                mol = from_fasta(seq)
                m3d, info = embed_and_optimize(
                    mol, force_field=FF, name=name, title=name,
                    conformation=pconf_key, align=align_axes,
                )
            _finalize(m3d, info, f"Peptide {seq.strip()[:30]} ({pep_conf})")
        except Exception as exc:
            _handle_error(exc)

else:  # Upload file
    up = st.sidebar.file_uploader(
        "Structure file", type=["pdb", "cif", "xyz", "mol", "sdf"]
    )
    reopt = st.sidebar.checkbox(
        "Re-optimize geometry after loading", value=False,
        help="PDB/XYZ coordinates are kept as-is by default (experimental "
        "geometry preserved). Enable to run a force-field cleanup.",
    )
    xyz_bonds = st.sidebar.checkbox(
        "XYZ: perceive bonds + charge (rdDetermineBonds)", value=True,
        help="For .xyz files with no connectivity, infer bonds and bond orders "
        "from the geometry and derive the spin multiplicity.",
    )
    cif_mode = st.sidebar.selectbox(
        "CIF handling", ["cluster (finite, gas-phase)", "unit cell", "supercell 2×2×2"],
        help="How to turn a crystal into a Gaussian-ready structure.",
    )
    if up is not None and st.sidebar.button("Load structure", type="primary"):
        try:
            raw = up.read().decode("utf-8", errors="replace")
            ext = up.name.rsplit(".", 1)[-1].lower()
            info = "Loaded from file"
            if ext == "pdb":
                m3d = from_pdb_block(raw, is_cif=False)
            elif ext == "cif":
                mode_map = {"cluster (finite, gas-phase)": ("cluster", (1, 1, 1)),
                            "unit cell": ("cell", (1, 1, 1)),
                            "supercell 2×2×2": ("supercell", (2, 2, 2))}
                cmode, scell = mode_map[cif_mode]
                m3d = from_cif(raw, mode=cmode, supercell=scell)
                info = f"Loaded CIF ({cmode})"
            elif ext == "xyz":
                if xyz_bonds:
                    try:
                        mol, radicals = xyz_perceive_bonds(raw, charge=0)
                        m3d = Molecule3D.from_rdkit(mol, name=up.name, title=up.name)
                        if radicals:
                            m3d.multiplicity = radicals + 1
                        info = "Loaded XYZ (bonds perceived)"
                    except Exception:
                        m3d = from_xyz_block(raw)
                        info = "Loaded XYZ (coordinates only; bond perception failed)"
                else:
                    m3d = from_xyz_block(raw)
            else:  # mol / sdf -> parse via RDKit then embed if needed
                from rdkit import Chem

                mol = Chem.MolFromMolBlock(raw, removeHs=False) if ext == "mol" \
                    else next(iter(Chem.ForwardSDMolSupplier(io.BytesIO(raw.encode()))), None)
                if mol is None:
                    raise InputError(f"Could not parse {up.name}")
                if mol.GetNumConformers() == 0:
                    m3d, info = embed_and_optimize(
                        mol, force_field=FF, name=up.name, align=align_axes
                    )
                else:
                    m3d = Molecule3D.from_rdkit(mol, name=up.name, title=up.name)
            if reopt:
                note = editor.cleanup_geometry(m3d, force_field=FF)
                info = note or info + " (cleanup unavailable)"
            m3d.name = up.name.rsplit(".", 1)[0]
            _finalize(m3d, info, f"Uploaded {up.name}")
        except Exception as exc:
            _handle_error(exc)


# --------------------------------------------------------------------------- #
# Main area
# --------------------------------------------------------------------------- #
st.title("MolGauss")
st.markdown(
    "Generate, edit and export **3D molecular structures ready for "
    "GaussView / Gaussian**. Pick an input type in the sidebar to begin."
)

mol: Molecule3D = st.session_state.mol

if mol is None:
    st.info(
        "No molecule yet. Enter a SMILES, name, polymer, peptide, or upload a "
        "structure file in the sidebar, then click **Generate 3D**."
    )
    with st.expander("Quick examples"):
        st.markdown(
            "- **SMILES**: `CC(=O)Oc1ccccc1C(=O)O` (aspirin)\n"
            "- **Name**: `caffeine`, `benzene`, `glucose`\n"
            "- **Polymer**: monomer `[*]CC(c1ccccc1)[*]`, n = 10 (polystyrene)\n"
            "- **Peptide**: `GIVEQCCTSICSLYQLENYCN` (insulin A-chain fragment)\n"
            "- **Metal**: SMILES `[Fe+2]` or upload a `.pdb` complex (UFF)"
        )
    st.stop()

if st.session_state.source_note:
    st.success(f"Loaded — {st.session_state.source_note}")

left, right = st.columns([1, 1])

with left:
    st.subheader("3D preview")
    pv1, pv2 = st.columns([2, 1])
    style = pv1.selectbox("Display style", ["stick", "ball-and-stick", "spacefill", "line"])
    show_idx = pv2.checkbox(
        "Atom indices", value=False,
        help="Overlay each atom's table index (and symbol) on the 3D view so "
        "you can tell exactly which atom a table row refers to.",
    )
    render_3d(mol, style=style, show_indices=show_idx)
    show_summary(mol)
    if any(a.frozen for a in mol.atoms):
        n_frozen = sum(a.frozen for a in mol.atoms)
        st.caption(f"❄ {n_frozen} atom(s) frozen for partial optimization (shown in blue).")

with right:
    st.subheader("Charge & spin")
    cc1, cc2 = st.columns(2)
    new_charge = cc1.number_input("Total charge", value=int(mol.charge), step=1, key="charge_edit")
    new_mult = cc2.number_input("Multiplicity", value=int(mol.multiplicity), min_value=1, step=1, key="mult_edit")
    if new_charge != mol.charge or new_mult != mol.multiplicity:
        mol.charge, mol.multiplicity = int(new_charge), int(new_mult)
    sug = mol.suggest_multiplicity()
    st.caption(
        f"Electrons: {mol.n_electrons} · suggested multiplicity for lowest "
        f"spin: {sug}" + ("" if sug == mol.multiplicity else "  ⚠️ differs from current")
    )

# --------------------------------------------------------------------------- #
# Atom editor (CRUD)
# --------------------------------------------------------------------------- #
st.markdown("---")
st.subheader("Atom editor")
st.caption(
    "Edit coordinates or element symbols directly in the table. Add rows at the "
    "bottom, delete rows with the row menu, then click **Apply table edits**. "
    "Substituting an element (e.g. C → Si) keeps the atom's position. Tick "
    "**Freeze** to hold an atom fixed during the Gaussian optimization "
    "(written as a `-1` freeze code)."
)

df = pd.DataFrame(
    [{"Index": i, "Symbol": a.symbol, "X": a.x, "Y": a.y, "Z": a.z, "Freeze": a.frozen}
     for i, a in enumerate(mol.atoms)]
)

edited = st.data_editor(
    df,
    num_rows="dynamic",
    use_container_width=True,
    hide_index=True,
    column_config={
        "Index": st.column_config.NumberColumn("Index", disabled=True),
        "Symbol": st.column_config.TextColumn("Symbol", help="Element symbol, e.g. C, N, Fe"),
        "X": st.column_config.NumberColumn("X (Å)", format="%.6f"),
        "Y": st.column_config.NumberColumn("Y (Å)", format="%.6f"),
        "Z": st.column_config.NumberColumn("Z (Å)", format="%.6f"),
        "Freeze": st.column_config.CheckboxColumn(
            "Freeze", help="Freeze this atom during optimization (Gaussian -1 code)."
        ),
    },
    key="atom_table",
)

ec1, ec2, ec3, ec4 = st.columns([1.1, 1.1, 1, 2])
if ec1.button("Apply table edits"):
    try:
        new_atoms = []
        for _, row in edited.iterrows():
            sym = str(row["Symbol"]).strip()
            if not sym or sym.lower() == "nan":
                continue
            symbol_to_atomic_number(sym)  # validate
            frozen = bool(row.get("Freeze", False))
            new_atoms.append(
                Atom(sym, row["X"] or 0.0, row["Y"] or 0.0, row["Z"] or 0.0, frozen)
            )
        if not new_atoms:
            st.error("The table has no valid atoms.")
        else:
            editor.replace_atoms(mol, new_atoms)
            n_frozen = sum(a.frozen for a in new_atoms)
            msg = f"Applied edits — {len(new_atoms)} atoms"
            msg += f", {n_frozen} frozen." if n_frozen else "."
            st.success(msg)
            st.rerun()
    except Exception as exc:
        st.error(f"Could not apply edits: {exc}")

if ec2.button("Freeze / free all"):
    all_frozen = all(a.frozen for a in mol.atoms) and len(mol.atoms) > 0
    editor.freeze_all(mol, not all_frozen)
    st.rerun()

if ec3.button("Force-field cleanup"):
    try:
        note = editor.cleanup_geometry(mol, force_field=FF)
        if note:
            st.session_state.method_info = note
            st.success(f"Geometry relaxed — {note}")
            st.rerun()
        else:
            st.warning("Cleanup unavailable (Open Babel force field could not run).")
    except Exception as exc:
        st.error(f"Cleanup failed: {exc}")

ec4.caption(
    "Use **Force-field cleanup** after manual edits to relax atomic clashes "
    "(rebuilds bonds from geometry and runs a short UFF/MMFF minimization). "
    "**Freeze / free all** toggles every atom at once."
)

# --------------------------------------------------------------------------- #
# Gaussian export
# --------------------------------------------------------------------------- #
st.markdown("---")
st.subheader("Export for GaussView / Gaussian")

# Export MODE: curated route sections with SCF-stabilizing keywords.
mode_labels = {k: v["label"] for k, v in EXPORT_MODES.items()}
mode_key = st.radio(
    "Export mode",
    list(EXPORT_MODES.keys()),
    format_func=lambda k: mode_labels[k],
    horizontal=False,
)
st.caption(EXPORT_MODES[mode_key]["description"])

# Dynamic resource suggestion scaled to the atom count.
sug_mem, sug_nproc = suggest_resources(len(mol))

gc1, gc2 = st.columns([1, 1])
with gc1:
    route_source = st.radio(
        "Route from", ["Export mode", "Preset", "Custom"], horizontal=True,
    )
    if route_source == "Export mode":
        default_route = mode_route(mode_key)
    elif route_source == "Preset":
        preset = st.selectbox("Route preset", list(ROUTE_PRESETS.keys()))
        default_route = ROUTE_PRESETS[preset]
    else:  # Custom — start from the current mode route, fully editable
        default_route = mode_route(mode_key)
    route = st.text_input("Route section", value=default_route)
    title = st.text_input("Job title", value=mol.title)
    chk = st.text_input("%chk file", value=f"{mol.name}.chk")
with gc2:
    st.caption(f"Suggested for {len(mol)} atoms: %mem={sug_mem}, %nprocshared={sug_nproc}")
    mem = st.text_input("%mem", value=sug_mem)
    nproc = st.number_input("%nprocshared", 1, 256, value=sug_nproc, step=1)
    fmt = st.radio("File format", [".gjf", ".com"], horizontal=True)

# --- genecp / mixed-basis helper (ECP for heavy metals) ---------------------
heavy_present = sorted({a.symbol for a in mol.atoms} & HEAVY_ELEMENTS_FOR_ECP)
is_genecp = "genecp" in route.lower() or "gen" == route.lower().split("/")[-1].strip()

auto_tail = ""
if is_genecp or heavy_present:
    with st.expander(
        "Mixed basis set / ECP (genecp) for heavy metals"
        + (f" — detected: {', '.join(heavy_present)}" if heavy_present else ""),
        expanded=is_genecp,
    ):
        st.caption(
            "Heavy metals (Pt, Au, Pd, Ru, …) have no all-electron Pople basis "
            "such as 6-31G(d), so they need an effective core potential. With a "
            "`genecp` route this builds the basis/ECP tail automatically: a "
            "regular basis for light atoms and an ECP basis for the heavy ones."
        )
        bc1, bc2, bc3 = st.columns(3)
        light_basis = bc1.text_input("Light-atom basis", value="6-31G(d)")
        heavy_basis = bc2.text_input("Heavy-atom basis", value="LANL2DZ")
        ecp_name = bc3.text_input("ECP (pseudopotential)", value="LANL2DZ")
        if heavy_present:
            auto_tail = build_genecp_tail(
                mol, light_basis=light_basis, heavy_basis=heavy_basis, ecp=ecp_name
            )
            st.markdown("**Auto-generated genecp tail**")
            st.code(auto_tail or "(no heavy atoms detected)", language="text")
        else:
            st.info("No heavy (ECP) elements detected in the current molecule.")

tail = st.text_area(
    "Tail section (custom basis / ECP, appended after coordinates)",
    value=auto_tail,
    height=140 if auto_tail else 80,
    help="Auto-filled for genecp routes with heavy atoms; edit freely. Leave "
    "empty for standard single-basis jobs.",
)

cfg = RouteConfig(route=route, chk=chk, mem=mem, nprocshared=int(nproc), tail=tail)
mol.title = title

gjf_text = write_gjf(mol, cfg)
xyz_text = write_xyz(mol)

st.markdown("**Preview of the Gaussian input file**")
st.code(gjf_text, language="text")

dc1, dc2 = st.columns(2)
dc1.download_button(
    f"⬇️ Download {fmt}",
    data=gjf_text,
    file_name=f"{mol.name}{fmt}",
    mime="chemical/x-gaussian-input",
    type="primary",
)
dc2.download_button(
    "⬇️ Download .xyz",
    data=xyz_text,
    file_name=f"{mol.name}.xyz",
    mime="chemical/x-xyz",
)

st.caption(
    "The .gjf is written with the strict Gaussian layout (Link0 → route → "
    "blank → title → blank → charge/multiplicity → coordinates → mandatory "
    "trailing blank line) so GaussView opens it directly."
)
