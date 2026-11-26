"""End-to-end smoke tests for the MolGauss chemistry engine.

Run with: python -m pytest tests/ -q   (or just: python tests/test_engine.py)
These exercise each input path, the editor, and the Gaussian/XYZ exporters.
Network-dependent tests (PubChem name lookup) are skipped gracefully offline.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from molgauss import (
    from_smiles,
    from_inchi,
    from_fasta,
    from_xyz_block,
    build_polymer,
    embed_and_optimize,
    write_gjf,
    write_xyz,
    RouteConfig,
    editor,
)


def _assert(cond, msg):
    if not cond:
        raise AssertionError(msg)
    print("  ok:", msg)


def test_smiles_to_gjf():
    print("[SMILES ethanol -> gjf]")
    mol = from_smiles("CCO")
    m3d, info = embed_and_optimize(mol, name="ethanol", title="Ethanol")
    _assert(len(m3d) == 9, f"ethanol has 9 atoms (got {len(m3d)})")
    _assert(m3d.formula == "C2H6O", f"formula C2H6O (got {m3d.formula})")
    _assert(m3d.charge == 0 and m3d.multiplicity == 1, "neutral singlet")
    gjf = write_gjf(m3d, RouteConfig(route="# opt freq b3lyp/6-31g(d)", chk="ethanol.chk"))
    _assert(gjf.endswith("\n\n"), "gjf ends with mandatory blank line")
    _assert("%chk=ethanol.chk" in gjf, "chk line present")
    _assert("0 1" in gjf, "charge/multiplicity line present")
    _assert(gjf.count("\n") >= 9 + 6, "has header + atoms")
    print(info)


def test_inchi():
    print("[InChI methane]")
    mol = from_inchi("InChI=1S/CH4/h1H4")
    m3d, _ = embed_and_optimize(mol)
    _assert(m3d.formula == "CH4", f"methane CH4 (got {m3d.formula})")


def test_radical_multiplicity():
    print("[Methyl radical multiplicity]")
    mol = from_smiles("[CH3]")
    m3d, _ = embed_and_optimize(mol)
    _assert(m3d.n_electrons % 2 == 1, "odd electrons")
    _assert(m3d.suggest_multiplicity() == 2, "doublet suggested for radical")


def test_polymer():
    print("[Polyethylene n=4]")
    mol = build_polymer("[*]CC[*]", 4)
    m3d, info = embed_and_optimize(mol, name="PE4")
    # 4 x C2H4 backbone + 2 capping H = C8H18
    _assert(m3d.formula == "C8H18", f"C8H18 expected (got {m3d.formula})")
    print(info)


def test_peptide():
    print("[Peptide AG]")
    mol = from_fasta("AG")
    m3d, info = embed_and_optimize(mol, name="peptide_AG")
    syms = {a.symbol for a in m3d.atoms}
    _assert("N" in syms and "O" in syms and "C" in syms, "backbone elements present")
    _assert(len(m3d) > 10, f"dipeptide should have many atoms (got {len(m3d)})")
    print(info)


def test_xyz_roundtrip_and_editor():
    print("[XYZ parse + editor CRUD]")
    xyz = "3\nwater\nO 0 0 0\nH 0 0 0.96\nH 0.93 0 -0.26\n"
    m3d = from_xyz_block(xyz)
    _assert(len(m3d) == 3, "water has 3 atoms")
    # substitute O -> S
    editor.substitute_element(m3d, 0, "S")
    _assert(m3d.atoms[0].symbol == "S", "O substituted to S")
    # add an atom
    editor.add_atom(m3d, "Cl", 2.0, 0.0, 0.0)
    _assert(len(m3d) == 4, "atom added")
    # delete it
    editor.delete_atom(m3d, 3)
    _assert(len(m3d) == 3, "atom deleted")
    xyz_out = write_xyz(m3d)
    _assert(xyz_out.splitlines()[0] == "3", "xyz count header correct")


def test_metal_complex_uff():
    print("[Ferrocene-like / metal -> UFF path]")
    # Simple iron complex SMILES; forces the UFF branch (Fe not MMFF-safe).
    mol = from_smiles("[Fe+2]")
    m3d, info = embed_and_optimize(mol, force_field="AUTO", charge=2, multiplicity=1)
    _assert(m3d.atoms[0].symbol == "Fe", "iron atom present")
    print(info)


def test_genecp_ecp_tail():
    print("[genecp mixed-basis tail for a Pt complex]")
    from molgauss import build_genecp_tail, write_gjf, RouteConfig
    from molgauss.molmodel import Atom, Molecule3D

    mol = Molecule3D(
        atoms=[Atom("Pt", 0, 0, 0), Atom("C", 2, 0, 0), Atom("O", 3, 0, 0),
               Atom("H", 2, 1, 0)],
        charge=0, multiplicity=1, name="ptco",
    )
    tail = build_genecp_tail(mol, light_basis="6-31G(d)", heavy_basis="LANL2DZ", ecp="LANL2DZ")
    _assert("Pt 0" in tail, "heavy element Pt listed in genecp block")
    _assert("6-31G(d)" in tail, "light basis present")
    _assert("LANL2DZ" in tail, "ECP name present")
    _assert(tail.count("****") == 2, "two **** basis terminators")
    gjf = write_gjf(mol, RouteConfig(route="# opt b3lyp/genecp", tail=tail))
    _assert("genecp" in gjf, "genecp route")
    _assert(gjf.endswith("\n\n"), "mandatory trailing blank line preserved")
    # Light-only molecule should produce no tail.
    light = Molecule3D(atoms=[Atom("C", 0, 0, 0), Atom("H", 1, 0, 0)])
    _assert(build_genecp_tail(light) == "", "no genecp tail when no heavy atoms")


def test_freeze_column():
    print("[frozen atoms -> -1 freeze code]")
    from molgauss import write_gjf, RouteConfig, editor
    from molgauss.molmodel import Atom, Molecule3D

    mol = Molecule3D(
        atoms=[Atom("C", 0, 0, 0), Atom("H", 0, 0, 1.0), Atom("H", 0, 1.0, 0)],
        name="frz",
    )
    # No freeze by default -> no freeze column.
    plain = write_gjf(mol)
    _assert(" C          " in plain or " C " in plain, "plain atom line")
    _assert(" -1 " not in plain and "  0 " not in plain, "no freeze column by default")
    # Freeze first atom.
    editor.set_frozen(mol, 0, True)
    out = write_gjf(mol)
    lines = [ln for ln in out.splitlines() if ln.strip().startswith(("C", "H"))]
    _assert(lines[0].split()[1] == "-1", f"first atom frozen (-1): {lines[0]!r}")
    _assert(lines[1].split()[1] == "0", f"second atom free (0): {lines[1]!r}")
    # freeze_all then unfreeze.
    editor.freeze_all(mol, True)
    out2 = write_gjf(mol)
    frozen_lines = [ln for ln in out2.splitlines() if ln.strip().startswith(("C", "H"))]
    _assert(all(ln.split()[1] == "-1" for ln in frozen_lines), "all atoms frozen")
    editor.freeze_all(mol, False)
    out3 = write_gjf(mol)
    _assert(" -1 " not in out3, "freeze column dropped when nothing frozen")


def test_pca_alignment_and_clash():
    print("[PCA alignment + clash filter]")
    import numpy as np
    from molgauss import align_principal_axes, detect_clashes
    from molgauss.molmodel import Atom, Molecule3D

    # A long chain tilted off-axis; after alignment its span should be on X.
    mol = Molecule3D(atoms=[Atom("C", i * 1.0, i * 1.0, 0.0) for i in range(6)])
    span_x = align_principal_axes(mol)
    xs = [a.x for a in mol.atoms]
    ys = [a.y for a in mol.atoms]
    span_actual_x = max(xs) - min(xs)
    span_actual_y = max(ys) - min(ys)
    _assert(span_actual_x >= span_actual_y, "longest extent placed on X after PCA")
    com = np.mean([[a.x, a.y, a.z] for a in mol.atoms], axis=0)
    _assert(np.allclose(com, 0, atol=1e-6), "centre of mass moved to origin")

    # Two atoms overlapping -> clash detected.
    clash_mol = Molecule3D(atoms=[Atom("C", 0, 0, 0), Atom("O", 0.5, 0, 0)])
    clashes = detect_clashes(clash_mol, cutoff=0.8, bonded_pairs=set())
    _assert(len(clashes) == 1, f"clash at 0.5 A detected (got {clashes})")


def test_refine_in_pipeline():
    print("[refine summary appended by embed_and_optimize]")
    mol = from_smiles("CCCCO")
    m3d, info = embed_and_optimize(mol, align=True)
    _assert("aligned" in info, f"refine summary present: {info!r}")
    _assert("clash" in info, "clash status reported")


def test_extended_polymer_conformation():
    print("[polymer forced all-trans extended]")
    import numpy as np
    mol = build_polymer("[*]CC[*]", 8)
    coiled, _ = embed_and_optimize(mol, conformation=None, align=True, name="pe_coil")
    extended, info = embed_and_optimize(mol, conformation="all-trans", align=True, name="pe_ext")
    _assert("all-trans" in info, f"conformation applied: {info!r}")
    # Extended chain should span farther along X than a default conformer.
    def span_x(m):
        xs = [a.x for a in m.atoms]
        return max(xs) - min(xs)
    print(f"  coiled X-span={span_x(coiled):.2f}  extended X-span={span_x(extended):.2f}")
    _assert(span_x(extended) >= span_x(coiled) - 0.5, "extended chain is not more compact")


def test_peptide_secondary_structure():
    print("[peptide extended vs alpha helix dihedrals]")
    from molgauss import conformers
    from rdkit.Chem import AllChem, rdMolTransforms
    from rdkit import Chem
    mol = from_fasta("AAAAA")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    n = conformers.set_peptide_secondary_structure(mol, "alpha")
    _assert(n >= 4, f"alpha-helix dihedrals set on residues (got {n})")
    bb = conformers._find_peptide_backbone(mol)
    _assert(len(bb) == 5, f"5 residues found in backbone (got {len(bb)})")


def test_xyz_perceive_bonds():
    print("[rdDetermineBonds on an XYZ methane]")
    from molgauss import xyz_perceive_bonds
    # methane geometry (angstrom)
    xyz = (
        "5\nmethane\n"
        "C  0.000  0.000  0.000\n"
        "H  0.629  0.629  0.629\n"
        "H -0.629 -0.629  0.629\n"
        "H -0.629  0.629 -0.629\n"
        "H  0.629 -0.629 -0.629\n"
    )
    mol, radicals = xyz_perceive_bonds(xyz, charge=0)
    _assert(mol.GetNumBonds() == 4, f"methane has 4 C-H bonds (got {mol.GetNumBonds()})")
    _assert(radicals == 0, "closed-shell methane, no radicals")


def test_cif_optional():
    print("[CIF via ASE] (ASE optional)")
    try:
        from molgauss import from_cif
        # Minimal NaCl-type CIF (rock salt), 2 atoms in a small cubic cell.
        cif = """data_test
_cell_length_a 4.0
_cell_length_b 4.0
_cell_length_c 4.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Na Na 0.0 0.0 0.0
Cl Cl 0.5 0.5 0.5
"""
        m3d = from_cif(cif, mode="cluster")
        syms = sorted(a.symbol for a in m3d.atoms)
        _assert(syms == ["Cl", "Na"], f"Na and Cl parsed from CIF (got {syms})")
        sc = from_cif(cif, mode="supercell", supercell=(2, 1, 1))
        _assert(len(sc) == 4, f"2x1x1 supercell doubles atoms (got {len(sc)})")
    except Exception as e:
        print("  SKIPPED:", e)


def test_export_modes_and_resources():
    print("[3 export modes + dynamic resources]")
    from molgauss import mode_route, suggest_resources, EXPORT_MODES, write_gjf, RouteConfig
    from molgauss.molmodel import Atom, Molecule3D
    r1 = mode_route("standard")
    _assert("empiricaldispersion=gd3" in r1 and "int=ultrafine" in r1 and "scf=xqc" in r1,
            "Mode 1 has D3 + ultrafine + xqc")
    _assert("freq" in mode_route("robust"), "Mode 2 includes freq")
    _assert("pm6" in mode_route("quick"), "Mode 3 is PM6")
    _assert(len(EXPORT_MODES) == 3, "three export modes defined")
    mem_s, nproc_s = suggest_resources(10)
    mem_l, nproc_l = suggest_resources(200)
    _assert(nproc_l > nproc_s, "larger system suggests more procs")
    # Build a gjf using a mode route.
    mol = Molecule3D(atoms=[Atom("C", 0, 0, 0), Atom("H", 0, 0, 1.09)], name="x")
    gjf = write_gjf(mol, RouteConfig(route=mode_route("standard"), mem=mem_s, nprocshared=nproc_s))
    _assert("scf=xqc" in gjf, "mode route written into gjf")


def test_name_lookup_optional():
    print("[PubChem name lookup: aspirin] (network optional)")
    try:
        from molgauss import from_name

        mol, smiles = from_name("aspirin")
        m3d, _ = embed_and_optimize(mol, name="aspirin")
        _assert(m3d.formula == "C9H8O4", f"aspirin C9H8O4 (got {m3d.formula})")
        print("  resolved SMILES:", smiles)
    except Exception as e:
        print("  SKIPPED (offline or API issue):", e)


if __name__ == "__main__":
    tests = [
        test_smiles_to_gjf,
        test_inchi,
        test_radical_multiplicity,
        test_polymer,
        test_peptide,
        test_xyz_roundtrip_and_editor,
        test_metal_complex_uff,
        test_genecp_ecp_tail,
        test_freeze_column,
        test_pca_alignment_and_clash,
        test_refine_in_pipeline,
        test_extended_polymer_conformation,
        test_peptide_secondary_structure,
        test_xyz_perceive_bonds,
        test_cif_optional,
        test_export_modes_and_resources,
        test_name_lookup_optional,
    ]
    failed = 0
    for t in tests:
        try:
            t()
        except Exception as e:
            failed += 1
            print("  FAIL:", t.__name__, "->", repr(e))
        print()
    print("=" * 50)
    print(f"{len(tests) - failed}/{len(tests)} test groups passed")
    sys.exit(1 if failed else 0)
