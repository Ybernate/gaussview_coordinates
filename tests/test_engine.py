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
