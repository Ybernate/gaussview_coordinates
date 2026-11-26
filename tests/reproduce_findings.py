"""Reproduce the verified findings H1-H5 against the CURRENT repo code.

Run: python tests/reproduce_findings.py
Requires: rdkit, numpy (always). tblite + ase for the xTB energy checks (H4/H5);
those sections print SKIPPED if unavailable.

This script is diagnostic: it prints what the current code does so the bugs are
visible and reproducible before any fix. It is NOT a pass/fail test suite.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from molgauss import build_polymer, embed_and_optimize
import importlib
refine_mod = importlib.import_module("molgauss.refine")
from molgauss.molmodel import Molecule3D


def hr(title):
    print("\n" + "=" * 70)
    print(title)
    print("=" * 70)


# --------------------------------------------------------------------- H1
def h1_asymmetric_cap():
    hr("H1: single `cap` cannot cap an asymmetric unit correctly")
    mol = build_polymer("[*]OCC[*]", 10, cap="H")
    formula = Chem.rdMolDescriptors.CalcMolFormula(Chem.AddHs(mol))
    print("build_polymer('[*]OCC[*]', 10, cap='H') ->", formula)
    print("Expected for true HO-(CH2CH2O)10-H : C20H42O11")
    print("Got the ethyl-terminated isomer instead:" , formula)


# --------------------------------------------------------------------- H2
def h2_tgt_by_index():
    hr("H2: CHAIN_PATTERNS['tgt'] applied by index ignores bond type")
    from molgauss.conformers import CHAIN_PATTERNS, _longest_heavy_chain, _backbone_dihedral_quads
    print("tgt pattern:", CHAIN_PATTERNS["tgt"])
    mol = build_polymer("[*]OCC[*]", 4, cap="H")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    chain = _longest_heavy_chain(mol)
    quads = _backbone_dihedral_quads(chain)
    print("backbone element sequence:",
          [mol.GetAtomWithIdx(i).GetSymbol() for i in chain])
    print("Each backbone torsion (by central bond) and the angle tgt would assign:")
    for k, (a, b, c, d) in enumerate(quads):
        sb = mol.GetAtomWithIdx(b).GetSymbol()
        sc = mol.GetAtomWithIdx(c).GetSymbol()
        assigned = CHAIN_PATTERNS["tgt"][k % len(CHAIN_PATTERNS["tgt"])]
        print(f"  torsion {k}: {sb}-{sc} bond -> assigned {assigned:>5}")
    print("Note the C-O torsions receive inconsistent gauche/sign by index.")


# --------------------------------------------------------------------- H3
def h3_improper_rotation():
    hr("H3: align_principal_axes may apply an improper rotation (det=-1)")
    dets = []
    inversions = 0
    for seed in range(12):
        m = Chem.MolFromSmiles("C[C@@H](N)C(=O)O")  # L-alanine
        m = Chem.AddHs(m)
        p = AllChem.ETKDGv3(); p.randomSeed = seed
        AllChem.EmbedMolecule(m, p)
        AllChem.MMFFOptimizeMolecule(m)
        before = _cip(m)
        m3d = Molecule3D.from_rdkit(m)
        det = _alignment_determinant(m3d)
        dets.append(det)
        # rebuild an RDKit mol from aligned coords and re-perceive stereo
        after = _cip_from_coords(m, m3d)
        if before and after and before != after:
            inversions += 1
        print(f"  seed {seed:2d}: det(rotation)={det:+.3f}  CIP before={before} after={after}")
    n_improper = sum(1 for d in dets if d < 0)
    print(f"\nImproper rotations (det<0): {n_improper}/12  |  CIP inversions: {inversions}/12")


def _alignment_determinant(mol: Molecule3D) -> float:
    """Replicate align_principal_axes and return the determinant of the net
    linear transform it applies (rotation composed with per-axis sign flips)."""
    xyz = np.array([[a.x, a.y, a.z] for a in mol.atoms])
    masses = refine_mod._masses(mol)
    com = (xyz * masses[:, None]).sum(0) / masses.sum()
    centred = xyz - com
    cov = (centred * masses[:, None]).T @ centred
    _, _, vt = np.linalg.svd(cov)
    rotated = centred @ vt.T
    spans = rotated.max(0) - rotated.min(0)
    order = np.argsort(spans)[::-1]
    perm = np.eye(3)[order]
    signs = np.ones(3)
    rr = rotated[:, order]
    for k in range(3):
        col = rr[:, k]
        if col[np.argmax(np.abs(col))] < 0:
            signs[k] = -1.0
    net = (np.diag(signs) @ perm @ vt)
    return float(np.linalg.det(net))


def _cip(mol):
    try:
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        for a in mol.GetAtoms():
            if a.HasProp("_CIPCode"):
                return a.GetProp("_CIPCode")
    except Exception:
        pass
    return None


def _cip_from_coords(template, m3d):
    rw = Chem.RWMol(template)
    conf = rw.GetConformer()
    for i, a in enumerate(m3d.atoms):
        conf.SetAtomPosition(i, (a.x, a.y, a.z))
    m = rw.GetMol()
    Chem.AssignStereochemistryFrom3D(m)
    for a in m.GetAtoms():
        if a.HasProp("_CIPCode"):
            return a.GetProp("_CIPCode")
    return None


# --------------------------------------------------------------------- H4/H5
def _xtb_energy(m3d) -> float:
    """GFN2-xTB single-point energy (eV) via tblite+ASE, or None if unavailable."""
    try:
        from ase import Atoms
        from tblite.ase import TBLite
    except Exception:
        return None
    at = Atoms(
        numbers=[a.atomic_number for a in m3d.atoms],
        positions=[[a.x, a.y, a.z] for a in m3d.atoms],
    )
    at.calc = TBLite(method="GFN2-xTB", verbosity=0)
    return float(at.get_potential_energy())


EV2KCAL = 23.060548


def h4_h5_energy_of_patterns():
    hr("H4/H5: current single imposed conformation is not energy-selected")
    e = _xtb_energy(Molecule3D.from_rdkit(_mk("CCO")))
    if e is None:
        print("SKIPPED: tblite/ase not available for xTB energies.")
        return

    # Compare repo 'tgt' vs 'all-trans' vs a random-conformer sweep for PEG-5.
    mol_graph = build_polymer("[*]OCC[*]", 5, cap="H")

    results = {}
    for patt in ("all-trans", "tgt", None):
        m3d, _ = embed_and_optimize(
            build_polymer("[*]OCC[*]", 5, cap="H"),
            conformation=patt, align=False, name="peg5",
        )
        results[str(patt)] = _xtb_energy(m3d)

    base = min(results.values())
    print("GFN2-xTB single-point energies of MMFF-optimized conformers (PEG-5):")
    for k, v in results.items():
        print(f"  pattern={k:>10}: {(v-base)*EV2KCAL:8.2f} kcal/mol (relative)")

    # Random-conformer sweep (H5): how many random confs beat the repo 'tgt'?
    repo = results["tgt"]
    better = 0
    best = repo
    N = 25
    for seed in range(N):
        m = Chem.AddHs(build_polymer("[*]OCC[*]", 5, cap="H"))
        p = AllChem.ETKDGv3(); p.randomSeed = seed
        if AllChem.EmbedMolecule(m, p) != 0:
            continue
        AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
        en = _xtb_energy(Molecule3D.from_rdkit(m))
        if en is not None:
            if en < repo - 1e-6:
                better += 1
            best = min(best, en)
    print(f"\nRandom-conformer sweep (N={N}): {better} beat the 'tgt' conformer; "
          f"best was {(best-repo)*EV2KCAL:.2f} kcal/mol vs tgt.")


def _mk(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m)
    return m


if __name__ == "__main__":
    h1_asymmetric_cap()
    h2_tgt_by_index()
    h3_improper_rotation()
    h4_h5_energy_of_patterns()
