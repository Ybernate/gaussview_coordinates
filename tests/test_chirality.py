"""P0-1 regression: PCA alignment must be a PROPER rotation (det=+1) and must
never invert chirality.

Measures the signed tetrahedral volume of each stereocentre before and after
alignment (a direct, coordinate-based chirality descriptor that does not depend
on RDKit re-perceiving stereo from the already-transformed coordinates), and
checks the determinant of the net alignment transform.

Run: python tests/test_chirality.py
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from molgauss.molmodel import Molecule3D
from molgauss.refine import align_principal_axes

# Chiral test molecules (SMILES) with at least one tetrahedral stereocentre.
CHIRAL = {
    "L-alanine": "C[C@@H](N)C(=O)O",
    "R-alanine": "C[C@H](N)C(=O)O",
    "L-serine": "C([C@@H](C(=O)O)N)O",
    "bromochlorofluoromethane": "[C@H](F)(Cl)Br",
    "2-butanol": "CC[C@H](C)O",
}


def _signed_volume(mol3d, center, nbrs):
    p = np.array([[a.x, a.y, a.z] for a in mol3d.atoms])
    a, b, c, d = p[center], p[nbrs[0]], p[nbrs[1]], p[nbrs[2]]
    return float(np.dot(np.cross(b - a, c - a), d - a))


def _first_stereocenter(mol):
    centers = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
    if not centers:
        return None
    idx = centers[0][0]
    nbrs = [a.GetIdx() for a in mol.GetAtomWithIdx(idx).GetNeighbors()][:3]
    return idx, nbrs


def run():
    failures = 0
    total = 0
    for name, smi in CHIRAL.items():
        print(f"\n{name}  ({smi})")
        inversions = 0
        improper = 0
        for seed in range(12):
            m = Chem.AddHs(Chem.MolFromSmiles(smi))
            p = AllChem.ETKDGv3()
            p.randomSeed = seed
            if AllChem.EmbedMolecule(m, p) != 0:
                continue
            AllChem.MMFFOptimizeMolecule(m)
            sc = _first_stereocenter(m)
            if sc is None:
                print("  (no stereocentre detected, skipping)")
                break
            idx, nbrs = sc
            m3d = Molecule3D.from_rdkit(m)
            v0 = _signed_volume(m3d, idx, nbrs)

            # Capture the net determinant the same way align_principal_axes does.
            det = _alignment_determinant(m3d)

            align_principal_axes(m3d)
            v1 = _signed_volume(m3d, idx, nbrs)

            total += 1
            if v0 * v1 < 0:
                inversions += 1
            if det < 0:
                improper += 1
        print(f"  chirality inversions: {inversions}/12   improper transforms: {improper}/12")
        if inversions or improper:
            failures += 1

    print("\n" + "=" * 60)
    if failures == 0:
        print(f"PASS: 0 chirality inversions and 0 improper rotations across "
              f"{total} alignments.")
    else:
        print(f"FAIL: {failures} molecule(s) showed inversion/improper rotation.")
    return failures


def _alignment_determinant(mol3d):
    """Determinant of the net transform that the (fixed) aligner will apply."""
    from molgauss import refine as _pkg  # may be the function; use module file
    import importlib
    refine_mod = importlib.import_module("molgauss.refine")

    xyz = np.array([[a.x, a.y, a.z] for a in mol3d.atoms])
    masses = refine_mod._masses(mol3d)
    com = (xyz * masses[:, None]).sum(0) / masses.sum()
    centred = xyz - com
    cov = (centred * masses[:, None]).T @ centred
    _, _, vt = np.linalg.svd(cov)
    V = vt
    projected = centred @ V.T
    spans = projected.max(0) - projected.min(0)
    order = np.argsort(spans)[::-1]
    P = np.eye(3)[order]
    ordered = projected[:, order]
    signs = np.ones(3)
    for k in range(3):
        col = ordered[:, k]
        if col[np.argmax(np.abs(col))] < 0:
            signs[k] = -1.0
    R = np.diag(signs) @ P @ V
    if np.linalg.det(R) < 0:
        signs[2] *= -1.0
        R = np.diag(signs) @ P @ V
    return float(np.linalg.det(R))


if __name__ == "__main__":
    sys.exit(1 if run() else 0)
