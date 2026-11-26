"""Geometry refinement for Gaussian-ready structures.

This is the mandatory pre-optimization pipeline that every generated structure
passes through before export. Its job is to hand Gaussian a geometry that
converges fast and avoids common failures:

1. Standard-bond sanity: report bonds far from their DFT equilibrium length so
   the force-field step starts from sensible distances.
2. Force-field minimization (handled upstream in ``geometry.py``; this module
   provides the clash filter that validates its result).
3. Steric-clash filter: flag any non-bonded pair closer than a cutoff (0.8 A by
   default) — such contacts make the SCF hard to converge.
4. Principal-axis (PCA) alignment: translate the centre of mass to the origin
   and rotate so the longest molecular axis lies on X. This makes GaussView
   open on a tidy orientation and helps Gaussian detect point-group symmetry
   (C2, Cs, C2v, ...), which cuts integral evaluation time.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np

from rdkit import Chem

from .molmodel import Molecule3D

# Equilibrium bond lengths (angstrom) used as reference targets. These are
# typical DFT/experimental single-bond values; they are references for the
# diagnostic, not hard constraints.
STANDARD_BONDS: Dict[frozenset, float] = {
    frozenset(("C", "C")): 1.52,
    frozenset(("C", "N")): 1.47,
    frozenset(("C", "O")): 1.42,
    frozenset(("C", "H")): 1.09,
    frozenset(("O", "H")): 0.96,
    frozenset(("N", "H")): 1.01,
    frozenset(("C", "S")): 1.82,
    frozenset(("C", "F")): 1.35,
    frozenset(("C", "Cl")): 1.77,
}

# Covalent radii (angstrom) for a light-weight bond test when no graph exists.
_COVALENT_RADII = {
    "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57, "P": 1.07,
    "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39, "Si": 1.11, "B": 0.84,
    "Fe": 1.32, "Pt": 1.36, "Au": 1.36, "Pd": 1.39, "Ru": 1.46, "Na": 1.66,
    "K": 2.03, "Mg": 1.41, "Ca": 1.76, "Zn": 1.22, "Cu": 1.32, "Ni": 1.24,
}


@dataclass
class RefineReport:
    """Diagnostics produced by the refinement pass."""

    clashes: List[Tuple[int, int, float]] = field(default_factory=list)
    off_standard_bonds: List[Tuple[int, int, str, float, float]] = field(
        default_factory=list
    )
    aligned: bool = False
    rms_long_axis: float = 0.0

    @property
    def has_clashes(self) -> bool:
        return bool(self.clashes)

    def summary(self) -> str:
        parts = []
        parts.append("aligned to principal axes" if self.aligned else "not aligned")
        if self.clashes:
            parts.append(f"{len(self.clashes)} steric clash(es) < cutoff")
        else:
            parts.append("no steric clashes")
        if self.off_standard_bonds:
            parts.append(f"{len(self.off_standard_bonds)} bond(s) off standard length")
        return "; ".join(parts)


def _coords(mol: Molecule3D) -> np.ndarray:
    return np.array([[a.x, a.y, a.z] for a in mol.atoms], dtype=float)


def _masses(mol: Molecule3D) -> np.ndarray:
    pt = Chem.GetPeriodicTable()
    return np.array(
        [pt.GetAtomicWeight(a.atomic_number) for a in mol.atoms], dtype=float
    )


def detect_clashes(
    mol: Molecule3D, cutoff: float = 0.8, bonded_pairs: set | None = None
) -> List[Tuple[int, int, float]]:
    """Return non-bonded atom pairs closer than ``cutoff`` angstrom.

    ``bonded_pairs`` is an optional set of frozenset({i, j}); pairs in it are
    skipped (they are legitimate bonds). Without it, any pair closer than the
    cutoff is reported, which is still a valid "something is overlapping" check
    because real bonds are >= ~0.95 A and the default cutoff is 0.8 A.
    """
    xyz = _coords(mol)
    n = len(xyz)
    out: List[Tuple[int, int, float]] = []
    for i in range(n):
        for j in range(i + 1, n):
            if bonded_pairs and frozenset((i, j)) in bonded_pairs:
                continue
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < cutoff:
                out.append((i, j, round(d, 3)))
    return out


def _bonded_pairs_from_distance(mol: Molecule3D, tol: float = 0.45) -> set:
    """Infer bonded pairs from covalent radii (used when no graph is present)."""
    xyz = _coords(mol)
    pairs = set()
    n = len(xyz)
    for i in range(n):
        ri = _COVALENT_RADII.get(mol.atoms[i].symbol, 0.77)
        for j in range(i + 1, n):
            rj = _COVALENT_RADII.get(mol.atoms[j].symbol, 0.77)
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < ri + rj + tol:
                pairs.add(frozenset((i, j)))
    return pairs


def check_standard_bonds(
    mol: Molecule3D, tolerance: float = 0.15
) -> List[Tuple[int, int, str, float, float]]:
    """Report bonded pairs whose length deviates from the standard reference.

    Returns tuples ``(i, j, "C-O", actual, expected)`` for bonds off by more
    than ``tolerance`` angstrom. Bonds are inferred from covalent radii.
    """
    xyz = _coords(mol)
    bonded = _bonded_pairs_from_distance(mol)
    out = []
    for pair in bonded:
        i, j = tuple(pair)
        key = frozenset((mol.atoms[i].symbol, mol.atoms[j].symbol))
        expected = STANDARD_BONDS.get(key)
        if expected is None:
            continue
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        if abs(d - expected) > tolerance:
            label = "-".join(sorted((mol.atoms[i].symbol, mol.atoms[j].symbol)))
            out.append((i, j, label, round(d, 3), expected))
    return out


def align_principal_axes(mol: Molecule3D) -> float:
    """Translate COM to the origin and rotate the longest axis onto X.

    Uses PCA (SVD of mass-weighted, centred coordinates). The principal axis
    with the largest spatial extent becomes X, the next Y, the last Z.

    The net linear transform is assembled explicitly as a single matrix
    ``R = S @ P @ V`` (principal-axis rotation ``V``, axis permutation ``P``,
    per-axis sign flips ``S``) and constrained to be a PROPER rotation
    (``det(R) = +1``). A naive independent sign flip per axis can yield a
    reflection (``det = -1``), which mirror-inverts the geometry and silently
    swaps the enantiomer of a chiral molecule. If the assembled transform is
    improper, the sign of the least-significant axis (Z) is flipped to restore
    ``det = +1`` while preserving the long-axis-on-X convention.

    Returns the extent (span) along the resulting X axis.
    """
    if len(mol.atoms) < 2:
        return 0.0
    xyz = _coords(mol)
    masses = _masses(mol)
    com = (xyz * masses[:, None]).sum(axis=0) / masses.sum()
    centred = xyz - com

    # Mass-weighted covariance -> principal directions (rows of V).
    cov = (centred * masses[:, None]).T @ centred
    _, _, vt = np.linalg.svd(cov)
    V = vt  # each row is a principal axis, largest-variance first

    # Spatial span along each principal axis, to order longest -> X.
    projected = centred @ V.T
    spans = projected.max(axis=0) - projected.min(axis=0)
    order = np.argsort(spans)[::-1]
    P = np.eye(3)[order]  # permutation matrix

    # Deterministic sign convention: along each (ordered) axis, make the atom
    # of largest magnitude positive. Removes the SVD sign ambiguity so repeated
    # runs are reproducible.
    ordered = projected[:, order]
    signs = np.ones(3)
    for k in range(3):
        col = ordered[:, k]
        if col[np.argmax(np.abs(col))] < 0:
            signs[k] = -1.0
    S = np.diag(signs)

    # Net transform and determinant check.
    R = S @ P @ V
    if np.linalg.det(R) < 0:
        # Improper (reflection). Flip the least-significant axis (Z, index 2)
        # so handedness is preserved and the X/Y long-axis convention stays.
        signs[2] *= -1.0
        S = np.diag(signs)
        R = S @ P @ V

    rotated = centred @ R.T
    for a, pos in zip(mol.atoms, rotated):
        a.x, a.y, a.z = float(pos[0]), float(pos[1]), float(pos[2])
    return float(rotated[:, 0].max() - rotated[:, 0].min())


def refine(
    mol: Molecule3D,
    clash_cutoff: float = 0.8,
    align: bool = True,
    check_bonds: bool = True,
) -> RefineReport:
    """Run the full refinement pass in place and return diagnostics.

    This does NOT run the force field itself (that happens in ``geometry.py``
    during embedding); it validates the result (clashes, bond sanity) and
    performs the PCA alignment that helps Gaussian find symmetry.
    """
    report = RefineReport()

    bonded = _bonded_pairs_from_distance(mol)
    report.clashes = detect_clashes(mol, cutoff=clash_cutoff, bonded_pairs=bonded)
    if check_bonds:
        report.off_standard_bonds = check_standard_bonds(mol)

    if align:
        report.rms_long_axis = align_principal_axes(mol)
        report.aligned = True
    return report
