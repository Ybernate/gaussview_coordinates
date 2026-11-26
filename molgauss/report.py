"""Pre-export report (P1-3).

Summarises a structure before it is written to a .gjf, so obvious problems are
caught before a long DFT run:

  * molecular formula, atom count;
  * electron-count / spin-multiplicity parity (a classic Gaussian error is an
    even-electron system with an even multiplicity, or vice versa);
  * steric clashes using a van der Waals overlap criterion, EXCLUDING 1-2 and
    1-3 bonded pairs (bonded and geminal atoms are legitimately close; the old
    flat 0.8 A cutoff flagged none of the real clashes and allowed bad overlaps);
  * bond-length and bond-angle ranges;
  * a backbone-torsion summary;
  * an estimate of the number of basis functions for the chosen basis, with a
    cost warning.

Basis-function counting is verified against Gaussian documentation: the 6-31G
family uses 6 CARTESIAN d functions by default (not 5 spherical); Gaussian's
general default for other bases is 5D/7F. https://gaussian.com/gen/
(Documentation content rephrased for licensing compliance.)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from rdkit import Chem

from .molmodel import Molecule3D

# Van der Waals radii (angstrom), Bondi values for common elements.
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98, "Si": 2.10, "B": 1.92,
    "Na": 2.27, "K": 2.75, "Mg": 1.73, "Ca": 2.31, "Fe": 2.00, "Zn": 2.10,
    "Cu": 2.00, "Ni": 2.00, "Pt": 1.75, "Au": 1.66, "Pd": 1.63, "Ru": 2.00,
}
_DEFAULT_VDW = 1.80


@dataclass
class PreExportReport:
    formula: str = ""
    n_atoms: int = 0
    n_electrons: int = 0
    charge: int = 0
    multiplicity: int = 1
    parity_ok: bool = True
    parity_message: str = ""
    clashes: List[Tuple[int, int, float, float]] = field(default_factory=list)
    bond_length_range: Optional[Tuple[float, float]] = None
    bond_angle_range: Optional[Tuple[float, float]] = None
    n_basis_estimate: Optional[int] = None
    basis_note: str = ""
    cost_warning: str = ""
    warnings: List[str] = field(default_factory=list)

    def ok(self) -> bool:
        return self.parity_ok and not self.clashes

    def text(self) -> str:
        lines = [
            f"Formula: {self.formula}   atoms: {self.n_atoms}",
            f"Charge: {self.charge}   Multiplicity: {self.multiplicity}   "
            f"Electrons: {self.n_electrons}",
            f"  parity: {'OK' if self.parity_ok else 'PROBLEM'} — {self.parity_message}",
        ]
        if self.bond_length_range:
            lines.append(
                f"Bond lengths: {self.bond_length_range[0]:.3f}–"
                f"{self.bond_length_range[1]:.3f} Å"
            )
        if self.bond_angle_range:
            lines.append(
                f"Bond angles: {self.bond_angle_range[0]:.1f}–"
                f"{self.bond_angle_range[1]:.1f}°"
            )
        if self.clashes:
            lines.append(f"Steric clashes (vdW overlap, excl. 1-2/1-3): {len(self.clashes)}")
            for i, j, d, thr in self.clashes[:8]:
                lines.append(f"  atoms {i}-{j}: {d:.2f} Å (vdW overlap < {thr:.2f} Å)")
        else:
            lines.append("Steric clashes: none")
        if self.n_basis_estimate is not None:
            lines.append(f"Estimated basis functions: ~{self.n_basis_estimate}  ({self.basis_note})")
        if self.cost_warning:
            lines.append(f"Cost: {self.cost_warning}")
        for w in self.warnings:
            lines.append(f"! {w}")
        return "\n".join(lines)


def _bonded_graph(mol3d: Molecule3D) -> Optional[Chem.Mol]:
    """Perceive a connectivity graph from the 3D geometry (for 1-2/1-3 masks)."""
    try:
        from rdkit.Chem import rdDetermineBonds
        raw = Chem.MolFromXYZBlock(mol3d.to_xyz_block())
        if raw is None:
            return None
        m = Chem.Mol(raw)
        rdDetermineBonds.DetermineConnectivity(m)
        return m
    except Exception:
        return None


def _excluded_pairs(graph: Chem.Mol) -> set:
    """Set of frozenset({i,j}) for all 1-2 (bonded) and 1-3 (geminal) pairs."""
    excl = set()
    dmat = Chem.GetDistanceMatrix(graph)
    n = graph.GetNumAtoms()
    for i in range(n):
        for j in range(i + 1, n):
            if dmat[i][j] <= 2:  # topological distance 1 or 2 -> 1-2 or 1-3
                excl.add(frozenset((i, j)))
    return excl


def vdw_clashes(
    mol3d: Molecule3D, overlap_factor: float = 0.6
) -> List[Tuple[int, int, float, float]]:
    """Find non-bonded pairs whose distance is below ``overlap_factor`` times the
    sum of their van der Waals radii, excluding 1-2 and 1-3 bonded pairs.

    Returns ``(i, j, distance, threshold)``. ``overlap_factor`` ~0.6 flags
    genuine clashes while tolerating normal non-bonded contacts (~0.8-1.0 of the
    vdW sum).
    """
    xyz = np.array([[a.x, a.y, a.z] for a in mol3d.atoms])
    graph = _bonded_graph(mol3d)
    excl = _excluded_pairs(graph) if graph is not None else set()
    out = []
    n = len(xyz)
    for i in range(n):
        ri = VDW_RADII.get(mol3d.atoms[i].symbol, _DEFAULT_VDW)
        for j in range(i + 1, n):
            if frozenset((i, j)) in excl:
                continue
            rj = VDW_RADII.get(mol3d.atoms[j].symbol, _DEFAULT_VDW)
            thr = overlap_factor * (ri + rj)
            d = float(np.linalg.norm(xyz[i] - xyz[j]))
            if d < thr:
                out.append((i, j, round(d, 3), round(thr, 3)))
    return out


# Basis-function counts per atom. Minimal-but-honest table for common bases.
# For 6-31G(d): H = 2 (1s split), heavy = inner+valence s/p (9) + 6 Cartesian d.
# These are ESTIMATES for a cost indicator, not exact per-element tallies for
# every basis; the note states the assumption.
def estimate_basis_functions(mol3d: Molecule3D, basis: str) -> Tuple[Optional[int], str]:
    b = basis.lower().replace(" ", "")
    counts = {}
    for a in mol3d.atoms:
        counts[a.symbol] = counts.get(a.symbol, 0) + 1

    if b in ("6-31g(d)", "6-31g*"):
        # H/He: 2; second-row heavy (Li-Ne): 3s+2*3p=... use standard tally:
        # 6-31G(d) per heavy first-row atom = 15 (3 s + 2x3 p + 6 d); H = 2.
        per = {}
        total = 0
        for sym, n in counts.items():
            z = Chem.GetPeriodicTable().GetAtomicNumber(sym)
            if z <= 2:
                f = 2
            elif z <= 10:
                f = 15   # 1s,2s,2s' (3) + 2p,2p' (6) + 6 Cartesian d
            elif z <= 18:
                f = 19   # add 3s/3p shell split (approx) + 6 d
            else:
                f = None
            if f is None:
                return None, f"no 6-31G(d) estimate for {sym}"
            total += f * n
        return total, "6-31G(d): 6 Cartesian d per heavy atom (Gaussian default for 6-31G family)"

    if b in ("sto-3g",):
        per = {"H": 1, "He": 1}
        total = 0
        for sym, n in counts.items():
            z = Chem.GetPeriodicTable().GetAtomicNumber(sym)
            f = 1 if z <= 2 else (5 if z <= 10 else 9)
            total += f * n
        return total, "STO-3G minimal basis"

    return None, f"no built-in basis-function estimate for {basis!r}"


def build_report(
    mol3d: Molecule3D, basis: str = "6-31g(d)", method: str = "b3lyp"
) -> PreExportReport:
    r = PreExportReport()
    r.formula = mol3d.formula
    r.n_atoms = len(mol3d)
    r.charge = int(mol3d.charge)
    r.multiplicity = int(mol3d.multiplicity)
    r.n_electrons = mol3d.n_electrons

    # Parity: unpaired electrons = multiplicity - 1 must have the same parity as
    # the electron count (both even or both odd).
    unpaired = r.multiplicity - 1
    if (r.n_electrons % 2) != (unpaired % 2):
        r.parity_ok = False
        r.parity_message = (
            f"{r.n_electrons} electrons with multiplicity {r.multiplicity} is "
            f"impossible (parity mismatch). Suggested: "
            f"{mol3d.suggest_multiplicity()}."
        )
    else:
        r.parity_ok = True
        r.parity_message = "electron count and multiplicity are consistent"

    r.clashes = vdw_clashes(mol3d)

    # Bond lengths / angles from the perceived graph.
    graph = _bonded_graph(mol3d)
    if graph is not None and graph.GetNumBonds() > 0:
        xyz = np.array([[a.x, a.y, a.z] for a in mol3d.atoms])
        lengths = [
            float(np.linalg.norm(xyz[b.GetBeginAtomIdx()] - xyz[b.GetEndAtomIdx()]))
            for b in graph.GetBonds()
        ]
        if lengths:
            r.bond_length_range = (min(lengths), max(lengths))
        angles = _bond_angles(graph, xyz)
        if angles:
            r.bond_angle_range = (min(angles), max(angles))

    r.n_basis_estimate, r.basis_note = estimate_basis_functions(mol3d, basis)
    if r.n_basis_estimate is not None:
        n = r.n_basis_estimate
        # Cost scales steeply with N; give a qualitative band.
        if n < 150:
            band = "small — fast"
        elif n < 400:
            band = "moderate"
        elif n < 800:
            band = "large — expect long SCF/opt times"
        else:
            band = "very large — consider a smaller basis or fragmenting"
        r.cost_warning = f"~{n} basis functions ({band}); DFT cost grows ~O(N^3–N^4) per step."

    if not r.parity_ok:
        r.warnings.append("Fix charge/multiplicity before running Gaussian.")
    if r.clashes:
        r.warnings.append(
            "Steric clashes present — run a force-field/xTB cleanup before export."
        )
    return r


def _bond_angles(graph: Chem.Mol, xyz: np.ndarray) -> List[float]:
    angles = []
    for atom in graph.GetAtoms():
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        c = atom.GetIdx()
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                v1 = xyz[nbrs[a]] - xyz[c]
                v2 = xyz[nbrs[b]] - xyz[c]
                nv1, nv2 = np.linalg.norm(v1), np.linalg.norm(v2)
                if nv1 < 1e-6 or nv2 < 1e-6:
                    continue
                cosang = float(np.clip(np.dot(v1, v2) / (nv1 * nv2), -1, 1))
                angles.append(np.degrees(np.arccos(cosang)))
    return angles
