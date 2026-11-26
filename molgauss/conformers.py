"""Backbone conformation control for chains (polymers and peptides).

Random ETKDG conformers of long chains come out as compact coils. Those are the
worst possible starting points for a DFT optimizer: many soft dihedrals, lots of
weak non-bonded contacts, and no exploitable symmetry. This module rewrites the
main-chain dihedrals to a regular, extended (or helical) geometry so Gaussian
starts near a sensible minimum and converges in far fewer steps.

Two use cases:

* **Polymers / linear molecules** – set every rotatable backbone dihedral to
  180 deg (all-trans, fully extended) or a repeating trans-gauche-trans pattern
  (a crude helix, useful for poly-ethers such as PEG).
* **Peptides** – set each residue's (phi, psi) to a secondary-structure target:
  extended beta strand (phi=-135, psi=135) or alpha helix (phi=-57, psi=-47).
"""

from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import rdMolTransforms

# Secondary-structure backbone dihedral targets (degrees).
SS_TARGETS = {
    "extended": (-135.0, 135.0),   # beta strand
    "alpha": (-57.0, -47.0),       # right-handed alpha helix
    "polypro": (-75.0, 150.0),     # polyproline-II-like
}

# Repeating torsion patterns for generic polymer backbones (degrees).
CHAIN_PATTERNS = {
    "all-trans": [180.0],
    "tgt": [180.0, 60.0, 180.0, -60.0],  # trans-gauche, PEG-like helix
}


def _longest_heavy_chain(mol: Chem.Mol) -> List[int]:
    """Return atom indices of the longest path through heavy atoms.

    Used to identify a polymer's main chain. Finds the two most distant heavy
    atoms (graph diameter endpoints) and returns the shortest path between them,
    which for a linear polymer is the backbone.
    """
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    if len(heavy) < 2:
        return heavy
    # Double BFS (diameter): farthest atom from an arbitrary start, then the
    # farthest from that one; the path between them is the longest chain.
    dmat = Chem.GetDistanceMatrix(mol)
    best_pair = (heavy[0], heavy[0])
    best = -1
    for i in heavy:
        for j in heavy:
            d = dmat[i][j]
            if d > best:
                best, best_pair = d, (i, j)
    path = Chem.GetShortestPath(mol, int(best_pair[0]), int(best_pair[1]))
    return [idx for idx in path if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1]


def _backbone_dihedral_quads(chain: Sequence[int]) -> List[Tuple[int, int, int, int]]:
    """Consecutive 4-atom tuples defining each backbone torsion along a chain."""
    quads = []
    for k in range(len(chain) - 3):
        quads.append((chain[k], chain[k + 1], chain[k + 2], chain[k + 3]))
    return quads


def set_chain_extended(
    mol: Chem.Mol, pattern: str = "all-trans", conf_id: int = -1
) -> int:
    """Set a linear molecule's backbone dihedrals to a regular pattern.

    ``pattern`` is a key of ``CHAIN_PATTERNS`` (``all-trans`` or ``tgt``).
    Returns the number of dihedrals rewritten. Requires an existing conformer.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no conformer; embed it before setting torsions")
    angles = CHAIN_PATTERNS.get(pattern)
    if angles is None:
        raise ValueError(f"Unknown chain pattern {pattern!r}; use {list(CHAIN_PATTERNS)}")

    chain = _longest_heavy_chain(mol)
    if len(chain) < 4:
        return 0
    conf = mol.GetConformer(conf_id)
    quads = _backbone_dihedral_quads(chain)
    count = 0
    for k, (a, b, c, d) in enumerate(quads):
        # Only rotate about rotatable single bonds (skip rings / double bonds).
        bond = mol.GetBondBetweenAtoms(b, c)
        if bond is None or bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        target = angles[k % len(angles)]
        try:
            rdMolTransforms.SetDihedralDeg(conf, a, b, c, d, target)
            count += 1
        except Exception:
            continue
    return count


def _find_peptide_backbone(mol: Chem.Mol) -> List[Tuple[int, int, int]]:
    """Return per-residue (N, CA, C) atom-index triples in sequence order.

    Relies on PDB residue info that RDKit attaches to peptides built with
    ``MolFromSequence``. Residues are ordered by residue number.
    """
    residues: dict = {}
    for atom in mol.GetAtoms():
        info = atom.GetPDBResidueInfo()
        if info is None:
            continue
        name = info.GetName().strip()
        resnum = info.GetResidueNumber()
        if name in ("N", "CA", "C"):
            residues.setdefault(resnum, {})[name] = atom.GetIdx()
    triples = []
    for resnum in sorted(residues):
        r = residues[resnum]
        if {"N", "CA", "C"} <= set(r):
            triples.append((r["N"], r["CA"], r["C"]))
    return triples


def set_peptide_secondary_structure(
    mol: Chem.Mol, structure: str = "extended", conf_id: int = -1
) -> int:
    """Set every residue's (phi, psi) backbone dihedrals to a target structure.

    phi = C(i-1)-N(i)-CA(i)-C(i); psi = N(i)-CA(i)-C(i)-N(i+1).
    ``structure`` is a key of ``SS_TARGETS``. Returns residues modified.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("Peptide has no conformer; embed it before setting torsions")
    target = SS_TARGETS.get(structure)
    if target is None:
        raise ValueError(f"Unknown structure {structure!r}; use {list(SS_TARGETS)}")
    phi_t, psi_t = target

    bb = _find_peptide_backbone(mol)
    if len(bb) < 2:
        return 0
    conf = mol.GetConformer(conf_id)
    modified = 0
    for i, (n, ca, c) in enumerate(bb):
        # phi needs the previous residue's carbonyl C.
        if i > 0:
            c_prev = bb[i - 1][2]
            try:
                rdMolTransforms.SetDihedralDeg(conf, c_prev, n, ca, c, phi_t)
            except Exception:
                pass
        # psi needs the next residue's amide N.
        if i < len(bb) - 1:
            n_next = bb[i + 1][0]
            try:
                rdMolTransforms.SetDihedralDeg(conf, n, ca, c, n_next, psi_t)
            except Exception:
                pass
        modified += 1
    return modified
