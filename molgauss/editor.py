"""Atom-level CRUD operations on a ``Molecule3D``.

These operate on the explicit atom list, which is exactly what the editable
table in the UI exposes. They never try to re-perceive bonds; after a structural
edit the caller may optionally re-run a force-field cleanup (see
``cleanup_geometry``) to relax clashes introduced by the change.
"""

from __future__ import annotations

from typing import List, Optional

from .molmodel import Atom, Molecule3D, symbol_to_atomic_number


def add_atom(
    mol: Molecule3D, symbol: str, x: float, y: float, z: float
) -> Molecule3D:
    """Append a new atom. Returns the same (mutated) molecule for chaining."""
    mol.atoms.append(Atom(symbol, x, y, z))
    return mol


def delete_atom(mol: Molecule3D, index: int) -> Molecule3D:
    """Delete the atom at ``index`` (0-based)."""
    if not 0 <= index < len(mol.atoms):
        raise IndexError(f"Atom index {index} out of range (0..{len(mol.atoms)-1})")
    del mol.atoms[index]
    return mol


def delete_atoms(mol: Molecule3D, indices: List[int]) -> Molecule3D:
    """Delete multiple atoms by index (e.g. remove a water or a group)."""
    for i in sorted(set(indices), reverse=True):
        delete_atom(mol, i)
    return mol


def substitute_element(mol: Molecule3D, index: int, new_symbol: str) -> Molecule3D:
    """Change the element of one atom in place (e.g. C -> Si), keeping its
    coordinates. Charge/multiplicity are NOT auto-adjusted because the correct
    values depend on chemical intent; the UI surfaces them for the user."""
    if not 0 <= index < len(mol.atoms):
        raise IndexError(f"Atom index {index} out of range (0..{len(mol.atoms)-1})")
    symbol_to_atomic_number(new_symbol)  # validate
    mol.atoms[index].symbol = new_symbol.strip().capitalize()
    return mol


def update_coordinates(
    mol: Molecule3D, index: int, x: float, y: float, z: float
) -> Molecule3D:
    if not 0 <= index < len(mol.atoms):
        raise IndexError(f"Atom index {index} out of range (0..{len(mol.atoms)-1})")
    a = mol.atoms[index]
    a.x, a.y, a.z = float(x), float(y), float(z)
    return mol


def replace_atoms(mol: Molecule3D, atoms: List[Atom]) -> Molecule3D:
    """Replace the entire atom list (used when the UI table is edited wholesale)."""
    mol.atoms = list(atoms)
    return mol


def set_frozen(mol: Molecule3D, index: int, frozen: bool = True) -> Molecule3D:
    """Freeze (``-1`` in Gaussian) or unfreeze a single atom for optimization."""
    if not 0 <= index < len(mol.atoms):
        raise IndexError(f"Atom index {index} out of range (0..{len(mol.atoms)-1})")
    mol.atoms[index].frozen = bool(frozen)
    return mol


def freeze_all(mol: Molecule3D, frozen: bool = True) -> Molecule3D:
    """Freeze or unfreeze every atom."""
    for a in mol.atoms:
        a.frozen = bool(frozen)
    return mol


def cleanup_geometry(
    mol: Molecule3D, force_field: str = "AUTO", max_iters: int = 500
) -> Optional[str]:
    """Relax the current coordinates with a force field after manual edits.

    Rebuilds connectivity from the 3D geometry (via Open Babel) and runs a short
    UFF/MMFF minimization. Returns a description string, or ``None`` if the
    geometry could not be cleaned (coordinates are then left untouched).
    """
    try:
        from openbabel import pybel
        from openbabel import openbabel as ob
    except ImportError:
        return None

    obmol = ob.OBMol()
    for a in mol.atoms:
        atom = obmol.NewAtom()
        atom.SetAtomicNum(a.atomic_number)
        atom.SetVector(a.x, a.y, a.z)
    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()

    pm = pybel.Molecule(obmol)
    ff = "uff" if force_field in ("UFF", "AUTO") else "mmff94"
    try:
        pm.localopt(forcefield=ff, steps=max_iters)
    except Exception:
        try:
            pm.localopt(forcefield="uff", steps=max_iters)
            ff = "uff"
        except Exception:
            return None

    for a, atom in zip(mol.atoms, pm.atoms):
        a.x, a.y, a.z = atom.coords
    return f"Open Babel {ff.upper()} cleanup ({max_iters} steps)"
