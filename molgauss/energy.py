"""Semi-empirical energy and optimization backend (optional).

GFN2-xTB is used to rank and relax conformers because it is far closer to DFT
than a classical force field while remaining cheap. It is an OPTIONAL dependency
(``tblite`` + ``ase``); when unavailable, callers fall back to MMFF energies.

IMPORTANT: GFN2-xTB is NOT DFT. Relative energies (and therefore the ordering of
conformers) can change at the DFT level. The UI and reports say so explicitly.

Backend choice, verified against the libraries actually installed:
  * ``tblite`` provides single-point energies/gradients only (no optimizer).
  * Geometry optimization is driven by ASE's optimizer (BFGS) using the
    ``tblite.ase.TBLite`` calculator. This is the documented, supported path.
"""

from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from .molmodel import Molecule3D

EV2KCAL = 23.060547830619026  # eV -> kcal/mol


def xtb_available() -> bool:
    """True if GFN2-xTB (tblite + ase) can be used."""
    try:
        import tblite.ase  # noqa: F401
        import ase  # noqa: F401
        return True
    except Exception:
        return False


def _to_ase(mol: Molecule3D):
    from ase import Atoms

    return Atoms(
        numbers=[a.atomic_number for a in mol.atoms],
        positions=[[a.x, a.y, a.z] for a in mol.atoms],
    )


def xtb_singlepoint_eV(mol: Molecule3D, charge: int = 0) -> Optional[float]:
    """GFN2-xTB single-point energy in eV, or None if the backend is missing."""
    if not xtb_available():
        return None
    from tblite.ase import TBLite

    at = _to_ase(mol)
    at.calc = TBLite(method="GFN2-xTB", charge=int(charge), verbosity=0)
    return float(at.get_potential_energy())


def xtb_optimize(
    mol: Molecule3D, charge: int = 0, fmax: float = 0.05, steps: int = 500
) -> Optional[Tuple[float, int]]:
    """Relax ``mol`` in place with GFN2-xTB (ASE BFGS).

    Returns ``(final_energy_eV, n_steps)`` or None if the backend is missing.
    ``fmax`` is the force convergence criterion in eV/A; ``n_steps`` is the
    optimizer step count, used as a reproducible proxy for optimization effort.
    """
    if not xtb_available():
        return None
    from tblite.ase import TBLite
    from ase.optimize import BFGS

    at = _to_ase(mol)
    at.calc = TBLite(method="GFN2-xTB", charge=int(charge), verbosity=0)
    opt = BFGS(at, logfile=None)
    opt.run(fmax=fmax, steps=steps)
    energy = float(at.get_potential_energy())
    nsteps = int(opt.get_number_of_steps())
    pos = at.get_positions()
    for a, p in zip(mol.atoms, pos):
        a.x, a.y, a.z = float(p[0]), float(p[1]), float(p[2])
    return energy, nsteps
