"""3D coordinate generation and force-field pre-optimization.

Primary engine is RDKit: explicit hydrogens, ETKDGv3 embedding, then MMFF94
(organics) or UFF (anything with metals / elements MMFF cannot handle). Open
Babel is used as a fallback when RDKit fails to embed (common for transition
metal complexes and some inorganics), and as an independent optimizer.
"""

from __future__ import annotations

from typing import Optional, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from .molmodel import Molecule3D

# Elements MMFF94 does not parameterize well; force UFF when any is present.
# (MMFF covers common organic main-group elements; metals need UFF.)
_MMFF_SAFE = set(
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Br I".split()
)


class GeometryError(RuntimeError):
    """Raised when a 3D geometry cannot be produced by any available engine."""


def _needs_uff(mol: Chem.Mol) -> bool:
    syms = {a.GetSymbol() for a in mol.GetAtoms()}
    return not syms.issubset(_MMFF_SAFE)


def _run_ff(mol: Chem.Mol, force_field: str, max_iters: int) -> str:
    """Run MMFF94 or UFF on an already-embedded mol. Returns an info string."""
    chosen = force_field
    if force_field == "AUTO":
        chosen = "UFF" if _needs_uff(mol) else "MMFF94"
    info = chosen
    try:
        if chosen == "MMFF94" and AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=max_iters)
        else:
            if chosen == "MMFF94":
                info = "UFF (MMFF params missing)"
            AllChem.UFFOptimizeMolecule(mol, maxIters=max_iters)
    except Exception as exc:  # pragma: no cover - defensive
        info += f" (optimization warning: {exc})"
    return info


def _apply_conformation(mol: Chem.Mol, conformation: str) -> str:
    """Rewrite backbone dihedrals to an extended/helical target.

    ``conformation`` is one of:
      chain:  "all-trans", "tgt"
      peptide: "extended", "alpha", "polypro"
    Returns a short description of what was applied (or "" if nothing).
    """
    from . import conformers

    if conformation in conformers.CHAIN_PATTERNS:
        n = conformers.set_chain_extended(mol, pattern=conformation)
        return f"chain {conformation} ({n} dihedrals)" if n else ""
    if conformation in conformers.SS_TARGETS:
        n = conformers.set_peptide_secondary_structure(mol, structure=conformation)
        return f"peptide {conformation} ({n} residues)" if n else ""
    return ""


def _rdkit_embed_optimize(
    mol: Chem.Mol,
    force_field: str,
    max_iters: int,
    random_seed: int,
    conformation: Optional[str] = None,
) -> Tuple[Chem.Mol, str]:
    """Embed with ETKDGv3, optionally impose a chain conformation, and optimize.

    For chains, the sequence is embed -> short FF cleanup -> set backbone
    dihedrals (extended/helix) -> final FF polish. Imposing the conformation
    after a first cleanup, then polishing, keeps the regular backbone while
    relaxing only local clashes, which is exactly the Gaussian-friendly start we
    want. Returns (mol, info-string).
    """
    mol = Chem.AddHs(mol, addCoords=True)

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.useSmallRingTorsions = True

    status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        # Retry with random coordinates (helps strained / macrocyclic systems).
        params.useRandomCoords = True
        status = AllChem.EmbedMolecule(mol, params)
    if status != 0:
        raise GeometryError("RDKit ETKDGv3 embedding failed")

    info = "RDKit ETKDGv3 + " + _run_ff(mol, force_field, max_iters)

    if conformation:
        desc = _apply_conformation(mol, conformation)
        if desc:
            # Polish with a few iterations so the imposed geometry relaxes
            # locally without collapsing the regular backbone.
            _run_ff(mol, force_field, max_iters=max(100, max_iters // 10))
            info += f" | {desc}"
    return mol, info


def _openbabel_embed_optimize(
    smiles: str,
    force_field: str,
    steps: int,
) -> Tuple[Molecule3D, str]:
    """Open Babel fallback: build 3D from SMILES and optimize.

    Used when RDKit cannot embed (e.g. many coordination complexes).
    """
    from openbabel import pybel

    mol = pybel.readstring("smi", smiles)
    mol.addh()
    ff = "uff" if force_field in ("UFF", "AUTO") else "mmff94"
    try:
        mol.make3D(forcefield=ff, steps=max(50, steps // 10))
        mol.localopt(forcefield=ff, steps=steps)
    except Exception:
        # Open Babel falls back to UFF internally if mmff94 is unavailable.
        mol.make3D(forcefield="uff", steps=max(50, steps // 10))
        mol.localopt(forcefield="uff", steps=steps)

    from .molmodel import Atom

    atoms = []
    for atom in mol.atoms:
        x, y, z = atom.coords
        sym = _obabel_symbol(atom.atomicnum)
        atoms.append(Atom(sym, x, y, z))
    charge = mol.charge
    out = Molecule3D(atoms=atoms, charge=int(charge), multiplicity=1)
    out.multiplicity = out.suggest_multiplicity()
    return out, f"Open Babel make3D + {ff.upper()}"


def _obabel_symbol(z: int) -> str:
    from .molmodel import atomic_number_to_symbol

    return atomic_number_to_symbol(z)


def embed_and_optimize(
    mol: Chem.Mol,
    force_field: str = "AUTO",
    max_iters: int = 2000,
    random_seed: int = 0xF00D,
    charge: Optional[int] = None,
    multiplicity: Optional[int] = None,
    title: str = "Molecule generated by MolGauss",
    name: str = "molecule",
    allow_openbabel_fallback: bool = True,
    conformation: Optional[str] = None,
    align: bool = True,
    clash_cutoff: float = 0.8,
) -> Tuple[Molecule3D, str]:
    """Generate a 3D, force-field-minimized ``Molecule3D`` from an RDKit mol.

    Parameters
    ----------
    force_field : "AUTO" | "MMFF94" | "UFF"
        AUTO picks MMFF94 for organics and UFF when metals/heavy elements appear.
    allow_openbabel_fallback :
        If RDKit embedding fails, retry through Open Babel using the molecule's
        canonical SMILES.
    conformation :
        Optional backbone conformation to impose on chains before the final
        polish: ``"all-trans"`` / ``"tgt"`` for polymers, ``"extended"`` /
        ``"alpha"`` / ``"polypro"`` for peptides.
    align :
        Run PCA alignment (COM to origin, longest axis on X) after optimization.
    clash_cutoff :
        Non-bonded distance (angstrom) below which a steric clash is reported.

    Returns ``(Molecule3D, method_description)``. The description includes the
    refinement summary (alignment, clashes, bond sanity).
    """
    from .refine import refine as _refine_fn

    def _finish(result: Molecule3D, info: str) -> Tuple[Molecule3D, str]:
        report = _refine_fn(result, clash_cutoff=clash_cutoff, align=align)
        return result, f"{info} | {report.summary()}"

    try:
        opt_mol, info = _rdkit_embed_optimize(
            mol, force_field, max_iters, random_seed, conformation=conformation
        )
        result = Molecule3D.from_rdkit(
            opt_mol, charge=charge, multiplicity=multiplicity, title=title, name=name
        )
        return _finish(result, info)
    except GeometryError:
        if not allow_openbabel_fallback:
            raise
        try:
            smiles = Chem.MolToSmiles(mol)
        except Exception as exc:
            raise GeometryError(
                "RDKit embedding failed and molecule could not be serialized "
                "to SMILES for the Open Babel fallback"
            ) from exc
        result, info = _openbabel_embed_optimize(smiles, force_field, max_iters)
        result.title, result.name = title, name
        if charge is not None:
            result.charge = int(charge)
        if multiplicity is not None:
            result.multiplicity = int(multiplicity)
        return _finish(result, info)
