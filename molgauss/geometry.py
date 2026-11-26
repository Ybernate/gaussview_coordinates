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


# ---------------------------------------------------------------------------
# Energy-selected conformer search (P1-1)
# ---------------------------------------------------------------------------
from dataclasses import dataclass, field  # noqa: E402


@dataclass
class ConformerCandidate:
    """One relaxed conformer with its energy and provenance."""

    source: str                      # how it was generated
    energy: float                    # kcal/mol relative to the best (filled later)
    abs_energy: float                # absolute energy (eV if xTB, FF units if MMFF)
    method: str                      # "GFN2-xTB" or "MMFF94"
    opt_steps: Optional[int]         # xTB/ASE optimizer steps (None for MMFF)
    mol: Molecule3D = field(repr=False, default=None)


@dataclass
class ConformerSearchResult:
    best: Molecule3D
    candidates: list                 # List[ConformerCandidate], energy-sorted
    method: str                      # ranking method actually used
    used_xtb: bool
    note: str = ""

    def table(self) -> str:
        lines = [f"{'rank':>4}  {'source':<22}  {'ΔE(kcal/mol)':>12}  {'steps':>6}"]
        for i, c in enumerate(self.candidates):
            steps = "" if c.opt_steps is None else str(c.opt_steps)
            lines.append(f"{i:>4}  {c.source:<22}  {c.energy:>12.2f}  {steps:>6}")
        return "\n".join(lines)


def _embed_variants(
    mol: Chem.Mol,
    chain_families: Optional[list],
    peptide_structures: Optional[list],
    n_random: int,
    force_field: str,
    max_iters: int,
    base_seed: int,
) -> list:
    """Generate labelled (source, Chem.Mol-with-conformer) candidates.

    Candidates = one ETKDG embed per type-aware chain family / peptide structure
    (each then FF-relaxed), plus ``n_random`` random ETKDG conformers. Nothing
    here assumes which is best; the caller ranks by energy.
    """
    from . import conformers as _conf

    variants = []

    def _fresh():
        m = Chem.AddHs(mol, addCoords=False)
        p = AllChem.ETKDGv3()
        p.randomSeed = base_seed
        if AllChem.EmbedMolecule(m, p) != 0:
            p.useRandomCoords = True
            if AllChem.EmbedMolecule(m, p) != 0:
                return None
        return m

    for fam in (chain_families or []):
        m = _fresh()
        if m is None:
            continue
        _run_ff(m, force_field, max_iters)          # relax embed
        try:
            _conf.set_chain_conformation(m, family=fam)
        except Exception:
            continue
        _run_ff(m, force_field, max(100, max_iters // 10))  # polish
        variants.append((f"chain:{fam}", m))

    for ss in (peptide_structures or []):
        m = _fresh()
        if m is None:
            continue
        _run_ff(m, force_field, max_iters)
        try:
            _conf.set_peptide_secondary_structure(m, structure=ss)
        except Exception:
            continue
        _run_ff(m, force_field, max(100, max_iters // 10))
        variants.append((f"peptide:{ss}", m))

    for s in range(n_random):
        m = Chem.AddHs(mol, addCoords=False)
        p = AllChem.ETKDGv3()
        p.randomSeed = base_seed + 1 + s
        if AllChem.EmbedMolecule(m, p) != 0:
            continue
        _run_ff(m, force_field, max_iters)
        variants.append((f"random:seed={base_seed + 1 + s}", m))

    return variants


def search_conformers(
    mol: Chem.Mol,
    chain_families: Optional[list] = None,
    peptide_structures: Optional[list] = None,
    n_random: int = 8,
    force_field: str = "AUTO",
    max_iters: int = 2000,
    base_seed: int = 0xF00D,
    charge: Optional[int] = None,
    use_xtb: bool = True,
    xtb_fmax: float = 0.05,
    name: str = "molecule",
    title: str = "Molecule generated by MolGauss",
) -> ConformerSearchResult:
    """Generate candidate conformers and return the lowest-energy one.

    Ranking uses GFN2-xTB (optimize + energy) when available, otherwise MMFF
    energies. Candidates come from the type-aware chain families / peptide
    structures requested, plus ``n_random`` random ETKDG conformers. The result
    carries the full ranked table and flags whether xTB was used.
    """
    from .energy import xtb_available, xtb_optimize, xtb_singlepoint_eV, EV2KCAL

    q = 0 if charge is None else int(charge)
    variants = _embed_variants(
        mol, chain_families, peptide_structures, n_random,
        force_field, max_iters, base_seed,
    )
    if not variants:
        raise GeometryError("No conformer could be embedded for the search")

    want_xtb = use_xtb and xtb_available()
    cands: list = []

    for source, m in variants:
        m3d = Molecule3D.from_rdkit(m, charge=charge, name=name, title=title)
        if want_xtb:
            res = xtb_optimize(m3d, charge=q, fmax=xtb_fmax)
            if res is not None:
                abs_e, nsteps = res
                cands.append(ConformerCandidate(source, 0.0, abs_e, "GFN2-xTB", nsteps, m3d))
                continue
        # MMFF fallback energy (relative, kcal/mol within MMFF).
        ff_e = _mmff_energy(m)
        cands.append(ConformerCandidate(source, 0.0, ff_e, "MMFF94", None, m3d))

    method = "GFN2-xTB" if want_xtb and any(c.method == "GFN2-xTB" for c in cands) else "MMFF94"
    cands.sort(key=lambda c: c.abs_energy)
    base = cands[0].abs_energy
    scale = EV2KCAL if method == "GFN2-xTB" else 1.0  # MMFF energies already kcal/mol
    for c in cands:
        c.energy = (c.abs_energy - base) * scale

    note = (
        "Ranked by GFN2-xTB (NOT DFT; relative order may change at DFT level)."
        if method == "GFN2-xTB"
        else "xTB unavailable -> ranked by MMFF energy (less reliable; install "
        "tblite+ase for GFN2-xTB)."
    )
    return ConformerSearchResult(
        best=cands[0].mol, candidates=cands, method=method,
        used_xtb=(method == "GFN2-xTB"), note=note,
    )


def _mmff_energy(mol: Chem.Mol) -> float:
    """MMFF energy (kcal/mol) of an embedded mol; large value if unavailable."""
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol)
        if props is None:
            return 1e6
        ff = AllChem.MMFFGetMoleculeForceField(mol, props)
        return float(ff.CalcEnergy()) if ff is not None else 1e6
    except Exception:
        return 1e6
