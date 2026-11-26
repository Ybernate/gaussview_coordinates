"""Multi-format molecular input parsing.

Every public function returns an RDKit ``Chem.Mol`` (without 3D coordinates,
unless a coordinate-bearing format such as PDB/XYZ was supplied) OR a finished
``Molecule3D`` when the source already carries geometry. The geometry module
then handles embedding/optimization for graph-only inputs.

Supported sources
-----------------
SMILES, InChI, common/IUPAC name (via PubChem), polymer from monomer SMILES,
FASTA peptide sequence, PDB/CIF block, and XYZ block.
"""

from __future__ import annotations

import re
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

from .molmodel import Atom, Molecule3D, symbol_to_atomic_number


class InputError(ValueError):
    """Raised for malformed or unresolvable molecular input."""


# --------------------------------------------------------------------- SMILES
def from_smiles(smiles: str) -> Chem.Mol:
    smiles = (smiles or "").strip()
    if not smiles:
        raise InputError("Empty SMILES string")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InputError(f"RDKit could not parse the SMILES: {smiles!r}")
    return mol


# ---------------------------------------------------------------------- InChI
def from_inchi(inchi: str) -> Chem.Mol:
    inchi = (inchi or "").strip()
    if not inchi:
        raise InputError("Empty InChI string")
    if not inchi.lower().startswith("inchi="):
        inchi = "InChI=" + inchi
    mol = Chem.MolFromInchi(inchi, sanitize=True, removeHs=False)
    if mol is None:
        raise InputError(f"RDKit could not parse the InChI: {inchi!r}")
    return mol


# ----------------------------------------------------- name / IUPAC (PubChem)
def from_name(name: str, timeout: int = 20) -> Tuple[Chem.Mol, str]:
    """Resolve a common or IUPAC name to a mol via PubChem.

    Returns ``(mol, resolved_smiles)``. Requires network access.
    """
    name = (name or "").strip()
    if not name:
        raise InputError("Empty compound name")
    try:
        import pubchempy as pcp
    except ImportError as exc:  # pragma: no cover
        raise InputError("pubchempy is not installed") from exc

    try:
        compounds = pcp.get_compounds(name, "name")
    except Exception as exc:
        raise InputError(
            f"PubChem lookup failed for {name!r} (network issue?): {exc}"
        ) from exc
    if not compounds:
        raise InputError(f"No PubChem compound found for name {name!r}")

    smiles = None
    for attr in ("smiles", "isomeric_smiles", "canonical_smiles", "connectivity_smiles"):
        try:
            smiles = getattr(compounds[0], attr, None)
        except Exception:
            smiles = None
        if smiles:
            break
    if not smiles:
        raise InputError(f"PubChem returned no SMILES for {name!r}")
    return from_smiles(smiles), smiles


# -------------------------------------------------------------------- polymer
def build_polymer(monomer_smiles: str, n_units: int, cap: str = "H") -> Chem.Mol:
    """Construct a linear polymer from a monomer with two attachment points.

    The monomer SMILES must mark its two connection points with ``[*]`` (dummy
    atoms), e.g. polyethylene is ``[*]CC[*]``. ``n_units`` copies are joined
    head-to-tail and the dangling ends are capped (default: hydrogen).
    """
    monomer_smiles = (monomer_smiles or "").strip()
    if n_units < 1:
        raise InputError("Degree of polymerization must be >= 1")

    unit = Chem.MolFromSmiles(monomer_smiles)
    if unit is None:
        raise InputError(f"Could not parse monomer SMILES: {monomer_smiles!r}")

    dummies = [a.GetIdx() for a in unit.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 2:
        raise InputError(
            "Monomer must contain exactly two attachment points marked as [*]. "
            f"Found {len(dummies)}. Example: [*]CC[*] for polyethylene."
        )

    # Build the chain by repeatedly joining copies with RDKit's RWMol, bonding
    # the second dummy's neighbor of one unit to the first dummy's neighbor of
    # the next, then stripping the consumed dummy atoms.
    chain = Chem.RWMol(Chem.MolFromSmiles(monomer_smiles))
    for _ in range(n_units - 1):
        chain = _join_unit(chain, monomer_smiles)

    # Replace remaining dummy atoms with the cap element.
    capped = _cap_dummies(chain, cap)
    Chem.SanitizeMol(capped)
    return capped


def _attachment_pairs(mol: Chem.Mol) -> List[Tuple[int, int]]:
    """Return (dummy_idx, neighbor_idx) pairs for every dummy atom."""
    pairs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            nbrs = atom.GetNeighbors()
            if not nbrs:
                raise InputError("Attachment point [*] has no neighbor atom")
            pairs.append((atom.GetIdx(), nbrs[0].GetIdx()))
    return pairs


def _join_unit(chain: Chem.RWMol, monomer_smiles: str) -> Chem.RWMol:
    """Append one monomer unit to the growing chain at its trailing dummy."""
    unit = Chem.MolFromSmiles(monomer_smiles)
    combined = Chem.RWMol(Chem.CombineMols(chain.GetMol(), unit))

    pairs = _attachment_pairs(combined)
    # Chain's trailing dummy = last dummy that belongs to the original chain;
    # with CombineMols, chain atoms keep their indices and unit atoms follow.
    n_chain = chain.GetNumAtoms()
    chain_dummies = [p for p in pairs if p[0] < n_chain]
    unit_dummies = [p for p in pairs if p[0] >= n_chain]

    # Join chain's last dummy to unit's first dummy.
    (d1, n1) = chain_dummies[-1]
    (d2, n2) = unit_dummies[0]
    combined.AddBond(n1, n2, Chem.BondType.SINGLE)

    # Remove the two consumed dummies (highest index first to keep indices valid).
    for d in sorted([d1, d2], reverse=True):
        combined.RemoveAtom(d)
    return combined


def _cap_dummies(mol: Chem.RWMol, cap: str) -> Chem.Mol:
    cap_z = symbol_to_atomic_number(cap)
    rw = Chem.RWMol(mol)
    dummies = [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0]
    for idx in dummies:
        rw.GetAtomWithIdx(idx).SetAtomicNum(cap_z)
    return rw.GetMol()


# ------------------------------------------------------------- FASTA / peptide
# One-letter code -> three-letter residue name for RDKit's sequence reader.
_AA1 = set("ACDEFGHIKLMNPQRSTVWY")


def from_fasta(sequence: str) -> Chem.Mol:
    """Build a peptide from a one-letter amino-acid sequence (FASTA body).

    Header lines beginning with ``>`` are ignored. RDKit builds the backbone
    with standard residues; geometry is embedded downstream.
    """
    seq = "".join(
        line.strip()
        for line in (sequence or "").splitlines()
        if not line.strip().startswith(">")
    ).upper().replace(" ", "")
    if not seq:
        raise InputError("Empty FASTA sequence")
    bad = sorted(set(seq) - _AA1)
    if bad:
        raise InputError(
            "Unsupported amino-acid one-letter codes: "
            + ", ".join(bad)
            + ". Use the 20 standard residues (ACDEFGHIKLMNPQRSTVWY)."
        )
    mol = Chem.MolFromSequence(seq)
    if mol is None:
        raise InputError(f"RDKit could not build a peptide from sequence {seq!r}")
    return mol


# ----------------------------------------------------------------- PDB / CIF
def from_pdb_block(block: str, is_cif: bool = False) -> Molecule3D:
    """Parse a PDB (or mmCIF) block that already contains 3D coordinates.

    Returns a finished ``Molecule3D`` (coordinates preserved as given); no
    re-embedding is done so experimental geometries are respected.
    """
    block = block or ""
    if is_cif:
        mol = Chem.MolFromMmCIFBlock(block, sanitize=False, removeHs=False) \
            if hasattr(Chem, "MolFromMmCIFBlock") else None
        if mol is None:
            raise InputError(
                "Could not parse the mmCIF block (RDKit mmCIF support may be "
                "unavailable in this build). Try converting to PDB first."
            )
    else:
        mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
        if mol is None:
            raise InputError("Could not parse the PDB block")

    if mol.GetNumConformers() == 0:
        raise InputError("PDB/CIF block contains no 3D coordinates")
    Chem.SanitizeMol(mol, catchErrors=True)
    return Molecule3D.from_rdkit(mol, title="Structure from PDB/CIF", name="structure")


# ----------------------------------------------------------------------- XYZ
_XYZ_LINE = re.compile(
    r"^\s*([A-Za-z]{1,2})\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)"
)


def from_xyz_block(block: str) -> Molecule3D:
    """Parse a standard XYZ block into a ``Molecule3D`` (coordinates preserved)."""
    lines = [ln for ln in (block or "").splitlines()]
    if not lines:
        raise InputError("Empty XYZ block")

    # A well-formed XYZ has: count, comment, then atom lines. Be lenient and
    # accept files that start directly with atom lines.
    atom_lines = lines
    if lines[0].strip().isdigit():
        atom_lines = lines[2:]

    atoms: List[Atom] = []
    for ln in atom_lines:
        if not ln.strip():
            continue
        m = _XYZ_LINE.match(ln)
        if not m:
            raise InputError(f"Malformed XYZ atom line: {ln!r}")
        sym, x, y, z = m.groups()
        atoms.append(Atom(sym, float(x), float(y), float(z)))
    if not atoms:
        raise InputError("No atoms found in XYZ block")
    mol = Molecule3D(atoms=atoms, title="Structure from XYZ", name="structure")
    mol.multiplicity = mol.suggest_multiplicity()
    return mol


def xyz_perceive_bonds(block: str, charge: int = 0) -> Tuple[Chem.Mol, int]:
    """Parse an XYZ block and perceive connectivity with ``rdDetermineBonds``.

    XYZ files carry no bond information. RDKit's ``rdDetermineBonds`` infers
    both connectivity and bond orders from the 3D geometry, which also lets us
    assign a sensible total charge and spin multiplicity. Returns
    ``(mol, radical_electron_count)``; a non-zero radical count means the system
    is open-shell and the multiplicity should reflect that.

    Falls back to a connectivity-only perception if bond-order assignment fails
    (common for radicals / unusual valences).
    """
    from rdkit.Chem import rdDetermineBonds

    # Build a bare mol from the coordinates via an XYZ block RDKit can read.
    raw = Chem.MolFromXYZBlock(block)
    if raw is None:
        # Normalise through our own parser then re-serialise to a clean block.
        m3d = from_xyz_block(block)
        raw = Chem.MolFromXYZBlock(m3d.to_xyz_block())
    if raw is None:
        raise InputError("RDKit could not read the XYZ coordinates")

    mol = Chem.Mol(raw)
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=int(charge))
    except Exception:
        try:
            rdDetermineBonds.DetermineConnectivity(mol)
            Chem.SanitizeMol(mol, catchErrors=True)
        except Exception as exc:
            raise InputError(
                f"Could not perceive bonds from the XYZ geometry: {exc}"
            ) from exc

    radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    return mol, radicals


# ----------------------------------------------------------- CIF crystal (ASE)
def from_cif(
    block: str, mode: str = "cell", supercell: Tuple[int, int, int] = (1, 1, 1)
) -> Molecule3D:
    """Load a crystallographic CIF via ASE and return a ``Molecule3D``.

    ``mode``:
        * ``"cell"``     – atoms of the (optionally repeated) unit cell.
        * ``"supercell"``– unit cell repeated by ``supercell`` factors.
        * ``"cluster"``  – same atoms but with cell/PBC info dropped, giving a
          finite molecular cluster suitable for a gas-phase Gaussian job.

    Requires ASE. The cell vectors are reported in the title for reference.
    """
    try:
        import ase.io
        from io import StringIO
    except ImportError as exc:  # pragma: no cover
        raise InputError("ASE is required to read CIF files (pip install ase)") from exc

    try:
        atoms_obj = ase.io.read(StringIO(block), format="cif")
    except Exception as exc:
        raise InputError(f"ASE could not parse the CIF: {exc}") from exc

    if mode == "supercell" and tuple(supercell) != (1, 1, 1):
        atoms_obj = atoms_obj.repeat(supercell)

    from .molmodel import Atom as _Atom

    cell = atoms_obj.get_cell()
    atoms = [
        _Atom(sym, float(p[0]), float(p[1]), float(p[2]))
        for sym, p in zip(atoms_obj.get_chemical_symbols(), atoms_obj.get_positions())
    ]
    if not atoms:
        raise InputError("CIF produced no atoms")

    title = f"Structure from CIF ({mode})"
    try:
        lengths = atoms_obj.cell.lengths()
        title += f" a,b,c = {lengths[0]:.3f},{lengths[1]:.3f},{lengths[2]:.3f} A"
    except Exception:
        pass

    mol = Molecule3D(atoms=atoms, title=title, name="crystal")
    mol.multiplicity = mol.suggest_multiplicity()
    return mol
