"""Gaussian / GaussView export (strict .gjf / .com formatting) and XYZ export.

The Gaussian input file layout is intentionally strict because GaussView is
unforgiving about it. The structure written here is:

    [Link 0 commands]   %chk=..., %mem=..., %nprocshared=...
    [Route section]     # <method/basis + keywords>
    [blank line]
    [Title line]
    [blank line]
    [Charge Multiplicity]
    [Cartesian coordinate block, one atom per line]
    [MANDATORY trailing blank line]

A trailing blank line after the molecule specification is required by Gaussian;
omitting it is the single most common cause of GaussView refusing the file.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

from .molmodel import Molecule3D

# Convenient presets for the route card.
ROUTE_PRESETS = {
    "Optimization (B3LYP/6-31G(d))": "# opt b3lyp/6-31g(d)",
    "Opt + Freq (B3LYP/6-31G(d))": "# opt freq b3lyp/6-31g(d)",
    "Single-point energy (B3LYP/6-31G(d))": "# b3lyp/6-31g(d)",
    "Opt (wB97XD/def2-SVP)": "# opt wb97xd/def2svp",
    "Opt + Freq (M06-2X/6-311+G(d,p))": "# opt freq m062x/6-311+g(d,p)",
    "Semi-empirical quick opt (PM6)": "# opt pm6",
    "Transition-metal opt (UB3LYP/LANL2DZ)": "# opt ub3lyp/lanl2dz",
    "Mixed basis / ECP for heavy metals (B3LYP/genecp)": "# opt b3lyp/genecp",
    "Mixed basis / ECP opt+freq (B3LYP/genecp)": "# opt freq b3lyp/genecp",
}

# Elements that genecp treats as "heavy": they need an effective core potential
# (ECP) because all-electron Pople bases such as 6-31G(d) are not defined for
# them. Light main-group atoms keep the regular basis. This covers the common
# transition / heavy metals (Pt, Au, Pd, Ru, ...) plus period-5+ p-block.
HEAVY_ELEMENTS_FOR_ECP = {
    # 4d / 5d / lanthanide-ish transition metals most users reach for
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "La", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    # heavy p-block often paired with LANL2DZ
    "Ga", "Ge", "As", "Se", "Br", "In", "Sn", "Sb", "Te", "I",
    "Tl", "Pb", "Bi", "Po", "At",
}


@dataclass
class RouteConfig:
    """User-configurable Gaussian job settings."""

    route: str = "# opt b3lyp/6-31g(d)"
    chk: str = "molecule.chk"
    mem: str = "4GB"
    nprocshared: int = 4
    extra_link0: List[str] = field(default_factory=list)
    # Optional trailing section (e.g. basis-set/ECP specs for heavy elements).
    tail: str = ""
    # Emit a per-atom freeze column (Symbol -1 x y z for frozen atoms). When
    # None the column is written only if at least one atom is flagged frozen.
    write_freeze_column: Optional[bool] = None

    def normalized_route(self) -> str:
        r = self.route.strip()
        if not r.startswith("#"):
            r = "# " + r
        return r


def _format_atom_line(
    symbol: str, x: float, y: float, z: float, freeze_code: Optional[int] = None
) -> str:
    """GaussView-friendly fixed-width Cartesian line.

    ``freeze_code`` is Gaussian's optional flag placed between the element and
    the coordinates: ``-1`` freezes the atom during optimization, ``0`` lets it
    move. ``None`` omits the column entirely (standard full optimization).
    """
    if freeze_code is None:
        return f" {symbol:<2} {x:>16.8f} {y:>16.8f} {z:>16.8f}"
    return f" {symbol:<2} {freeze_code:>2d} {x:>16.8f} {y:>16.8f} {z:>16.8f}"


def build_genecp_tail(
    mol: Molecule3D,
    light_basis: str = "6-31G(d)",
    heavy_basis: str = "LANL2DZ",
    ecp: str = "LANL2DZ",
) -> str:
    """Build the genecp trailing section for a molecule with heavy metals.

    Partitions the molecule's elements into "light" (regular all-electron basis)
    and "heavy" (ECP basis). Returns a block ready to drop in as ``RouteConfig.
    tail`` for a ``genecp`` route. Empty string if no heavy element is present.

    Layout produced (blank line between basis and ECP sections)::

        C H N O 0
        6-31G(d)
        ****
        Pt 0
        LANL2DZ
        ****

        Pt 0
        LANL2DZ
    """
    present = []
    for a in mol.atoms:
        if a.symbol not in present:
            present.append(a.symbol)
    heavy = [s for s in present if s in HEAVY_ELEMENTS_FOR_ECP]
    light = [s for s in present if s not in HEAVY_ELEMENTS_FOR_ECP]
    if not heavy:
        return ""

    lines: List[str] = []
    # --- basis-set section ---
    if light:
        lines.append(" ".join(light) + " 0")
        lines.append(light_basis)
        lines.append("****")
    lines.append(" ".join(heavy) + " 0")
    lines.append(heavy_basis)
    lines.append("****")
    # --- blank line then ECP (pseudopotential) section ---
    lines.append("")
    lines.append(" ".join(heavy) + " 0")
    lines.append(ecp)
    return "\n".join(lines)


def write_gjf(mol: Molecule3D, config: Optional[RouteConfig] = None) -> str:
    """Return the full text of a Gaussian .gjf / .com input file."""
    cfg = config or RouteConfig()
    lines: List[str] = []

    # --- Link 0 section ---
    if cfg.chk:
        lines.append(f"%chk={cfg.chk}")
    if cfg.mem:
        lines.append(f"%mem={cfg.mem}")
    if cfg.nprocshared and int(cfg.nprocshared) > 0:
        lines.append(f"%nprocshared={int(cfg.nprocshared)}")
    for extra in cfg.extra_link0:
        extra = extra.strip()
        if extra:
            lines.append(extra if extra.startswith("%") else f"%{extra}")

    # --- Route section ---
    lines.append(cfg.normalized_route())

    # --- Blank line + title + blank line ---
    lines.append("")
    title = (mol.title or "Molecule generated by MolGauss").strip() or "Molecule"
    lines.append(title)
    lines.append("")

    # --- Charge / multiplicity ---
    lines.append(f"{int(mol.charge)} {int(mol.multiplicity)}")

    # --- Cartesian coordinates (optionally with a freeze column) ---
    any_frozen = any(getattr(a, "frozen", False) for a in mol.atoms)
    if cfg.write_freeze_column is None:
        emit_freeze = any_frozen
    else:
        emit_freeze = bool(cfg.write_freeze_column)
    for a in mol.atoms:
        if emit_freeze:
            code = -1 if getattr(a, "frozen", False) else 0
            lines.append(_format_atom_line(a.symbol, a.x, a.y, a.z, freeze_code=code))
        else:
            lines.append(_format_atom_line(a.symbol, a.x, a.y, a.z))

    # --- Optional tail section (ECPs, custom basis, etc.) ---
    tail = (cfg.tail or "").strip("\n")
    if tail:
        lines.append("")
        lines.extend(tail.splitlines())

    # --- MANDATORY trailing blank line ---
    text = "\n".join(lines)
    if not text.endswith("\n"):
        text += "\n"
    text += "\n"  # the required blank line
    return text


def write_xyz(mol: Molecule3D) -> str:
    """Return a standard XYZ file (count, comment/title, atom lines)."""
    return mol.to_xyz_block()
