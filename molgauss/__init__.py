"""MolGauss: generate, edit and export 3D molecular structures for GaussView/Gaussian.

The package is split into a pure-Python chemistry engine (independent of the UI)
and a Streamlit front-end (``app.py`` at the project root).

Engine modules
--------------
molmodel  : lightweight atom-list model used by the editor and exporters.
inputs    : multi-format parsing (SMILES, InChI, name/PubChem, polymer, FASTA, PDB, XYZ).
geometry  : 3D embedding (ETKDGv3) and force-field pre-optimization (MMFF94 / UFF).
editor     : atom CRUD operations (add / delete / substitute element).
gaussian  : strict Gaussian .gjf/.com writer and .xyz export.
"""

from .molmodel import Atom, Molecule3D
from .inputs import (
    InputError,
    from_smiles,
    from_inchi,
    from_name,
    from_fasta,
    from_pdb_block,
    from_xyz_block,
    xyz_perceive_bonds,
    from_cif,
    build_polymer,
    END_GROUP_PRESETS,
)
from .geometry import (
    embed_and_optimize,
    GeometryError,
    search_conformers,
    ConformerSearchResult,
    ConformerCandidate,
)
from . import energy
from .gaussian import (
    write_gjf,
    write_xyz,
    RouteConfig,
    build_genecp_tail,
    suggest_resources,
    mode_route,
    ROUTE_PRESETS,
    EXPORT_MODES,
    HEAVY_ELEMENTS_FOR_ECP,
)
from .refine import refine, align_principal_axes, detect_clashes, RefineReport
from .report import build_report, PreExportReport, vdw_clashes
from .gaussian import build_route
from . import editor
from . import conformers

__all__ = [
    "Atom",
    "Molecule3D",
    "InputError",
    "GeometryError",
    "from_smiles",
    "from_inchi",
    "from_name",
    "from_fasta",
    "from_pdb_block",
    "from_xyz_block",
    "xyz_perceive_bonds",
    "from_cif",
    "build_polymer",
    "END_GROUP_PRESETS",
    "embed_and_optimize",
    "search_conformers",
    "ConformerSearchResult",
    "ConformerCandidate",
    "energy",
    "write_gjf",
    "write_xyz",
    "RouteConfig",
    "build_genecp_tail",
    "suggest_resources",
    "mode_route",
    "build_route",
    "build_report",
    "PreExportReport",
    "vdw_clashes",
    "ROUTE_PRESETS",
    "EXPORT_MODES",
    "HEAVY_ELEMENTS_FOR_ECP",
    "refine",
    "align_principal_axes",
    "detect_clashes",
    "RefineReport",
    "editor",
    "conformers",
]

__version__ = "1.0.0"
