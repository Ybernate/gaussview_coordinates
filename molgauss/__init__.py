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
    build_polymer,
)
from .geometry import embed_and_optimize, GeometryError
from .gaussian import (
    write_gjf,
    write_xyz,
    RouteConfig,
    build_genecp_tail,
    ROUTE_PRESETS,
    HEAVY_ELEMENTS_FOR_ECP,
)
from . import editor

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
    "build_polymer",
    "embed_and_optimize",
    "write_gjf",
    "write_xyz",
    "RouteConfig",
    "build_genecp_tail",
    "ROUTE_PRESETS",
    "HEAVY_ELEMENTS_FOR_ECP",
    "editor",
]

__version__ = "1.0.0"
