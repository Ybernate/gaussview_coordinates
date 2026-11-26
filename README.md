# MolGauss

Generate, edit and export **3D molecular structures ready for GaussView / Gaussian**.

MolGauss takes almost any chemical input — a SMILES/InChI string, a common or
IUPAC name, a polymer monomer, a peptide sequence, or an uploaded PDB/CIF/XYZ
file — builds a force-field-optimized 3D geometry, lets you edit atoms directly,
and exports a strictly formatted `.gjf` / `.com` (or `.xyz`) that GaussView
opens without complaint.

## Features

| Area | What it does |
|------|--------------|
| **Input** | SMILES, InChI, name/IUPAC (resolved via PubChem), polymer from a monomer + degree of polymerization, peptide from a FASTA sequence, and PDB / CIF / XYZ / MOL / SDF file upload. |
| **3D engine** | Explicit hydrogens (`AddHs`), ETKDGv3 embedding, and MMFF94 / UFF minimization. `AUTO` picks MMFF94 for organics and UFF when metals or heavy elements appear. Open Babel is a fallback when RDKit cannot embed (common for coordination complexes). |
| **Gaussian-ready refinement** | Every structure passes a refinement pass before export: a steric-clash filter (flags non-bonded contacts < 0.8 Å), a standard-bond-length sanity check, and **PCA alignment** — centre of mass moved to the origin and the longest molecular axis rotated onto X, so GaussView opens on a tidy orientation and Gaussian can detect point-group symmetry (faster integrals). |
| **Chain conformations** | Polymers and peptides are built in a regular, extended backbone instead of a random coil — a much better DFT starting point. Polymers: all-trans (fully extended) or trans-gauche (PEG-like helix). Peptides: β-strand (φ≈−135, ψ≈135), α-helix (φ≈−57, ψ≈−47), or polyproline-II, applied via `rdMolTransforms.SetDihedralDeg`. |
| **XYZ bond perception** | `.xyz` files with no connectivity are run through RDKit `rdDetermineBonds` to infer bonds, bond orders, and the spin multiplicity. |
| **CIF crystals** | Crystallographic `.cif` files are read with ASE; export the unit cell, a supercell, or a finite gas-phase cluster. |
| **Atom editor** | An editable table of `Index, Symbol, X, Y, Z`. Add atoms, delete atoms (e.g. strip a water), substitute elements (e.g. C → Si) keeping position, and optionally relax clashes with a force-field cleanup. |
| **Charge / spin** | Charge and multiplicity are auto-suggested (closed-shell singlet by default; doublet for odd-electron radicals) and fully editable — essential for ions, radicals, and metal complexes. |
| **Gaussian export** | Strict layout: Link 0 (`%chk / %mem / %nprocshared`) → route card → blank → title → blank → `charge multiplicity` → Cartesian block → **mandatory trailing blank line**. Route presets plus **three curated export modes**: (1) enhanced standard DFT opt with D3 dispersion + `int=ultrafine` + `scf=xqc`, (2) robust opt+freq with extra optimizer/SCF cycles and `calcfc`, (3) fast PM6 pre-screen. `%mem` / `%nprocshared` are auto-sized to the atom count (editable). |
| **Mixed basis / ECP (genecp)** | For heavy metals (Pt, Au, Pd, Ru, …) that have no all-electron Pople basis, the app auto-builds the `genecp` tail: a regular basis (default 6-31G(d)) for light atoms and an ECP basis (default LANL2DZ) for the heavy ones, with the pseudopotential block. Heavy elements are detected automatically; bases are editable. |
| **Freeze atoms (partial optimization)** | A **Freeze** checkbox per atom in the editor table writes the Gaussian `-1` freeze code (`C  -1  x  y  z`) so you can optimize an active site while holding the rest fixed. Frozen atoms are tinted in the 3D view. |
| **Preview** | Interactive 3D view (py3Dmol / stmol) with stick, ball-and-stick, spacefill, and line styles, plus an optional **atom-index overlay** (`#:symbol`) so a table row maps unambiguously to an atom in 3D. |

## Install & run

### Option A — local Python

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

streamlit run app.py
```

Then open the URL Streamlit prints (default http://localhost:8501).

### Option B — Docker (no local Python / RDKit needed)

```bash
docker build -t molgauss .
docker run --rm -p 8501:8501 molgauss
```

Then open http://localhost:8501. The image is based on `python:3.11-slim` and
installs RDKit and Open Babel from their manylinux wheels, so there is nothing
to configure on the host.

## Project layout

```
app.py                 Streamlit front-end (UI only)
molgauss/
  __init__.py          Public API
  molmodel.py          Molecule3D / Atom model + RDKit <-> model conversion
  inputs.py            SMILES / InChI / name / polymer / FASTA / PDB / XYZ / CIF parsers
  geometry.py          ETKDGv3 embedding + MMFF94/UFF optimization (+ Open Babel fallback)
  conformers.py        Backbone conformations (all-trans/helix, β-strand/α-helix)
  refine.py            PCA alignment, steric-clash filter, standard-bond check
  editor.py            Atom CRUD (add / delete / substitute / freeze) + geometry cleanup
  gaussian.py          Strict .gjf/.com writer, export modes, ECP/genecp, route presets
tests/
  test_engine.py       End-to-end engine smoke tests
requirements.txt
Dockerfile             Container image (python:3.11-slim + pip wheels)
.dockerignore
run.sh                 Local venv bootstrap + launcher
```

## Using the engine without the UI

The chemistry engine is a plain Python package, usable in scripts:

```python
from molgauss import from_smiles, embed_and_optimize, write_gjf, RouteConfig

mol = from_smiles("CC(=O)Oc1ccccc1C(=O)O")          # aspirin
m3d, method = embed_and_optimize(mol, name="aspirin", title="Aspirin opt")
gjf = write_gjf(m3d, RouteConfig(route="# opt freq b3lyp/6-31g(d)",
                                 chk="aspirin.chk", mem="4GB", nprocshared=4))
open("aspirin.gjf", "w").write(gjf)
```

Build a polymer in an extended conformation, PCA-aligned:

```python
from molgauss import build_polymer, embed_and_optimize
mol = build_polymer("[*]CC(c1ccccc1)[*]", n_units=10)   # polystyrene, 10 units
m3d, info = embed_and_optimize(mol, force_field="AUTO",
                               conformation="all-trans", align=True)
# info includes the refinement summary (alignment, clashes, bond sanity)
```

Build a peptide as an α-helix and export with the enhanced DFT mode:

```python
from molgauss import from_fasta, embed_and_optimize, write_gjf, RouteConfig
from molgauss import mode_route, suggest_resources
mol = from_fasta("ACDEFG")
m3d, _ = embed_and_optimize(mol, conformation="alpha", align=True, name="hexapep")
mem, nproc = suggest_resources(len(m3d))
gjf = write_gjf(m3d, RouteConfig(route=mode_route("standard"), mem=mem, nprocshared=nproc))
```

Read an XYZ with bond perception, or a CIF crystal:

```python
from molgauss import xyz_perceive_bonds, from_cif, Molecule3D
mol, radicals = xyz_perceive_bonds(open("frag.xyz").read(), charge=0)
cluster = from_cif(open("crystal.cif").read(), mode="cluster")       # finite cluster
supercell = from_cif(open("crystal.cif").read(), mode="supercell", supercell=(2, 2, 2))
```

Edit atoms:

```python
from molgauss import from_xyz_block, editor
m3d = from_xyz_block(open("water.xyz").read())
editor.substitute_element(m3d, 0, "S")   # O -> S
editor.add_atom(m3d, "Cl", 2.0, 0.0, 0.0)
editor.delete_atom(m3d, 3)
```

Freeze atoms for a partial optimization (writes the `-1` code):

```python
from molgauss import editor, write_gjf, RouteConfig
editor.set_frozen(m3d, 0, True)          # freeze atom 0
editor.set_frozen(m3d, 1, True)          # freeze atom 1
gjf = write_gjf(m3d)                      # freeze column appears automatically
# -> "O  -1  x  y  z" for frozen atoms, "H   0  x  y  z" for the rest
```

Heavy metals with a mixed basis / ECP (`genecp`):

```python
from molgauss import build_genecp_tail, write_gjf, RouteConfig
tail = build_genecp_tail(m3d, light_basis="6-31G(d)",
                         heavy_basis="LANL2DZ", ecp="LANL2DZ")
gjf = write_gjf(m3d, RouteConfig(route="# opt b3lyp/genecp", tail=tail))
```

## Tests

```bash
source .venv/bin/activate
python tests/test_engine.py
```

Covers each input path, the editor CRUD, the metal/UFF branch, and the
Gaussian/XYZ exporters. The PubChem name-lookup test is skipped gracefully when
offline.

## Notes on correctness

- **Charge & multiplicity** are suggestions for neutral closed-shell molecules.
  For radicals, ions, and open-shell metal complexes you must set them yourself
  (the UI exposes both; the suggested lowest-spin value is shown for reference).
- **Force field is a pre-optimization**, not a substitute for the quantum
  calculation. It removes atomic clashes so Gaussian starts from a sane
  geometry. Run the actual DFT/ab-initio optimization in Gaussian.
- **Uploaded PDB/XYZ geometries are preserved as-is** by default so experimental
  coordinates are respected; enable "Re-optimize" to relax them.
