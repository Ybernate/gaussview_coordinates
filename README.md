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
| **Atom editor** | An editable table of `Index, Symbol, X, Y, Z`. Add atoms, delete atoms (e.g. strip a water), substitute elements (e.g. C → Si) keeping position, and optionally relax clashes with a force-field cleanup. |
| **Charge / spin** | Charge and multiplicity are auto-suggested (closed-shell singlet by default; doublet for odd-electron radicals) and fully editable — essential for ions, radicals, and metal complexes. |
| **Gaussian export** | Strict layout: Link 0 (`%chk / %mem / %nprocshared`) → route card → blank → title → blank → `charge multiplicity` → Cartesian block → **mandatory trailing blank line**. Route presets for opt, opt+freq, single-point, dispersion-corrected functionals, semi-empirical, transition-metal, and mixed-basis `genecp` jobs. |
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
  inputs.py            SMILES / InChI / name / polymer / FASTA / PDB / XYZ parsers
  geometry.py          ETKDGv3 embedding + MMFF94/UFF optimization (+ Open Babel fallback)
  editor.py            Atom CRUD (add / delete / substitute) + geometry cleanup
  gaussian.py          Strict .gjf/.com writer, .xyz writer, route presets
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

Build a polymer:

```python
from molgauss import build_polymer, embed_and_optimize
mol = build_polymer("[*]CC(c1ccccc1)[*]", n_units=10)   # polystyrene, 10 units
m3d, _ = embed_and_optimize(mol, force_field="AUTO")
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
