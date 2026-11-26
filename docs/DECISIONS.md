# Decision log — Gaussian-ready geometry work

Every entry: the decision, the hypothesis behind it, how it was measured, the
result, and (where relevant) what did **not** work. All energies are **GFN2-xTB
as a proxy for DFT** unless stated. GFN2-xTB is not DFT — relative orderings can
change at the DFT level. No Gaussian binary was available, so DFT step counts /
energies were **not measured**.

Reproduce everything with:

```bash
python tests/reproduce_findings.py   # H1–H5 against the pre-fix code
python tests/test_chirality.py       # P0-1 regression (chirality)
python tests/measure.py              # P0-3 / P1-1 / global acceptance
python tests/test_engine.py          # full engine suite (20 groups)
```

Measurement environment: Python 3.11, RDKit 2026.09.2, ASE 3.29.0, NumPy 2.4.6,
tblite (GFN2-xTB) via `tblite.ase.TBLite` + ASE `BFGS`, `fmax=0.05 eV/Å`.

---

## P0-1 — Proper-rotation PCA alignment (chirality)

**Decision.** Constrain the PCA alignment transform to a proper rotation
(`det = +1`). Assemble the net transform `R = S·P·V` explicitly (principal-axis
rotation `V`, axis permutation `P`, sign flips `S`); if `det(R) < 0`, flip the
sign of the least-significant axis (Z).

**Hypothesis.** The old code flipped each axis sign independently with no
determinant control, so ~half the time it produced a reflection that
mirror-inverts a chiral molecule.

**Measurement.** Signed tetrahedral volume of the stereocentre before/after
alignment (a direct, coordinate-based chirality descriptor), plus the transform
determinant. L-alanine, 12 seeds, and 5 chiral molecules × 12 seeds for the
regression.

**Result.**
- Before fix: det = −1 in **6/12** seeds; signed volume inverted in exactly the
  same **6/12** (physical enantiomer swap).
- After fix: **0/60** inversions, **0/60** improper rotations
  (`tests/test_chirality.py`).

**What did not work.** My first instrument — re-perceiving CIP with
`AssignStereochemistryFrom3D` — reported 0/12 inversions and hid the bug,
because it perceives stereo from the already-transformed coordinates. The signed
volume is the correct, transform-independent instrument.

---

## P0-2 — Independent end groups (`cap_left` / `cap_right`)

**Decision.** Replace the single `cap` with per-end `cap_left` / `cap_right`
presets (H, OH, OCH₃, CH₃, F, Cl, NH₂, COOH, vinyl) or any `[*]`-SMILES
fragment. Keep `cap` as a backward-compatible fallback. Show the resulting
formula in the UI.

**Hypothesis (H1).** A single symmetric cap cannot terminate an asymmetric
repeat unit correctly.

**Measurement.** `CalcMolFormula` of PEG decamers.

**Result.**
- `build_polymer("[*]OCC[*]", 10, cap="H")` → **C20H42O10** (H1 reproduced: the
  ethyl-terminated isomer, not canonical PEG).
- `build_polymer("[*]OCC[*]", 10, cap_left="H", cap_right="OH")` → **C20H42O11**,
  SMILES `OCCOCC…OCCO`, matching the explicit reference `HO-(CH2CH2O)10-H`.

**What did not work / chemistry caught.** Capping the O-terminus of `[*]OCC[*]`
(or using `[*]CCO[*]` with OH on both ends) gives **C20H42O12** — a peroxide
(O–O). This is chemically correct behaviour for that cap choice, not a bug; the
fix is to surface the formula so the user sees it. "Both OH" is ambiguous for an
ether repeat unit; the canonical PEG needs `H` on the O-end and `OH` on the
C-end.

---

## P0-3 — Type-aware backbone torsions

**Decision.** Assign each backbone torsion from its central-bond element pair
`(elem_b, elem_c)` rather than by index. Provide candidate families (all-trans,
single-sign gauche on heteroatom bonds, alternating gauche) and let the energy
search choose. Report final torsions after MMFF to show drift.

**Hypothesis (H2).** The old `tgt = [180, 60, 180, −60]` applied by `k % 4` lands
gauche values on the wrong bonds of a PEG backbone (bond types repeat with
period 3), scrambling C–O torsion signs.

**Measurement.** `tests/measure.py` prints, per backbone bond, the type and the
angle the old index scheme assigns vs. the new type scheme; and the
post-MMFF drift.

**Result.** H2 reproduced: the old scheme assigned 60/180/−60 to C–O bonds by
index. The new scheme assigns **C–C → 180, C–O/O–C → 60 consistently**. Post-MMFF
drift on C–O torsions was 15–26° (MMFF pulls away from the imposed gauche), which
is itself useful diagnostic output and motivates P1-1 (don't trust one pattern).

---

## P1-1 — Energy-selected conformer search

**Decision.** Generate candidates (type-aware families + N random ETKDG),
relax each, rank by energy, export the lowest. Use GFN2-xTB (`tblite`+`ase`)
when available, else MMFF. Show the ranked table and a clear "not DFT" caveat.

**Hypothesis (H4/H5).** A single hand-fixed pattern is not the lowest-energy
conformer; sampling several and selecting by energy is better.

**Measurement.** `tests/measure.py`, PEG-5 and PEG-10, GFN2-xTB single-point /
optimize, `fmax=0.05 eV/Å`.

**Result.**
- PEG-5: the winner was a **random** conformer at ΔE = 0, beating the type-aware
  families by ~3.5 kcal/mol — confirming no pattern should be assumed best.
- PEG-10 global acceptance: the new energy-selected best was **10.58 kcal/mol
  lower** (12 random conformers) than the old single-`tgt` conformer, at a
  comparable xTB/BFGS step count (29 vs 26).

**What did not reproduce.** The specific magnitudes in H4 (a by-type pattern
≈12.7 kcal/mol below the repo pattern) and H5 (39/60 random better, best ≈8
kcal/mol). My protocol (different builder, seeds, single start, gas phase) gave
the same **direction** but different numbers: H5 → 13/25 random better, best
−4.32 kcal/mol vs the legacy pattern. Documented as not numerically reproduced;
the qualitative conclusion stands and drives the design.

---

## P1-2 — Export modes, version-aware grid, resources

**Decision.** Build the route programmatically. Add a **G09/G16 selector**: emit
`int=ultrafine` only for G09. Make `scf=xqc` opt-in. Offer `opt freq` as one
job. Size `%nprocshared` from physical cores and `%mem` from RAM (with
headroom), not from atom count alone.

**Hypothesis (H6-adjacent).** `int=ultrafine` is redundant on G16.

**Verification (documentation, not measurement).** UltraFine is the **default
DFT grid in Gaussian 16** ([gaussian.com/integral](https://gaussian.com/integral/),
[gaussian.com/dft](https://gaussian.com/dft/)); G09's default is the coarser
FineGrid. So adding `int=ultrafine` on G16 is redundant; on G09 it is meaningful.
(Content rephrased for licensing compliance.)

**Result.** `mode_route("standard", "G16")` omits `int=ultrafine`;
`mode_route("standard", "G09")` includes it. `scf=xqc` appears only in the robust
mode / when opted in. Verified in `tests/test_engine.py`.

---

## P1-3 — Pre-export report

**Decision.** Before writing the `.gjf`, report formula, atom count,
electron/multiplicity parity, vdW clashes (excluding 1-2/1-3), bond
length/angle ranges, and a basis-function estimate with a cost band.

**Hypothesis.** The old flat 0.8 Å clash cutoff is too permissive — it is below
any real bond length, so it flags essentially nothing and misses genuine
non-bonded overlaps.

**Verification + measurement.**
- vdW criterion: flag non-bonded pairs below `0.6–0.75 × (r_vdW,i + r_vdW,j)`,
  excluding topological 1-2/1-3 pairs (from perceived connectivity). Verified on
  a folded butane: the 1-4 C···C pair at 2.088 Å is flagged; all 1-2/1-3 pairs
  are correctly excluded.
- Basis count: 6-31G(d) uses **6 Cartesian d functions** per heavy atom
  (Gaussian default for the 6-31G family; general default elsewhere is 5D/7F,
  [gaussian.com/gen](https://gaussian.com/gen/)). Water/6-31G(d) → 19 functions
  (O = 15, 2×H = 4), verified in the test suite. (Content rephrased for
  licensing compliance.)

---

## P2 / H6 — PCA does **not** help Gaussian symmetry (claim removed)

**Decision.** Remove the README claim that PCA alignment helps Gaussian detect
symmetry and speeds up integrals. Keep alignment only as a tidy default
viewing orientation.

**Verification (documentation).** When symmetry is used, Gaussian rotates the
molecule into its own **standard orientation** before computing and detects the
point group itself; the input orientation does not drive symmetry detection or
integral cost ([gaussian.com/symmetry](https://gaussian.com/symmetry/)).
(Content rephrased for licensing compliance.)

**Result.** README corrected; the benefit of alignment is stated honestly
(reproducible viewing orientation and a proper, chirality-preserving transform),
not a performance claim.

---

## Additional bug found (not in the brief)

**Name collision in `molgauss/__init__.py`.** The package exported both a
submodule `refine` and a function `refine`; the function shadowed the submodule
as a package attribute, so `import molgauss.refine` / `from . import refine`
resolved to the function and broke access to module-level helpers. Worked around
in callers with `importlib.import_module("molgauss.refine")`; the
`embed_and_optimize` internal import was changed to
`from .refine import refine as _refine_fn`. A cleaner long-term fix is to not
export a bare `refine` function name alongside the module.

---

## Open items / not measured

- **DFT step count and energy were not measured** (no Gaussian binary). All
  "fewer steps / lower energy" claims here are GFN2-xTB proxies and may differ at
  the DFT level.
- GFN2-xTB relative orderings can reorder at DFT; treat the "best" conformer as a
  strong starting guess, not a guaranteed global/DFT minimum.
- The basis-function estimate covers 6-31G(d) and STO-3G explicitly; other bases
  report "no estimate" rather than guessing.
