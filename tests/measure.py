"""Reproducible measurements for P0-3 and P1-1, and the global acceptance test.

Everything here uses GFN2-xTB as a PROXY for the DFT optimization that this app
ultimately feeds. GFN2-xTB is NOT DFT; relative energies and the conformer
ordering can change at the DFT level. The .gjf is produced correctly, but the
DFT step count / energy are NOT measured here (no Gaussian binary available).

Metrics:
  * Energy: GFN2-xTB final energy (eV), reported as ΔE in kcal/mol vs the best.
  * Optimization effort: ASE-BFGS step count to fmax=0.05 eV/A, a reproducible
    stand-in for "how far from a minimum did we start".

Run:
  python tests/measure.py            # full run (needs tblite+ase for xTB)
  python tests/measure.py --quick    # fewer random conformers
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

from molgauss import build_polymer, search_conformers
from molgauss.energy import xtb_available, xtb_optimize, EV2KCAL
from molgauss.molmodel import Molecule3D
from molgauss import conformers as conf

QUICK = "--quick" in sys.argv
N_RANDOM = 4 if QUICK else 12


def hr(t):
    print("\n" + "=" * 72 + f"\n{t}\n" + "=" * 72)


# --------------------------------------------------------------------- P0-3
def p0_3_type_aware_and_drift():
    hr("P0-3: type-aware torsions + drift after MMFF (PEG-5, [*]OCC[*])")
    m = Chem.AddHs(build_polymer("[*]OCC[*]", 5, cap_left="H", cap_right="OH"))
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m)
    applied = conf.set_chain_conformation(m, family="gauche-1sign")
    print("Assigned by BOND TYPE (central bond -> target):")
    for (a, b, c, d, eb, ec, tgt) in applied:
        print(f"  {eb}-{ec} bond: target {tgt:+6.1f}")
    AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    final = conf.measure_backbone_torsions(m)
    print("\nFinal torsions after MMFF polish (drift from target):")
    amap = {(a, b, c, d): tgt for (a, b, c, d, _, _, tgt) in applied}
    for (a, b, c, d, eb, ec, ang) in final:
        tgt = amap.get((a, b, c, d))
        drift = "" if tgt is None else f"  Δ={_angdiff(ang, tgt):+6.1f}"
        print(f"  {eb}-{ec}: {ang:+7.1f}{drift}")


def _angdiff(a, b):
    d = (a - b + 180) % 360 - 180
    return d


# --------------------------------------------------------------------- P1-1
def p1_1_energy_search():
    hr("P1-1: energy-selected conformer search (PEG-5)")
    if not xtb_available():
        print("xTB unavailable -> this will rank by MMFF. Install tblite+ase.")
    mol = build_polymer("[*]OCC[*]", 5, cap_left="H", cap_right="OH")
    res = search_conformers(
        mol,
        chain_families=["all-trans", "gauche-1sign", "gauche-alt"],
        n_random=N_RANDOM,
        use_xtb=True,
        name="peg5",
    )
    print(f"Ranking method: {res.method}  ({res.note})")
    print(res.table())
    print(f"\nWinner: {res.candidates[0].source}")


# --------------------------------------------------------- global acceptance
def acceptance_peg10():
    hr("GLOBAL ACCEPTANCE: PEG-10, canonical HO-(CH2CH2O)10-H")
    # (1) formula
    mol = build_polymer("[*]OCC[*]", 10, cap_left="H", cap_right="OH")
    formula = rdMolDescriptors.CalcMolFormula(Chem.AddHs(mol))
    print(f"(1) Formula: {formula}  (expected C20H42O11)  "
          f"{'OK' if formula == 'C20H42O11' else 'FAIL'}")

    if not xtb_available():
        print("(2) SKIPPED energy/steps comparison: xTB backend unavailable.")
        print("    DFT step count / energy NOT measured (no Gaussian).")
        return

    # (2a) 'repo-style' baseline: the OLD index-based 'tgt' pattern on a single
    # conformer (what the previous code would export).
    base = _single_pattern_energy(mol, legacy_pattern="tgt")
    # (2b) new: energy-selected best from the search.
    res = search_conformers(
        mol, chain_families=["all-trans", "gauche-1sign", "gauche-alt"],
        n_random=N_RANDOM, use_xtb=True, name="peg10",
    )
    best = res.candidates[0]
    print(f"(2) GFN2-xTB final energies (proxy for DFT; NOT DFT):")
    print(f"    OLD single 'tgt' conformer : {base[0]:.6f} eV, "
          f"{base[1]} BFGS steps")
    print(f"    NEW energy-selected best   : {best.abs_energy:.6f} eV, "
          f"{best.opt_steps} BFGS steps  (source={best.source})")
    dE = (best.abs_energy - base[0]) * EV2KCAL
    print(f"    ΔE(new - old) = {dE:+.2f} kcal/mol "
          f"({'new is lower' if dE < 0 else 'old is lower'})")
    print("    NOTE: DFT optimization step count NOT measured (no Gaussian binary).")


def _single_pattern_energy(graph_mol, legacy_pattern="tgt"):
    """Replicate the OLD behaviour: one conformer, legacy index pattern, then
    GFN2-xTB relax; return (final_energy_eV, n_steps)."""
    m = Chem.AddHs(graph_mol)
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(m, maxIters=2000)
    try:
        conf.set_chain_extended(m, pattern=legacy_pattern)
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    except Exception:
        pass
    m3d = Molecule3D.from_rdkit(m, name="peg10_old")
    r = xtb_optimize(m3d, charge=0, fmax=0.05)
    return r if r is not None else (1e9, -1)


if __name__ == "__main__":
    p0_3_type_aware_and_drift()
    p1_1_energy_search()
    acceptance_peg10()
