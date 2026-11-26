# rank1: exact rank bounds for the seven conditional ratios (run 5, F2)

**Experiment, not counted; tier D of `LEDGER.md`.** Scripts and outputs written by an AI model (Claude, the "Fable"
session of run 5); reviewed by a separate AI agent (`reviews/review_x3.md`: holds). Not read by a person.

Exact rank bounds from the repository's `scripts/g1_full/rank_indep.py` and `rank_second.py` (24d6ef6, unchanged) for
the five curves behind the seven conditional ratios: 3 (E1 of (3540,2581)), 5 (E1 of (4440,2231)), 3 (E1 of
(6480,4961)), 3 (E1 of (8249,3720); was 5: Θ₂ matrix of F₂-rank 1), 3 (D1 of (12007,1799); was 4: no ℚ₂-point on any
class of the open coset). None reaches 1; the seven ratios stay conditional on rank 1. Computed (exact).

- `f2_curves.py` (+ `f2_curves_hints.json`): the five curves and their hint points (checked on the curve).
- `f2_pairing.py k` (+ `f2_pairing_0.log` … `_2.log`): the descent and the pairing ⟨y0,y1⟩ of `rank_second.analyse`;
  curves 0 and 1 have no open coset on one side, curve 2 pairs to +1 everywhere (no drop). Curves 3 and 4 have more than
  one open coset and stop on the script's own assertion; those logs are not kept.
- `f2_general.py k` (+ `f2_general_3.log`, `_4.log`): `rank_second.analyse_general`; on curve 3 (8249) the Θ₂ matrix has
  F₂-rank 1, bound 5 → 3; on curve 4 (D1) no class is certified and the bound stays 4.
- `f2_general_neg.py` (+ `f2_general_negative.txt`): the repository's `rank_second.negative()` and point-class controls.
- `f2_d1_local2.py` (+ `.txt`, `_negative.txt`): the four classes of D1's open E′ coset have no ℚ₂-point
  (`rank_indep.solvable` at depth 40), bound 4 → 3. The review confirmed this with its own 2-adic solver.
- The session's original `sympy_pari_shim.py` and `run_shim.py` are not committed: the scripts use sympy.
