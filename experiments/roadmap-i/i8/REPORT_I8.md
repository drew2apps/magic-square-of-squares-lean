# REPORT_I8 (checkpoint, 28 Nov 2026, about 00:10 UTC): t10's pairing test on the ratios that need a good prime

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review z17) holds.


Written by an AI model; not read by a person. Experiment, not counted; tier D. Nothing here changes the census. **Checkpoint:** this report covers lines 1 to 360 of `i2/i2_needgood_list.txt` (all 360 are from H2's exhaustive list). The run goes on over all 3,876 ratios, in list order; G2's 937 are lines 2,940 to 3,876 and are not reached yet.

## Bottom line

- **3,160 of 3,166 residual tuples are removed by the Cassels–Tate pairing** on Sel² of a pair curve y² = F_j·F_k, in 357 of the 360 ratios. *(computed, exact, given t10's Sel² and pairing code as reviewed in z11)*
- **6 tuples in 3 ratios are "neither".** The pairing is trivial on their classes on every pair curve, and no pair curve has rank 0. They are (18313,120) and its partner (18433,18193), 2 tuples each, killed at the good prime 61, and (20056,19825), 2 tuples killed at 17. For each of them, an E± curve has rank exactly 1 by PARI, and Demjanenko–Manin through z ↦ −z finds no four-corner point. So each is closed, conditionally (see "The three 'neither' ratios").
- **11 ratios are also excluded by a pair curve of rank 0** (dim Sel² = 2, torsion complete, every point over t ∈ {0, ±1, ∞}): 7 by E_LU and 4 by E_08. Their tuples are pairing-removed as well.
- **The mixed curves do most of the work.** Kills per curve (a tuple counts on every curve that kills it): 0L, 0U, 8L and 8U 3,058 each, 08 1,836, LU 1,794. 462 tuples are killed only on mixed curves (0L, 0U, 8L, 8U) and on neither 08 nor LU.
- **The residual file and t10 agree:** recomputing 50 ratios (seed 1, from the first 2,900 of the list) with t10's own `residual()` gives the same tuples and killing primes for 50 of 50.

## The three "neither" ratios

| Ratio | Tuples | Killing good prime | Per pair curve (08, 0L, 0U, 8L, 8U, LU) | E± (PARI 2.19.1 and 2.17.2 agree) |
|---|---|---|---|---|
| (18313,120) | [1,2,1,335351569], [1,2,670703138,2] | 61 | kernel ×5, LU: class of a point found | E+: root number −1, ellrank [1,1]; E−: [1,3] |
| (18433,18193), its partner | [2,670703138,2,1], [335351569,1,2,1] | 61 | 08: class of a point found, kernel ×5 | E+: [1,3]; E−: [2,2] |
| (20056,19825) | [2,1,1,1], [2,1,2,2] | 17 | kernel ×6 | E+: root number +1, [2,2]; E−: root number −1, [1,1] |

- **E± closes them, conditionally.** The repository's `experiments/census2000/dm_pari.py` (Demjanenko–Manin with z ↦ −z, PARI only):
  - E+ of (18313,120): bound |a| ≤ 58, 119 candidates, no four-corner point, so no square at this ratio. That covers the partner too.
  - E− of (20056,19825): |a| ≤ 55, 113 candidates, no point.
  - (20056,19825)'s partner (39881,231) is not in the list. Its curve E also has rank 1 by PARI and gives no point (|a| ≤ 28).
- **What it rests on:**
  - The rank 1 comes from PARI's `ellrank` alone. The repository's exact 2-isogeny descent with a second descent (`scripts/g1_full/rank_indep.py`, its `rank_bound` called on these curves) bounds the rank only by 3 on both curves, like 8 of the 12 census curves.
  - The generator is saturated only at primes ≤ 1000.
  - The height bounds are in floating point.
  - So this is *conditional*, not one of the census's exact closures.
- **What it means:** on these 360 ratios, the pairing does not reproduce every good-prime kill. The 6 tuples above are killed by a good prime (17 or 61) but by no pairing on the six pair curves. Their ratios fall to a rank-1 curve instead.

## How it was checked

- `i8_pairct.py` runs t10_pairct.analyse_ratio on each ratio with the tuples from `i2_residual.jsonl`. Before use, each record is checked: every tuple is allowed at every bad prime and refused at its own killing prime, and that prime is good. Per tuple it records the verdict, the curves that kill it (`killed_by`), and t10's per-curve label.
- **Rank 0** on a pair curve needs:
  - dim Sel² = 2;
  - the torsion found by PARI, closed under the group law and checked on the curve, equal in size to the gcd of #E(F_p) over 25 good odd primes, counted exactly in Python;
  - every point over t ∈ {0, ±1, ∞}.
- **Controls** (`--negative`, `i8_negative.txt`):
  - the true record of (15001,2081) is accepted;
  - a forged tuple and a forged killing prime (the good prime 23, below the true one) are refused;
  - (5623,1017) is excluded by E_08 of rank 0 (torsion 8 of bound 8), as in t10 and z11;
  - a torsion list with a point left out, and a rank-0 claim on E_0L (dim 6), are refused;
  - (4457,287) gives t10's 0L, 0U, 8L, 8U.
- **PARI and CI.** PARI 2.19.1 produced every output here; CI has 2.17.4.
  - **Basis-independent outputs:** the diffed outputs (`i8_negative.txt`, `i8_compare50.txt`, `i2_residual_negative.txt`) print only facts that don't depend on PARI's choice of Sel² basis or of points: membership of a class in Sel², the set of curves on which the pairing with some element of Sel² is −1 (the pairing is linear and trivial on point images, so this doesn't depend on the basis), dim Sel², the complete torsion group, and killing primes.
  - **Checked under PARI 2.17.2:** the diffed outputs print exactly the committed text under PARI 2.17.2 (cypari2 2.2.3 in a venv), and so does the analysis of (18313,120). 2.17.4 itself was not available here.
  - **Not checked under 2.17.4:** all of the above, the checkpoint files, and the E± outputs in `neither/`.
  - **PARI-dependent labels:** in the jsonl, a per-curve label "image" (the class of a point found) versus "kernel" depends on which points PARI finds. The verdicts do not.
- PARI's `nbthreads` is set to 1 (its threads can hang under `ulimit -v`). Each worker runs under `ulimit -v 7000000` and is restarted every 40 ratios, because t10's memory grows.

## Files

- `i8_pairct.py` (the run, `--merge`, `--compare50`, `--negative`), `run_i8.sh` (3 workers).
- `i8_checkpoint.py` takes the first K ratios in list order from the workers' files, and refuses if any is missing.
- `i8_first360.jsonl` (lines 1 to 360 of the list, in list order), `i8_first360_summary.txt`, `i8_negative.txt`, `i8_compare50.txt`, `times.txt`.
- `neither/`:
  - `epm.py`: E± ranks by PARI, a finder;
  - `ri.py`: the exact rank bound with rank_indep;
  - `dm_*.jsonl`: dm_pari on each ratio and its partner;
  - their outputs.
- Each record: `{m, n, source, tuples: [{tau, kill, pairs: {08, 0L, 0U, 8L, 8U, LU}, killed_by, verdict}], rank0}`.

## Proposed CI lines

```
python3 experiments/roadmap-i/i8/i8_pairct.py --negative      # 13 s (12 s under PARI 2.17.2)
python3 experiments/roadmap-i/i8/i8_pairct.py --compare50     # 5 s
python3 experiments/roadmap-i/i8/i8_checkpoint.py 360 3       # under 1 s, prints i8_first360_summary.txt (needs the workers' files)
```

The full run took about 25 to 70 s per ratio per worker, depending on the load.

## Proposed LEDGER/CHANGES text (draft)

> **Experiment, not counted (tier D).** I8: on the first 360 of the 3,876 ratios that need a good prime, the Cassels–Tate pairing on Sel² of the six pair curves removes 3,160 of the 3,166 residual tuples. The other 6, in (18313,120), its partner (18433,18193) and (20056,19825), are in the pairing's kernel on every pair curve. Each of those ratios has an E± curve of rank 1 by PARI, on which Demjanenko–Manin finds no point (conditional on PARI's rank and generator). Computed. Written by an AI model; not read by a person.
