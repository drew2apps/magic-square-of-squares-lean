# Experiments beyond the census (26 November 2026)

Written by an AI model (Claude) in a separate session, and committed here unchanged except for the corrections listed below; no person has
read any of it. **Nothing in this folder is counted in the census** (which stops at `m ≤ 1000`, `README.md`), nothing is
in Lean, and every result is tier D of `LEDGER.md` (computed, with the labels each write-up gives). It is kept here so
that the results can be cited and rerun.

## What is here

| Folder | Write-up | What it claims (labels as in the write-up) |
|---|---|---|
| `census2000/` | `CENSUS3000.md` | the census pipeline for `1000 < m ≤ 3000` |
| `census-wide/` | `CENSUS5000.md`, `CENSUS7000.md`, `CENSUS9000.md` | the same for `3000 < m ≤ 9000`. Together: every coprime ratio with `1000 < m ≤ 9000` is excluded except **(7735,4968)**, which run 5 (`census-wide/g4k/`) excludes only **conditionally on rank C(ℚ) = 1** of the quotient C: y² = g4·K (that rank is PARI's bound only; the exact descent gives 2); **seven** exclusions are **conditional on a rank being 1** ((3540,2581), (4440,2231), (6121,959), (6671,2209), (6480,4961), (6903,5104), (8249,3720); the exact bounds for their five curves are 3, 5, 3, 3, 3, `census-wide/rank1/`); the nine rank-1-route exclusions are conditional on (a) the generator `G` generating `E1(ℚ)` modulo torsion (saturation checked at primes `≤ 2000`) and (b) an AI reading of the constants of Checcoli, Veneziano and Viada; every Demjanenko–Manin exclusion is conditional on PARI's generators (saturation checked only at primes `≤ 1000`); every rank bound is PARI's; (6557,324) and (6881,6233) rest on a **numerical** `L(E,1) ≈ 1.44409` (no rigorous tail bound); `census-wide/g4k/` gives them a second route, also conditional on rank C(ℚ) = 1 of their g4·K (PARI's bound; the exact descent gives 3 and 5) |
| `open-ratio/` | `WRITEUP.md` | (627,364) over `ℚ(i)`: 2-descent and φ-descents (rank `≤ 2`), the Sha ratios across the isogeny class, a cyclic 4-isogeny descent (bounds 3 to 5), and a new PARI implementation of the Cassels–Tate pairing on Sel², which vanishes on Sel² of `E₁` (as `scripts/second_descent/kctp*.py` found) and, in run 5, on the deciding value of Sel² of `E₄`, so Sel² with its pairing cannot prove rank 0 for `E₄` (a second implementation, unreviewed) |
| `HANDOFF_26NOV.md` | | the session's one-page summary |

The pipeline is: congruences (`census2000/gen_pairs.py`) → the repository's class criterion → second descent
(`census-wide/wide_pari.py`, `sd_pari.py`) → rank-0 quotients (the repository's `scripts/g1_full/global.py`, unchanged) →
Demjanenko–Manin (`census2000/dm_pari.py`) → the rank-1 route (`census-wide/route/route_pari.py`). Each stage's inputs and
outputs are the `.json`/`.jsonl` files beside the scripts; the logs are kept as they were written.

## Rerunning

- Tools: Python 3 with cypari2 (the runs used PARI 2.19.1; this repository's CI pins cypari2 2.2.9, which bundles an
  older PARI), numpy, and `gp`.
- Paths: the scripts were run from the session's own folders and still name them. From the repository root, a copy that
  runs here is made by
  ```
  cp -r experiments /tmp/exp && cd /tmp/exp
  grep -rl /tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo . | xargs sed -i "s|/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo|$OLDPWD|g"
  grep -rl /mnt/project-files/mss-experiments . | xargs sed -i "s|/mnt/project-files/mss-experiments|/tmp/exp|g"
  ```
  (the first path is a clone of this repository at `391e693`; the census programs it uses — `scripts/g1_full/global.py`,
  `scripts/gauss_descent/`, `scripts/second_descent/` — are identical at `391e693` and at this commit, apart from
  `scripts/second_descent/lball_exact.py`, which these experiments do not import, and new files).
- Run times are hours for the full census stages; no CI job reruns this folder (`OUTPUTS.md`, class B).

## Not included

The run for `9000 < m ≤ 11000` was still going when the bundle was made; its partial files are not included.

## Review and corrections

A separate AI agent reviewed this folder before it was committed (26 November 2026, UTC). It recomputed every count of
every stage table from the `.json`/`.jsonl` files (each ratio that passes the congruences appears exactly once; the
partner map is right, also for partners outside the range), re-ran `gen_pairs.py` on four sub-ranges (identical
listings), and re-ran a sample of closures of every stage under PARI 2.17.4 (10 second-descent records, 9 rank-0
ratios, 18 Demjanenko–Manin curves, 5 route runs, the conditional Demjanenko–Manin run): all reproduce. It read the
Hilbert-reciprocity shortcut of `open-ratio/ctp/ctp.gp` and found it correct under the hypotheses now stated in
`WRITEUP.md`, and recomputed the deciding pairing value on `E₄` (0). Its verdict: include after fixes. Applied here:

- `CENSUS3000.md`: "every ratio with m ≤ 3000" corrected to `1000 < m ≤ 3000` ((627,364) and (991,263) are not excluded).
- The headline sentences of `CENSUS3000/5000/7000/9000.md`, the running total of `CENSUS9000.md` and `HANDOFF_26NOV.md`
  now list every condition: PARI's rank bounds, the Demjanenko–Manin generators, the nine route exclusions and their
  conditions (a), (b), and the numerical `L`-value of (6557,324), (6881,6233).
- `CENSUS7000.md`: (6903,5104) rests on rank 1 and the Demjanenko–Manin generator condition (not (a), (b)); "rank 1 or 3"
  and "0 or 2" now say "if Ш[2^∞] is finite"; the tail of the `L`-series has an elementary bound (below 0.01); the BSD
  ratio (15.9999999997) is a consistency check, not evidence.
- `CENSUS9000.md`: PARI's lower bound 2 for two curves may come from parity (one point), so rank 2 there needs Ш finite;
  no exclusion depends on it.
- `open-ratio/`: after the review, the session sent an updated bundle with the drivers and logs of the pairing's
  validation (`ctp/validation/`: 9 curves over ℚ, six of them with pairing rank 2, the shortcut equal to the direct
  method on all 9; 8 curves over `ℚ(i)`, two with pairing rank 2) and a second logged run on `E₄` (`ctp/run627b.log`,
  the deciding entries 0, each twice); `WRITEUP.md` and `HANDOFF_26NOV.md` now cite these files, and `run627.out` for the
  first run. The check used over `ℚ(i)` holds only without elements of order 4 in Ш (noted there). Still without
  drivers or logs here: the 4-isogeny descent's six test curves and `l1_sieve.py` against `ellL1` (marked so).
  `ctp/validation/val_Q.err` holds one line, `INCONSISTENT at inf [1, 1, 0]`, which the write-up does not explain.

Not changed: `open-ratio/ctp/torsion_E4.gp` still tests "is a square" at 300 split primes (heuristic in that direction;
the review made it exact in its own copy, with the same result; so did review agent x3 of run 5,
`referee_x3/torsion_exact.py`).

## Run 5 (26 November 2026)

The same session sent a second bundle (run 5, part 1), reviewed by a separate AI agent (`reviews/review_x3.md`) and
committed with these changes: its PARI stand-in for sympy (`sympy_pari_shim.py`, `run_shim.py`) is left out and the
scripts import sympy (with the stand-in, 7 of 11 of its checks fail under PARI 2.17.4; with sympy all reproduce the
committed outputs); the scripts' fallback paths to the session's folders are removed; logs that only show a script
stopping on its own assertion are left out. Everything is tier D, written by an AI model and not read by a person.

- `census-wide/g4k/` (F4, `NOTES.md`): Demjanenko–Manin on C: y² = g4·K in exact arithmetic. **If rank C(ℚ) = 1**, every
  solution of (7735,4968), (6557,324) and (6881,6233) is degenerate. Rank 1 is PARI's `ellrank` bound only; the exact
  descent gives 2, 3, 5. Conditional.
- `census-wide/rank1/` (F2, `NOTES.md`): exact rank bounds 3, 5, 3, 3, 3 for the five curves behind the seven
  conditional ratios ((8249,3720) lowered from 5, D1 of (12007,1799) from 4). None reaches 1.
- `open-ratio/` (F1): the real-place fix in `ctp/ctp.gp` and the logs that explain `val_Q.err`
  (`ctp/validation/val_Q_run1.err`, `val_Q_run2.err`; the rerun's `val_Q.err` is empty and `val_Q.log` unchanged); the
  driver and log of the 4-isogeny descent's six test curves (`four/validation/`); `WRITEUP.md` updated to match.
  `ctp/run627b.log` is kept as committed (the session's copy only adds a hand-written stop note, now in `WRITEUP.md`).
