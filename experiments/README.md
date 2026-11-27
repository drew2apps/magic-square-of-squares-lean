# Experiments beyond the census (26 November 2026)

Written by an AI model (Claude) in a separate session, and committed here unchanged except for the corrections listed below; no person has
read any of it. **Nothing in this folder is counted in the census** (which stops at `m ≤ 1000`, `README.md`), nothing is
in Lean, and every result is tier D of `LEDGER.md` (computed, with the labels each write-up gives). It is kept here so
that the results can be cited and rerun.

## What is here

| Folder | Write-up | What it claims (labels as in the write-up) |
|---|---|---|
| `census2000/` | `CENSUS3000.md` | the census pipeline for `1000 < m ≤ 3000` |
| `census-wide/` | `CENSUS5000.md`, `CENSUS7000.md`, `CENSUS9000.md`, `CENSUS11000.md` | the same for `3000 < m ≤ 9000`, and (run 5, part 2) for `9000 < m ≤ 11000`: there every ratio is excluded under the same conditions, (10512,145) and its partner (10657,10367) only **conditionally on a rank being 1** (rank-1 route; PARI bounds [1,3], root number −1), (10871,4751) through its route-closed partner (7811,3060). Together, for `1000 < m ≤ 9000`: every coprime ratio with `1000 < m ≤ 9000` is excluded except **(7735,4968)**, which run 5 (`census-wide/g4k/`) excludes only **conditionally on rank C(ℚ) = 1** of the quotient C: y² = g4·K (that rank is PARI's bound only; the exact descent gives 2); **seven** exclusions are **conditional on a rank being 1** ((3540,2581), (4440,2231), (6121,959), (6671,2209), (6480,4961), (6903,5104), (8249,3720); the exact bounds for their five curves are 3, 5, 3, 3, 3, `census-wide/rank1/`); the nine rank-1-route exclusions are conditional on (a) the generator `G` generating `E1(ℚ)` modulo torsion (saturation checked at primes `≤ 2000`) and (b) an AI reading of the constants of Checcoli, Veneziano and Viada; every Demjanenko–Manin exclusion is conditional on PARI's generators (saturation checked only at primes `≤ 1000`); every rank bound is PARI's; (6557,324) and (6881,6233) rest on a **numerical** `L(E,1) ≈ 1.44409` (no rigorous tail bound); `census-wide/g4k/` gives them a second route, also conditional on rank C(ℚ) = 1 of their g4·K (PARI's bound; the exact descent gives 3 and 5), and `census-wide/lvalue/` an exact interval L(E,1) ∈ [1.444085955623, 1.444085956081] **conditional on PARI's `ellap`** at the good primes (not used by the census table) |
| `census-wide/table/` | `NOTES.md` | one table of all 20,682 ratios `1000 < m ≤ 11000` that pass the congruences, with the first step that closes each and that step's label (bookkeeping over the committed records): class 17,636, second descent 2,582, repository L-value 23, rank 0 231, Demjanenko–Manin 175, repository route 9, special 26; (7735,4968) has no closing record (OPEN there); 9 rows conditional on a rank being 1; 25 on a numerical L-value (the 2 above and 23 repository values for `m ≤ 3000`, PARI's with a tail bound) |
| `census-wide/lvalue/` | `NOTES.md` | the L-value above for w·f0·f8 of (6557,324) (the repository's `lball_exact` arithmetic, streamed over 1.21·10⁹ terms; conditional on PARI's `ellap`; the six sums from one run on the session's machine), and a rerun of `census_lvalue_big.py` for (700,529) that reproduces the committed record |
| `open-ratio/` | `WRITEUP.md` | (627,364) over `ℚ(i)`: 2-descent and φ-descents (rank `≤ 2`), the Sha ratios across the isogeny class, a cyclic 4-isogeny descent (bounds 3 to 5), and a new PARI implementation of the Cassels–Tate pairing on Sel², which vanishes on Sel² of `E₁` (as `scripts/second_descent/kctp*.py` found) and, in run 5, on the deciding value of Sel² of `E₄`, so Sel² with its pairing cannot prove rank 0 for `E₄` (a second implementation, unreviewed) |
| `roadmap-g/` | `g1/` … `g4/NOTES.md` (with the session's `REPORT_G1.md` … `REPORT_G4.md`) | run 6, Roadmap G: a stress test of the four-class reduction on 340,000 random ratios to `m = 10⁹` (no survivor) with an exact local test at odd primes (`g1/`); which primes empty the four-class set, test 6 of `COMPLETE.md` §7 answered **no** (`g2/`); the census of the two degenerate types for `11000 < m ≤ 13000` (4 open, 3 conditional; `g3/`); a point search on the square-type fibres to `m = 20,000`, a search, not a certificate (`g4/`) |
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

The prime and `a_p` arrays of the (6557,324) and (700,529) L-value runs (about 465 MB for (6557,324); `lvalue/r6557/run.sh`
rebuilds them in about 2 hours on 4 cores). (Run 5's first bundle said the `9000 < m ≤ 11000` run was not included; it
is now, see "Run 5, part 2".)

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

## Run 5, part 2 (committed in run 6)

The session's second bundle (F0, F5), reviewed by a separate AI agent (`reviews/review_x4.md`; its own checks are in
`referee_x4/`) and committed with these changes: no PARI stand-in for sympy (`table/class_step.sh` now calls the
repository's programs with sympy); the census table reads the class step from the repository's committed outputs (the
session's byte-identical rerun copies, 2.1 MB, are not committed); the label of (6903,5104) in the table corrected (rank 1
and the Demjanenko–Manin generator condition, not (a), (b)); `CENSUS11000.md`'s headline now lists the conditions; the
binaries are left out. Everything is tier D, written by an AI model and not read by a person.

- `census-wide/CENSUS11000.md` and its records (`s11_*`, `g11k*`, `dm11k*`, `pipe11k.*`, `route/*10512*`): the range
  `9000 < m ≤ 11000`, 6536 ratios. The review recomputed the stage counts and re-ran a sample of every stage under
  PARI 2.17.4 (10 second-descent ratios, 10 rank-0 inputs, 10 Demjanenko–Manin inputs, the route and effort-8 rank bounds
  of (10512,145), the class step for `10000 < m ≤ 11000` with sympy): all reproduce. (10512,145), (10657,10367) are
  correctly conditional.
- `census-wide/table/` (`NOTES.md`): the table above. The review recounted it with separate code and regenerated the
  congruence lists (the same 20,682 ratios), and checked 25 rows of every step kind and every hand-entered "special" row
  against its record. The 26 special rows are typed in from the write-ups; the script checks no record for them.
- `census-wide/lvalue/` (`NOTES.md`): (6557,324): L(E,1) ∈ [1.444085955623, 1.444085956081], so L(E,1) ≠ 0, **conditional
  on PARI's `ellap`** at the good primes ≤ 1.21·10⁹ (sampled against the repository's exact `a_p` certificate at
  50,600 primes, no disagreement) and on one run of the six sums outside CI. With the torsion step, which the review
  made exact (`referee_x4/tors6557.py`), this would exclude (6557,324) and (6881,6233) with PARI used only through `ellap` (and the six sums
  from that one run). The census table keeps their numerical label. (700,529): the rerun equals the committed record.

## Roadmap G (committed in run 6)

A bundle from a separate AI session (27 November 2026, UTC; run with PARI 2.19.1 and real sympy 1.14), reviewed by a
separate AI agent (`reviews/review_x6.md`; its own checks are in `referee_x6/`) and committed with the changes listed in
each folder's `NOTES.md`. Everything is **tier D** (computed; experiments, not counted), written by an AI model and not
read by a person. Nothing changes the census, which stops at `m ≤ 1000`.

- `roadmap-g/g1/` (`NOTES.md`): **no survivor among 340,000 random coprime ratios** (20,000 per dyadic range,
  `10⁴ ≤ m ≤ 10⁹`; 3 of them of the degenerate types): each of the other 339,997 is excluded by the class criterion of `COMPLETE.md` §7, or (63 of them) by the four-class descent
  at the primes of the supports and the good primes `≤ 150`. The review emptied the 63 again with the repository's own
  local test only. New: `exact_local.py`, the local test of `tuple_descent.Local` at an odd prime without a depth cap
  (residue symbols on the residues of a disc, Weil's bound for their sign patterns, roots shared in `ℤ[i]` closed in
  form); tested against the repository's test and against explicit local points, not proved. One degenerate control,
  (1109503,597671), keeps four tuples besides those of its degenerate points; the review found them **locally soluble
  everywhere** (explicit points at 2 and at every bad prime; the good primes by `exact_local.py` and a short argument);
  whether a rational point has one of these classes is open.
- `roadmap-g/g2/` (`NOTES.md`): **the bad primes alone do not always empty the four-class set**: they leave tuples for 101
  of the 997 ratios of test 4 and 836 of the 8,739 of test 7, each tuple with explicit local points at 2 and (found by
  the review) at every odd bad prime; good primes up to 181, all `≡ 1 (mod 4)`, empty all 9,736; irredundant kill sets;
  no residue-symbol law found and no conjecture stated.
- `roadmap-g/g3/` (`NOTES.md`): the 611 ratios of the two degenerate types with `11000 < m ≤ 13000` that pass the
  congruences, each with its closing step and label (`CENSUS_TABLE_G3.tsv`, the format of `census-wide/table/`): **open**
  (11165,7548), (12201,7360), (12204,9853), (12703,2767); **conditional on a rank being 1** (11441,1519), (11969,4529),
  (12007,1799); (12416,5313) by the rank-1 route, conditional on (a), (b). PARI rank bounds throughout; other types with
  `m > 11000` not covered. One row's label regenerated from the corrected table.
- `roadmap-g/g4/` (`NOTES.md`): on the 3,573 square-type ratios with `m ≤ 20,000`, none of the points listed (height up
  to 64 times the largest generator's, on the subgroup of PARI's points, wherever PARI's rank bounds meet) gives a
  solution or a non-degenerate `t` with four of the five forms squares. **A search, not a certificate** (PARI's ranks,
  saturation at primes `≤ 100`, floating-point heights; 1,741 ratios with a full search on E+ or E−); its outputs are
  PARI 2.19.1's and are not reproduced bit for bit under 2.17.4.
- Not committed: nothing of the bundle except the two input lists of `g2/` (copies of repository files) and four
  per-worker files of `g2/` (their lines are the `_all` files). The scripts find the repository through `MSS_REPO`, or
  else `../../..`; `g3/`'s pipeline calls `census-wide/` scripts that still name the session's folders ("Rerunning" above).
