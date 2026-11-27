# G4: a point search on the square-type fibres (Roadmap G, run 6)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (x6, `reviews/review_x6.md`, its
own checks in `../../referee_x6/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`, as a search:
computed (numerical), not a certificate.** `REPORT_G4.md` is the session's report, unchanged; the notes below take
precedence over it.

## What it claims (as corrected)

On the 3,573 coprime ratios with `m² + n²` a square and `m ≤ 20,000` (`square_20000.json`, no congruence filter), on E±
and the 12 quotient curves of `scripts/g1_full/global.py`, wherever PARI's rank bounds meet (rank `≤ 3`): the points of
canonical height up to 64 times the largest generator's, on the subgroup that PARI's points generate (saturated at primes
`≤ 100`), were listed, mapped to `t` and tested exactly. **None of the points listed gives a solution, or a
non-degenerate `t` with four of the five forms squares** (every recorded hit is a degenerate `t` with five squares).
Coverage: 1,741 ratios have a full search on E+ or E−, 1,473 only a partial one (fewer points than the rank bound), 359
neither, 23 no curve at all.

It is **a search, not a certificate**: PARI's rank bounds (where the lower bound exceeds the points found, the rank is
conditional and the search partial), generators saturated only at primes `≤ 100`, floating-point heights and a
floating-point Fincke–Pohst bound, and no bound tying the height of a solution's point to the search bound.

## Corrections to the session's text

- The hand-off's bottom line ("no solution, and no non-degenerate t with four of the five forms squares, on the 3,573
  ratios") and `REPORT_G4.md`'s proposed LEDGER sentence ("No solution, and no non-degenerate t ...") are about the
  points listed, not about the ratios: read "none of the points listed gives ...". `g4_summary.txt`'s line
  "SOLUTIONS: none" is the same statement.
- The planted test (the point of `t = ±1` found on every curve that carries it; 15,724 of 15,724) exercises the whole
  search path but weakly, since that point is among, or is appended to, the generators. The review checked the
  enumeration itself: on the curves searched for 7 ratios, every combination `Σ aᵢgᵢ`, `|aᵢ| ≤ 4`, of height below the
  bound (9,495 points) is listed (`../../referee_x6/g4_enum.txt`).

## PARI versions

The committed search outputs (`g4_*_3.json`, `pass2_*_2.jsonl`) are PARI 2.19.1's and are not reproduced bit for bit
under PARI 2.17.4 (CI's): of 5 ratios rerun, 2 are identical and 3 differ on one quotient curve each (PARI's points;
on (195,28) the curve w·f8·K is searched under 2.19.1 and not under 2.17.4, so two degenerate hits fewer), with no
solution in either (`../../referee_x6/samples/compare.txt`). Only the bookkeeping (`g4_summary.py`), the negative mode
and the input (`gen_square.py`) are checked by CI.

## Rerunning

From this folder: `python3 g4_summary.py` (under a second), `python3 g4_search.py --negative` (about 4 s),
`python3 gen_square.py 20000` (about 9 s), each asserting its committed output. `bash run_g4.sh i 3` (first pass, about
4.6 h per worker on three workers) and `bash run_g4_pass2.sh i 2` (second pass on E±) remake the search under the PARI
at hand. `g4_search.py` finds the repository through `MSS_REPO`, or else `../../..`.
