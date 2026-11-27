# G1: stress test of the reduction on random ratios to m = 10⁹ (Roadmap G, run 6)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (x6, `reviews/review_x6.md`, its
own checks in `../../referee_x6/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.** Nothing here
changes the census, which stops at `m ≤ 1000`. `REPORT_G1.md` is the session's report, unchanged; the corrections
below take precedence over it.

## What it claims (as corrected)

- **No survivor among 340,000 random coprime ratios** (20,000 per dyadic range, `10⁴ ≤ m ≤ 10⁹`, seed `1000 + k`). Each
  is excluded either by the repository's class criterion (`class_selmer.empty`, `class_selmer_linalg`, unchanged;
  `COMPLETE.md` §7) or, for the **63** that keep a class, by the four-class descent of `LAST_PAIR.md` §12 at the primes of
  the supports and the good primes `≤ 150` (2-adic test at depth 12). Computed; tier D.
- The review ran the 63 again with **the repository's own `tuple_descent.Local` only** (primes up to `2·10⁵`, larger
  support primes left out): nothing is left for 63 of 63 (`../../referee_x6/kept63_repo_only.txt`). So this result
  does not depend on the new local test `exact_local.py`.
- **The control (1109503,597671)**: twice-square type, `m² + n² = 2·891125²` with `891125 = 5³·7129` (the report's
  "2·(5·7129)²" is wrong). The descent run to good primes 20,000 keeps six tuples: the two of its degenerate points
  `t = 0, ∞` and four more, `(dL, dU) = (2, 959)` or `(1918, 1)` with either pair `(d0, d8)`. The review
  (`../../referee_x6/control_1109503.txt`) found explicit local points for each of the six at 2 and at every odd bad prime
  (18 bad primes, up to `1.19·10¹⁵`), and ran `exact_local.py` at the 23,397 good primes `≡ 1 (mod 4)` between 20,000 and
  600,000: no tuple fails. Above 595,717 every tuple passes at every good prime: there the eight roots of the components
  are distinct modulo `p`, and the Weil step of `exact_local.py` puts every sign pattern on the other residues; at a good
  prime `≡ 3 (mod 4)` every candidate tuple has the same labels as the degenerate ones. So **the four extra tuples are
  locally soluble everywhere** (computed, with that short argument; tier D). Whether a rational point has one of these classes is **open**. (Run 7: none does — `y² = w·f0·f8` has rank 0 and the pair curves `E_0L`, `E_0U` rank exactly 3; `../../roadmap-h/h1/NOTES.md`, review z8.) The report's "the descent can leave classes that no point explains" goes beyond this: read "the
  four-class descent, with every local condition, keeps four classes besides those of the degenerate points; whether a
  rational point has one of them is open".

## Corrections to `REPORT_G1.md`

1. "m² + n² = 2·(5·7129)²": it is `2·(5³·7129)² = 2·891125²`.
2. The proposed LEDGER text says the descent used "exact_local.py, no prime left out". In G1 only the primes of the
   supports (of `2mn` and `2(m² − n²)`) and the primes up to `pmax` are tested; the primes of `m² + n²` and of the octic
   above `pmax` are left out (that can only keep more tuples, so "nothing is left" is unaffected). "No prime left out"
   is true of G2, not of G1.
3. "Patterns ... settled by Weil's bound for p > 5·10⁶ and by listing below that" (also in the docstring of
   `exact_local.py`): in the code the Weil bound is tried first at every prime and decides whenever its inequality holds
   (for eight roots from `p = 595,717`, for one root from `p = 3`); listing is used only below that, and `5·10⁶` is the
   cap asserted for listing, which can never be reached (`../../referee_x6/weil_check.txt`). `MixedLocal` is unused.
4. "The descent can leave classes that no point explains": see above.

## `exact_local.py`, as reviewed

The local test of `tuple_descent.Local` at an odd prime, computed without a depth cap: each form is linear on a disc;
a constant class where `v(a) < v(b)`; otherwise residue symbols of `x − r` on the residues that are not roots, whose
sign patterns are listed or come from Weil's bound `|Σ χ(f)| ≤ (d − 1)√p` (this step can only add patterns, so an error
in it would make a set too large, never exclude a ratio); roots shared exactly (tested in `ℤ[i]`) are closed in form,
other roots refined. The review read it against `LAST_PAIR.md` §14 (R1, R3, R4), checked the Weil inequality and the
count of characters, and compared it on 107 (ratio, prime) cases (primes from 3 to `2⁶¹ − 1`, every odd bad prime of
six ratios) with the repository's test (equal on all 90 cases where that can run) and with labels of explicit points
computed by separate code (the exact set is exactly the set of labels reached, in all 107), and with listing forced in
place of the Weil step on 593 cases (no difference): `../../referee_x6/local_compare.txt`, `weil_check.txt`. Tested,
not proved; at `p = 2` the repository's one-sided test (depth 12) is kept.

## Rerunning (PARI 2.17.4, as in CI, reproduces the committed outputs)

From this folder (`F = grep -v '^PARI stack size'`):

| Command | Output | Time (2-core container) |
|---|---|---|
| `python3 test_exact_local.py` | `test_exact_local.txt` (asserted) | 3 s |
| `python3 g1_stress.py 13 20000 200` | `g1_k13.txt` (asserted) | 17 s |
| `python3 g1_stress.py 20 20000 200` | `g1_k20.txt` (holds the control) | see `OUTPUTS.md` |
| `python3 g1_stress.py 29 20000 200` | `g1_k29.txt` (rerun by the review: identical) | 13.5 min |
| `python3 g1_stress.py --negative` | `g1_negative.txt` | a few s |
| `python3 g1_summary.py` | `g1_summary.txt` (compare with `diff`) | under 1 s |
| `bash run_g1.sh 13 … 29` | all ranges, `--write` | about 1.5 h on two cores (`times.txt`) |

`exact_local.py` finds the repository through `MSS_REPO`, or else `../../..`. Factorisations in `g1_stress.py` are
PARI's with every factor proved prime by `isprime`.
