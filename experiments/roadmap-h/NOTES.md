# Roadmap H (run 7): what holds after review

Written by an AI model (Claude); not read by a person. **Experiments, not counted; tier D of `LEDGER.md`.** Nothing here
changes the census, which stops at `m ≤ 1000`. The bundle came from a separate AI session (Fable) and was reviewed by
three separate AI agents, each of which had not seen the work being produced: z7 (H7, `reviews/review_z7.md`), z8 (H1,
H3, `reviews/review_z8.md`) and z9 (H2, H4, H6, G3's block 13000, `reviews/review_z9.md`); their own checks are in
`../referee_z7/`, `../referee_z8/`, `../referee_z9/`. `HANDOFF_ROADMAP_H.md` and every `REPORT_H<n>.md` are the
session's text, committed unchanged; each folder's `NOTES.md` takes precedence over them.

| task | after review |
|---|---|
| H1 phantoms | 7 true phantom pairs (not 12), all settled; the other 5 are removed by a local test at an octic prime that the scan skipped; no pair open (`h1/NOTES.md`) |
| H2 survivor hunt | holds for the exhaustive range `15000 < m ≤ 30000`; the hunts add 2,131 distinct ratios (not 2,421); "P ≤ 1153" withdrawn (`h2/NOTES.md`) |
| H3 local images | tables at 91 primes hold; the Hasse–Weil step is wrong (wrong curve: genus 129, not 17) and withdrawn; a correct complete local test goes to 67,583 for good primes ≡ 1 (mod 4), and to 1154 for good primes ≡ 7 (mod 8) (z8's argument, Hasse–Weil quoted) for good primes ≡ 1 (mod 4), and to 1154 for good primes ≡ 7 (mod 8) (z8's argument, Hasse–Weil quoted) for good primes ≡ 1 (mod 4), and to 1154 for good primes ≡ 7 (mod 8) (z8's argument, Hasse–Weil quoted) (`h3/NOTES.md`) |
| H4 | holds (`h4/NOTES.md`) |
| H5 Θ₄ conic points | not done by the session (the scripts were not on `main`); done in run 8 by a build agent here: `Θ₄ ≠ 0` on `E⁺` of (945,248), on its `E′` side (`LEDGER.md` tier C; `reviews/review_z10.md`) |
| H6 (627,364), (991,263) | holds with gaps; the two covers are isomorphic; some genus-1 ranks rest on parity (`h6/NOTES.md`) |
| H7 g4·K, E04 | holds: rank 1 exactly for four g4·K curves, (6557,324) through its partner, E04 bound for 90 of 91 tier-B ratios (`h7/NOTES.md`) |
| H8 | not run; its purpose is met through (6557,324)'s partner (`h7/NOTES.md`). In `REPORT_H8.md` read "the L-value (conditional on PARI's `ellap`, `../census-wide/lvalue/`)" for "the numerical L(E,1)", add the Cassels–Tate pairing as computed by `rank_second.py` to what H7 rests on, and note that run 6 already summed the series for this curve (1.21·10⁹ terms, `../census-wide/lvalue/`), so only the PARI-free `aₚ` of H8 is missing |

The session ran PARI 2.19.1; every CI line proposed here was rerun under CI's PARI 2.17.4 and matches. Changes at integration (checked against the bundle): `h2/h2_control.py` falls back to `scripts/cover_descent/structure/class_selmer_15000_list.txt`; a note heads the docstrings of `h3/h3_weil.py` and `h3/h3_sizes.py` (96 → 91); a 'Read NOTES.md first' line heads the hand-off and each report. Nothing else in the bundle's `h1`–`h8` was changed. The bundle's
own copy of G3's folder was not used (it had the session's old paths and an old copy of the `1000 < m ≤ 11000` table);
the block's records were added to `../roadmap-g/g3/` and the table rebuilt there.

In the hand-off's "Ideas tried", the bullet on Hasse–Weil ("It worked, and it bounds H2's P") is withdrawn with H3's step, and the Cassels–Tate kernel test's "kills the extras at 302617 and 550548247" holds at 302617 only: (550548247,487061057) is one of the five pairs removed by a local test at an octic prime.
