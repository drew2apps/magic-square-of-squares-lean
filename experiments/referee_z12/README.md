# referee_z12: checks of review z12 (run 10) of the closed pairing formula on `E_08`, `E_LU`

Written by an AI model (Claude), review agent z12 of `../research-pairing/` (a separate AI research session's
`PAIRING_FORMULA.md`); not read by a person. Experiment, not counted; tier D; nothing here changes the census. The
review itself: `reviews/review_z12.md`. The scripts were written in the reviewer's scratch folder and committed with
one change: their absolute paths are replaced by paths relative to this folder (`../conductor-t10` for the pilot's
files and `t10_pairct.py`, `../research-pairing` for `ct08.py` and `parse.py`), through a first line `HERE = ...`.
PARI 2.17.4 (cypari2 2.2.9), sympy 1.14.

| Script | What it does | Output | CI |
|---|---|---|---|
| `fulltest.py` | the full test: on a curve `E_08` or `E_LU` of a pilot ratio, rank_second's `full_selmer2` and `DFull.pairing` (through `t10_pairct.PairData`, `DFullX`) against `ct08.pairing`, for `a` = the class of every residual tuple (both signs) and every `(1, δ, δ)` in `Sel²`, `b` over the basis of `Sel²`; seeds 1 and 7; one JSON line per curve | `full0.jsonl` (146 ratios on `E_08`), `full1.jsonl` (133 on `E_LU`); job lists `jobs0.txt`, `jobs0b.txt`, `jobs1.txt`; logs `*.log` | no (about 90 minutes); `rerun_compare.py` reruns 6 curves |
| `fulltest_lowmem.py` | the same with a memory cap, for the restart after the out-of-memory stop | appended to `full0.jsonl` | no |
| `rerun_compare.py` (the conductor's) | reruns `fulltest.run` from scratch on 6 curves (4 `E_08`, 2 `E_LU`) and compares every row with the committed records | `rerun_compare.txt` | `descents` |
| `analyse.py` | the counts of the review's bottom line from `full[01].jsonl`: 28,918 distinct triples, 0 disagreements; seeds 1 and 7 equal; linearity over all of `Sel²` | `analyse.txt` | `descents` |
| `classcheck.py` | the class formula of §1 on every residual tuple of the sample | `classcheck.txt` | `descents` |
| `zct.py`, `zlocal.py` | an independent evaluator of Cassels' local terms (own local points, all 8 sign choices, the real place included, local constancy checked); the odd-`p` table of §2 against it on 66 coprime pairs, 15 with `m < n`, and non-coprime pairs (where it fails) | `zlocal.txt` | `descents` |
| `brone.py` | the `b_r = 1` product formula against rank_second (on the complete `full[01].jsonl`: 14,110 + 7,456 triples; the review's text quotes the counts at the time of its check) | `brone.txt` | `descents` |
| `claims.py` | the counts of `PAIRING_FORMULA.md` per kill line and per distinct (ratio, class, `b`), and the other assignment (conv 1) on the kill lines | `claims.txt` | `descents` |
| `convtest.py` | all six assignments of `(b_r, b_m, b_*)` to the three tangent forms, by `zct.py`, on the first 100 curves (`part0.jsonl`, `part1.jsonl`) | `conv_part.txt` | no |
| `seedtest.py` | six more seeds on 2,826 triples | in the review | no |
| `zparse.py` | reads the pilot's files and rebuilds the residual tuples | — | — |

Rerun a check: `python3 -I X.py 2>/dev/null | grep -v 'PARI stack' | diff - X.txt` from any folder.
