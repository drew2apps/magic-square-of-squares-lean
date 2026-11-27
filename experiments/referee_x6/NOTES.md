# referee_x6: the review's own checks of Roadmap G (run 6)

Written by an AI model (review agent x6, run 6); not read by a person. The review is `reviews/review_x6.md`. Every
script here only reads committed files and prints; none writes a file (`samples/samples.sh` writes into a scratch
folder it is given). Run from anywhere with `python3 -I`; the repository root is found from the script's folder (or
`MSS_REPO`). `mylocal.py` holds the review's own local labels at explicit points (its own square root of −1 modulo
`p⁶⁰`, valuations and residue symbols; points at random in both charts and close to every `p`-adic root; at 2 the
repository's `key2`); it does not use `exact_local.py` or `tuple_descent.Local` for labels. Times on a 2-core container.

| Script | What it checks | Output | Time |
|---|---|---|---|
| `local_compare.py` | `g1/exact_local.py` against the repository's `tuple_descent.Local` (where it can run) and against labels of explicit points, on 107 (ratio, prime) cases, primes 3 to `2⁶¹ − 1` | `local_compare.txt` | 1 min |
| `local_mutation.py` | five mutations of `exact_local.py` (in memory) that would make a set too small, each caught by explicit points on some of 24 cases; the unmutated code on none | `local_mutation.txt` | 7 s |
| `weil_check.py` | the thresholds of the Weil step, and 593 cases with listing forced in its place | `weil_check.txt` | 15 s |
| `kept63_repo_only.py` | G1's 63 kept ratios, emptied by the four-class descent with the repository's local test only | `kept63_repo_only.txt` | about 10 min |
| `control_1109503.py` | the control (1109503,597671): arithmetic, its six tuples, explicit points at 2 and at every bad prime, the good primes `≡ 1 (mod 4)` up to 600,000 | `control_1109503.txt` | about 14 min |
| `killsets_witness.py` | 20 kill sets of G2 (empty; irredundant, each drop proved by explicit points); 10 "bad primes alone fail" ratios witnessed at every bad prime | `killsets_witness.txt` | 1 min |
| `badalone_witness_all.py` | every one of G2's 937 "bad primes alone fail" ratios: every tuple left has explicit points at 2 and at every odd bad prime | `badalone_witness_all.txt` | about 30 min |
| `g3_family.py` | G3's family for `11000 < m ≤ 13000` by brute force over all pairs | `g3_family.txt` | 1 min |
| `g3_sample.py` | 10 random rows and the 11 open or special rows of `CENSUS_TABLE_G3.tsv` against the raw records | `g3_sample.txt` | under 1 s |
| `g4_enum.py` | G4's enumeration lists every small combination of the generators, not only the generators | `g4_enum.txt` | 1 min (PARI-dependent) |
| `samples/samples.sh`, `samples/compare.py` | reruns under PARI 2.17.4 of G3 stages (4 second descents, 2 rank-0 inputs, 2 Demjanenko–Manin inputs, the route on (12416,5313)) and of G4's search on 5 ratios, compared with the committed records | `samples/*`, `samples/compare.txt` | 5 min (reruns); under 1 s (comparison) |

Outputs with times in them (`kept63_repo_only.txt`, `killsets_witness.txt`, `control_1109503.txt`) are logs of the runs
made for the review; rerun, only the times change.
