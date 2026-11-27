# referee_z8: review z8 of roadmap H (H1, H3)

Written by an AI model (Claude), as the hostile review z8 of `experiments/roadmap-h/h1` and `h3`; not read by a person.
Experiment, not counted; nothing here changes the census. The report is `reviews/review_z8.md`; "claim n" below is the numbering of that report.

Layout: meant to sit at `experiments/referee_z8/`. The scripts find the repository through `MSS_REPO`, or else `../..`;
G1's folder as `<repo>/experiments/roadmap-g/g1` (or `G1_DIR`), H1's and H3's as `<repo>/experiments/roadmap-h/h1`, `h3`
(or `H1_DIR`, `H3_DIR`). Each script prints exactly its `.txt` when run without arguments (stdout only; PARI writes
nothing to stdout), and `--write` rewrites it. Check with `python3 X.py 2>/dev/null | diff - X.txt`.
Times: one core of a shared 2-core container, PARI 2.17.4 with cypari2, sympy 1.14, numpy. Deterministic (fixed seeds).
CI: "yes" means under 2 minutes.

| Script | Output | Time | CI | Backs |
|---|---|---|---|---|
| `my_local.py` | `my_local.txt` | 2 s | yes | claim 3: my own local test at good p ≡ 1 mod 4 (from LAST_PAIR.md §11–12) equals H3's `h3_weil.json` on every good residue at 15 primes from 17 to 1153; library for the next scripts |
| `scan_losses.py` (driver; `scan_losses.py LO HI` for other ranges) | `scan_losses.txt` (every good residue, 1153 ≤ p ≤ 5000) | 731 s | no | claim 1 fails: good primes above 1153 remove vectors (last loss to 5000: 3853; 90 primes above 1153 with a loss) |
| `brute_p2.py` | `brute_p2.txt` | 12 s | yes | claim 1 fails, third method: brute force over ℙ¹(ℤ/p²), (5,1) and (17,1) at p = 1181 |
| `methA_above.py` | `methA_above.txt` | 2 s | yes | claim 1 fails with the bundle's own method A (G1's `exact_local`) at 1181, 1201, 1249 |
| `p3mod4.py` | `p3mod4.txt` | 2 s | yes | claim 2 fails: good p ≡ 3 mod 4 can empty the set (p = 7, 23, 31, 47, 79, 103, 127, 151; none from 167 to 1151) |
| `p3_methA.py` | `p3_methA.txt` | 2 s | yes | claim 2 with G1's `exact_local`: empty sets at good p ≡ 7 mod 8 |
| `h3_counts.py` | `h3_counts.txt` | 2 s | yes | claim 3: H3's numbers recounted from its json files (91 primes, not 96; the rest as stated) |
| `badprimes_left.py` | `badprimes_left.txt` | 6 s | yes | claim 5: the primes of m² + n² and of the octic above 20,000 that H1's scan skipped; `exact_local` there cuts every extra tuple of 10 of the 24 ratios (5 pairs, the "open" pair among them) |
| `octic_cut.py` | `octic_cut.txt` | 1 s | yes | claims 5, 8: a short exact argument at an octic prime (two root residues meet with valuation 1, so e_ja + e_jb = χ(κ)); the "open" pair (11875080,6497849) at p = 5462461 |
| `octic_cut_all.py` (reads `badprimes_left.txt`) | `octic_cut_all.txt` | 2 s | yes | claim 5: the argument cuts exactly where `exact_local` cuts, at all 90 octic primes |
| `full_descent.py` | `full_descent.txt` | 405 s | no | claim 5: G1's descent with every bad prime on H1's 24 ratios: 14 keep H1's extra tuples (7 pairs), 10 keep only the known ones |
| `phantom_mid.py` | `phantom_mid.txt` | 373 s | no | "local everywhere" at good primes for H1's 24 ratios: no extra tuple cut at any good p ≡ 1 mod 4 < 1200 or in (20000, 67583); 67,583 is the Hasse–Weil threshold for the genus-129 cover |
| `witness2.py` | `witness2.txt` | 1 s | yes | "local everywhere" at 2: explicit 2-adic points for all 108 known and extra tuples of the 24 ratios; three forged tuples get none |
| `witness_odd.py` | `witness_odd.txt` | 16 s | yes | at the odd primes of m² + n² and of the octic (outside the supports), explicit points for the extra tuples of the 7 true phantom pairs (all but p = 5) |
| `witness5.py` | `witness5.txt` | 7 s | yes | the same at p = 5 (mod 5⁸) |
| `witness_supp.py` | `witness_supp.txt` | 1 s | yes | at the odd support primes (p ∣ mn(m² − n²)), explicit points for every known and extra tuple of the 7 true pairs |
| `rank0_check.py` | `rank0_check.txt` | 265 s | no | claim 6, cross-check with PARI (not a proof): rank 0 on the isogeny class of y² = w·f0·f8 and t over its torsion ∈ {−1, 0, 1, ∞} at the 7 ratios of H1 part 6 |
| `z8common.py` | — | — | — | paths, output, the 24 ratios read from H1's committed files |

Not packaged (reruns of the bundle, reported in the review): H3's and H1's cheap checks, `h1_pairs.py`, and
`h1_rank0_exact.analyse` / `h1_pairs.analyse` on single ratios, all matching their committed outputs on PARI 2.17.4.

What these rest on: the definitions of LAST_PAIR.md §11–12; Legendre symbols and exact integers; PARI's `factor` with
`isprime` for factorisations; G1's `exact_local.py` where named (methA_above, p3_methA, badprimes_left, full_descent);
Hasse–Weil (quoted) for the threshold 67,583; PARI's `ellrank`/`elltors` in `rank0_check.py` only.
