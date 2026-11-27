# referee_z9: checks of review z9 (run 7, roadmap H)

Written by an AI model (Claude), review z9 of roadmap H (H2, H4, H6, and G3's block 13000 < m ≤ 14000); not read by a
person. Experiment, not counted; nothing here changes the census. The review itself: `reviews/review_z9.md`. Paths are relative to this folder (repository root `../..`, or `MSS_REPO`). Every script prints to
stdout exactly its `.txt` and then, on stderr, `# matches NAME.txt` (exit 0) or `# DIFFERS from NAME.txt` (exit 1);
`--write` rewrites the `.txt`. PARI's own messages are kept off stdout. Outputs carry no timings. Run with PARI 2.17.4,
sympy 1.14 (CI's versions). Times are CPU times on a shared 2-core container. **CI** marks those under 2 minutes.

| Script | Command | Time | CI | What it checks (claim) |
|---|---|---|---|---|
| `recount.py` | `python3 -I recount.py` | 8 s | CI | H2's counts from its outputs, own parser: 15000 < m ≤ 30000 every m once, coprime = φ(m), 205,172,879 non-degenerate, 30,434 class, 2,807 need a good prime, largest P 229, the per-block table; hunts: 2,251 lines = 2,234 distinct ratios (17 repeats), 273 with m ≤ 30000, so **2,131** class-keeping ratios with m > 30000 (1,961 + 170), not 2,421; 3 class timeouts at k = 36; ranges reached; 31,322 core-seconds (0.58 core-h per block) |
| `degcount.py` | `python3 -I degcount.py` | 2 s | CI | H2's degenerate counts (4,579; per block and per m) from Pythagorean parameters |
| `sample_h2.py` | `python3 -I sample_h2.py` | 35 s | CI | 23 class-keeping ratios of H2 (20 random, the 3 largest P) and 20 uniform ratios: class status by the repository's criterion with sympy factorisations, P from G2's `g2_killsets.run`: all agree |
| `sample_h2b.py` | `python3 -I sample_h2b.py` | 123 s | | 20 ratios with points everywhere (class or not) and the whole count lines of m = 25789 and 29999, recomputed: all agree |
| `h4_one.py` | `ulimit -v 6500000; python3 -I h4_one.py` | 350 s (about 2 min per ratio; 15 min wall on the loaded container) | | H4: G1's three out-of-memory ratios rerun through `h4_split.line`: lines identical to `h4_split.txt` |
| `h6_ident.py` | `python3 -I h6_ident.py` | 1 s | CI | H6: h6_split's forms against `G1_FULL_SYSTEM.md`; σ, τ and the four genus-2 pieces F0F4L, F0F4U, F8F4L, F8F4U as one curve; the partner map t′ = (t + 1)/(t − 1) carries the forms of (m + n, m − n) to L, U, 4F4, 16F0, 16F8 of (m, n): the covers of (627,364) and (991,263) are isomorphic |
| `c5.py` | `python3 -I c5.py` | 1 s | CI | H6: y² = L·U of (627,364) has the Jacobian of y² = F0·F8 of (991,263), torsion ℤ/2 × ℤ/4, root number −1, and exactly its 8 obvious points at t = 0, ∞, ±1, so rank 0 is not excluded by counting; the PARI classes whose lower bound exceeds the points returned (parity only) |
| `g3_family_z9.py` | `python3 -I g3_family_z9.py` | 11 s | CI | G3 11000 < m ≤ 14000: the family by brute force over all pairs (916, identical to `fam_11000_30000.json` and to the table's rows), counts by step and type, block 13000 against `stage_13000.json`, the row (12007,1799) with run 6's label, every record of the six open ratios and their partners, the route row (13103,4271) through (8687,4416) |
| `g3_samples.sh` + `g3_samples_compare.py` | `bash g3_samples.sh [SCRATCH]`, then (done by the script) `python3 -I g3_samples_compare.py SCRATCH` | 193 s (8.5 min wall on the loaded container) | | G3 block 13000 reruns with the repository's programs: 3 second descents, 3 rank-0 inputs, 2 Demjanenko–Manin inputs with partners, against the committed records (identical apart from timings, the last digit of one floating height and one PARI rank bound [1,3] against [3,3], same status) |

Scripts of the bundle that this review ran as they are (PARI 2.17.4), all reproducing the committed outputs: from
`../roadmap-h/h2`, `python3 -I h2_exhaust.py 14900 15000 0 1` (44 s, `# matches`), `h2_control.py`,
`h2_control.py --negative`, `h2_summary.py`; from `../roadmap-h/h4`, `h4_split.py --negative`; from
`../roadmap-h/h6`, `h6_split.py 627 364`, `h6_split.py 991 263`, `h6_ranks.py --negative`, `h6_richelot.py --negative`,
`h6_richelot.py --controls` (about 70 s), `h6_richelot.py 627 364` (155 s), `h6_ranks.py 627 364` (82 s; differs only in the
line that names the PARI version); from `../roadmap-g/g3`, `g3_table.py`, `g3_table.py --negative`.

`_common.py` holds the paths and the output helper. `g3_samples_in/` holds the inputs of `g3_samples.sh`.
