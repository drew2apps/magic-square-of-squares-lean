# referee_z11: checks of review z11 (run 9) of t10's pair-curve Cassels–Tate test

Written by an AI model (Claude), review agent z11 of `../conductor-t10/` (t10, commit 714b45c, notes cc8e409); not
read by a person. Experiment, not counted; tier D; nothing here changes the census. The review itself:
`reviews/review_z11.md`. Paths are relative to this folder: repository root `../..` (or `MSS_REPO`), t10's folder
`../conductor-t10` (or `T10_DIR`), z7's folder `../referee_z7` (or `Z7_DIR`). Every script prints to stdout exactly
its `.txt` and then, on stderr, `# matches NAME.txt` (exit 0) or `# DIFFERS from NAME.txt` (exit 1); `--write`
rewrites the `.txt`. Outputs carry no timings and no PARI messages (check with
`python3 -I X.py 2>/dev/null | grep -v '^PARI stack' | diff - X.txt`). Runs are deterministic (fixed seeds;
`setrand(1)` before every PARI point search). PARI 2.17.4, sympy 1.14. Times are CPU times on a shared 2-core
container. **CI** marks those under 2 minutes.

**Own code** (`z11_pair.py`): forms and base points from `LAST_PAIR.md` §11, Washington's map (`wash_check.py`),
the roles of the three roots of `E_jk` found by a finite-field test (not t10's resolvent formula), and the class of a
tuple by the divisor form of the Kummer map, `c(τ) = (F_j(t0), ±d_j d_k D(t0), product)`, with **both signs tested**.
**Reused, read first:** z7's `my2desc.Curve` (complete 2-descent) and `ct_mine.Cov` (Cassels–Tate pairing). Two
changes are made in memory: the primes come from a PARI factorisation, with primality checked; and the real local
points are drawn in `(e1, e2)` and `(e3, e3 + 3(e3 − e1))`, because z7's draw misses them when the roots cluster.
G1's `exact_local.py` gives the residual tuples, as it does for G2 and t10.
**Where PARI is used:**
- integer and cubic factorisations, with the roots checked exactly;
- square roots mod p;
- z7's p-adic square roots, `hilbert` and `qfsolve`, with the points checked (z7's real place uses mpmath: a
  cross-check, not a certificate);
- `ellrank` and `elltors` only as finders of points, every point checked: in `forge_mine.py`, in `forge_t10.py`
  through t10's finder, and in `z11_pair.analyse(points=True)`;
- `rank0.py` uses no `elltors`, `ellap` or `ellrank`: the torsion bound comes from exact point counts.

| Script | Command | Time | CI | What it checks (claim) |
|---|---|---|---|---|
| `rank0.py` | `python3 -I rank0.py` | 1 s | CI | **The five ratios of the 201 with a pair curve of dim Sel² = 2 have no four-form point**: (5623,1017) by `E_08`; (3320,2303), (9191,7688), (9532,4445) and (14885,6076) by `E_LU`. In each case dim Sel² = 2, so the rank is 0. The gcd of the exact counts #E(F_p) over the good primes 3 ≤ p ≤ 199 is 8, so the torsion has order at most 8. C(Q) is therefore exactly its 8 obvious points, over t = 0, ±1, ∞, and at none of them are the four forms squares together. So the 4 tuples that t10 could not pairing-kill are removed by rank 0 |
| `z11_pair.py` | `python3 -I z11_pair.py` | 42 s | CI | Degenerate controls (4,3), (12,5), (15,8), (23,7), (21,20) and (31,17), with the pairing always computed: every known tuple survives every pair curve, with both signs. (12,5)'s two extra tuples (good prime 29) are killed on 0L, 0U, 8L and 8U. Also a library: `python3 -I z11_pair.py m,n` prints any ratio |
| `z11_compare.py` | `python3 -I z11_compare.py` | 40 s | CI | My code against t10's committed lines for (465,416), (5623,1017), (3320,2303) and (14885,6076): 108 (tuple, curve) verdicts agree, and every class and killing b that t10 prints is equal; both signs agree |
| `z11_compare.py --set mid` | `python3 -I z11_compare.py --set mid` | 66 s | CI | The same for (2999,1081), (3073,412) and (3143,2567), from the middle of the recovered `t10_997_a.txt` (`z11_compare_mid.txt`) |
| `z11_compare.py --set big` | `python3 -I z11_compare.py --set big` | 228 s | | The same for (9751,5084), which has 216 tuples: 1,296 verdicts (`z11_compare_big.txt`) |
| `z11_compare.py --set batch30` | `python3 -I z11_compare.py --set batch30` | 351 s | | The same for 30 ratios drawn by `random.Random(11)` among those with ≤ 40 tuples: 960 verdicts (`z11_compare_batch30.txt`) |
| `recount.py` | `python3 -I recount.py` | 0.1 s | CI | t10's counts from its outputs, own parser: 201 ratios and 1,766 tuples (598 + 1,168), equal to G2's "N left" for each ratio; 1,762 pairing-killed; per-curve counts and first killing primes; **8,930 kill lines carrying 18,222 local −1 terms** (an odd number on each line); 4,976 distinct (ratio, curve, class) kills; the largest first killing prime equals the largest good prime of G2's kill set for all 201; t10's 997 selection and seed-1 sample reproduced; the five curves of dim Sel² = 2 |
| `wash_check.py` | `python3 -I wash_check.py` | 0.5 s | CI | Washington's Thm 2.17 and its inverse, symbolically modulo v² = g(s): residues 0 |
| `pmat.py` | `python3 -I pmat.py` | 4 s | CI | Full pairing matrices for (5623,1017): `E_08` dim 2; `E_0L` dim 6, matrix identically +1, so rank ≤ 4; `E_LU` dim 5, matrix rank 2, so rank ≤ 1 (t10's "pairing trivial on Sel²") |
| `forge_mine.py` | `python3 -I forge_mine.py` | 70 s | CI | My method never kills the class of a rational point. Tested on every subset sum of the points found, on all six curves of (7,1), (4,3) and (465,416) (all 32 classes of Sel² on the mixed curves of (7,1) and (4,3)), with the same classes recomputed by the divisor formula on C, both signs |
| `forge_t10.py` | `python3 -I forge_t10.py` | 16 s | CI | t10's committed method, imported and not changed, with PARI's points withheld so that point classes and the known tuples of degenerate ratios go through its pairing: nothing killed on (7,1), (4,3), (12,5), (627,364), (465,416) or (14885,6076) |
| `h1_mine.py` | `python3 -I h1_mine.py` | 10 s | CI | H1's control with my code and other base points (m/n, −m/n): (1109503,597671). The known tuples survive. The four extra tuples are killed on exactly the curve that the VERDICT lines of `../roadmap-h/h1/h1_pairs.txt` name, with both signs |
| `oom_check.py` | `python3 -I oom_check.py` | 335 s | | t10's committed script (it only prints) rerun on indices 35–40 and 74–78 of `g2_997_all.txt`: every block is identical to the one in `t10_997_a.txt`, the file recovered from the run that was killed for lack of memory |

`_z11common.py` holds the paths and the output helper. t10's own controls were also rerun as they are, from
`../conductor-t10`: `python3 -I t10_pairct.py --h1control` and `--negative` print exactly `t10_h1control.txt` and
`t10_negative.txt`, apart from PARI's stack message.
