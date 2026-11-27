# referee_z7: independent checks of roadmap H's H7 and of `scripts/second_descent/rank0_ct.py`

Written by an AI model (review agent z7, 27 November 2026); not read by a person. Experiment, not counted; nothing here
changes the census. The review itself is `reviews/review_z7.md` (run 7) (verdicts: H7's rank step holds; rank0_ct holds).

The code here was written by the reviewer: `my2desc.py` from the textbook before reading the repository's descent,
and `ct_mine.py` from the reviewer's own derivation of the Cassels–Tate pairing, following the *approach* of
`reviews/review_t2.md` (not its code, and not `scripts/g1_full/rank_second.py`'s). The drivers compare their results
with the committed outputs of H7 (`../roadmap-h/h7`) and of `rank0_ct.py`. Paths: repository root `../..`, H7 at
`../roadmap-h/h7`, H1 at `../roadmap-h/h1` (override with `MSS_REPO`, `H7_DIR`, `H1_DIR`).

Every script prints exactly its `.txt`. Check with `python3 -I X.py 2>&1 | grep -v '^PARI stack' | diff - X.txt`.
Runs are deterministic: fixed seeds, and PARI only where each header says. PARI's random point search is used only by
the `--find` modes, which wrote `z7_points.json`; the points are checked exactly on use.

Times were measured on a shared 2-core machine with `nice` (PARI 2.17.4, sympy 1.14). **CI** marks the runs under 2 minutes.

| script | output | time | CI | what it backs |
|---|---|---|---|---|
| `my2desc.py` | `my2desc.txt` | 1–3 s | CI | Self-test of the complete 2-descent: dim Sel² on 31 curves y² = x³ − n²x equals R + s + 2 of PARI's ellrank, recorded at review time (no PARI here) |
| `ct_mine.py` | `ct_mine.txt` | 8–24 s | CI | Self-test of the Cassels–Tate pairing: the pairing rank equals PARI's s, and the bound equals R, on 7 congruent-number curves (n = 17, 73, 205 with Ш[2] ≠ 0) |
| `velu_check.py` | `velu_check.txt` | 2–3 s | CI | H7 claim 1: Vélu's 2-isogeny formula as a polynomial identity; the shift of g4·K; 16b a square for all m, n; E′ = X(X − 16Smn(m−n)²)(X + 16Smn(m+n)²) |
| `partner_cells.py` | `partner_cells.txt` | 5–9 s | CI | H7 claim 3: the partner map permutes all nine cells with one common factor; (6881,6233) → (6557,324) |
| `run_ct_g4k.py` | `run_ct_g4k.txt` | 45–61 s | CI | H7 claims 1, 5: on E′ of g4·K for 8 ratios, my dim Sel², pairing rank and bound equal H7's; rank 1 exactly at (12201,7360), (7735,4968), (6881,6233), (13348,3795), because the kernel is spanned by the images of T1, T2 and one point |
| `ct_rank0.py` | `ct_rank0.txt` | 27–39 s | CI | rank0_ct.py: the same curve for (609,200) and (809,409) (exact isomorphisms); the shift and E′; my pairing gives dim 6, rank 4, so rank 0, on all three full-2-torsion members (two seeds on E′); controls (861,620) and (440,279) agree |
| `classes.py` | `classes.txt` | about 70 s | CI | H7 claim 4: the classes of E04 (320,289), (1015,961), (1015,713) each have 2 full-2-torsion curves (bounds 1/1, 3/1, 3/3); the walk's "6 curve(s)" are 2 curves up to isomorphism; no D1 (16) or σ (12) curve has full 2-torsion in its class |
| `run_ct_e04.py` | `run_ct_e04.txt` | 3.5–4 min | no | H7 claim 4: 91 E04 lines, 66 curves; my results on E04 itself match H7's first curve on every curve except (320,289), where the bundle gave up and I get bound 1; tally (3,0,1)×46, (5,2,1)×16, (7,4,1)×2, (5,0,3)×2 |
| `forge_e04.py` | `forge_e04.txt` | 3.5–5 min | no | H7 claim 4: the bundle's least_bound is never below an exact lower bound from checked points (24 random E04, 9 congruent-number curves; 0 violations); an order-4 point is refused as a generator |

`z7_points.json` holds the points used by `run_ct_g4k.py` (on E′) and by `forge_e04.py`. They were found by PARI's
`ellrank` (`--find`; `setrand` fixed) and are verified on the curve on every run.

What all of this rests on (quoted, tier E): the definition of Sel² and the order of E(Q_v)/2E(Q_v); the Cassels–Tate
pairing and its Brauer–Manin form (Cassels 1998; Skorobogatov 6.2); Vélu's formula; isogeny invariance of the rank.
The real place in `ct_mine.py` is evaluated with mpmath at 120 digits: a cross-check, not a certificate.
