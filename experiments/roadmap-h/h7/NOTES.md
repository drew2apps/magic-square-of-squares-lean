# H7: F4's g4·K route, rank 1 by a full-2-torsion isogenous curve, and E04 of tier B (Roadmap H, run 7)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (z7, `reviews/review_z7.md`, its
own complete 2-descent and Cassels–Tate pairing in `../../referee_z7/`); not read by a person. **Experiment, not
counted; tier D of `LEDGER.md`.** Nothing here changes the census, which stops at `m ≤ 1000`. `REPORT_H7.md` is the
session's report, unchanged; the notes below take precedence over it.

## What it claims (as corrected)

- **Rank exactly 1 for F4's quotient `C: y² = g4·K`** (the curve of `../../census-wide/g4k/dm_g4k.py`) at four ratios,
  (12201,7360), (7735,4968), (6881,6233) and (13348,3795): C has one rational 2-torsion point, its 2-isogenous curve
  `E′` (Vélu) has full rational 2-torsion (for every `m, n`), and the complete 2-descent with the Cassels–Tate pairing of
  `scripts/g1_full/rank_second.py` gives `rank ≤ dim Sel² − 2 − (rank of the pairing)`: 3/0, 3/0, 7/4, 5/2, so 1; a
  point of infinite order gives 1 from below. Review z7's own descent and pairing (its own tangent and local points)
  give the same dimensions and pairing ranks, the kernel of the pairing equals the image of the known points, and PARI
  2.17.4's `ellrank` on `E′` agrees.
- **With F4's exact Demjanenko–Manin chain** (`g4k_exact.py`, run 5, `reviews/review_x3.md`), whose one condition is
  exactly "rank C(ℚ) = 1", every solution of these four ratios is degenerate. This **no longer depends on PARI's rank**
  (the report's "unconditional" is too strong: it rests on the quoted statements below). PARI still finds the root, the
  point, the conic points of the pairing and the factorisations, each checked exactly.
- **(6557,324)** is excluded **through its partner** (6881,6233), as the census already moves exclusions between partners
  (`CENSUS9000.md`; the forms permute with one common factor, `partner_forms` in `G1Relations.lean`). Its own curve gets
  only rank ≤ 3 (5/0). So the L-value of (6557,324) (conditional on PARI's `ellap`, `../../census-wide/lvalue/`) is no
  longer needed for this exclusion; it would be a second, independent one. The same move excludes (12703,2767), open in G3's table, through its partner (7735,4968) (the session's "G3's other four open ratios … need another route" overlooks it; noted by the pre-publication referee, not separately reviewed).
- (7735,4968) was the one ratio open in `CENSUS9000.md`'s range and in the `1000 < m ≤ 11000` table (conditional on rank 1
  by this route); (12201,7360) and (13348,3795) are open in G3's table (`../../roadmap-g/g3/`).
- **E04 of tier B:** for 90 of the 91 ratios that tier B closes on `E04: y² = x(x + 2Sm²)(x + 2Sn²)` (66 curves, a pair
  and its partner sharing one), the same descent and pairing give rank ≤ 1 on 65 of the 66 curves (63 on `E04` itself,
  2 on a 2-isogenous curve with full 2-torsion), and each line's generator is checked to have infinite order. So their
  upper rank bound rests on no PARI, Sage or Magma output. The exception is (864,151), whose curve (run at (1015,713))
  gets 3: its isogeny class has two curves with full rational 2-torsion (the walk visits each three times, as scaled
  copies); `E04` gives 3 and the pairing's local search gives up on the other (z7: dim Sel² = 7 and a pairing of rank 2
  there, so 3 again). Only the rank bound changes: the generators, heights and candidates of those closures stay as
  `scripts/g1_full/dm_evidence.txt` records them. Review z7 reproduced the 91/66/65/90 counts and every curve's result.
- **Not covered:** the 25 `D1` lines and 6 `E0L` lines of tier B (no `D1` curve has full 2-torsion in its isogeny class;
  checked by z7), the 12 σ curves of `G1_FULL_SYSTEM.md` Result 4 (none either), and the seven exclusions of the
  `1000 < m ≤ 9000` census conditional on a rank being 1 (the report's "rank-1-route exclusions" is the wrong label).
- **G3 block `13000 < m ≤ 14000`:** finished; see `../../roadmap-g/g3/NOTES.md`.
- Computed (exact), resting on Vélu's formula, isogeny invariance of the rank, Mazur's torsion theorem, the Cassels–Tate
  pairing in its Brauer–Manin form (Cassels 1998; Skorobogatov, *Torsors and rational points*, 6.2) as computed by
  `rank_second.py` (`reviews/review_t2.md`), and Silverman and Cremona–Prickett–Siksek Theorems 2 and 15 (F4's quotes,
  tier E). Tier D.

## The same step elsewhere in run 7
`scripts/second_descent/rank0_ct.py` (written by the conductor, reviewed by z7: holds) applies it to the curve of
(609,200) and (809,409): dim Sel² = 6 and a pairing of rank 4 on the full-2-torsion curve of its class, so rank 0, where
`SECOND_DESCENT_CHAIN.md` had rank 0 only from PARI's pairing and Magma.

## Rerunning
```
python3 h7_g4k_rank1.py             # 7 s;  --negative 11 s
python3 h7_e04_rank.py              # 57 s; --negative 1 s
```
Each prints `# matches` under PARI 2.17.4 (CI's version; the session used 2.19.1).
