# H7: G3 toward 30,000, and F4's g4·K route

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Machine: PARI 2.19.1 with cypari2, not CI's PARI 2.17.4. Real sympy 1.14.0. PARI is a finder (roots, points, generators, minimal models), and everything found is checked exactly, except where a line says "PARI rank".

## Bottom line

- **F4's g4·K route applies to one of G3's four open ratios, (12201,7360), and there it is now unconditional.** *(computed, exact; rests on the statements quoted below)*
  - PARI's rank of C: y² = g4·K is [1,1] at (12201,7360), and [2,2], [2,2], [2,4] at (11165,7548), (12204,9853), (12703,2767). The route needs rank 1, so it does not apply to those three (h7_g4k.jsonl, *PARI rank*).
  - F4's exact chain (g4k_exact.py, unchanged) at (12201,7360): torsion ℤ/2, ĥ(G) enclosed in [9.870024, 9.884853], G saturated at every prime ≤ 199, an exact lower height bound on the minimal model, then the Demjanenko–Manin bound a ≤ 18 and 37 candidates. Only t = ±1 come back, both degenerate. So **if rank C(ℚ) = 1, every solution t of (12201,7360) is degenerate** (h7_g4k_exact.txt).
  - F4's own exact rank bound (g4k_rank_exact.py, descent via 2-isogeny on C) gives only 3 (h7_g4k_rank_exact.txt).
  - **New: rank C(ℚ) = 1 exactly** (h7_g4k_rank1.py). C has one rational 2-torsion point. Its 2-isogenous curve E′ (Vélu) has full rational 2-torsion. The repository's complete 2-descent with the Cassels–Tate pairing (rank_second.py, through h1_rank0_exact.bound_full, unchanged) gives dim Sel² = 3 on E′, so rank ≤ 1. A point of infinite order gives rank ≥ 1. With F4's chain, every solution of (12201,7360) is degenerate.
- **The same step makes two of F4's three conditional ratios unconditional.** (7735,4968) has dim Sel² = 3, so rank ≤ 1. (6881,6233) has dim Sel² = 7 with a pairing of rank 4, so rank ≤ 1. Both are now rank 1 exactly, so F4's g4k_exact.txt applies without its condition. **(6557,324)'s own curve stays at bound 3** (it follows anyway as the partner of (6881,6233), below): dim Sel² = 5 with a pairing of rank 0 gives only rank ≤ 3, and descents via 2-isogeny on C and on the other two curves of its isogeny class also give 3 (scratch, below). *(computed, exact)*
  - (7735,4968) is the one ratio that CENSUS9000.md leaves open in 1000 < m ≤ 9000. (6557,324) and its partner (6881,6233) are excluded there only on a numerical L-value. (6557,324) = ((6881 + 6233)/2, (6881 − 6233)/2) is the partner of (6881,6233). The census already moves exclusions between partners, because their forms are the same (`partner_forms`, G1Relations.lean). So (6557,324) follows from (6881,6233) without its own rank. That is what H8's L-value was for. If Drew accepts this step, none of these conditions is needed. This is not a census closure; the census stops at m ≤ 1000.
- **Beyond H7: an exact upper rank bound for the curves E04 of the census's Demjanenko–Manin tier.** Tier B of LEDGER.md counts 122 ratios closed by Demjanenko–Manin on curves of rank 1, and their upper rank bound rests on PARI, Sage and Magma. 91 of them use E04: y² = x(x + 2Sm²)(x + 2Sn²), which has full 2-torsion. The repository's complete 2-descent with the Cassels–Tate pairing gives rank ≤ 1 on 65 of the 66 E04 curves (on E04 itself or on a curve 2-isogenous to it), and each line's generator is checked to be a point of infinite order. **So 90 of the 91 E04 ratios have rank exactly 1 without PARI, Sage or Magma** (h7_e04_rank.py). The exception is (864,151), whose curve, run at (1015,713), gets only 3. Only the rank bound changes. The generator, heights and candidates of those closures stay as dm_evidence.txt records them, on PARI, Sage and eclib. *(computed, exact; experiment, not counted)*
- **G3 block 13000 < m ≤ 14000: finished at 12:23 UTC** (g3/times.txt; second descent 7829 s, stages 443 s). g3_table.py now covers 11000 < m ≤ 14000, 916 ratios. The block leaves two new open ratios. **(13348,3795) is closed** the same way as (12201,7360): F4's route gives PARI rank [1,1] and "IF rank 1 then degenerate" (g4k_exact.py), and H7's step gives rank exactly 1 (dim Sel² = 5, pairing of rank 2 on the 2-isogenous curve). *(computed, exact)* **(13899,1180) stays open:** g4·K has PARI rank [2,2] and the exact bound is 2. One more ratio is closed by the rank-1 route, (13103,4271), under that route's conditions (a), (b). Files: run_h7_b13000.sh, h7_b13000_g4k.jsonl, h7_b13000_g4k_exact.txt, h7_b13000_rank1.txt, and g3/g3_table.txt and CENSUS_TABLE_G3.tsv rewritten (both check with `# matches`; `--negative` 4 of 4). Block 14000 was started and stopped at the deadline; nothing from it is claimed.

## Part 1: the g4·K route on G3's open ratios (run_h7_g4k.sh, h7_g4k.jsonl)

dm_g4k.py (F4, unchanged) on each open ratio. The step uses PARI's ranks and generators and floating-point heights, so its "closed" is a finder's answer. It is *computed (numerical)*. Part 2 is the exact check.

| ratio | PARI rank of g4·K | root number | dm_g4k |
|---|---|---|---|
| (11165,7548) | [2,2] | +1 | rank not 1: does not apply |
| (12201,7360) | [1,1] | −1 | closed (a ≤ 39, 79 candidates, t = ±1 degenerate) |
| (12204,9853) | [2,2] | +1 | rank not 1: does not apply |
| (12703,2767) | [2,4] | +1 | rank not 1: does not apply |

## Part 2: F4's exact chain on (12201,7360) (run_h7_g4k_exact.sh, h7_g4k_exact.txt, h7_g4k_rank_exact.txt)

g4k_exact.py and g4k_rank_exact.py, unchanged, in their print-only mode on the one ratio. The output is quoted in the bottom line. g4k_exact's own `--negative` (forged generator 3G, a curve with a small point, a point not on the curve, a classifier forged to accept t = 1) is in the repository (`experiments/census-wide/g4k/g4k_exact_negative.txt`) and was not rerun.

## Part 3: rank 1 by a complete 2-descent on the isogenous curve (h7_g4k_rank1.py, h7_g4k_rank1.txt)

1. C is `dm_g4k.g4k_curve(m, n)`, the curve of F4's scripts. A rational root e of its cubic is found by PARI and checked. The shift X = x − e gives y² = X(X² + aX + b) over ℚ.
2. E′: y² = X(X² − 2aX + a² − 4b) is 2-isogenous to it (Vélu). The script also checks that both curves have the same number of points mod every good prime below 1000, as isogenous curves must. Isogenous curves have the same rank.
3. When 16b is a square, E′ has full rational 2-torsion. The complete 2-descent and the Cassels–Tate pairing of the repository's rank_second.py (`full_selmer2`, `DFull.pairing`, as H1 uses them) give rank E′(ℚ) ≤ dim Sel² − 2 − (rank of the pairing).
4. The point G that PARI finds is checked on the curve, and checked to have infinite order: nG ≠ O for n = 1, …, 12 (Mazur).

| ratio | dim Sel²(E′) | pairing rank | rank bound | rank |
|---|---|---|---|---|
| (12201,7360) | 3 | 0 | 1 | **1** |
| (7735,4968) | 3 | 0 | 1 | **1** |
| (6557,324) | 5 | 0 | 3 | 1 to 3 |
| (6881,6233) | 7 | 4 | 1 | **1** |

The F4 bound on C itself, by descent via 2-isogeny, was 3, 2, 3, 5 for these four curves. The complete 2-descent on E′ is sharper because E′ has full 2-torsion and the pairing applies.

**Controls (h7_g4k_rank1_negative.txt):**
- a forged "isogenous" curve, with b′ changed, is refused by the point counts;
- a forged point is refused;
- positive control: at (11165,7548), where PARI finds two independent points, the exact bound is 2, as it must be at least 2;
- at (12204,9853), the bound is 2 and PARI finds no point at this effort, which is consistent.

## Part 4: E04 of the census's Demjanenko–Manin tier (h7_e04_rank.py, h7_e04_rank.txt)

- For each E04 line of `scripts/g1_full/dm_evidence.txt` (91 ratios, 66 curves, since a pair and its partner share one), the model is dm.py's [0, 2S², 0, 4S²m²n², 0].
- The generator abscissa in that line is checked to be the x of a point, and the point is checked to have infinite order (nG ≠ O for n ≤ 12).
- Then the complete 2-descent and the pairing run on E04. If that does not give 1, they run on the full-2-torsion curves reached by 2-isogenies (at most 6), and the least bound is kept.
- **Result:** rank 1 on 65 curves. On E04 itself: 63 curves (46 with dim Sel² = 3; 16 with dim Sel² = 5 and a pairing of rank 2; 1 with dim Sel² = 7 and a pairing of rank 4). On a 2-isogenous curve: the curves of (320,289), where the pairing's local search gives up on E04, and (1015,961), where E04 gives 3. On the curve of (1015,713) all 6 curves give 3, or the pairing's local search gives up. That curve carries (864,151).
- **Controls (h7_e04_rank_negative.txt):**
  - a forged model is refused;
  - a forged generator (x + 1) is refused;
  - the 2-torsion point as generator is refused;
  - positive control: the bound on the full-2-torsion curve of (11165,7548) (PARI: two independent points) is 2.
- **D1 and E0L (31 ratios), scratch, not in the bundle:**
  - No D1 curve has full 2-torsion anywhere in its isogeny class.
  - The descent via 2-isogeny with second descent (rank_indep.py) gives rank 1 only for the D1 curves of (388,287) and (920,371), and 3 to 5 for the other 14. That includes the two that rest on Magma alone, (631,71) and (1389,59), with bounds 3 and 4.
  - E0L was not tried.

## Checks

- The DM step and its negatives are F4's, reviewed by a separate AI agent (`reviews/review_x3.md`), which re-derived the bounds with its own code. Nothing in it was changed.
- The complete 2-descent and the pairing are the repository's, as used for H1's rank bounds, with H1's negatives (h1/h1_rank0_exact_negative.txt).
- Every output matches its committed file without `--write`.

## What it rests on

- **Exact:** the shift, the point counts, the Selmer group and the pairing, the point's order, and all of F4's chain (torsion, the height enclosure, saturation at primes ≤ 199, the minimal model, the height search, every DM candidate).
- **Quoted:** Vélu's formula for the 2-isogeny; Mazur's torsion theorem; the soundness of the complete 2-descent and of the Cassels–Tate pairing as computed by `rank_second.py`; Silverman's height bounds and Cremona–Prickett–Siksek Theorems 2 and 15 (F4's quotes, tier E).
- **Not used:** PARI's rank. PARI only found the root, the point, G's saturated form and the minimal model, each checked.

## What failed

- **(6557,324) on its own curve:** the full 2-descent on E′ gives 3 (not needed, since it follows from its partner). Descents via 2-isogeny on C and on the other two curves of the isogeny class (rank_indep.py; it does not take E′, which has full 2-torsion) all give 3 too, though PARI gives [1,1] on one of them (scratch, not in the bundle). A second 2-descent on E′, or a 4-descent, would give a second route.
- **The 12 σ curves of G1_FULL_SYSTEM.md Result 4** (eight of them rank 1 only by PARI and Magma): none has a 2-isogenous curve with full 2-torsion (16b is not a square), so this step does not apply to them (scratch).
- **The D1 curves resting on Magma alone, (631,71) and (1389,59):** no full 2-torsion in their isogeny classes, so this step does not apply; the descent via 2-isogeny gives 3 and 4 (scratch).
- **The seven exclusions of CENSUS5000 to 9000 conditional on rank 1:** PARI already gives [1,3] on every curve of their isogeny classes, with full 2-descents. So this step, which is the same descent done exactly, cannot do better (not run).
- **H1's open pair:** dm_g4k reports "closed" on both sides (PARI rank 1). But F4's g4 = 4F₄ uses F₄, which is not one of the four forms of H1's curve, so it says nothing about H1's extra tuples. The quotient y² = w·K that H1's curve does have also has PARI rank 1 on both sides. There, though, t → 1/t acts as translation by the 2-torsion point (0,0) (x x′ = b), so Demjanenko–Manin gives nothing (scratch).
- **G3 block 13000:** see the bottom line. Two earlier runs were lost to a container restart and to a stop by hand that freed cores for H2 (g3/times.txt).

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** Rank 1 for F4's quotient y² = g4·K at (12201,7360), (7735,4968) and (6881,6233), by the repository's complete 2-descent and Cassels–Tate pairing on the 2-isogenous curve with full 2-torsion (experiments/roadmap-h/h7). With F4's exact Demjanenko–Manin chain (g4k_exact.py), every solution of these three ratios is degenerate, no longer conditional on PARI's rank. (12201,7360) was one of G3's four open ratios in 11000 < m ≤ 13000. The same holds at (13348,3795), one of the two open ratios of G3's block 13000 < m ≤ 14000; the other, (13899,1180), stays open. (6557,324) follows as the partner of (6881,6233), so the numerical L-value behind that pair is no longer needed. Separately, the same descent gives rank exactly 1 for the curve E04 of 90 of the 91 ratios of tier B closed on E04, so their upper rank bound no longer needs PARI, Sage or Magma (generators and heights unchanged); (864,151) is not covered. Computed (exact), resting on Vélu, Mazur, Silverman and Cremona–Prickett–Siksek (quoted). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-h/h7        # F = grep -v '^PARI stack size'; MSS_REPO, H1_DIR as in the script
python3 h7_g4k_rank1.py 2>&1 | $F | tail -1 | grep -q '^# matches'              # about 9 s
python3 h7_g4k_rank1.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'   # about 14 s
python3 h7_e04_rank.py 2>&1 | $F | tail -1 | grep -q '^# matches'               # about 46 s
python3 h7_e04_rank.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'    # about 1 s
# print-only runs of F4's scripts (outputs committed here): g4k_exact.py 12201,7360 (2 s), g4k_rank_exact.py 12201,7360 (20 s), dm_g4k.py on the four ratios (2 to 9 s each)
```
