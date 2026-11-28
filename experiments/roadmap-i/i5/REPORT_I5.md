# REPORT_I5: exact rank bounds for the rows conditional on a rank being 1

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review z17) holds.


Written by an AI model; not read by a person. Experiment, not counted. Nothing here changes the census.

## Bottom line

- **No row closes by a new rank computation.** Each of the 9 curves has an isogeny class of 2 curves, and each curve has exactly one rational 2-torsion point. H7's step needs full rational 2-torsion somewhere in the class, so it does not apply anywhere. rank_indep.py's exact descent (via 2-isogeny, with the second descent) gives rank ≤ 3 on 6 curves and ≤ 5 on 3. Each curve has a point of infinite order, so rank ≥ 1. *(computed, exact)*
- PARI gives root number −1 for all 9. Under the parity conjecture each rank is odd, so 1 or 3 (or 5). That is *conditional (on the parity conjecture)* and does not reach 1. The conductors have 49 to 60 digits, so an analytic rank is out of reach.
- **One G3 open ratio follows without a rank: (12703,2767).** It is the partner of (7735,4968): 12703 = 7735 + 4968, 2767 = 7735 − 4968, and 12703² + 2767² = 2·9193², 7735² + 4968² = 9193². H7 made (7735,4968) unconditional, exactly. The census already moves exclusions between partners (`partner_forms`, `partner_forms_half`, G1Relations.lean), as H7 did for (6557,324). This rests on Drew accepting H7's step. Neither REPORT_H7 nor HANDOFF_ROADMAP_H states it. *(proved, given H7's result and the partner map)*
- **Two more G3 open ratios have their partners in I4's range.** (11165,7548) has partner (18713,3617) (block 18000). (13899,1180) has partner (15079,12719) (block 15000). If I4 closes either partner, the open ratio follows the same way. (12204,9853)'s partner (22057,2351) is outside 20000.

## The table

| Curve | Rows resting on its rank | Class | Exact bound | Root number |
|---|---|---|---|---|
| E1(3540,2581), d = −11740078 | (3540,2581) [CENSUS5000], (6121,959) [CENSUS7000, partner] | 2 curves, one 2-torsion point each | 3 | −1 |
| E1(4440,2231), d = −13342 | (4440,2231) [CENSUS5000], (6671,2209) [CENSUS7000, partner] | same | 5 | −1 |
| E1(6480,4961), d = 709342 | (6480,4961) [CENSUS7000] | same | 3 | −1 |
| E1(8249,3720), d = −108415202 | (8249,3720) [CENSUS9000] | same | 5 | −1 |
| E1(10512,145), d = +220962238 | (10512,145), (10657,10367) [CENSUS11000] | same | 3 | −1 |
| E1(10512,145), d = −220962238 | the same two rows | same | 3 | −1 |
| D1(12007,1799) | (6903,5104) [CENSUS7000], (12007,1799) [G3] | same | 3 | −1 |
| Sig+(12204,9853) | (12204,9853) [G3 open]: its one curve with PARI bounds [1,3] | same | 3 | −1 |
| Sig+(22057,2351) | (12204,9853) through its partner | same | 5 | −1 |

That covers the 9 conditional rows of the 1000–11000 table and G3's (12204,9853). G3's other three:
- (12703,2767) follows from its partner, as above.
- (11165,7548) and (13899,1180): PARI's lower rank bound is 2 or more on every usable curve (REPORT_G3, REPORT_H7), so a rank-1 step cannot help. Their partners are the way in (I4).
- (7735,4968), (6881,6233), (12201,7360) and (13348,3795) are H7's and were not redone.

## How it was checked

- The model is rebuilt from m, n (and d) by the census's own code. E1 comes from census-wide/sd_pari.py (`family(m,n)['+13']`, `weierstrass(d·Q, r)`). D1 and Sig+ come from census2000/dm_pari.py's `all_curves(m, n)`.
- PARI is a finder for the minimal model, the isogeny class and the points. The minimal model is checked by c4, c6 up to u⁴, u⁶. Each isogenous curve is checked by equal point counts mod every good prime < 1000. A point is checked on the curve and checked to be of infinite order (n·G ≠ O for n ≤ 12, Mazur). The bound is rank_indep.py's `rank_bound`, unchanged.
- The root number is PARI's `ellrootno`. It is reported, not used for any bound.
- Rests on: Vélu, Mazur, the 2-isogeny descent and second descent as rank_indep.py implements them (quoted), and the repository's code.
- PARI: every output comes from PARI 2.19.1 (cypari2) on this machine. CI has 2.17.4, which I can't run here. The class is sorted before printing, so PARI's order doesn't show. The descent lines print how many classes the known points fill. That count depends on which points PARI's search finds, so check it under 2.17.4.

## What failed

- H7's step: no curve with full rational 2-torsion in any of the 9 classes.
- The exact descent stops at 3 or 5. Closing a row would need a 4-descent, or a 2-descent over the 2-torsion field with the Cassels–Tate pairing. Those are not in the repository, and I did not start one.

## Files

- i5_rank.py, i5_rank.txt (`--write`; without it the script prints the committed output and `# matches`).
- i5_rank_negative.txt (`--negative`): a forged isogenous curve (a4 + 4·10⁶) is refused, a forged point (y + 1) is refused, and the true 2-isogenous curve and the hint point pass.
- times.txt: 175 s for the 9 curves; under 10 s for `--negative`.

## Proposed LEDGER/CHANGES text

> **Experiment, not counted.** I5: rank_indep.py's exact descent on the 9 curves behind the rows conditional on a rank being 1 (CENSUS5000 to CENSUS11000, G3's (12007,1799) and (12204,9853)). Every class is 2 curves with one rational 2-torsion point each, so H7's full-2-torsion step does not apply. The bound is 3 or 5, with rank ≥ 1, so no row closes. (12703,2767), open in G3, is the partner of (7735,4968) and follows from H7's exact result through `partner_forms`. Computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
python3 experiments/roadmap-i/i5/i5_rank.py            # 175 s here, over 2 minutes: the 9 descents each take 5 to 40 s, so put it in a nightly job, not the per-push CI
python3 experiments/roadmap-i/i5/i5_rank.py --negative # under 10 s
```
