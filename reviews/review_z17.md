# Referee z17: roadmap I, parts 1 and 2

Written by an AI model (reviewer z17), not read by a person.

> Committed as written apart from this note. In this repository the reviewer's folder is `experiments/referee_z17/; the bundle reviewed is `experiments/roadmap-i/``. Paths under `/home/claude/` are the reviewer's scratch folders.

Worktree: /home/claude/conductor/z17 (branch run11_z17 from run10 = 315ad1e). Bundle copied to experiments/roadmap-i; my scripts and outputs in experiments/referee_z17/. Container: cypari2 2.2.9 / PARI 2.17.4 (to be confirmed below).

Findings are appended as I go.

## 1. CI lines under PARI 2.17.4 (cypari2 2.2.9, sympy 1.14.0, Python 3.13.18)

Confirmed: `version()` = [2, 17, 4]. Machine shared with other reviewers (2 cores, load about 11), so wall times are inflated 2-4x over an idle machine. Outputs in experiments/referee_z17/ci/ (stdout .out, stderr .err, times.txt).

| Line | Exit | Byte for byte | Wall here | Claimed |
|---|---|---|---|---|
| `python3 -I i2/i2_residual.py --negative` | 0 | `# matches` | 7.9 s | 3 s |
| `python3 -I i8/i8_pairct.py --negative` | 0 | `# matches` | 47.6 s | 13 s |
| `python3 -I i8/i8_pairct.py --compare50` | 0 | `# matches` (50 of 50) | 32.6 s | 5 s |
| `python3 -I i4/i4_exchange.py --negative` | 0 | `# matches` | 2.8 s | 3 s |
| `python3 -I i4/i4_exchange.py` (the report itself) | 0 | `# matches` | 0.4 s | - |
| `python3 -I i5/i5_rank.py --negative` | 0 | `# matches` | 5.9 s | <10 s |
| `python3 -I i1/i1_laws.py` | 0 | `# matches` | 5.0 s | 1 s |
| `python3 -I i1/i1_laws.py --negative` | 0 | `# matches` | 3.5 s | 1 s |
| `python3 -I i1/i1_negative.py` | **1** | ModuleNotFoundError: i1_hunt | 0.4 s | 2 s |
| `python3 i1/i1_negative.py` (no -I) | 0 | `# matches` | 17.7 s | 2 s |

- [note] i1_negative.py does `import i1_hunt` from its own folder, so it fails under `python3 -I` (Python 3.11+ -I implies -P). The CI's `hm` helper runs without -I from the script's folder, so it would pass there; either add `sys.path.insert(0, HERE)` or keep it out of any `-I` loop. i1/times.txt itself records "i1_negative (rewrite) exit=1 3 s", unexplained in REPORT_I1.
- Later runs: `python3 -I i5/i5_rank.py` (9 curves) exit 0, `# matches`, output byte-identical to i5_rank.txt, 500.0 s here (175 s claimed; nightly). `python3 -I i3/i3_scan.py --negative` exit 0 `# matches` 48.9 s (27 s claimed); `python3 -I i3/i3_scan.py 13` exit 0 `# matches` 294.0 s (88 s claimed). `i8_checkpoint.py 360 3`: exit 1 (see 2). neither/epm.py, ri.py and dm_pari.py on in_*.json: identical outputs (section 3).
- (Sections below were appended in the order I did them: 2, 4, 5a, 3, 5b, 5c, 6, 7, 8.)

## 2. Spot checks with the repository's own code

### I2 (experiments/referee_z17/z17_i2check.py -> z17_i2check.txt, 29 s)
- Recounted from `i2_residual.jsonl` with my own code: 3,876 records, same ratios and order as `i2_needgood_list.txt`; fields exactly m, n, residual, kill, killing_prime in that order; **32,474 tuples; |S| from 2 to 834, every |S| even; no null kill; largest killing prime 229**; killing_prime = max(kill) everywhere; residual sorted, no repeats. Killing-prime counts 29: 1,592, 37: 713, 41: 491, 53: 260, 61: 250 as stated. All reproduce.
- Every killing prime is good (divides none of 2mn(m^2-n^2), m^2+n^2, the octic; checked by exact division, no factoring) and prime: 0 exceptions.
- **Seeded sample (seed 17), 50 ratios: 30 from the H2 lines and 20 from the G2 lines (which I8's --compare50 never draws from: it samples only the first 2,900), recomputed with conductor-t10's `residual()` (exact_local + X.join; I2 used h4_split.chunked): 50 of 50 same tuples and killing primes (326 tuples).** Both runs share exact_local's local test, so this checks the join and the prime loop, not the local test itself.
- [fix] **"G2's 937 are lines 2,940 to 3,876"** (hand-off part 2 item 1, REPORT_I8 checkpoint line) is wrong by 9: by the list's own labels, G2 lines are 2,949-3,876 (928 lines: 827 `g2-15000` + 101 `g2-997`); the other 9 of G2's 836 (checked against roadmap-g/g2/g2_15000_need_good.txt: 836 = 827 + 9) carry the label `h2-hunt` and sit among lines 2,808-2,948. The |S| distribution in REPORT_I2 (1,547 / 893 / 393 / 266) is for the shape run's 3,838; over all 3,876 it is 1,559 / 904 / 398 / 268 (say which).

### I8 checkpoint (z17_i8count.py -> z17_i8count.txt)
- Recounted from `i8_first360.jsonl`: 360 records = lines 1-360 (all `h2-exhaustive`); tuples and killing primes equal those of i2_residual lines 1-360 (3,166 tuples); verdicts recomputed from the per-curve labels: pairing 3,160, neither 6, 0 inconsistent; kills per curve 08:1836 0L/0U/8L/8U:3058 LU:1794; mixed-only 462; 11 ratios excluded by a rank-0 pair curve (7 LU, 4 08, torsion 8 of bound 8 each). **All headline numbers reproduce from the file.** No kill is "outside Sel2": all 15,862 kill labels are pairing -1.
- [fix] The proposed CI line `i8_checkpoint.py 360 3` **cannot run on the bundle**: it reads the workers' files i8_<i>of3.jsonl, which are not in it; here it exits 1 with "not done yet: 360 of the first 360 ratios". Either commit a recount from i8_first360.jsonl (like z17_i8count.py) or drop the line.

### 2b. I8 from scratch (z17_i8rerun.py, run_i8rerun.sh -> z17_i8rerun.jsonl; z17_i8compare.py -> z17_i8compare.txt)
- conductor-t10's `analyse_ratio(m, n)` with **its own residual()** (not I2's file), one fresh process per ratio, PARI 2.17.4, on the 3 "neither" ratios and 15 others drawn with seed 17 from lines 1-360 (1,253 s in all).
- **18 of 18 ratios (108 tuples): same tuples, same killing primes, same killing curves (verdicts), and even the same PARI-dependent 'image'/'kernel' labels as i8_first360.jsonl.**
- **The 6 "neither" tuples are really in the kernel or the image on every pair curve**, and no pair curve of the three ratios has dim Sel^2 = 2, so I8's rank-0 test cannot apply: (18313,120) dims 6,6,6,6,6,5 (08,0L,0U,8L,8U,LU); (18433,18193) 5,6,6,6,6,6; (20056,19825) 7,7,7,7,7,3.
- [note] On (20056,19825), E_LU has dim Sel^2 = 3, so its rank is <= 1 exactly; if its rank is 1 (PARI says [1,1]; I did not check a point of infinite order on E_LU exactly), Sel^2 is the image of E_LU(Q), the two tuples' LU class is the class of a rational point, and no pairing on LU can ever remove them. Their 'kernel' label there means only that t10's point images span just the 2-torsion classes (ri = 2) although PARI reports [1,1]; the same holds on LU of (18313,120). Labels depend on PARI's points, as REPORT_I8 says; verdicts do not.

## 4. I5's partner claim: (12703,2767) through (7735,4968)

- **Partner map holds.** 7735 + 4968 = 12703 (odd), 7735 - 4968 = 2767, gcd(12703,2767) = 1; 12703^2 + 2767^2 = 169,022,498 = 2 * 9193^2 and 7735^2 + 4968^2 = 9193^2. `partner_forms` (MagicSquareOfSquares/G1Relations.lean:103) states that the five forms of (m+n, m-n) at t' = (t+1)/(t-1) are the five forms of (m,n) at t in another order times squares; `partner_coprime` covers m+n odd, as here. This is the step CENSUS9000.md and H7 already use ((6557,324) from (6881,6233)).
- **H7's result for (7735,4968) is exact in rank and reproduces here** (`h7_g4k_rank1.py`, `# matches` under 2.17.4, 29 s): on the 2-isogenous curve with full 2-torsion, dim Sel^2 = 3 and the Cassels-Tate pairing has rank 0, so rank <= 1; a point of infinite order gives 1. The rest of the exclusion is F4's exact chain `census-wide/g4k/g4k_exact.txt` line 1-8 (exact rational height enclosure, saturation at primes <= 199 plus a height lower bound, so G generates given rank 1; DM a <= 18: only t = +-1, degenerate), which is in CI.
- So the claim holds, **with the status word fixed**: REPORT_I5 says "(proved, given H7's result and the partner map)" and "H7 made (7735,4968) unconditional, exactly". h7/NOTES.md already corrected that word: "the report's 'unconditional' is too strong: it rests on the quoted statements". [fix] Write "computed (exact), resting on H7's and F4's quoted statements (Velu, Mazur, Cassels-Tate in Brauer-Manin form, Silverman, Cremona-Prickett-Siksek) and the census's partner step"; not "proved", not "unconditional".
- [note] REPORT_I5's "Neither REPORT_H7 nor HANDOFF_ROADMAP_H states it" is true of those two files but h7/NOTES.md line 24 and roadmap-g/g3/NOTES.md line 58 already state it (hand-off part 1 item 6 says so too). Nothing new is needed in main.
- [note] Hand-off part 1 item 7 (a sentence repeated three times in roadmap-h/NOTES.md's H3 row and in h7/NOTES.md): true at the session's base 97c8564, **already fixed in main by run 9 (d073fd7, "repeated text of run 7 collapsed")**; run10 has each once. Drop the item.

## 5a. I4: open rows and the exchange test

- Chunk bookkeeping reproduces: chunk_14000 has 301 ratios and sd_14000_{0,1} 301 unique records, no duplicates left; chunk_16000 305/305; chunk_18000 306 with 109 done (REPORT_I4 says 111: [note]). stage_14000.json / stage_16000.json say 301 -> 45 -> 11 -> 2 and 305 -> 47 -> 22 -> 2, open (14017,2856), (14640,10679) and (16065,2032), (16873,11161), as reported. I did not rerun the second descent or DM (hours); the stage counts are taken from the files.
- (14017,2856) <-> (16873,11161) is a partner pair ((14017+2856, 14017-2856) = (16873,11161)), so "that pair stays open" is right as a pipeline statement. (14640,10679) is still open too: the task brief's list omits it, the hand-off part 1 has it.
- **What the exchange test tests.** The five forms are symmetric under (m,n,u,v) <-> (u,v,m,n); I checked the algebra: with t = v/u and t' = n/m each of the five forms and the four linear cells at (u:v, t') equals the one at (m:n, t) times u^2/m^2, so the nine-cell verdict must agree. For a nontrivial SOLUTION it then requires u:v to be a valid class-keeping ratio. On the data it is **vacuous on its substantive branch**: all 1,927 points are 'no' (828) or trivial (1,099, t in {0, +-1, inf}); no nontrivial point at all, so only the exact re-classification and the trivial branch ran. REPORT_I4 says this ("can only bite on a SOLUTION, and none has come up"); the hand-off (part 2 item 4: "0 failures over 1,927 reported points") does not. [fix] add "all 1,927 are 'no' or trivial; the ratio check never ran on data".
- [fix] REPORT_I4: "The same test on G3's own outputs (roadmap-g/g3) gives 767 'no', 1,098 trivial and 0 failures". On main's experiments/roadmap-g/g3 (run10) I get **no 1,064, trivial 1,501, 0 failures** (12 files, blocks 11000-13000); no subset of those files gives 767/1,098. The verdict (0 failures) holds; the counts do not reproduce.
- Files: sd_14000_1.err (355 KB, 5,639 lines of PARI "not enough memory, new thread stack" warnings), the empty sd_*.err, done_14000 (a marker), and the partial chunk_18000/sd_18000_* should not be committed as results; i4_exchange.txt depends on the partial sd_18000 files, so if they are not committed, regenerate it without them.

## 3. The three "neither" closures (z17_neither.py parts a-e -> z17_neither_{a..e}.txt; neither/ reruns in referee_z17/neither/)

- **Curves and models.** The DM runs are on census2000/dm_pari.py's `Sig+` of (18313,120), `Sig-` of (20056,19825) and `E` of (39881,231). I8 calls the first two "E+" and "E-" (sigma_exact.model(m,n,+-2), neither/epm.py). Checked: the a-invariants are **identical** (so isomorphic), and Sig+/Sig- are not isomorphic to the other sign's model. The naming is consistent; [note] REPORT_I8 should say once that "E+" = dm_pari's Sig+, "E-" = Sig-, since the dm_*.jsonl files use the Sig names.
- **Under PARI 2.17.4 everything in neither/ reproduces**: epm_*.txt and ri_*.txt byte for byte (ri with s = 2 for (18313,120), s = -2 for (20056,19825); the committed ri files do not record s: [note]); dm_pari.py rerun on in_*.json gives the same 8 + 8 records as dm_*.jsonl, field by field without `sec`, including the floating-point hG, B, amax.
- **What each closure rests on (all three are dm_pari closures):** (i) rank exactly 1 from PARI's `ellrank` bounds over the isogeny class (the upper bound is PARI's own descent, taken on trust; it may rest on class-group computations that PARI does not certify); (ii) the generator is PARI's point, passed through PARI's `ellsaturation(.., 1000)`: no index bound, so "G generates E(Q)/tors" is unproved (a missed index would have only prime factors > 1000); (iii) BUP, U and hG are floats (1.0001 and 1.02 margins), not certified enclosures. The relation and candidate test are the repository's. I checked exactly: PARI's point on each curve is on the curve and kG != O for k <= 12 (part e), so rank >= 1.
- **Can the repository's exact descents bound the rank?** (part d)
  - (18313,120) Sig+ (and so the partner (18433,18193)): isogeny class of 2 curves, one rational 2-torsion point each, no full-2-torsion member (so no H7 / rank0_ct step). rank_indep.rank_bound gives **rank <= 3** on both members (Sel(phi) 64 / Sel(phi^) 2, three cosets not excluded). Stays on PARI's rank.
  - (20056,19825) Sig-: class of 2, one 2-torsion point each; rank_indep gives **rank <= 3** on both. Stays on PARI's rank.
  - **(39881,231) E = y^2 = x(x + 2Sm^2)(x + 2Sn^2) (tier B's E04 shape): class of 6 curves, 2 with full rational 2-torsion. H7's step (h1_rank0_exact.bound_full: complete 2-descent + Cassels-Tate pairing of rank_second.py) gives dim Sel^2 = 7 with a pairing of rank 4 on one, dim 3 with rank 0 on the other: rank <= 1. With the point of infinite order: rank E(39881,231) = 1, computed (exact).** (rank_indep on the one-torsion members gives 5, 5, 1, 1.) So **(20056,19825) closes through its partner (39881,231) with its rank exact**, leaving only (ii) and (iii). REPORT_I8 mentions this partner closure only in passing and calls it conditional on PARI's rank like the others: [fix].
- **The right labels.**
  - (18313,120) and its partner (18433,18193): "closed by Demjanenko-Manin on Sig+ (= E+) of (18313,120); conditional: rank 1 from PARI's bounds only (exact descent: <= 3), generator saturated by PARI at primes <= 1000 only, heights in floating point". In the census's tables this is the footing of every dm_pari row (PARI [1,1]); in I5's sense ("rests on a rank being 1 that no exact descent gives") it is **conditional on a rank being 1**. REPORT_I8's word "conditional" is right; "closed" alone is not.
  - (20056,19825) and (39881,231): "closed by DM on E of (39881,231), rank exactly 1 (H7's step, exact); conditional only on PARI's saturation at primes <= 1000 and floating-point heights". Not "conditional on a rank being 1" any more.
  - The pairing result itself stands: 6 tuples in 3 ratios are in the kernel or the image on all six pair curves (reproduced from scratch, section 2b), so "the pairing does not reproduce every good-prime kill" is right as stated for these 360 ratios.
- To remove (ii) and (iii) the pattern exists in the repository: census-wide/g4k/g4k_exact.py does exact rational height enclosures, exact saturation to 199 and a height lower bound that makes G generate given rank 1. For E(39881,231) that would make (20056,19825) exact apart from quoted theorems. Not attempted here.

## 5b. I1 designs

- Recounted from the seven hunt files: taken 1,569; no local point 589; no class 908; class kept 50 (29 distinct: A 3, B 16, C 10, no overlaps); P = 0 for 49 lines, P = 29 for (651432,64265) only; i2/i1_classkept_list.txt = the same 29. **"No survivor in 1,569 designed ratios" holds as computed by the repository's pipeline** (SURVIVORS: none in every file).
- [fix] 21 of the 1,569 "taken" lines are rounds where nothing passed the congruences (`- K none`: 6 + 7 + 7 in design B, 1 in design A seed 2), not ratios. REPORT_I1 says "6 rounds in B" (it is 20). So 1,548 ratios were tested; B's table row 519 = 0 + 463 + 36 + 20 empty rounds.
- [fix] REPORT_I1's killing-set section is stale: it says 25 of 29 done, 6 over the memory cap and 4 "still running"; the final i2/i2_i1_classkept_list_0of1.txt has **19 done and 10 over the cap** (the 4 were marked MEMORY). The .part beside it is superseded. i1/times.txt records the shape run on that list as exit=137 (killed).
- [note] i1_negative.py's positive control uses H1's (204724,97893), which I3 confirms is one of the 7 true phantom pairs, so "not emptied to 5000" is the right expectation.
- Design C "works backwards from H3's tables" as described (draw_C: CRT on residues at 29, 37, 53, 61 whose H3 mask allows the target pattern); I did not rerun the hunts (270-4,772 s each).

## 5c. I3 pairs and rescan vs H1 and review z8

- i3_pairs.txt: 14 of 24 ratios still phantoms (7 pairs), 10 not (5 pairs), each pair losing its extra tuples at one octic prime: 5462461, 218077, 10057297, 21062406847889269, 1650997632721. **These are exactly the first octic primes of reviews/review_z8.md lines 83-87** (z8 also lists further octic primes that kill). Agrees with main's h1/NOTES.md ("H1's 12 pairs are 7").
- Per-range counts recounted from i3_scan_k13..k29.txt: 0 0 0 0 1 0 3 1 0 0 0 0 0 0 2 0 0, 17 x 1177 = 20,009 ratios: as stated.
- `i3_scan.py --negative`: `# matches` under 2.17.4, 48.9 s here (27 s claimed). `i3_scan.py 13`: see section 1.
- [note] "This agrees with review z8, and was found independently": the session's base already carried z8's correction (hand-off part 1: "H1's 12 pairs are 7" applied), so "independently" means "recomputed", not "found without knowing the answer".
- The 50,014 new ratios: not finished, no result; nothing to check. i3_scan_k*.txt are H1's ratios, not the new ones.

## 6. Claims editor (one line per headline number; "repro" = I recomputed it, "file" = recounted from the committed file, "not checked" = taken on trust)

Hand-off part 1
- I1 "no survivor in 1,569 designed ratios (three designs)": file; 21 of the 1,569 are empty rounds, so 1,548 ratios. Status "computed, exact" fine.
- I1 "29 distinct ratios keep a class; bad primes empty all but (651432,64265), emptied at good 29": file. Holds.
- I1 "Design C: 10 class-keeping ratios, every one emptied with P = 0": file. Holds.
- I2 "2,929 of the 3,876" (item 2) vs "2,932 of 3,876" (Per task): inconsistent; i2_residual_partial.jsonl has 2,930 records. Superseded by part 2 (all 3,876).
- I2 "69% of the sets are cosets": superseded; 2,693 of 3,838 = 70% (file).
- I2 "killing prime at most 229": repro (all 32,474 tuples).
- I3 "7 of H1's 12 pairs are true phantoms; 5 are not (each killed at an octic prime)": file; agrees with z8 line by line.
- I3 "rescan of all 20,009 differs from H1 only by those false pairs": file (per-range counts). Holds.
- I4 "two open rows so far, (14017,2856) and (14640,10679)": file (stage_14000.json). Pipeline status, "consistent with PARI" at rank 0 as G3.
- I4 "0 failures over every point reported so far": repro, but vacuous (no nontrivial point).
- I5 "no conditional row closes by a rank step; 9 classes of 2 curves with one 2-torsion point": repro under 2.17.4 (i5_rank.py # matches). Bounds 3 (6 curves), 5 (3 curves): file. Holds.
- I5/(12703,2767) via (7735,4968): holds (section 4); "proved"/"unconditional" too strong.
- Item 7 (repeated sentences in main's notes): already fixed in run 9. Drop.
- Base: "67,583 for good p = 1 mod 4, below 1154 for p = 7 mod 8; H1's 12 pairs are 7": agrees with main's roadmap-h/NOTES.md.

Hand-off part 2
- "6 residual tuples in 3 ratios are neither": repro from scratch (t10 rerun) and from the file.
- Per-curve table for the 3 ratios (kernel/image pattern): repro, labels identical under 2.17.4.
- "E+ rank 1 by PARI; DM 119 candidates; no four-corner point" (18313,120): repro under 2.17.4 (dm_pari rerun identical). E+ = dm_pari's Sig+.
- "(18433,18193): no curve of rank 1 (E+ [1,3], E- [2,2])": repro (dm records Sig+-partner [1,3], Sig--partner [2,2]).
- "(20056,19825) E- rank 1 by PARI; DM 113 candidates": repro. "Its partner (39881,231) closes the same way on its curve E": repro, and **its rank is exactly 1** (new here, section 3), so "conditional on PARI's rank" is too weak a description for this one.
- "rank_indep.rank_bound gives only rank <= 3 on both curves": repro (and on both members of each class).
- "generator saturated only at primes <= 1000; heights in floating point": correct description of dm_pari.
- "same ranks under PARI 2.17.2": and under 2.17.4 (repro).
- I8 "3,160 of 3,166 pairing-removed, 6 neither, 11 ratios rank 0 (7 LU, 4 08), 462 mixed-only": file; and 18 ratios repro from scratch (section 2b).
- "G2's 937 are lines 2,940 to 3,876": **wrong**: G2-labelled lines are 2,949-3,876 (928); 9 more are labelled h2-hunt at earlier lines.
- I2 "32,474 tuples, every |S| even, every tuple killed by a good prime <= 229": repro.
- I2 shape run "3,838 (38 over its memory cap); 70% cosets; killing sets of 3 to 5 primes for 2,557; separating quotient at P never 1 coordinate (2,477 / 474 / 887)": file, all reproduce; |S| and P of the 3,838 agree with the residual file (0 of 3,838 disagree).
- PARI: "2.17.4 was not available here, so nothing was run under 2.17.4": now done here; every diffed output matches under 2.17.4.
- I4 "chunk 16000 done, 2 of 305 open: (16065,2032) and (16873,11161)": file. "exchange test 0 failures over 1,927 points": repro, vacuous.
- I4 "chunk 18000 paused part way" (REPORT: 111 ratios): the bundle's sd_18000 files have 109.
- I3 "50,014 new ratios started, neither range finished": nothing to check.
- Times table: 3 / 13 / 5 / 3 s claimed; here (loaded, 2.17.4) 7.9 / 47.6 / 32.6 / 2.8 s.

Status words
- REPORT_I5 "(proved, given H7's result and the partner map)" -> "computed (exact), resting on H7's and F4's quoted statements and the partner step".
- REPORT_I5 "H7 made (7735,4968) unconditional, exactly" -> "no longer conditional on PARI's rank" (h7/NOTES.md's own correction).
- REPORT_I8/hand-off "So each is closed, conditionally": right for (18313,120)/(18433,18193); for (20056,19825)/(39881,231) say "rank exact (H7's step), generator and heights by PARI".
- REPORT_I3 "found independently": "recomputed".
- REPORT_I8 "(computed, exact, given t10's Sel^2 and pairing code as reviewed in z11)": fine for the 3,160; t10's labels 'image'/'kernel' are PARI-dependent, the verdicts are not (my rerun: identical anyway).

## 7. CI lines: what passes under 2.17.4, and how they must be wired

- Every committed-output line passes under PARI 2.17.4 (section 1 table), plus `i3_scan.py --negative` (48.9 s), `i3_scan.py 13` (294 s here, 88 s claimed) and `i5_rank.py` (500 s here, 175 s claimed): all `# matches`.
- [fix] **The proposed lines are bare `python3 ... --negative`. These scripts print `# matches` or `# DIFFERS` and exit 0 either way**, so as written a mismatch would pass CI. Run them through the existing `hm` helper (ci.yml line 646: greps that the last line is `# matches` and no line says DIFFERS), which also does the `cd` that i1_negative.py needs.
- [fix] Drop `i8_checkpoint.py 360 3` (needs uncommitted worker files; exits 1).
- Proposed block (about 3-4 minutes on an idle runner by the session's times; about 4 minutes here under load without i5_rank/i3 k13):
```
          R=experiments/roadmap-i
          hm $R/i1 i1_laws.py; hm $R/i1 i1_laws.py --negative; hm $R/i1 i1_negative.py
          hm $R/i2 i2_residual.py --negative
          hm $R/i3 i3_scan.py --negative
          hm $R/i4 i4_exchange.py --negative
          hm $R/i5 i5_rank.py --negative
          hm $R/i8 i8_pairct.py --negative; hm $R/i8 i8_pairct.py --compare50
          # review z17 (experiments/referee_z17), under a minute
          for s in z17_i8count z17_i8compare; do (cd experiments/referee_z17 && python3 -I $s.py | tail -1 | grep -qx '# matches'); done
          (cd experiments/referee_z17 && python3 -I z17_i2check.py 2>/dev/null | grep -v 'PARI stack' | tail -1 | grep -qx '# matches')
          for p in a b c e; do (cd experiments/referee_z17 && python3 -I z17_neither.py $p 2>/dev/null | grep -v 'PARI stack' | tail -1 | grep -qx '# matches'); done
          # nightly / [full ci] only: hm $R/i5 i5_rank.py (175-500 s); hm $R/i3 i3_scan.py 13 (88-294 s); z17_neither.py d (218 s)
```
- `i4_exchange.py` without arguments matches only with the partial sd_18000_*.jsonl committed; leave it out until chunk 18000 is final.
- Sizes: the bundle is 5.27 MB in 127 files. **Do not commit** (2.29 MB): i2_residual_0of1.jsonl.part (822 KB) and i2_residual_partial.jsonl (820 KB) (both superseded; I checked every record equals i2_residual.jsonl, but i8_pairct.residual_records() reads the .part *after* the final file, so a stale .part would silently override it: [note]), the three shape .part files (290 KB), i4/sd_14000_1.err (355 KB of PARI thread-stack warnings) and the three empty .err files, i4/done_14000. Probably not as results: the partial i4/chunk_18000.json and sd_18000_*.jsonl (22 KB). Keep: i2_residual.jsonl (1.06 MB, the deliverable), i8_first360.jsonl (707 KB; the checkpoint's evidence), the shape outputs (376 KB). Without the listed files the bundle is about 2.95 MB. Also exclude any __pycache__ (created in i1/ and i2/ by running the scripts; .gitignore already covers it).
- Nothing outside experiments/roadmap-i and experiments/referee_z17 changed (`git status` clean apart from those), and the copied bundle is byte-identical to the untrusted original after all runs.

## 8. Proposed LEDGER/CHANGES text

> **Experiment, not counted (tier D).** Roadmap I (experiments/roadmap-i), parts 1 and 2, reviewed by z17 (experiments/referee_z17). I1: three designed hunts, 1,548 ratios (1,569 rounds), no survivor; 29 distinct ratios keep a class and the bad primes alone empty the four-class set of all but (651432,64265), emptied at the good prime 29. I2: for the 3,876 ratios of H2 and G2 that need a good prime, `i2_residual.jsonl` lists the 32,474 tuples the bad primes leave and the good prime that kills each; every set has even size and every killing prime is at most 229 (50 ratios recomputed by review z17 with conductor-t10's own code, G2's included). I3: 7 of H1's 12 pairs are true phantoms, as review z8 found; the rescan of H1's 20,009 ratios agrees. I5: the exact descent gives rank 3 or 5 on the 9 curves behind the rows conditional on a rank being 1, so none closes; (12703,2767) follows from (7735,4968) (H7's exact rank and F4's exact chain) through the census's partner step. I8 (checkpoint, lines 1-360 of 3,876): the Cassels-Tate pairing on the six pair curves removes 3,160 of 3,166 residual tuples; 6 tuples in (18313,120), its partner (18433,18193) and (20056,19825) are in the pairing's kernel on every pair curve (rerun from scratch by z17). Those ratios are closed by Demjanenko-Manin: (18313,120) and its partner on E+ (dm_pari's Sig+), conditional on a rank being 1 (PARI's bounds; the exact descent gives 3), on PARI's saturation at primes <= 1000 and on floating-point heights; (20056,19825) through its partner (39881,231) on E, whose rank is exactly 1 (H7's step, z17), conditional only on the saturation and the heights. I4 (interim): G3's pipeline leaves (14017,2856) with its partner (16873,11161), (14640,10679) and (16065,2032) open in 14000-15000 and 16000-17000. Every committed output reproduces under PARI 2.17.4. Computed. Written by AI models; not read by a person.

## Bottom line per item (short)
- I1: holds (counts from files; 1,548 ratios, not 1,569; killing-set section stale).
- I2: holds; reproduced on 50 ratios (20 from G2's lines) with t10's code; "lines 2,940-3,876" wrong.
- I3: holds; agrees with z8 prime by prime; "independently" overstated.
- I4: counts from files hold; exchange test is vacuous on the data; G3 comparison counts (767/1,098) do not reproduce (1,064/1,501, still 0 failures).
- I5: reproduces byte for byte under 2.17.4; partner claim holds with "computed (exact)", not "proved".
- I8: all checkpoint numbers reproduce; 18 ratios rerun from scratch identical; the 6 "neither" tuples are genuine; closure of (20056,19825) upgrades to exact rank; (18313,120) pair stays conditional on a rank being 1.
- CI: everything with a committed output passes under 2.17.4, but the lines must go through `hm` (scripts exit 0 on DIFFERS), `i8_checkpoint.py 360 3` must go, and i1_negative.py cannot run under -I.

## 9. Commit and loose ends
- Committed on branch run11_z17 in /home/claude/conductor/z17 as **4150673** (not pushed): experiments/roadmap-i (the bundle minus the 10 do-not-commit files and __pycache__: 117 files, 2.98 MB) and experiments/referee_z17 (scripts, outputs, ci/ logs, a short README in the referee folders' style).
- On that tree, through CI's `hm` helper (ci/hm_clean_tree.txt, machine less loaded by then): i1_laws 1 s, i1_laws --negative 1 s, i1_negative 3 s, i2_residual --negative 4 s, i4_exchange --negative 1 s, i4_exchange 0 s, i5_rank --negative 2 s, i8_pairct --negative 21 s, i8_pairct --compare50 9 s, i3_scan --negative 27 s: **all PASS**. Control on the wiring: with a tampered expected output the bare line prints `# DIFFERS` and **exits 0**; `hm` fails, as it must.
- [fix] `scripts/verify_manifest.py` fails on the branch (1,677 listed, 1,839 tracked): the conductor's manifest commit is needed, as in run 10.
- [note] I4: DM on (16065,2032)'s partner (18097,14033) already ran inside chunk 16000's DM (all 8 curves "rank not 1"), so chunk 18000 can only change its second-descent and rank-0 stages; REPORT_I4's "the partner ... in chunk 18000" reads as if it were untouched.
- My own reruns: z17_i2check, z17_i8count, z17_i8compare, z17_neither a-e all print `# matches` when rerun (part d: 92 s on the second run).
