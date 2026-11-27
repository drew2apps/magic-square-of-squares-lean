# z9 — hostile referee and claims editor (run 7): H2, H4, H6, G3 block 13000

Start: Fri Nov 27 13:19:28 UTC 2026. Written by an AI model; not read by a person.

Scratch: /tmp/z9/. int7 not edited.

## Log

### H2 (13:20-13:45 UTC)
- Recount of h2_ex_15000_{0,1,2}of3.txt with my own parser (/tmp/z9/recount.py): 15,000 values of m (15001..30000, none
  missing, none twice); coprime count = phi(m) for every m; per-m sums consistent; total coprime 205,177,458, degenerate
  4,579, non-degenerate 205,172,879; no point 203,267,359; no class 1,875,086; class 30,434 (= number of ratio lines);
  need a good prime 2,807; largest P 229 ((25789,12060)); no SURVIVOR. Per-block table of REPORT_H2 matches exactly.
- Degenerate count 4,579 (and per block 304, 305, ...) recomputed from Pythagorean parameters (/tmp/z9/degcount.py): equal.
- Definition: every coprime 0 < n < m, minus the two degenerate types, then CL.system/CL.solve (system calls C.empty).
  Same as class_selmer_linalg --list, which made G2's class_selmer_15000_list.txt. No congruence pre-filter in the
  exhaustive run (the hunt pre-filters by modp_classes.txt, which are obstructed classes: only a speed-up).
- h2_exhaust.py 14900 15000 0 1 rerun here (PARI 2.17.4, 44 s): output identical to h2_ex_14900_0of1.txt ("# matches").
  h2_control.py and --negative: "# matches" (128 = 128; P = 0 iff G2 bad-alone EMPTY; P = G2's largest kill-set good prime on 13 of 13).
  h2_summary.py: "# matches". Every line of REPORT_H2's four tables is a line of h2_summary.txt.
- Sample (/tmp/z9/sample_h2.py, output sample_h2.txt): 20 random class-keeping ratios + the 3 with largest P (229, 197, 197),
  class status by the repository's criterion with sympy factorisations (not PARI), P from G2's g2_killsets.run (irredundant
  kill set; P = its largest good prime): 23/23 agree on dim, log2 space and P. 20 uniformly random ratios: all "no point",
  none listed: agree. Further sample (sample_h2b.py) below.
- Hunt rows inside 15000 < m <= 30000 (185 distinct) agree exactly (dim, log2, P) with the exhaustive output; the 88 hunt
  rows with m <= 15000 are all in G2's class_selmer_15000_list.
- **FAILS (count):** the hunt has 2,251 class-keeping LINES but only **2,234 distinct ratios** (sampling with replacement:
  16 ratios occur twice or three times, 17 repeats; e.g. (15785,1879) three times in k13/k13_r1). And 273 of the 2,234 have
  m <= 30000 (88 with m <= 15000, G2's range; 185 inside the exhaustive range). So "2,421 ratios keep a class outside the
  exhaustive range" is false: distinct class-keeping ratios with m > 30000 are 1,961 (hunt) + 170 (design) = **2,131**
  (still above the roadmap's 2,000). Distinct over both hunts: 2,404.
- **FAILS (statement):** "Time limits 600 s (class) and 1800 s (descent); no ratio hit either." h2_hunt_k36.txt: "timed out:
  class 3". Three ratios at k = 36 hit the class limit (class status unknown; not counted as kept).
- Wording: REPORT_H2 part 1 says "the odd part of m^2 + n^2 at most 2"; the code keeps ratios whose odd part of m^2+n^2 has
  at most 2 distinct primes.
- Range: the hunt sampled k = 13..30, 33, 36 only: m < 2^37 (about 1.4e11), class-keeping ratios only for m < 2^34 (k = 33
  is the last with any). "m from 10^4 to 10^12" is the design, not what ran. Design k = 16..28 is 2^16 <= m < 2^29 (report: "2^16 to 2^28").
- Design (excess-guided): 3,849 taken, 171 class, 1 timeout ((117702451,107508765), k = 26, dim 2), 170 finished, 1 needs a good
  prime (P = 41), none of the 170 in the hunt. 3.4x: 4.44% vs 1.29%: holds.
- HANDOFF "about 2.5 core-hours per block of 1000": times.txt gives 8,274 + 11,480 + 11,568 = 31,322 core-s for 15 blocks =
  0.58 core-h per block on average (the last blocks have about 1.9x the ratios of the first). 2.5 is not supported.
- REPORT "the largest P per block rises from 101 to 229": block maxima are 149, 181, 137, 149, 137, 101, 109, 137, 173, 181,
  229, 193, 197, 193, 181: they range from 101 to 229, not monotone.
- "P <= 1153 by H3": REPORT_H2 states the dependence (Hasse-Weil, quoted, through H3); the proposed LEDGER text ("By H3, a
  cutoff above 1153 is never needed") does not. Also "emptied" rests on G1's exact_local.py (run 6: "tested, not proved")
  and the repository's one-sided 2-adic test; REPORT_H2 "What it rests on" omits exact_local's status.

### H4 (13:45-13:55 UTC)
- Chunking argument checked against the code: exact_local.join builds tuples coordinate by coordinate and keeps s + (d,)
  iff the label prefix is in the projection of each allowed set; at the last coordinate the projection is the allowed set
  itself, so join(divs, loc) = {tuples of the product whose label vector is allowed at every prime}, a per-tuple filter;
  dlabel(d) depends on d and the prime only. So join([[d0]] + divs[1:], loc) is the part with first coordinate d0 and the
  union over d0 is the whole join (disjoint, so sorted(concatenation) = sorted(join)). stop_first is used only in
  emptiness tests (left_bad uses first=False). Holds. (The hunt's need() also uses H4.chunked, for emptiness only.)
- h4_split.py --negative here: "# matches" (6.8 s). It is a consistency check (split = unsplit on 4 ratios; a degenerate
  ratio keeps its known tuples), not a forged-input test: no forged kill set is refused. Call it "consistency checks".
- Rerun of (88373348,70002375) through h4_split.line() under ulimit -v 6500000 (other agents' jobs share the 8 GB):
  identical line, 126 s. The other two (log2 34, 36) not rerun (memory shared with other agents).
- The table and LEDGER text agree with h4_split.txt. "All 63 of G1's kept ratios are now done: 60 by bad primes alone,
  3 with good primes, the largest 61": 57 (G2, g2_g1_kept_0.txt) + 3 = 60 holds (I did not recount G2's 57; it is run 6's reviewed file).

### H6 (13:55-14:20 UTC)
- Partner: 627 + 364 = 991, 627 - 364 = 263, gcd 1. Stronger than the report says: with (m', n') = (m + n, m - n) and
  t' = (t + 1)/(t - 1), (t - 1)^2 F'_j(t') equals L, U, 4 F4, 16 F0, 16 F8 for j = F0, F8, F4, L, U (sympy, generic m, n;
  /tmp/z9/h6_ident.py). The constants are squares, so the two (Z/2)^5-covers ARE isomorphic over Q, by an identity.
  "That the two covers are isomorphic is not proved here" can be replaced by that identity; the "mirror" is explained, and
  the matching of bounds across the pair checks model-independence of the code, not new data. All 10 Richelot bounds and
  the genus-1 ranks match under this dictionary (e.g. 627 F0F8F4L-wP 9 <-> 991 F0F4LU-wP 9; 627 c4 F4L rank 3 <-> 991 c2 F0F4 rank 3).
- sigma: t -> -t and tau: t -> 1/t (times t^2) identities, and t^6 F0F4L(1/t) = F8F4L, F0F4L(-t) = F0F4U,
  t^6 F0F4L(-1/t) = F8F4U: checked with sympy for generic m, n. The bundle has no script for these identities (only the
  Frobenius polynomials agreeing at p = 31). J(cover) ~ prod_S J(C_S) is a quoted decomposition (Kani-Rosen type); the four
  C_S are distinct factors, so J(cover) ~ A^4 x B holds given that theorem. DM with 4 independent maps needs rank J(A) <= 3: holds.
  Status "proved: ... written out here, not machine-checked": per the brief, write "a written argument resting on a quoted
  decomposition theorem; the identities checked by sympy (review z9)".
- h6_richelot.py read in full. Kummer map (G1(D), G2(D)) on J(C) is the connecting map of the dual isogeny; sampled points
  give differences delta(P - P0), quadratic pairs give delta(Q1 + Q2 - K) (Res with the monic h), points at infinity
  lc(G_i): all lie in the true image. Certification |W_p||W'_p| = |H^1(Q_p, mu_2^2)| = 16 / 64 / 4 is right (exact
  annihilators under Tate local duality), so found subgroups of product exactly that size are the images (and the code
  asserts the product never exceeds it). S = {2} u primes of disc(f) lc(f) contains every bad prime of J and of J'
  (isogenous). Rank formula r <= s1 + s2 - 2 - max(t1, t2): re-derived from |X/2X| = |X/phi'Y||Y/phi X| /
  [ker phi'(Q) : ker phi'(Q) n phi X] with dim X[2] >= 2 and t = dim delta(rational kernel): holds (and is not sharp).
  Gaps (minor): the docstring says the twist of C' is fixed by counts over F_p and F_{p^2}, the code counts over F_p only
  (at <= 12 primes; p = 1 mod 4 cannot tell Delta from -Delta); quadratic-pair images rest on PARI's nfislocalpower.
- Ran here (PARI 2.17.4): h6_richelot.py --negative "# matches" (3.6 s); --controls "# matches" (10/10 equal to PARI's
  sum of ranks); h6_ranks.py --negative "# matches"; h6_split.py 627 364 and 991 263: 39 pieces, "# matches".
  Negative (2) ("dropped generator ... product 8") is arithmetic on sizes, not a run of the code path.
- Genus-1 table of REPORT_H6 against h6_ranks_627_364.txt: every class, map count and rank matches (29 genus-1 pieces;
  10 genus-2 pieces, all "simple over Q" at the primes stated). 991's sentences (c7 4, c8 2, c9/c10 [1,3]; F0F4 3, F4L 2) match.
- **Gap in "no route opens" (status of the genus-1 ranks):** for several classes PARI's lower bound exceeds the number of
  points it returned, so it rests on parity (Sha[2^inf] finite), not on points: (627,364) c3 (F0L, F0U, F8L, F8U; rank 4,
  3 points), c5 (L*U; rank 1, **no point**), c6, c9, c16; (991,263) c1 (F0*F8, no point), c3, c7, c11, c12. Two matter:
  c3 ("one short of DM": unconditionally 3 <= rank <= 4, and rank 3 would give DM with 4 maps) and c5 of (627,364)
  = c1 of (991,263) (the same curve under the partner map): unconditionally rank 0 or 1 (root number -1, conductor
  34,500,574,027,920, my PARI 2.17.4 run; ellrank effort 4 on the first curve of the class: [1,3], no point). The Jacobian of y^2 = L*U
  (627,364) is the minimal model of y^2 = F0*F8 (991,263) (same model: the partner map again), torsion Z/2 x Z/4, and the
  quartic has exactly the 8 obvious points t = 0, inf, +-1, so rank 0 is not excluded by counting; rank 0 would close
  the ratio outright (only degenerate t). It is excluded only under finiteness of Sha (2-parity, root number -1). So "no Demjanenko-Manin or Chabauty route opens" is conditional on these PARI parity bounds.
- h6_richelot.py 627 364 rerun here (PARI 2.17.4, 155 s): output identical to h6_richelot_627_364.txt.

### G3 block 13000 < m <= 14000 (14:20-14:45 UTC)
- g3_table.py and --negative in my copy: "# matches g3_table.txt, CENSUS_TABLE_G3.tsv"; "# 4 of 4 as they must be".
  916 = 540 square (476 + 34 + 24 + 1 + 5) + 376 twice (327 + 28 + 15 + 5 + 1). Block 13000: 305 ratios = 270 second
  descent + 23 rank 0 + 9 DM + 1 special + 2 open, the same as stage_13000.json (after second descent 35, after rank 0 12,
  after DM 3: (13103,4271), (13348,3795), (13899,1180)).
- Family by my own brute force over all pairs 11000 < m <= 14000 (/tmp/z9/g3_family_z9.py, from gen_pairs.py's prose:
  coprime, m^2+n^2 a square or twice a square, no prime = 5 mod 8 to an odd power, not in the 484 classes): 916 (540 square,
  376 twice; none removed by the congruences in this range), identical to fam_11000_30000.json restricted, and to the
  916 rows (m, n, type) of CENSUS_TABLE_G3.tsv; block 13000: 305 (176 square, 129 twice).
- Open list: for each of the 6 and its partner, every sd/dm/r0 record read: all sd "open", all DM "rank not 1", no
  rank-0 curve "closes"; (12703,2767)'s partner (7735,4968) is OPEN in the 1000-11000 table. Holds.
- Route row (13103,4271): twice type, partner ((13103+4271)/2, (13103-4271)/2) = (8687,4416), square type (9745^2), whose
  row in CENSUS_TABLE_1000_11000.tsv is "special, rank-1 route, d = 111925826", label conditional on (a), (b); copied
  with its label. Holds.
- run 6's 55 partial block-13000 second-descent records (git show run6:...sd_13000_{0,1}.jsonl) are reproduced
  identically (apart from "sec") by the restarted block.
- Bundle copy vs int7: the bundle's CENSUS_TABLE_G3.tsv differs from the rebuilt one ONLY in the (12007,1799) row
  (label "... root number -1, PARI bounds [1,3]); plus (a), (b)" against run 6's "... (D1 of the partner: root number -1,
  PARI bounds [1,3]); plus the DM generator condition"). Cause: the bundle's g3_table.py reads the F0 table from
  ../run5/f0/ (the session's pre-run-6 copy). The bundle's g3_table.txt equals int7's. The bundle's g3_table.py,
  gen_family.py, pipe_g3.py, run_g3.sh, route/run_route.sh differ from run 6's only in relative paths (int7 keeps run 6's).
  The bundle also carries chunk_*.json, dm/r0 *_in*.json, run_g3.log, stages.log (inputs and logs, not integrated; run 6
  left such files out too).
- **Regression in int7:** experiments/roadmap-g/g3/REPORT_G3.md is the bundle's copy: it drops run 6's banner
  "Read NOTES.md first ..." and run 6's path fixes (../../census2000, ../../census-wide -> ../census2000, ../census-wide).
  It has no block-13000 content (still "11000 < m <= 13000", 611). Restore run 6's file: `git checkout run6 --
  experiments/roadmap-g/g3/REPORT_G3.md`. NOTES.md (the authority) still describes 611 ratios only: it needs a
  paragraph for block 13000 (text in the Final section).
- H7/HANDOFF on the block: (13348,3795) is OPEN in the G3 table and closed only by H7's separate g4*K step (PARI [1,1] in
  h7_b13000_g4k.jsonl, while G3's own rank-0 stage recorded g4*K at PARI [1,3] for the same pair; PARI bounds vary with
  effort, the exact rank 1 is H7's); (13899,1180): h7 gives PARI [2,2], exact bound 2 (dim Sel^2 = 4, pairing rank 0): holds.
  REPORT_H7 "finished at 12:23 UTC": chunk 10:04 + 7829 s + stages 443 s = about 12:22: holds.
- H4 (continued): the other two ratios rerun too, (120756965,14608109) under ulimit -v 5000000 (176 s) and
  (711343101,552021980) under 5500000 (177 s): both lines identical to h4_split.txt. All three reproduce on PARI 2.17.4.
- G3 sample rerun under PARI 2.17.4 (patched copies of census-wide's sd_pari.py and wide_pari.py in /tmp/z9/samp, as in
  referee_x6/samples/samples.sh): second descent on (13001,9679) closed, (13348,3795) open, (13447,1303) open: identical
  apart from "sec"; global.py on (13015,4872), (13049,5543), (13348,3795): identical (closes by w*f0*f8, w*K; open);
  dm_pari.py on (13439,7440), (13447,1303) with partners: 16 records, 14 identical, (13439,7440) Sig- differs in the last
  digit of the floating hG (6.659894683951019 vs ...018), (20879,5999) Sig+-partner has PARI bounds [1,3] here against
  [3,3] (status "rank not 1" either way). No closure changes.
- h6_ranks.py 627 364 rerun here (PARI 2.17.4, 82 s): identical to h6_ranks_627_364.txt except the line naming the PARI
  version (so "# DIFFERS" by that line only; not a CI line). Every class, bound and point count is the same.

## Final (13:55 UTC)

Written by an AI model; not read by a person. Every rerun here was on PARI 2.17.4 (CI's); the bundle used 2.19.1.

### Verdicts
| Claim | Verdict |
|---|---|
| H2 exhaustive 15000 < m <= 30000: 205,172,879 ratios, 30,434 keep a class, all emptied, largest P 229, 2,807 need a good prime | **holds** (recounted; degenerate count recomputed; 63 sample ratios (23 class-keeping incl. the 3 largest P, 20 uniform, 20 with points everywhere) and 2 whole values of m recomputed with sympy factorisations and G2's code; control block rerun identical). Gap: "exact" rests on G1's exact_local.py (tested, not proved) and the one-sided 2-adic test |
| H2 definition of the ratios vs G2 | **holds** (same criterion and degeneracy test as class_selmer_linalg --list; no congruence pre-filter needed) |
| H2 hunt + design "2,421 ratios keep a class outside the exhaustive range" | **fails**: 2,251 hunt lines are 2,234 distinct ratios; 273 of them have m <= 30000. Beyond 30000: 1,961 + 170 = **2,131** (still > 2,000) |
| H2 "no ratio hit either [time limit]" | **fails**: 3 class-criterion timeouts at k = 36 |
| H2 "m from 10^4 to 10^12" | **overstated**: sampled m < 2^37; class-keeping ratios m < 2^34; design 2^16 <= m < 2^29 |
| H2 unfinished (117702451,107508765) | holds (reported as unfinished; not in the 170) |
| H2 "P <= 1153 by H3" | **withdrawn** (conductor, 14:00: review z8 found H3's Hasse–Weil step wrong): delete every such sentence; my replacement texts below no longer contain it |
| H2 "largest P per block rises from 101 to 229"; HANDOFF "about 2.5 core-hours per block" | wording fails (range, not monotone); 2.5 not supported (0.58 core-h per block, times.txt) |
| H2 h2_summary.py, h2_control.py (patched), --negative | **holds** ("# matches" each) |
| H4 three ratios emptied by bad primes alone | **holds** (all three rerun identical; chunking = join checked against exact_local.join) |
| H6 partners | **holds, stronger**: t -> (t+1)/(t-1) maps the five forms up to square constants: the covers are isomorphic |
| H6 four pieces one curve, J ~ A^4 x B, DM needs rank <= 3, Richelot bound 8 | **holds with gaps**: identities checked (sympy); decomposition quoted; Richelot logic and certification re-derived, rerun identical; "proved" too strong; twist fixed by F_p counts only |
| H6 genus-1 table vs h6_ranks_*.txt | **holds** (and h6_ranks 627 364 identical on 2.17.4) |
| H6 "no Demjanenko-Manin or Chabauty route opens" | **holds with gaps**: conditional on PARI parity bounds (L*U of (627,364): no point, rank 0 or 1; c3: 3 points, rank 3 or 4) |
| G3 block 13000: 916 ratios, open list, route row, (12007,1799) correction kept | **holds** (brute-force family identical; records read; samples rerun identical up to one last-digit hG and one PARI bound with the same status) |
| G3 docs in int7 | **fail**: REPORT_G3.md is the bundle's copy (run 6's banner and paths lost); NOTES.md and LEDGER row still say 611 |

Tier: everything here is **LEDGER.md tier D** (experiment, computed), not released. Nothing is tier C.

### Exact replacement text

**1. REPORT_H2.md, "Proposed LEDGER / CHANGES text" (replace the quoted paragraph):**
> **Experiment, not counted.** Survivor hunt for the four-class descent (experiments/roadmap-h/h2). Of the 205,172,879 coprime ratios m : n with 15000 < m ≤ 30000, 0 < n < m, not of the degenerate types, 30,434 keep a class under the repository's criterion, and the four-class descent with their bad primes and the good primes ≤ 229 empties every one (2,807 need a good prime). Two random hunts beyond that range, a biased one (10⁴ ≤ m < 2³⁷) and one steered by H3's tables (2¹⁶ ≤ m < 2²⁹), add 2,131 distinct class-keeping ratios with m > 30,000, all emptied, the largest cutoff 137 (the descent of one more ratio did not finish). No survivor. Computed (exact arithmetic; the local test is G1's exact_local.py, tested in run 6, not proved; at 2 the repository's test). Tier D. Written by an AI model; not read by a person.

**2. REPORT_H2.md, local edits:**
- "Biased hunt, m from 10⁴ to 10¹²:" -> "Biased hunt, m from 10⁴ to below 2³⁷ (ranges k ≥ 37 not started; the last class-keeping ratios are at k = 33):"; same in the title ("m from 10⁴ to 2³⁷").
- "174,580 ratios sampled over 30 ranges, 2,251 keeping a class" -> "174,580 ratios sampled over 30 ranges (with replacement), 2,251 class-keeping draws, 2,234 distinct ratios".
- "m from 2¹⁶ to 2²⁸" -> "2¹⁶ ≤ m < 2²⁹".
- "Together 2,421 ratios keep a class outside the exhaustive range, above the 2,000 the roadmap asks for." -> "Together 2,404 distinct ratios keep a class; 2,131 of them have m > 30,000 (outside the exhaustive range and G2's), above the 2,000 the roadmap asks for."
- "In the exhaustive run the largest P per block rises from 101 to 229" -> "In the exhaustive run the largest P per block lies between 101 and 229 (149 in the first block, 229 in 25000–26000, 181 in the last)".
- "the odd part of m² + n² at most 2" -> "the odd part of m² + n² with at most 2 distinct primes".
- "Time limits 600 s (class) and 1800 s (descent); no ratio hit either." -> "Time limits 600 s (class) and 1800 s (descent); three ratios at k = 36 hit the class limit (class status unknown, not counted as kept); no descent of the hunt hit its limit."
- Delete every sentence resting on H3's bound: in the bottom line the sentences "By H3, good primes above 1153 never cut anything, so P ≤ 1153 ..." and "A survivor would have to keep a tuple through all 96 good primes p ≡ 1 mod 4 up to 1153" (keep "the largest P per block lies between 101 and 229 (G2: 181 for all m ≤ 15000)", with the heading "P grows slowly with m" -> "P per block"), and in "What it rests on" the line "Quoted: Hasse–Weil, through H3, for 'P ≤ 1153'". The excess-guided hunt used H3's tables only to choose ratios; nothing rests on them.
- "What it rests on", Exact: add "the local test of G1's exact_local.py (tested against the repository's test and explicit points in run 6, not proved) and the repository's one-sided test at 2".

**3. HANDOFF_ROADMAP_H.md, H2:**
- Bottom line -> "no survivor. Of the 205,172,879 coprime non-degenerate ratios with 15000 < m ≤ 30000, 30,434 keep a class, and each is emptied by its bad primes and good primes ≤ 229. A biased hunt (10⁴ ≤ m < 2³⁷; the ranges up to 10¹² were not reached) and an excess-guided hunt built on H3's tables (2¹⁶ ≤ m < 2²⁹) add 2,131 distinct class-keeping ratios with m > 30,000, all emptied, the largest cutoff 137; the descent of one more did not finish. The largest P per block of 1000 lies between 101 and 229. A rerun of 14900 < m ≤ 15000 agrees with G2 exactly."
- Next step -> "the exhaustive run past 30,000 cost about 0.6 core-hours per block of 1000 on average here (31,322 core-seconds for the 15 blocks, times.txt), growing with m."

**4. H4 (REPORT_H4.md and HANDOFF):** in the LEDGER text "Computed (exact)." -> "Computed (exact arithmetic; G1's exact_local.py, tested in run 6, not proved; at 2 the repository's test). Tier D." In "Checks", "**Negative mode (h4_negative.txt):**" -> "**Consistency checks (h4_negative.txt; no forged input is refused here):**".

**5. REPORT_H6.md, "Proposed LEDGER / CHANGES text" (replace):**
> **Experiment, not counted.** The genus-1 and genus-2 quotients of the (ℤ/2)⁵-cover for (627,364) and (991,263) (experiments/roadmap-h/h6). The two ratios are partners, and t ↦ (t + 1)/(t − 1) carries the five forms of (991,263) to those of (627,364) up to square constants, so the two covers are isomorphic. Every genus-2 piece is simple over ℚ, and a descent via a Richelot isogeny, with local images certified by the product formula of local duality, bounds its rank by 6 to 9. The four pieces F0·F4·L, F0·F4·U, F8·F4·L, F8·F4·U are one curve A, so J(cover) ~ A⁴ × B (by the quoted decomposition of the Jacobian of the cover), and Demjanenko–Manin through them would need rank J(A)(ℚ) ≤ 3; the bound is 8. Genus-1 ranks are PARI's bounds, and some of PARI's lower bounds rest on parity (Ш finite), not on points: y² = L·U of (627,364) (the same curve as y² = F0·F8 of (991,263)) has no known point of infinite order, so it has rank 1 only if Ш is finite (rank 0 would exclude the ratio); the class of y² = F0·L, which has four maps from the cover, has three known independent points and rank 4 only if Ш is finite. With these ranks, no Demjanenko–Manin or Chabauty route opens. Computed (genus-2 bounds in exact arithmetic on quoted theory: the Richelot construction, the Kummer map, Tate local duality, unramified images at good odd primes; PARI's nfislocalpower decides local squares for quadratic points). Tier D. Written by an AI model; not read by a person.

**6. REPORT_H6.md, local edits:**
- "This is an internal consistency check. That the two covers are isomorphic is not proved here." -> "The two covers are isomorphic: with (m′, n′) = (m + n, m − n) and t′ = (t + 1)/(t − 1), (t − 1)²·F′(t′) is L, U, 4F4, 16F0, 16F8 for F′ = F0, F8, F4, L, U (checked with sympy in review z9). So the equal bounds check that the code does not depend on the model; they are not new data."
- "*(proved: the argument is standard, written out here, not machine-checked)*" -> "*(a written argument resting on the quoted decomposition of the Jacobian of a (ℤ/2)⁵-cover into those of its hyperelliptic quotients; the identities for σ and τ checked with sympy in review z9)*".
- "The twist of C′ is fixed by equal point counts over F_p and F_{p²}." -> "The twist of C′ (Δ or −Δ) is fixed by equal point counts over F_p at 8 to 12 good primes." (same fix in the h6_richelot.py docstring).
- Under the genus-1 table add: "For c3, c5, c6, c9, c16 PARI's lower bound exceeds the number of points it returned (it uses parity, so Ш finite); c5 (L·U) has no known point of infinite order. For (991,263): c1, c3, c7, c11, c12."

**7. HANDOFF_ROADMAP_H.md, H6:** Bottom line -> "no route opens, with PARI's genus-1 ranks. The four genus-2 pieces F0·F4·L, F0·F4·U, F8·F4·L, F8·F4·U are one curve A up to isomorphism, so Demjanenko–Manin needs rank J(A) ≤ 3; the Richelot descent gives ≤ 8. Every other genus-2 piece has bound 6 to 9. The two ratios are partners (991 = 627 + 364, 263 = 627 − 364), and t ↦ (t + 1)/(t − 1) is an isomorphism of their covers. Two of PARI's genus-1 ranks rest on parity only: y² = L·U of (627,364) (rank 1, no point found; rank 0 would exclude the ratio) and the four-map class of y² = F0·L (rank 4 with three points)." Label -> "computed (genus-2 bounds exact on quoted theory; genus-1 ranks by PARI, some resting on Ш being finite)".

**8. HANDOFF_ROADMAP_H.md, H7, G3 sentence:** "G3 block 13000 < m ≤ 14000 finished (G3's table now runs to 14000, 916 ratios): of its two new open ratios, (13348,3795) closes the same way, exactly; (13899,1180) stays open (PARI rank 2)." -> "G3 block 13000 < m ≤ 14000 finished (the table now runs to 14000: 916 ratios, 6 open in the table, 3 conditional). Of the block's two open ratios, (13348,3795) is excluded the same way by H7's step (outside the table); (13899,1180) stays open (g4·K: PARI rank [2,2], exact bound 2). (13103,4271) closes through its partner (8687,4416), closed by the rank-1 route under (a), (b)."

**9. LEDGER.md, the G3 row (replace the claim cell):**
> **Degenerate types, `11000 < m ≤ 14000`**: of 916 ratios, 6 open in the table ((11165,7548), (12201,7360), (12204,9853), (12703,2767), (13348,3795), (13899,1180); (12201,7360) and (13348,3795) are excluded separately by H7's g4·K step, `experiments/roadmap-h/h7`), 3 conditional on a rank being 1 ((11441,1519), (11969,4529), (12007,1799)), the rest closed with the labels of `CENSUS_TABLE_G3.tsv` (PARI rank bounds; DM on PARI generators; (12416,5313) by the route, (12487,4241) and (13103,4271) through partners closed by the route, conditional on (a), (b)); other types not covered (run 6, after release 0.6.0, not released; reviewed: holds with gaps; block `13000 < m ≤ 14000` added in run 7 by a separate AI session; reviewed: holds)

**10. experiments/roadmap-g/g3:** `git checkout run6 -- experiments/roadmap-g/g3/REPORT_G3.md` (the bundle's copy only loses run 6's banner and path fixes). Append to NOTES.md:
> ## Block 13000 < m ≤ 14000 (run 7)
> Run by a separate AI session (roadmap H, H7) with the same programs; `g3_table.py --write` rebuilt the table here, so correction 1 above is kept (the bundle's own copy of the table had lost it, and differed in no other row). 305 ratios (176 square, 129 twice): 270 closed by the second descent, 23 by rank 0, 9 by Demjanenko–Manin, (13103,4271) through its partner (8687,4416), closed by the route in the 1000 < m ≤ 11000 table (conditions (a), (b)); open: (13348,3795), (13899,1180). The table now covers 11000 < m ≤ 14000: 916 ratios, 6 open, 3 conditional. (13348,3795) is excluded separately by H7's g4·K step (`../../roadmap-h/h7`), not in this table. Review z9 (`reviews/` or `conductor/reports/z9.md`): the family regenerated by brute force over all pairs (916, identical), every record of the six open ratios and their partners read, a sample rerun under PARI 2.17.4 (3 second descents, 3 rank-0 inputs, 2 Demjanenko–Manin inputs with partners: identical apart from timings, the last digit of one floating height and one PARI rank bound [1,3] against [3,3] with the same status), and run 6's 55 partial records of this block reproduced by the restarted block.

### Could not check
H3's Hasse–Weil reduction (H2's 1153 rests on it); H7's rank step for (13348,3795); the full exhaustive, hunt and design runs (only samples, the control block and the rows inside the exhaustive range); h6_richelot 991 263; PARI's ellrank, nfislocalpower and factor internals.

### Memory
Nothing from memory was used.

## Packaging (14:00-, conductor's request)
Conductor's note (14:00): review z8 found H3's Hasse–Weil step wrong, so every "P ≤ 1153 by H3" goes; my replacement
texts above were edited to drop it (H2 LEDGER text, HANDOFF bottom line; REPORT_H2 edit list now says to delete those
sentences). Checks copied to /home/claude/conductor/r7/referee_z9/ (to be committed as experiments/referee_z9/):
recount.py, degcount.py, sample_h2.py, sample_h2b.py, h4_one.py, h6_ident.py, c5.py (c5 and c5b merged; no ellrank),
g3_family_z9.py, g3_samples.sh + g3_samples_compare.py (+ g3_samples_in/), _common.py, README.md, each with its .txt.
Tested in a copy of int7 at /tmp/z9/int7/experiments/referee_z9 (paths ../.. and ../roadmap-h, ../roadmap-g).
Done 14:40 UTC. Each script was run with --write and then again without it ("# matches" for all but h4_one.py and
g3_samples*, which were run once with --write: 350 s and 193 s of CPU). CPU times: recount 8 s, degcount 2 s,
h6_ident 1 s, c5 1 s, g3_family_z9 11 s, sample_h2 35 s (CI); sample_h2b 123 s, h4_one 350 s, g3_samples 193 s (not CI).
g3_samples_compare.txt records the two differing DM fields (one hG last digit, one PARI bound [1,3] vs [3,3]).
