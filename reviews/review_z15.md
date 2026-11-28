# Review z15: hostile referee of EXPLORATION_4 (2026-11-27), sections 1-4 and 6

Written by an AI model (reviewer z15), not read by a person.

> Committed as written apart from this note. In this repository the reviewer's folder is `experiments/referee_z15/` (with `session_data/`, a copy of the session's records); the note reviewed is `experiments/research-pairing/explore4/`. Paths under `/home/claude/` are the reviewer's scratch folders.

Worktree: /home/claude/conductor/z15 (branch run11_z15). Scripts and outputs: experiments/referee_z15/.
Note reviewed: /mnt/user-data/uploads/2b6e4a44-EXPLORATION_4_2026-11-27.md (it is NOT inside
/home/claude/untrusted/explore4/explore4/; that directory holds only the scripts and outputs).

## Log (appended as I go)

### Task 1 (five survivors with the repository's code)

Ran the repository's `t10_pairct.analyse_ratio` (byte-identical copy `experiments/referee_z15/t10_pairct_copy.py`,
no hints, so PARI found every point anew) through the driver `z15_t10run.py` -> `t10_five.jsonl`.
Residual tuples are t10's own (`TD.magic` + `exact_local.allowed_at` + `join` at the bad primes).

| ratio | residual | survive all six (t10) | least killing good prime of survivors | dim Sel2 08/LU (t10) |
|---|---|---|---|---|
| (3320,2303) | 2 | (2,1,1,635399), (2,1,1270798,2) | 29, 29 | 5 / 2 |
| (5623,1017) | 2 | (2,1270798,1,2), (635399,1,1,2) | 29, 29 | 2 / 5 |
| (9113,3200) | 2 | (1,2,511,1759), (1,2,3518,1022) | 37, 37 | 4 / 3 |
| (12313,5913) | 2 | (1022,3518,2,1), (1759,511,2,1) | 37, 37 | 3 / 4 |
| (14903,9696) | 8 | (1,2,10414,49198), (1,2,24599,5207); the other 6 killed | 37, 37 | 4 / 6 |

Exactly the note's 10 tuples survive; Sel2 dimensions equal the session's `dims` for 08 and LU at all five
(the session's mixed-curve dims are 6/6/7/7/9 as t10's). No 'undecided', no 'outside Sel2'. CPU 5-28 s per ratio.
The note's table says "good primes 17, 29" for (3320,2303): 17 is in G2's kill set of the ratio, but it kills
neither survivor (t10's least killing prime is 29 for both) -- [note] wording.
In every case the survivors' verdict is 'image' on the curve where the ghost shape makes the class torsion
(LU for S1, 08 for S2) and 'kernel' (in the kernel, not a class of a point found) on the other five.

### Task 1, rank bounds -- `z15_rankbound.py` -> `z15_rankbound.txt`

Full Cassels-Tate matrix on a basis of Sel2 (t10's model, `rank_second.full_selmer2` and `DFull.pairing`; every
entry decided, matrices symmetric and alternating). rank <= dim Sel2 - 2 - rank(matrix):

| ratio | E_08: dim, pairing rank, bound | E_LU: dim, pairing rank, bound | PARI ellrank (effort 2) 08 / LU | w(08), w(LU) |
|---|---|---|---|---|
| (3320,2303) | 5, 2, <=1 | 2, 0, 0 | [1,1] / [0,0] | -1, +1 |
| (5623,1017) | 2, 0, 0 | 5, 2, <=1 | [0,0] / [1,1] | +1, -1 |
| (9113,3200) | 4, 0, <=2 | 3, 0, <=1 | [0,2] / [1,1] | +1, -1 |
| (12313,5913) | 3, 0, <=1 | 4, 0, <=2 | [1,1] / [0,2] | -1, +1 |
| (14903,9696) | 4, 0, <=2 | 6, 4, 0 | [0,2] / [0,0] | +1, +1 |

All ten bounds equal the note's table, and the note's "PARI agrees" holds. Every surviving tuple's class pairs
trivially with all of Sel2 on E_08 and E_LU (not only with a basis modulo point images).

### Point search -- `z15_findpoint.py` -> `z15_findpoint.txt`

At (3320,2303), PARI's `ellrank` at effort 5 (a finder) found a point of E_08 (t10's model, checked on the curve);
with it the surviving S1 class (1, 19505) **is the class of a rational point of E_08** (in the span of the torsion
images and that point's image). So t10's verdict 'kernel' there meant "not a class of a point found", and no
obstruction on E_08 at any descent level can remove these tuples: E_08 has rank >= 1 and the class is in
delta(E_08(Q)). (z15_findpoint2.py, running, repeats this and carries the point to E_LU of the partner (5623,1017).)

### Tasks 3-4 -- `z15_residual_all.py`, `z15_shapes.py`, `z15_shapes2.py` (own code; PARI `core`)

- Residual sets recomputed with the repository's code alone (TD.magic + exact_local at all bad primes + join) for
  all 937 ratios: 6,904 tuples, **identical to the session's out_*.jsonl at 937 of 937 ratios**; no known tuples.
- Ghost shape: 2,042 tuples (1,378 S1 + 664 S2) at 499 ratios; 2,032 killed (1,556 on E_08 or E_LU, 476 only on
  mixed curves); 10 survive, all S1 (6) or S2 (4). All numbers of section 2 reproduced.
- Consistency: no tuple with a = 1 is killed on E_08 and no tuple with b = 1 on E_LU (0 of 6,904): as it must be
  for a trivial / torsion class. At every ratio all S1 tuples have one E_08 verdict, all S2 one E_LU verdict.
- Cosets, exact test (2^rank(V - v0) = |V|) on all 937: coset at 682, **not a coset at 255** (the note's 247 plus its
  8 untested large sets, all non-cosets). 189 of the 255 fail already because |V| is not a power of 2.
- Sizes: the note's size list is reproduced; no odd size. 876 of 937 sets are invariant under a translation
  (283 by (1,1,c,c), c = sqf(2(m^2-n^2)); 141 by (c,c,1,1), c = sqf(2mn); 452 by other vectors); 61 even sets have no
  translation symmetry. The note gives no reason for evenness; it only states it.
- Rates: 101/997 = 10.13 %, 836/8739 = 9.57 %, 1/997 = 0.100 %, 4/8739 = 0.046 % -- reproduced.

### Task 5, root numbers -- `z15_models.py`, `z15_rootno.py` -> `z15_models.txt`, `z15_rootno.txt`

- The session's closed-form models of E_08, E_LU, E_mix are isomorphic over Q to t10's six curves (equal minimal
  models) at 120 of 120 (ratio, curve) checks (the 5 survivor ratios + 15 drawn with seed 5).
- The session's five populations, redrawn exactly as rootno.py draws them: every count of rootno.txt reproduced.
  The 937 line of rootno2.txt reproduced (249/216/230/242; w_mix 493/444, z = +1.60).
- Complete populations instead of samples: all 8,739 (5000<m<=15000): chi2(3) = 6.95; all 9,736: chi2 = 8.50
  (p ~ 0.04); the 8,799 without residual tuples: chi2 = 10.93 (p ~ 0.012), phi(w08,wLU) = -0.026 (z = -2.45); the
  excess is in the pattern (w08,wLU) = (-,+) (2,550 of 9,736 against 2,434 expected). So "about a quarter" is right to
  within ~1-2 %, but the complete populations show a small, statistically visible departure; "no correlation"
  holds for the 937 (chi2 2.68, phi z = +1.48) and the session's samples, not exactly for the full lists.
- Subpopulations the note does not mention: w_mix = +1 at 65 of the 101 ratios with m <= 5000 (z = +2.89); random
  ratios "with a 2-adic point" w_mix z = -2.08. With ~15 tests a single 2.9 sigma is borderline (p ~ 0.06 after a
  Bonferroni correction), not a finding either way.
- 2-Selmer parity (Monsky; Dokchitser-Dokchitser): (-1)^(dim Sel2 - 2) = w at 937 x 6 of the session's dims and at
  all t10 dims of my reruns. So the root numbers carry no information beyond the parity of the Sel2 dimensions the
  pipeline already has; and this is an independent consistency check of the session's Sel2 dims (passed).
- At the 5 survivor ratios only N = 5: no test of "root number vs survival" is possible. On the curve that carries
  the non-torsion ghost class, w = -1 at (3320,2303), (5623,1017) and w = +1 at the other three.
- Verdict: the note's numbers are right; "the four patterns about equally" is supported; "Parity gives no handle on
  (R)" is a conclusion the data cannot carry (5 survivors; the parity is a function of dim Sel2 anyway).

### Task 3, my own derivation -- `z15_degenerate.py` -> `z15_degenerate.txt`; `z15_findpoint2.txt`

- By hand (t10's Gaussian forms; i ~ 2 and 2i ~ 1, -1 ~ 1 modulo Q(i)*^2):
  t = 1: a0 a8 = iC, aL aU = (2i)^2 (m^2-n^2); t = -1: a0 a8 = -iC, aL aU = 4(m^2-n^2). With C a square:
  d0d8 ~ 2, dLdU ~ m^2-n^2, so (a,b) = (sqf(2mn), 1) = S1. t = 0: a0 a8 = mn, aL aU = -2iC; t = inf: a0 a8 = -mn,
  aL aU = 2iC. With C = 2c^2: d0d8 ~ mn, dLdU ~ 2, so (a,b) = (1, sqf(2(m^2-n^2))) = S2. The note is right.
- Machine check of the same (own class function sqf(2(Re z + |z|)), independent of tuple_descent): all 910
  degenerate ratios m <= 3000: 1,068 solutions S1 (square type), 752 S2 (twice-square type), 0 other; the products
  (d0d8, dLdU) equal those of tuple_descent.magic's known tuples at 910 of 910. (ghost.py checked 13.)
- Sign: products of square norm have a rational class only up to -1 (a square in Q(i)); on E_08 and E_LU the two
  signs differ by (1,-1,-1), the class of a 2-torsion point (computed on the models (0,A^2,C^2), (0,B^2,C^2)), so the
  sign never matters for torsion-ness or for the pairing.
- The E_08 claim: on y^2 = x(x-A^2)(x-C^2), T3 = (C^2,0) = 2P4 with x(P4) = C(C+B), delta(P4) = (C, 2mn, 2mnC); so
  (1,2mn,2mn) = delta(P4)(C,1,C); torsion classes are 1, (1,-1,-1), delta(P4), delta(P4+T1) (if the torsion is
  Z/2 x Z/4), so it is torsion iff C is a square. On E_LU, delta(P4') = (2C, 2A, AC), which is the "4-torsion
  shift" turning (2C, 2dLdU) into (1, A dLdU, A dLdU). Both verified.
- The split: a and b depend only on d0d8 and dLdU, so the ghost shape cannot depend on the (dL,dU) split; the
  mixed-curve classes (d0dL etc.) do. **At t10's degenerate controls every residual tuple, known or extra, has the
  ghost shape**: (12,5), (17,7), (627,364), (991,263) each have 2 extra tuples, S1 resp. S2, with a non-degenerate
  split. The pilot already found that the extras at (12,5), (17,7) are pairing-killed and those at (627,364),
  (991,263) are not (they die at 41). So the 10 survivors are the non-degenerate analogue of a phenomenon the pilot
  had at (627,364)/(991,263); the note does not connect the two.
- **New, unconditional:** at (3320,2303) PARI found (effort 4-5) a point of E_08 of infinite order (height 39.75),
  and the S1 class (1,19505) of both survivors is in the span of it and the torsion. Carried through the common
  minimal model to E_LU of the partner (5623,1017), the survivors' S2 class (16326209,1) is likewise the class of a
  rational point. So at the pilot's two exceptions the surviving tuples are classes of rational points on BOTH E_08
  and E_LU: no descent of any level on those two curves can remove them.

### The note's "Note" under the table of section 1 is wrong at (14903,9696) -- `z15_findpoint3.txt`, `z15_rankbound.txt`

The note says "a pair curve of rank 0 does not describe the exceptions; the description that fits all five is
section 2". But:
- (14903,9696): E_LU has dim Sel2 6, pairing rank 4, so rank 0 (rank_second; PARI [0,0]); its torsion is Z/4 x Z/2
  (elltors), so E_LU(Q) has 8 points and C_LU(Q) (which has the point t = 0) has exactly the 8 obvious points over
  t in {0, 1, -1, inf}. There F_0(+-1) = C and F_L(0) = 2C are not squares (C = 5^2 * 12644473). So this ratio has no
  four-form point by rank 0 of E_LU -- the pilot's mechanism for (5623,1017)/(3320,2303) (their E_08 resp. E_LU:
  rank 0, torsion [4,2], checked again here).
- (9113,3200) and (12313,5913): the curve that carries the non-torsion ghost class (E_08 resp. E_LU, root number +1,
  bound 2) has rank 0 or 2, not decided; PARI's ellrank at effort 3-5 finds no point (finder). The other curve has
  w = -1 and bound 1, so rank 1 if Sha is finite. The note's "neither E_08 nor E_LU has rank 0" is not established:
  if E_08 at (9113,3200) has rank 0 (torsion [4,2], so again 8 obvious points), the rank-0 mechanism closes this
  partner pair too.
- So "rank 0 of a pair curve" describes 3 of the 5 ratios for certain (modulo the repository's pairing code) and
  possibly all 5. The ghost shape fits all five, but it is not the only description, and the note's sentence
  should be replaced. [fix]

### Task 3: is "no tuple of any other shape survives" a tautology? -- `z15_shapes3.py` -> `z15_shapes3.txt`

Not a tautology of the definitions (survival = kernel on six curves; shape = a condition on d0d8, dLdU only), and
half of the shape is automatic for one curve: a = 1 means a trivial class on E_08, b = 1 a torsion class on E_LU,
never killed there (0 violations in 6,904). But **on this data the claim carries almost no evidence**:
- At 4 of the 5 survivor ratios the residual set has 2 tuples, both of ghost shape (228 ratios have such a set).
  Any survivor there is ghost-shaped whatever the pairing does.
- At (14903,9696) all 8 residual tuples have a = sqf(2mn) (d0d8 ~ 2); the 4 with b = 24599 are killed on E_LU,
  which has rank 0 (pairing kernel = torsion, so every non-torsion class is killed); the 4 S1 tuples reach the
  mixed curves and 2 survive. Again the shape of the survivors is forced before the mixed curves act.
- Ratio level, the only real test: tuples survive E_08 and E_LU at 217 ratios -- only ghost tuples at 119, only
  non-ghost at 88, both at 10. The 5 survivor ratios are all in the first group. Under "survival independent of
  shape" the chance is (129/217)^5 = 0.07, or (129/217)^3 = 0.21 if partner ratios count once.
- Tuple level: of tuples reaching the mixed curves, 10 of 486 ghost and 0 of 554 non-ghost survive, but the 10
  sit in 5 symmetric pairs at 5 ratios (3 up to partners), so the tuple count overstates the evidence.
Plainly: a real but weak regularity (p ~ 0.07-0.2), observed in 3 independent instances; the heading word
"answer" and "No tuple of any other shape survives" in bold are stronger than this. [fix]

### Task 6, claims editor (one line per number or status word of sections 1-4 and 6)

Legend: OK = reproduced by me (how); FIX = wrong or stronger than the evidence; NOTE = right, but wording/caveat.

Section 1
- "997 with m <= 5000 and 8,739 with 5000 < m <= 15000 (9,736)": OK (line counts of g2_997_all.txt, g2_15000_all.txt).
- "896 / 7,903 / 8,799 leave no tuple; 101 / 836 / 937 leave tuples": OK ('left' lines; arithmetic).
- "100 / 832 / 932 all killed; 1 / 4 / 5 survive": OK (session data; the 5 confirmed with t10_pairct).
- "pipeline.py takes the residual tuples from the repo's own code": OK (my repository-only recomputation equals
  out_*.jsonl at 937 of 937 ratios, 6,904 tuples).
- "The pilot had run 201 of these 937": OK (pilot set rebuilt from t10's seed equals the ratios in its files).
- "1,766 tuples, 10,596 decisions, 8,930 killed by both, 1,666 by neither": OK (my own parser of the pilot files).
- "All 6,904 map into Sel2 of all six curves (0 outside)": OK as data; NOTE: sel.py's `insel` tests only the
  primes of S, which is enough here because every class map only involves primes of 2mnC(m^2-n^2) (checked by
  hand); t10 tests it properly and finds no 'outside' in my reruns either.
- "Killed per curve 4,580 / 3,122 / 6,587 each; on E_08 or E_LU 5,864 of 6,904": OK (recounted). NOTE: "alone"
  in "Killed on E_08 or E_LU alone" reads as "only there"; the 5,864 are tuples killed on E_08 or E_LU (at least).
- Table, surviving tuples: OK (t10_pairct). m^2+n^2 factorisations: OK (PARI).
- Table, rank bounds <=1/0, 0/<=1, <=2/<=1, <=1/<=2, <=2/0 and "PARI agrees": OK (rank_second full matrix).
- Table, "repo's kill: good primes 17, 29" for (3320,2303): NOTE: 17 is in G2's kill set but kills neither
  survivor; t10's least killing prime is 29 (and 37 at the other three, as stated).
- "The first two are the pilot's exceptions; the other three are new": OK. Partners and (24599,5207): OK.
- "The census already closes all five with good primes": OK (29, 29, 37, 37, 37).
- "at (9113,3200) and (12313,5913) neither E_08 nor E_LU has rank 0": FIX (E_08 resp. E_LU: w = +1, bound 2,
  rank 0 or 2 undecided).
- "So 'a pair curve of rank 0' does not describe the exceptions": FIX (it describes (14903,9696): E_LU rank 0,
  torsion 8, only the obvious points -- as for the pilot's two; undecided for the partner pair).
- "The description that fits all five is section 2": NOTE (it fits, but at every one of the five it is forced by
  the residual set and rank 0 of E_LU at (14903,9696); see Task 3).

Section 2
- heading "(answer to the question in the task, on this data)": FIX -> "an observation on this data".
- "class on E_08 is (1,a,a); on E_LU after the 4-torsion shift (1,b,b)": OK (derived: delta(P4') = (2C,2A,AC)).
- "all 10 surviving tuples have (a,b) = S1 or S2": OK. "No tuple of any other shape survives" (bold): OK as a
  count, FIX as emphasis (weak evidence, p ~ 0.07-0.2 at ratio level).
- "exactly the E_08 and E_LU classes of the degenerate solutions (13 degenerate ratios)": OK (derived by hand;
  910 degenerate ratios m <= 3000 by machine).
- "(1,2mn,2mn) = class of the rational 4-torsion point times (C,1,C); torsion exactly when C is a square": OK
  (for C not a square the torsion is exactly Z/2 x Z/4, so the list of torsion classes is complete).
- "a surviving tuple is a 'ghost' of a degenerate solution": NOTE (a name; the mixed-curve classes differ).
- "(dL,dU) split not always the degenerate one (a factor 7 at (9113,3200))": OK (m+n = 7*1759, dL = 511 = 7*73).
  The shape does not depend on the split (only products); the mixed-curve kills do. The same non-degenerate
  splits occur in the pilot's degenerate controls (12,5), (17,7), (627,364), (991,263) (not mentioned).
- "2,042 at 499; 2,032 killed (1,556 / 476); 10 survive": OK (recomputed).
- "the only classes that get past bad primes and all six pairings are the classes of the degenerate solutions
  themselves": FIX (only their E_08/E_LU classes coincide; at these ratios there is no degenerate solution; and at
  (3320,2303)/(5623,1017) those classes are classes of rational points of E_08 and E_LU -- not ghosts there).
- "would reduce stage 2 to one explicit class per ratio": FIX (two shapes, so up to two classes, one on E_08 and
  one on E_LU; the ghost-shaped tuples per ratio number 2-32 with different mixed classes; and the explicit class
  can be a rational point's class, so the reduction would not by itself decide anything).
- "tested on 937 ratios, not proved, and I have no proof route yet": OK.

Section 3
- "coset at 682 of 929 (8 larger not tested), not at 247": NOTE: exact on all 937: 682 coset, 255 not (the 8 are
  all non-cosets); 189 of the 255 fail already by size (not a power of 2).
- "Sizes 6, 10, 12, 14, 18 occur": OK (also 20, 22, 24, 28, 30, 34, 36, 40, ...).
- "Reason: the local image of the genus-5 curve is a set, not a subgroup": OK as far as tested by me (odd primes:
  4,514 of 13,003 non-trivial local sets are not cosets; the note itself exhibits none).
- "Linearity holds after mapping to a pair curve; the pair-curve stage can be a matrix computation": OK (Sel2
  is a group; the pairing is bilinear) -- but the kill test is a pairing against Sel2, not a linear condition on
  tuples, so NOTE.
- "Every residual set has even size": OK (876 of 937 are translation-invariant; not explained in the note).
- "101/997 = 10.1 %, 836/8,739 = 9.6 %, 1/997 = 0.10 %, 4/8,739 = 0.046 %": OK.

Section 4
- All counts of rootno.txt and rootno2.txt: OK (redrawn draw for draw).
- "(w_mix) 493 / 444 (1.6 standard deviations)": OK (z = +1.60).
- "each of the four sign patterns about a quarter of the time in every population tested": OK for the tested
  samples; NOTE: on the complete lists the departure is small but visible (chi2 8.50 on 9,736, 10.93 on 8,799).
- "The local conditions for a solution do not select a parity": NOTE (supported for "having residual tuples";
  nothing can be said about survivors, N = 5).
- "Parity gives no handle on (R)": FIX (too strong; and parity = (-1)^(dim Sel2 - 2), already in the pipeline).

Section 6 (what was not done)
- The list is honest as far as it goes. Add: no rerun of the repository's own pairing code outside the pilot's
  201 ratios (done here for a sample); the survivors' classes were not tested against rational points (they are
  point classes at (3320,2303)/(5623,1017)); the rank-0 mechanism at (14903,9696) was missed; no timings.
- "A 4-descent on E_08 at the ghost ratios": NOTE: pointless at (3320,2303)/(5623,1017) (point classes); useful
  only at the three ratios whose ghost-carrying curve has w = +1 and bound 2, and only if that rank is 0.

### Task 4, the note's "reason" for non-cosets -- `z15_localcoset.py` -> `z15_localcoset.txt`

At the odd bad primes, the allowed label sets of `exact_local.allowed_at` (labels in (F_2^2)^4, dlabel a
homomorphism) are not cosets at 4,514 of 13,003 (ratio, prime) pairs with a non-trivial set, and every one of the
937 ratios has at least one non-coset odd local set -- including the 682 whose residual set is nevertheless a
coset (the intersection with the divisor supports can restore it). So the note's reason ("the local image is a
set, not a subgroup") is confirmed at the level of the local sets; p = 2 was not tested.

### Partners -- `z15_partners.py` -> `z15_partners.txt`

For both partner pairs, E_08 of one ratio and E_LU of the other have the same minimal model, and the mixed curves
agree too. So the five survivor ratios involve only three sets of curves: the 10 survivors are 3 independent
instances (each a symmetric pair of tuples), which is what the significance estimates above assume.

### Task 2, interim (sample 80 = seed 20261127 from the 736 non-pilot ratios; then ghost ratios, seed 7/8 order)

Interim at ~50 sample ratios + ~20 ghost ratios: 0 disagreements, residual sets and Sel2 dims identical, no
'undecided'. Final numbers below.

### Coordinator's question: tier A and section 2's "why this matters for (R)" -- `z15_tierA.py` -> `z15_tierA.txt`

(The note in the bundle is byte-identical to the copy I reviewed; the claims editor above uses its wording.)
Computed membership (the Lean status of the two criteria is the repository's, not re-checked here): of the 937,
844 lie in one of the 484 residue classes of `scripts/g1_full/modp_classes.txt` (m = r n mod p), 318 have a prime
= 5 mod 8 to an odd power in m^2+n^2 not dividing n, 866 either, 71 neither -- z16's numbers. All five survivor
ratios are excluded: (3320,2303) and (5623,1017) by the class mod 151 and by 2621, 6229 = 5 mod 8; (9113,3200),
(12313,5913) by classes mod 47, 151, 271; (14903,9696) by classes mod 89, 103, 151. (I count 5 ratios with a
survivor, 10 tuples; z16's "nine ghost ratios" is not a count I can reproduce from the note -- maybe a different
definition.) The degenerate controls (627,364), (991,263) are in neither.
Bearing on section 2: for magic squares the ghost tuples are irrelevant on this range -- they occur only at
ratios already excluded at tier A, and at the 71 ratios not excluded every residual tuple is pairing-killed (by
the session; by the repository's t10 at 14 (pilot) + 8 (my reruns) of them before I started a run on the other
49, below). So "why this matters for (R)" is about the four-form system (R) itself, not about any open case of
the magic square: the ghost observation neither threatens nor advances the census at m <= 15000.

### Task 2, final -- `z15_t10run.py` (driver around the repository's t10_pairct), `z15_compare.py` -> `z15_compare.txt`

Ratios rerun with t10_pairct (no hints; PARI finds every point anew): 174 distinct, 172 not in the pilot:
- all 80 of the fixed-seed sample (seed 20261127 from the 736 non-pilot ratios; `sample80_seed20261127.txt`);
- 43 ghost-shape ratios (the 84 at which a ghost class survives E_08/E_LU first, order seed 7; run stopped after 46);
- all 57 of the 71 ratios that are NOT excluded at tier A and that the pilot had not run (`notA_order.txt`, seed 71
  order) -- so with the pilot's 14, the repository's own pairing code has now decided every one of the 71;
- the 5 survivor ratios. (Overlaps were run twice: 20 ratios, identical output each time.)
Result: 1,330 residual tuples, 7,980 (tuple, curve) decisions: (session, t10) = 6,450 killed/killed, 1,530
not/not, **0 disagreements**; residual sets identical at all 174; Sel2 dimensions identical on every curve; no
'undecided', no 'outside Sel2'. Distinct decisions: 4,931 distinct (ratio, curve, class) [4,461 K/K, 470 N/N];
1,541 with the four mixed curves merged (one curve in the session's model) [1,245 K/K, 296 N/N]; a class never
gets two verdicts on the four isomorphic mixed models. 112 of the 174 ratios carry ghost-shape tuples (470 tuples).
The only tuples surviving all six curves are the note's 10. Least killing good prime of the 1,330 tuples
<= 157 (17: 218, 29: 550, 37: 270, ...). Also (A): my own parse of the pilot's files gives the note's 10,596 =
8,930 + 1,666 agreement, 0 disagreements, residual sets equal at 201 of 201.
Time: t10 CPU per ratio median 15 s, mean 28 s, max 426 s; 194 runs = 5,518 s CPU (wall 14,908 s on a 2-core
machine at load ~7). The session's cassels.py is ~100 times faster (its records say 0.1-0.3 s per ratio).
2-Selmer parity holds on all 194 x 6 t10 curves too.

## Bottom line

The note's computations stand: residual sets (937/937 equal to the repository's code), the 10 surviving tuples at
5 ratios (confirmed with t10_pairct), the rank bounds, every count in sections 1-4, the root-number tables. The
repository's own pairing (t10_pairct) agrees with the session's pairing on all 7,980 decisions at 174 more
ratios (172 outside the pilot), including every one of the 71 ratios not excluded at tier A. What does not stand
is the interpretation: (1) "a pair curve of rank 0 does not describe the exceptions" is false at (14903,9696)
(E_LU rank 0, torsion 8) and undecided at the partner pair (9113,3200)/(12313,5913); (2) the ghost shape of the
survivors is forced at every survivor ratio before the mixed curves act, so the evidence for "no other shape
survives" is weak (p ~ 0.07-0.2 at ratio level, 3 independent instances); "the classes of the degenerate solutions
themselves" and "reduce stage 2 to one explicit class per ratio" overstate it; at (3320,2303)/(5623,1017) the
survivors' classes are classes of rational points on both E_08 and E_LU (new, a PARI-found point checked);
(3) "Parity gives no handle on (R)" is too strong; the root number is the parity of dim Sel2. All five survivor
ratios are excluded at tier A (classes mod 47/89/103/151/271; primes = 5 mod 8), so nothing here touches an open
magic-square case at m <= 15000.
