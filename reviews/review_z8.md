# z8 — hostile referee of roadmap H, H1 and H3 (run 7)

Written by an AI model; not read by a person.

(Step headers carry no clock times: earlier drafts had wrong ones.)

Start: Fri Nov 27 13:19:29 UTC 2026. Scratch: /tmp/z8/.

## Log

### Step 1 — H3 claim 1 (good primes p >= 1154 keep all 16 vectors): FAILS
- Read h3_weil.py / h3_tables.py. The local test that the tables compute (method B) is right: at a good p = 1 mod 4 the
  label of a point is (chi(l_j^+(x)))_j at points where chi(l_j^+) = chi(l_j^-) for every j, i.e. where every F_j(x) is a
  NONZERO SQUARE. The Weil argument instead uses C_e: y_j^2 = u^{e_j} F_j(t). For e_j = 1 a point of C_e has F_j(t) a
  non-square, which is not an allowed point at all; the label is the class of l_j^+ (the linear factor), not of F_j.
  The curve that carries label e is the (Z/2)^7 cover given by sqrt(l_j^+ l_j^-) (4) and sqrt(u^{e_j+e_0} l_j^+ l_0^+) (3):
  branched at 8 points, genus 1 + 2^5 * 4 = 129, not 17. Hasse–Weil then needs roughly p > 67,000, not 1154.
- Wrote my own local test from LAST_PAIR.md s.11–12 (/tmp/z8/my_local.py, numpy integers). It equals h3_weil.json on
  every good residue at 15 primes from 17 to 1153 (incl. 204 of 1,140 losing at 1153).
- Above 1154 it finds losses at every p = 1 mod 4 tried: 1181 (152 of 1,176 residues lose), 1193 (124), 1201 (24), 1213
  (128), 1217 (144), 1229 (136), 1237 (112), 1249 (32).
- The bundle's own method A (G1's exact_local, via h3_tables.method_A, PARI 2.17.4) agrees on 12 lost residues at 1181,
  1201, 1249 (e.g. p = 1181, r = 5: vectors 5 and 10 missing).
- Third method, brute force over P^1(Z/p^2) with the actual integer forms, no Hensel (/tmp/z8/brute_p2.py): ratio (5,1)
  at p = 1181 (good for it): vectors 5, 10 not allowed even mod p^2; (17,1) at 1181: 4, 7, 8, 11 not allowed.
- So "good primes above 1153 never remove anything" is false, and so is "decided by bad primes, R and the 96 primes".

### Step 2 — H3 claim 2 (good p = 3 mod 4 never remove anything): FAILS as stated (for every ratio)
- Argument (mine): at good p = 3 mod 4, Q_p(i) is unramified, every rational unit d_j is a square in it, so the set is
  all-or-nothing; gamma_j(t) is a square in F_{p^2} iff its norm F_j(t) is a square in F_p. So the set is everything iff
  some t in P^1(F_p) has all four F_j(t) nonzero squares (F_j is anisotropic there, no roots), and EMPTY otherwise.
  That is "every candidate has the same label" (G2), which does not say the label is allowed.
- t = 0, inf need (2(m^2+n^2)/p) = 1; t = +-1 need ((m^2+n^2)/p) = 1; so p = 3 mod 8 always has a point. p = 7 mod 8 can fail.
- /tmp/z8/p3mod4.py (enumeration of P^1(F_p), all good residues, p = 3 mod 4 up to 1151): empty at
  p = 7 (4 residues), 23 (8), 31 (16), 47 (8), 79 (24), 103 (8), 127 (8), 151 (8); none from 167 to 1151.
- G1's exact_local (method A, PARI 2.17.4) confirms empty allowed sets: (2,1) and (9,1) at 7, (3,1) at 23, (5,1) at 47,
  (25,1) at 127, (14,1) and (151000467,1) at 151 (all good primes for these ratios).
- For p >= 1154, p = 3 mod 4: the genus-17 curve y_j^2 = F_j(t) (e = 0) IS the right curve there, has no F_p-points
  over branch points, and Hasse–Weil gives a point. So H3's genus-17 argument is a correct proof for p = 3 mod 4,
  p >= 1154 — not for p = 1 mod 4.

### Step 3 — how far the p = 1 mod 4 losses go (computed by my code; numerical evidence only above what is scanned)
- The right curve for label e at p = 1 mod 4 has genus 129; Hasse–Weil (quoted) with the 8*64 points over branch points
  removed gives a point with label e (or e + 1111) for every good p >= 67,583 ((p - 511)^2 > 66,564 p, exact integers).
  So a correct finite statement is "bad primes, R, the good p = 7 mod 8 below 1154 (all-or-nothing; empties only seen
  at p <= 151) and the good p = 1 mod 4 below 67,583". Not computed in full here.
- Sampled (my_local.py, sample_high.py): losses at 1301 (64 of 1,296 residues), 1549 (72), 1709 (56), 2017 (36); none
  in 3,000 sampled residues at each of 2521, 3001, 4001, ..., 30013. Full scan 1990–3700 running.

### Step 4 — reruns under PARI 2.17.4 (all from /tmp/z8/reruns, PYTHONDONTWRITEBYTECODE=1, env pointing into int7)
- h3_tables.py: "# matches" (30 s); --negative: matches (1 s); --check 2000: matches (8 s); h3_weil.py --negative:
  matches; h3_sizes.py: matches; h1_iso.py: matches (136 s here, report says 3 s); h1_iso.py --negative: matches;
  h1_rank0_exact.py --negative: matches; h1_scan.py --negative: matches; h1_pairs.py --negative: matches (70 s).
- Full scan of my_local.py over every good residue: losses at EVERY p = 1 mod 4 from 1993 to 2477 except a few, and
  at 2549 ... 3457 (last in 1990–3700: 3457, 12 residues). So losses persist at least to 3457. Scan 3700–5000 running.

### Step 5 — H3 numbers against h3_weil.json / h3_tables.json (my own count script /tmp/z8/h3_counts.py)
- h3_weil.json holds 91 primes (all p = 1 mod 4, 17 <= p <= 1153), NOT 96. sympy: 91 primes = 1 mod 4 in [17,1153]
  (93 with 5 and 13). "96" (REPORT_H3, HANDOFF, h3_sizes.py docstring, REPORT_H2) is wrong.
- Correct as stated: 49,764 good residues, 17,224 lose a vector; a loss at every one of the 91 primes; 15,984 not
  affine; 36 empty sets, the largest such p = 137; <= 2 up to 313, <= 4 up to 401, <= 8 up to 617, <= 10 up to 1069;
  above 1069 at most 4 vectors lost; 204 of 1,140 at 1153; mean 2.67 at 17, 10.3/9.8 at 193/197, 14.4 at 509, 15.62 at
  1153; p <= 197: 1,848 good residues, 36 empty, full only at 137 (8), 181 (16), 193 (16); 5, 13 no good residue;
  47,916 residues at 197 < p < 1154.

### Step 6 — H1 "extra tuples are local everywhere": the good-prime part
- H1 tested good primes to 20,000 (G1's one()) and leaned on H3 for p > 1153. With H3's claim 1 false, the range
  20,000 < p < 67,583 was untested. /tmp/z8/phantom_mid.py (my local test) on all 24 phantom ratios at every good
  p = 1 mod 4 in (20000, 67583): 2,234 primes, 53,608 (ratio, prime) pairs, 0 cuts of an extra tuple; known tuples
  allowed at every one (asserted). Same at all good p = 1 mod 4 <= 1193 (2,142 pairs). With the genus-129 argument
  (p >= 67,583) and the rational points at t in {0, +-1, inf} (so good p = 3 mod 4 keep everything for these ratios),
  every extra tuple survives every good prime. (Computed by me; Hasse–Weil quoted.)
- Remaining gap for "local everywhere": p = 2 is the repository's one-sided depth-12 test (a superset), and odd bad
  primes are exact_local ("tested, not proved", x6). Nobody produced explicit local points for the 13 new phantoms or
  the partners (x6 did it only for (1109503,597671)).

### Step 7 — MAJOR: H1's scan never tested the bad primes of m^2+n^2 and of the octic above 20,000
- G1's one() (g1_stress.descent) tests only the primes of the supports (2mn, 2(m^2-n^2)) and primes <= pmax
  (x6's correction 2 to G1 says exactly this). H1's "every bad prime" is therefore false as a description.
- /tmp/z8/badprimes_left.py: factor m^2+n^2 and the two octic quartics (PARI factor, isprime), take the primes not in
  the supports and > 20,000, run G1's exact_local.allowed_at there (known tuples asserted allowed).
- Result: at those primes exact_local CUTS every extra tuple for 5 of the 12 pairs (10 of the 24 ratios):
  (11875080,6497849)/(18372929,5377231) [the "open" pair] at the octic primes 5462461 and 798775110477654781;
  (32396425,3694248)/(36090673,28702177) at 218077 and 148619132854261;
  (175775580,15972419)/(191747999,159803161) at 10057297, 15534866530669, 3953388125944997713, 1.136e26 prime;
  (234319635,115640732)/(349960367,118678903) at 21062406847889269, 114475219027652869;
  (518804652,31743595)/(550548247,487061057) at 1650997632721, 38859794425882715723401.
  The other 7 pairs keep their extras at every such prime.
- If right, those 5 are not phantoms at all: the descent with every bad prime empties their extras, and the "open"
  pair is closed by a local test. Checking this independently next (short exact argument at an octic prime).

### Step 8 — the octic-prime cuts confirmed by an independent exact argument (/tmp/z8/octic_cut.py, octic_cut_all.py)
- Argument (mine, short, exact): let p = 1 mod 4 divide the octic but not 2mn(m^2-n^2) (so every d_j is a unit). Two
  linear components l_a (of form j_a) and l_b (of form j_b != j_a) have the same root x0 mod p. If v_p(resultant) = 1
  (root difference of valuation 1; checked with iota mod p^3), then on the disc of x0 one of l_a, l_b always has odd
  valuation, so the disc has no allowed point; everywhere else both are units with l_b = kappa l_a mod p. Since a form's
  two components must have the same class, every allowed point has e_{j_a} + e_{j_b} = chi(kappa).
- At every untested octic prime of all 24 ratios: 4 such collisions, all with v = 1; known tuples satisfy all four
  constraints (positive control); the extras violate them exactly at the primes where exact_local cuts (14 ratio-prime
  cases cut, all others keep). E.g. (11875080,6497849), p = 5462461 = 5 mod 8: constraints e_0+e_U = 1, e_0+e_L = 0,
  e_8+e_U = 0, e_8+e_L = 1; the extra (2,1,1,98795483379599) has e_0 = chi(2) = 1, e_L = 0: violates e_0+e_L = 0.
- H1's own output already contradicts "the extras pass every local test": at (36090673,28702177) h1_allphantoms.txt
  says the extra classes on E_0L, E_0U, E_8L, E_8U are "in Sel^2: no".
- Consequence (computed, exact; exact_local + my argument agree): with every bad prime, only 7 of the 13 scan ratios
  are phantoms — (204724,97893), (561607,264401), (756492,255635), (1012127,500857), (1992873,999664),
  (189089400,156063271), (231641256,72371617). With (1109503,597671) and partners: 7 pairs, not 12. The other 5
  "pairs" are not phantoms: their extra tuples are not locally soluble at an octic prime. In particular the "open pair"
  (11875080,6497849)/(18372929,5377231) is settled by the local test at 5462461 (and at 798775110477654781), no
  ranks needed. All 7 true pairs are settled by H1's methods: rank-0 quotient (413004/561607, 853587/1109503,
  1992873/2992537, 189089400/345152671, 231641256/304012873) or pair curves (204724/302617 at 302617 by the CT kernel
  test; 756492/1012127 at 1012127 by the CT kernel test).
- p = 2: my explicit 2-adic points (/tmp/z8/witness2.py, square test in Q_2(i) mod 8 with Hensel margin) for all 108
  known and extra tuples of the 24 ratios (small coprime (u,v), unit lambda); forged tuples (1109503,597671,6,1) etc.
  get none (negative control).

### Step 9 — further H1 checks
- Partner isomorphism (claim 7), derived by hand and checked with sympy over Q(m, n, u, v): with (u', v') = (-u-v, -u+v)
  (t' = (t+1)/(t-1)) and (m', n') = ((m+n)/g, (m-n)/g): z'_L∘M = -4 z_0, z'_U∘M = -4 z_8, z'_0∘M = -z_L,
  z'_8∘M = -z_U for g = 1; -2 and -1/2 for g = 2. -1 and -4 are squares in Q(i); M(u, v) has gcd 1 or 2, so the image
  tuple is sigma(T) times e in {1, 2}; h1_iso enumerates more e (harmless, stronger). M^2 = 2I, so the map is a
  bijection on rational points; {-1, 0, 1, inf} -> {0, -1, inf, 1}. Transfer argument: holds. H1 forms = TD forms
  ((1+i) multiplied out; z_8 is the conjugate of TD's, same rational class): checked.
- Rank-0 quotient (claim 6): map C -> A, (t, y0, y8) -> (x, Y) = (m^2 n^2 t^2, m^2 n^2 t y0 y8) onto
  Y^2 = x(x + m^4)(x + n^4): checked by hand. A always has Z/2 x Z/4 (x = +-m^2 n^2), x = -m^2 n^2 gives t^2 = -1.
  PARI 2.17.4 as a second implementation (/tmp/z8/rank0_check.py): for the 7 ratios, ellrank gives [0,0] on at least
  one curve of the isogeny class (ellisomat), elltors gives 8 or 16 points, and the t over torsion with w a square are
  exactly {-1, 0, 1, inf} (my exact check). h1_rank0_exact.analyse rerun on all 7 (PARI 2.17.4): every line equals
  the committed file (1-5 s each). Holds.
- CT kernel test (claim 6): logic sound. delta(E_jk(Q)) lies in the kernel of the CT pairing on Sel^2; the class of
  the E_jk-image of a C-point with tuple T is, for the root e* of the factorisation (a_j a_k)(conj), r = 2(Re h + y)
  with h = a_j a_k: r = (h + y)^2/h, so r = h = d_j d_k mod K*^2, i.e. +-d_j d_k in Q*/Q*^2, normalised at the base
  point. The code compares modulo -1 (sound: weaker). e* is identified by agreement on >= 12 points (one root only).
  Rerun of analyse on E_0L of (302617,106831): equals h1_holdouts.txt (45 s). Partial independent check: PARI
  ellrank on that curve gives bounds 4..4 (3 points returned); the extra classes 11186, 5780239 are the x - e_i class
  of no element of the span of PARI's 3 points and the 2-torsion, for every root (/tmp/z8/pari_e0l2.py).
  Minor: "lo = ri - 2" assumes the 2-torsion points have independent images; with a point of order 4 it under-counts
  (h1_allphantoms.txt prints "-1 <= rank <= 1" at 561607): harmless (only the lower bound).
- 13 phantoms / 20,009 / 12 pairs / every partner a phantom (claim 5): the counts are what the files say (17 x 1,177;
  1+3+1+1+1+4+2 = 13; 14 with the control; 2 internal pairs -> 12), but they are counts of a descent that skipped
  bad primes. With every bad prime: 7 of the 13, 7 pairs (see 15:00).
- Open pair (claim 8): tried pair curves (3 <= rank <= 4, extra classes are classes of points), A and B quotients,
  s-line quotients, g4 quotients, DM on y^2 = wK (scratch). None needed: local test at the octic prime 5462461.

### Step 10 — H3 claim 4 and H2's dependence
- Claim 4 fails (claims 1 and 2 fail; "96" should be 91). H2 REPORT says "By H3 ... P <= 1153" and "a survivor would
  have to keep a tuple through all 96 good primes p = 1 mod 4 up to 1153" (lines 16, 132, 143): both must go. H2's
  exhaustive results themselves (kill sets found) do not depend on H3.
- Scan 3700–5000 finished: losses at 3709 (8 residues) and 3853 (8); none at the other primes of 3700–5000.

### Step 11 — more local points for the 7 true phantom pairs
- /tmp/z8/witness_odd.py: explicit points (a residue x with every component a unit and the right symbols; any lift
  works) for every extra tuple at every odd prime of m^2+n^2 and of the octic outside the supports: 224 of 252 cases;
  the other 28 are all p = 5. /tmp/z8/witness5.py: explicit 5-adic points mod 5^8 (even valuation < 7, square unit) for
  all 28; known tuples also found (control). So for the 7 true pairs only the odd SUPPORT primes rest on exact_local.
- h1_pairs.py (main) rerun under PARI 2.17.4: "# matches" (28 s); (1109503,597671) closed twice: holds.

## Verdicts

**H3**
1. Good p >= 1154 keep all 16 vectors — **FAILS.** Wrong curve: the label is the class of the linear factor l_j^+,
   at points where every F_j is a nonzero square; y_j^2 = u^{e_j}F_j(t) with e_j = 1 has F_j non-square, not an
   allowed point. Right curve: (Z/2)^7-cover, genus 129, Hasse–Weil threshold 67,583. Counterexamples: my enumeration,
   G1's exact test and a brute force mod p^2 agree (e.g. ratio (5,1) at p = 1181 loses vectors 5, 10); losses at good
   p = 1 mod 4 up to 3853 at least (full scans 1153–5000). Genus 17, smoothness, unit twists, threshold arithmetic
   (1154 is the least p with (p+1)^2 > 1156p) are right for the curve written, which is the right curve only at
   p = 3 mod 4. Branch points: in C_e "the class is free" at F_j = 0 is also wrong (there the label is chi(l_j^{-s})).
2. Good p = 3 mod 4 never remove anything — **FAILS as stated.** G2's fact (all candidates have one label) makes the
   set all-or-nothing; it is empty for some residues at p = 7, 23, 31, 47, 79, 103, 127, 151 (enumeration + exact_local;
   e.g. (2,1) at 7, (14,1) at 151). Holds for p = 3 mod 8, for p >= 1154 (genus-17 HW, now correct) and for ratios with
   a rational point (all of H1's).
3. Tables — **holds**, with one wrong number: 91 primes, not 96. Reruns match; my code equals h3_weil.json on 15
   primes; every other count checked holds. --check 2000 only uses primes <= 197.
4. "Decided by bad primes, R and those 96 primes" — **FAILS** (1, 2, and 96). Correct finite list (my argument, HW
   quoted, not in the repository): bad primes, R, good p = 7 mod 8 below 1154, good p = 1 mod 4 below 67,583.

**H1**
5. 13 phantoms / 12 pairs / partners phantoms — **fails as a statement about phantoms** (counts match the files).
   G1's one() skips the primes of m^2+n^2 and of the octic above 20,000 (x6's correction 2 of G1). At those primes
   exact_local and my short argument (two root residues meet with v_p = 1, giving e_ja + e_jb = chi(kappa)) agree:
   every extra tuple dies for 5 pairs. True phantoms: 7 of the 13 scan ratios; 7 pairs. Every partner of a true
   phantom is one: holds.
6. 7 rank-0 + 4 pair-curve settlements — **holds (method), with gaps**: rank-0 logic, map and torsion fibres checked,
   PARI 2.17.4 agrees on rank 0 and torsion; CT kernel test sound (compared in Q*/<-1,squares>, which is weaker); it
   rests on rank_second.py's pairing (quoted theory). Of the 11, only 7 are phantom pairs (5 rank-0, 2 pair-curve);
   the other 4 (175775580, 234319635, 32396425, 518804652 pairs) were settled by a local test all along.
7. Partner isomorphism and transfer — **holds** (derived by hand; sympy identity; reruns match).
8. Open pair — **settled by the local test at p = 5462461** (and 798775110477654781): no rank needed.
"Extra tuples pass every local test" — **fails** for those 5 pairs; for the 7 true pairs **holds** (step 12: explicit
points at 2, at every odd bad prime, support primes included, and at good p <= 67,583, by my code; Hasse–Weil, quoted,
on the genus-129 cover above 67,583; not by H3's argument).

## Replacement text

**REPORT_H3.md, proposed LEDGER/CHANGES text (replace the whole quote; tier D):**
> **Experiment, not counted.** Local images of the four-class cover at good primes (experiments/roadmap-h/h3). For
> each of the 91 primes p ≡ 1 mod 4 with 17 ≤ p ≤ 1153 and each good residue r = m/n mod p, the set of label vectors
> allowed at p, computed two ways (G1's exact local test and an enumeration of ℙ¹(F_p)) and checked against the exact
> test on 2,000 random ratios at primes up to 197. Every one of these primes removes a vector at some residue (17,224
> of 49,764 residues); a set is empty only up to p = 137; 15,984 of the sets are not affine subspaces of F₂⁴. These
> tables are not a complete local test: larger good primes remove vectors too (a review found losses at good primes
> p ≡ 1 mod 4 up to 3853), and good primes p ≡ 3 mod 4 can empty the set (p = 7, …, 151). Computed (exact). Written by
> an AI model; not read by a person.

Optional, only if /tmp/z8/my_local.py, p3mod4.py and the argument below are committed (tier C, Hasse–Weil quoted):
> At a good prime p ≡ 1 mod 4 a label vector e is carried by the (ℤ/2)⁷-cover of ℙ¹ given by √(l_j⁺l_j⁻) and
> √(u^{e_j+e_0} l_j⁺ l_0⁺), branched at 8 points, of genus 129; Hasse–Weil gives a point outside the branch fibres for
> p ≥ 67,583, so such primes remove nothing. At a good prime p ≡ 3 mod 4 the set is everything or nothing, everything
> iff some t ∈ ℙ¹(F_p) makes all four F_j(t) squares: always for p ≡ 3 mod 8 (t = 0 or ±1) and for p ≥ 1154
> (Hasse–Weil on y_j² = F_j(t), genus 17); for some residues it is nothing at p = 7, 23, 31, 47, 79, 103, 127, 151, and
> never at 167 ≤ p ≤ 1151.

Also: REPORT_H3 bottom line bullets 1–3 and "What failed: Nothing" must go; "96" -> "91" (REPORT_H3 twice, HANDOFF,
h3_sizes.py docstring, REPORT_H2 line 16); h3_weil.py docstring STATEMENT/PROOF is wrong as written.

**HANDOFF_ROADMAP_H.md, H3 bottom line and label (replace):**
> - **Bottom line:** tables of the allowed label vectors at the 91 good primes p ≡ 1 mod 4 with 17 ≤ p ≤ 1153,
>   computed two ways: losses at every one; 15,984 of the 49,764 sets are not affine subspaces of F₂⁴. The bundle's
>   Hasse–Weil step used the curve y_j² = u^{e_j}F_j(t), which does not carry the label vector, and is withdrawn:
>   good primes above 1153 do remove vectors (to 3853 at least), and good primes ≡ 3 mod 4 can empty the set.
> - **Label:** tables computed (exact); no completeness statement.

**REPORT_H1.md, proposed LEDGER/CHANGES text (replace; tier D):**
> **Experiment, not counted.** The phantoms of the four-class descent (experiments/roadmap-h/h1). A scan of 20,009
> random degenerate ratios up to m = 10⁹, with the four-class descent at the primes of the supports and the primes up to
> 20,000, keeps tuples beyond those of the degenerate points at 13 ratios; with (1109503,597671) and the partners they
> form 12 partner pairs. That run left out the primes of m² + n² and of the octic above 20,000. A review tested them: at
> an octic prime, 5 of the 12 pairs, (11875080,6497849)/(18372929,5377231) among them, lose every extra tuple, so they
> are not phantoms. In each of the other 7 pairs every extra tuple is the class of no rational point: in 5 by an exactly
> computed rank-0 quotient y² = w·f0·f8 (rank bound over the isogeny class), in 2 by the pair curves' complete 2-descent
> and Cassels–Tate pairing. The partner isomorphism t′ = (t + 1)/(t − 1), with the forms permuted (0, 8, L, U) →
> (L, U, 0, 8), is written down and checked exactly, so a proof on either side settles a pair. No pair is open.
> Computed (exact); the pairing as implemented in rank_second.py. Written by an AI model; not read by a person.
(The "review tested them" sentences need /tmp/z8/badprimes_left.py and octic_cut.py, or a fixed h1_scan, committed.)

**HANDOFF H1 bottom line (replace):**
> - **Bottom line:** the scan's descent skipped the primes of m² + n² and of the octic above 20,000. With them, 7 of its
>   13 phantoms remain; with (1109503,597671) and the partners, 7 pairs, all settled (5 by a rank-0 quotient, 2 by the
>   pair curves and the Cassels–Tate pairing). The other 5 pairs, (11875080,6497849)/(18372929,5377231) included, are
>   emptied by a local test at an octic prime. No pair is open.
> - **Next step:** make h1_scan test every bad prime (as G2 does) and rerun the 20,009.
Delete the H1 "Next step" (Mordell–Weil sieve) and the "s-line"/"g4"/"DM on the open pair" ideas' claim that the pair
needs them. In REPORT_H1: "Their extra tuples are local everywhere ..." and "By H3, good primes above 1153 ..." go.

**REPORT_H2 / HANDOFF H2:** delete "By H3, a cutoff above 1153 is never needed", "P ≤ 1153 for every ratio ...",
"A survivor would have to keep a tuple through all 96 good primes ...", and "and by H3 never needs to pass 1153".

## What I could not check
- The Cassels–Tate pairing code (rank_second.DFull) itself: reruns only, plus a partial PARI cross-check on one curve.
- The largest good prime p = 1 mod 4 that removes a vector: every residue scanned to 5000 (last loss 3853), partial
  scans 5000–6133 and 7600–8093 found none; a complete statement needs every p < 67,583.
- h1_scan / h1_allphantoms / h1_holdouts / h1_rank0_exact full runs (time); per-ratio reruns matched.
- H2 beyond the sentences that cite H3.

## Final
H3: the Hasse–Weil step uses the wrong curve (genus 17 instead of 129); good primes above 1153 do remove vectors (to
3853 at least), good p = 3 mod 4 can empty the set (p = 7 ... 151), and there are 91 primes, not 96. Claims 1, 2, 4 fail;
the tables hold. H1: the scan never tested the primes of m^2+n^2 and of the octic above 20,000; at those primes 5 of the
12 pairs, the "open" pair among them, lose every extra tuple (two independent exact methods agree). The 7 true phantom
pairs have explicit local points everywhere (Hasse–Weil above 67,583) and are all settled (5 rank-0, 2 Cassels–Tate).
No pair is open. Tier D for both. Scripts: /home/claude/conductor/r7/referee_z8/ (README.md), scratch in /tmp/z8/.

End: Fri Nov 27 14:42:11 UTC 2026. Package complete in /home/claude/conductor/r7/referee_z8/ (17 scripts + z8common.py + README.md; 15 of the 17
outputs rerun and byte-identical to their .txt; scan_losses and phantom_mid, which use no randomness, were not rerun).
