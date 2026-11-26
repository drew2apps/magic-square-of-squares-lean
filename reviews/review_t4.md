# Hostile review of sess_t4 (SELMER_SQUARE.md, scripts/g1_full/selmer_square*.py)

Reviewer: separate AI agent (Claude), 2026-11-25. Scratch in /tmp/rev_t4/. Worktree not edited.
Diff reviewed: 8214be3..149c3cb (4 commits, 14 files, all new; no existing file changed).

## Log
- Read SELMER_SQUARE.md, reports/t4.md, UNIFORM_RANK.md sections 1-4 (model, Theorem 1, Lemma D and its remark).
- Independent re-implementation, written from the prose before reading selmer_square*.py: /tmp/rev_t4/indep.py
  (curve from the prose of UNIFORM_RANK.md, checked equal to sigma_exact.model(m, n, +-2) on 61 ratios; my own exact
  p-adic solubility search: generic subdivision for p < 1500 with Hensel/decided tests, root-branching for larger p,
  the two cross-checked on 3000 random quartics; real place by the closed max test; Selmer = signed squarefree d | b
  with class in every local image; my own implementation of the rule of Prop. S2 from its prose).
  Run on ALL 890 square-type ratios m <= 5000 (not only 40; 90 s): 1780 curves, 0 mismatches with
  selmer_square_5000.tsv (s, Sel(E), Sel(E') as sets, undecided = 0), 53784 local images = selmer_square_5000.json,
  rule of S2 = exact on 1780/1780, distribution 2:110 3:326 4:510 5:484 6:250 7:76 8:19 9:5 (= committed),
  2 in Sel(E) and -1 in Sel(E') on every curve; my product formula |im_p(E)||im_p(E')| = 4/8 and real 1*2 hold.
- Lemmas against my solver on random curves outside the family (/tmp/rev_t4/lemmas.py): L3 1500 random
  (a0 = 1 mod 8 of both signs, v2(kappa) in 4..12): 0 bad; L1 734, L2 2769 (hypotheses holding), L4 2000: 0 bad.
  L3 hypotheses are sharp: v2(kappa) = 3 gives im2(E) = {1, 10}, a0 = 5 mod 8 gives {1, 10}, a0 = 3 mod 8 gives {1};
  so 4 | mn (v2(kappa) = 2 + v2(mn) >= 4) is load-bearing; it holds at square type (even leg 2rs, r, s opposite parity).
- Scripts, run in an isolated copy (/tmp/rev_t4/wt, git archive of sess_t4 + git init), git status empty before/after
  every run without --write: selmer_square.py --negative and selmer_square_rule.py --negative: stdout = 
  selmer_square_negative.txt byte for byte (10/10 refused, exit 0); selmer_square_rule.py --max 5000 (2.6 min): stdout
  = selmer_square_rule_5000.txt; selmer_square_rank2.py: = selmer_square_rank2.txt, --negative 2/2; selmer_square_pari.py:
  = selmer_square_pari.txt except one extra first stdout line (cypari2 allocatemem banner) - harmless.
- Rank 2 side result re-checked with my own model: for E+ (40,9) and (144,17) the second hint point is on the curve and
  (x - e) times each of 1, HH-, 2H, 2H- is not a rational square; with my s = 2 this gives rank 2 (first hint = R).

## 1. Theorem S1, line by line (my own derivation, then compared with section 3)
- Groups: Sel(E) = classes d | b of y^2 = x(x^2+ax+b) (image of alpha, alpha(T) = b), Sel(E') = d | b'. Known images:
  alpha(T) = 4S^2 H H- = HH-, alpha(R) = 2SH = 2H (S = k^2), alpha(R+T) = 2H-, alpha'(T') = b' = -1. Correct. I
  re-checked independently: R = (-4mnS(m^2-n^2), -8mn(m+n)SH) is on the prose model and x(R) - e = 2SH (30 random
  (m, n)); b = 4S^2 HH-, b' = -(8mn(m-n)(m+n)^3)^2 on every curve.
- 2 in Sel(E), uniformly: real place (d = 2 > 0, e = 0); p = 2: 2SH/2 = k^2 H with H = 1 mod 8 (needs 4 | mn, true at
  square type) so 2 = alpha(R) in Q_2*/Q_2*^2 - this is the cleanest route and needs no L3; p in P_b: 2 is a unit and
  every unit class is in im_p (conic count, needs p >= 5: primes of k are 1 mod 4, of H, H- are 1 mod 12 since
  3 | mn gives u + v = m^2 - n^2 + 2mn != 0 mod 3); p in P_b': j = v_p(b') even, L1 says 2 fails only if (2/p) = -1 and
  (-2a/p) = +1; p | mn gives (-2a/p) = (2/p), p | m +- n gives (2/p) = 1 ((r +- s)^2 - 2s^2). Good primes: unit d.
  So S1 is proved for every square-type ratio, not only on the data. HOLDS.
- "Outside {1, HH-, 2H, 2H-}": 2 = 2H iff H square, 2 = 2H- iff H- square, HH- odd. Independence of 2, H, H- needs
  H, H-, HH- non-squares: H, H- by Lemma D through I9 (I re-verified I9 by sympy, and S H(m, n) = n^6 f(m/n), so the
  four points of y^2 = f, x in {0, inf}, exclude mn != 0); HH- by gcd(H, H-) = 1 (H - H- = -4mn(m^2-n^2),
  H(+-n, n) = 4n^4) - Lemma B is not needed. HOLDS.
- What S1 rests on: Lemma D (the only non-elementary input) and the point R (identity A5 of UNIFORM_RANK.md, needed
  for 2H in Sel(E)); plus textbook descent and elementary local facts. "Rests only on Lemma D" is accurate if read
  as "the only input beyond UNIFORM_RANK's identities and Silverman-Tate".
- L1, L2, L4, L0.1-L0.7 and the proof of S2 in section 3: re-derived, no gap found (L2's unit-class count, the
  j even/(A/p) = -1 exclusion by valuations 1 + 4mu, 2mu, j - 1; L1's conic t^2 - eps c^2 = delta1 with p + 1 points).

## 2. Lemma L3 (p = 2), re-derived
- C_{A,B}: x((x - 2a0)^2 + kappa^2); r = v(x) != 1 gives class 1 (r >= 2: 4((x/2 - a0)^2 + (kappa/2)^2) in 4(1+8Z2);
  r <= -1: v(kappa^2/X1^2) >= 10; r = 0: odd^2 + 0 mod 256); r = 1 needs v(X1) = v(kappa) = v, x1 = a0 mod 2^{v-1},
  class 2; existence via t = s0 + 2^{v-2} eta: eta(s0 + 2^{v-3} eta) covers all odd residues mod 8 for v = 4, 5 (shift
  by 2, 4) and v >= 6. Correct.
- C' (scaled by (4x, 8y), class-preserving): all nine sub-cases checked, including r + 1 > 2w forcing r >= 6 so that
  a0' = a0 + 2^{r-1} x1 = 1 mod 8, and v(delta) = r + 1 - 2w odd; the witness x = 2(2 - a0) has x + 2a0 = 4. Correct.
- Edge cases: a0 < 0 happens (e.g. (21,20)) and is fine (a0 = 1 mod 8 is a 2-adic square); v2(mn) = 2 vs >= 3 both
  covered (v2(kappa) = 4 and >= 5, tested); m +- n odd; v2(a) = 2, v2(b) = 2, v2(b') = 2(3 + v2(mn)). No change of
  model beyond the shift x = X - e (integral) and the internal (4x, 8y). All 1780 curves have exactly
  im2(E) = {1, 2}, im2(E') = {+-1, +-2} in my run. Hypotheses sharp (see log). HOLDS.

## 3-4. Re-implementation and scripts (details in the log above)
- Independent code on all 890 ratios / 1780 curves: Sel(E), Sel(E'), s and all 53784 local images equal to the
  committed TSV/JSON; my implementation of the S2 prose rule equals the exact groups on 1780/1780; the S3 example
  (k, H, H- prime: Sel(E) = <2, H, H-> iff k = 5 mod 8 or (k/p) = -1 for a split p) holds on all 8 such curves.
- Full regeneration in the isolated copy: selmer_square.py --max 5000 --jobs 2 --write (873 s, not 26 min):
  selmer_square_5000.tsv, .json, .txt byte-identical to the committed files; stdout = .txt. Real worktree:
  --negative of both scripts run with python3 -B, git status empty before and after, no file newer than a marker.

## 5. Prior art (short search)
- Legendre-symbol / F2-matrix descriptions of 2-isogeny Selmer groups of y^2 = x(x^2+ax+b) are classical: Monsky's
  matrix (appendix to Heath-Brown) and Aoki for congruent-number curves, T. Goto's explicit formulae for curves with a
  rational 2-torsion point (RIMS Kokyuroku 1319 and thesis), Kane's phi-Selmer matrices for j = 1728, the twist
  matrices of arXiv 2510.23291. Lower bounds on Selmer ranks forced in families by local conditions are also known
  (Klagsbrun arXiv 1201.5407 via Cassels' Tamagawa ratio; Lemmermeyer math/9903208, families with forced Sha classes).
  Zywina (arXiv 2502.01957) gets infinitely many rank-2 curves from a 2-descent with explicit Legendre conditions plus
  Tao-Ziegler (linear forms); here H has degree 4, so the doc's "Schinzel type" remark is the right caution.
- I found no source for this family or for S1 itself. Safe wording: "a classical mechanism, applied to a new family";
  neither "new phenomenon" nor "known result".

## Verdict
- FIX-FIRST (wording only): SELMER_SQUARE.md:3 says S1/S2 "are tier C of LEDGER.md", but the diff adds no LEDGER row
  (nor README/CHANGES); say "proposed for tier C" or add the row.
- NIT: SELMER_SQUARE.md:67 and t4.md "rests on Lemma D": also on the point R (UNIFORM_RANK A5) for 2H in Sel(E); the
  2-adic step is cleanest via alpha(R) = 2k^2 H = 2 (H = 1 mod 8), which makes L3 unnecessary for S1.
- NIT: selmer_square.py:240 prints "both curves of a ratio with bound 1", a weaker statistic than the claim; print
  "curves with bound 1: 0 of 1780".
- NIT: SELMER_SQUARE.md:71 "6000 random inputs": selmer_square.py:74 tests p in (2, 3, 5) but only p = 2 uses the new
  code, so 2000 effective comparisons.
- NIT: selmer_square_pari.py:18 allocatemem banner goes to stdout, so stdout != selmer_square_pari.txt by one line
  (do not diff it in CI as is). Run time is ~15 min on 2 cores here, not 26 (harmless).
- Optional addition: S1 gives unconditionally that if rank E = 1 then the class 2 is a non-trivial element of
  Sel \ alpha(E(Q)), i.e. a non-zero element of Sha[psi] (psi: E' -> E).
- HOLDS: S1 (uniform proof, both curves), S2, S3, L1-L4, L0, all data, negative tests, no writes without --write.
