# t11: Fisher's higher descents on the nine I5 curves

Written by an AI model (build agent t11), not read by a person. Experiment, not counted. Nothing here changes the census. Not reviewed by a separate agent; the pre-release maths referee (panel P1, `reviews/panel_p1.md`) read it: holds, with gaps (Fisher's paper read through summaries; the bounds hold up to level 4 only: a level-5 pairing could still lower them). Nothing here changes any row's status.

## Bottom line

- **None of the 9 curves reaches the bound 1. No row becomes unconditional.**
- Two bounds go down from 5 to 3 (tier C):
  - E1(8249,3720): through `Θ₂`, which has rank 1.
  - Sig+(22057,2351): through `Θ₃`, which has rank 2.
- On three curves the pairing that could have closed the curve was computed and is 0 (tier C):
  - E1(10512,145)−: `Θ'₃ = 0`.
  - E1(8249,3720): `Θ₄(113, 755496193) = 0`.
  - E1(4440,2231): `Θ'₃(2, ·) = 0`, so the class 2 is in the kernel and no bound below 2 is possible.

  On these curves Fisher's method, as far as Section 7 takes it, cannot give 1. The open classes survive to level 4.
  If the rank is 1, that means Sha[2^∞] has elements of order 4 or more. With conductors of 49 to 60 digits, that
  is plausible.
- Not finished:
  - D1(12007,1799) failed: no conic point was found over a real quadratic field with a large regulator.
  - E1(3540,2581) was not attempted: every field it needs has a regulator of 396790 or more.
  - E1(6480,4961), E1(10512,145)+ and Sig+(12204,9853) (type B, `Θ₄`) are unfinished: the time ran out while the
    conic finders searched over the real or large-class-group fields. See "What failed".

## The structure of the nine curves (rank_indep, exact)

Notation follows Fisher. rank_indep's "E side" (classes `d | b`) is `S'_1`, and its "E′ side" (classes `d' | b'`)
is `S_1`. "Open" means not excluded by the second descent and not the image of a known point. Every curve has
`a² − 4b = −□`, so `ℚ(E[2]) = ℚ(i)`, as the `Θ₄` code assumes.

| curve | E side: dim S'_1 / S'_2 / known | E′ side: dim S_1 / S_2 / known | type | exact bound (I5) |
|---|---|---|---|---|
| E1(3540,2581) | 3 / 3 / 1 | 2 / 2 / 2 | A: 2 open on E | 3 |
| E1(4440,2231) | 5 / 5 / 1 | 2 / 2 / 2 | A′: 4 open on E | 5 |
| E1(6480,4961) | 4 / 2 / 1 | 3 / 3 / 2 | B: 1 + 1 | 3 |
| E1(8249,3720) | 3 / 3 / 1 | 6 / 4 / 2 | C: 2 + 2 | 5 |
| E1(10512,145)+ | 2 / 2 / 1 | 5 / 3 / 2 | B | 3 |
| E1(10512,145)− | 3 / 3 / 1 | 4 / 2 / 2 | A | 3 |
| D1(12007,1799) | 3 / 3 / 1 | 2 / 2 / 2 | A | 3 |
| Sig+(12204,9853) | 3 / 3 / 2 | 4 / 2 / 1 | B | 3 |
| Sig+(22057,2351) | 2 / 2 / 2 | 5 / 5 / 1 | A′: 4 open on E′ | 5 |

What each type needs:

- **Type B** (one open class on each side, as on the five curves of rank_fisher.py) needs `Θ₄`, which is what the
  repository's code computes.
- **Type A** has all its open classes on one side. Then `Θ₂ = 0`, because its other argument ranges over the images
  of the rational points. What is needed is `Θ₃` (or `Θ'₃`), the alternating pairing of that side.
- **Type C** needs `Θ₂` first.

So for five of the nine curves `Θ₄` is not the pairing that matters. `Θ₃` is computed by the same Section 3 formula
on the level-3 curve `C3'`, with one Lemma 7.1 pushout:

`Θ₃(y, z) = Σ_v (f(P_v), z)_v`, for `z` of the same side as `y`.

## Per curve

| curve | pairing computed | matrix rank | bound now | where | status |
|---|---|---|---|---|---|
| E1(10512,145)− | `Θ'₃(2, 11878571404153441) = 0` (two lifts: 2 and 2·b0) | 0 | 3 | E side, pushouts over ℚ(√2) | computed; does not close |
| E1(8249,3720) | `Θ₂` on {13278577, 755496193} × {113, 2991753233} = [[0,1],[0,0]], both ways; then `Θ₄(113, 755496193) = 0` | 1, then 0 | 3 (from 5) | `Θ₂` rational; `Θ₄` on E′, ℚ(√±113) | 5 → 3; does not close |
| Sig+(22057,2351) | `Θ₃(14, ·) = (0,0,0,1,0)` and `Θ₃(274, ·) = 0` on (−1, 14, 274, 339, 108146) | 2 | 3 (from 5) | E′ side, ℚ(√−14), ℚ(√−274) | 5 → 3; does not close |
| E1(4440,2231) | `Θ'₃(2, ·) = 0` on a basis of S'_3 | ≤ 2 | 5 (3 at best) | E side, ℚ(√2) | does not close |
| D1(12007,1799) | lift of 354384937 built; no conic point over ℚ(√354384937) | – | 3 | E side | failed: finder |
| E1(3540,2581) | – | – | 3 | E side | not attempted: fields too large |
| E1(6480,4961) | E′ lift of 34129302 built | – | 3 | E′ side | not attempted: time |
| E1(10512,145)+ | E′ lift of 23868431349 built (λ = −13, ξ = 730); finder over ℚ(√−23868431349) (class group [4064,4,2,2,2]) did not finish | – | 3 | E′ side | unfinished: time |
| Sig+(12204,9853) | E′ lift of −959 (λ = −1); conic point over ℚ(√959) found in 13 s; the one over ℚ(√−959) did not finish | – | 3 | E′ side | unfinished: time |

How the bounds follow from the pairings:

- **E1(8249,3720).** After `Θ₂`, `S_3 = <known, 113>` (dim 3) and `S'_3 = <known, 755496193>` (dim 2). Each has one
  open class, so `Θ₃ = Θ'₃ = 0`, and Lemma 2.2 with n = 2 gives 3 + 2 − 2 = 3.
  - For `Θ₄`, the first level-3 twist ξ = 7054917115 had a non-zero level-4 obstruction. It is non-zero only on the
    classes of `S_1` that the second descent excludes: on (−1, 113, 1533, 9058, 11969, 380058) it is
    (+,+,−,+,+,−), and it is 0 on `S_2`, as Theorem 2.1 requires.
  - The twist ξ·11969 removes it.
- **Sig+(22057,2351).** 274 is in the kernel of `Θ₃`, and `Θ₃(14, 339) ≠ 0`. So `Θ₃` has rank exactly 2 on the
  4-dimensional open part, and dim `S_4` = 3. Lemma 2.2 with the primes swapped gives
  dim `S'_3` + dim `S_4` − 2 = 2 + 3 − 2 = 3.
- **E1(4440,2231).** The basis is b0, 2, 970682689, 6846322457 and 119601394217; the fifth class paired,
  1620825056929, is a product of three of them. The other rows need ℚ(√1941365378) (regulator about 895) or larger
  fields.

## What is checked exactly, and what is a written argument

The same split as the repository's `Θ₄` (LEDGER, tier C).

**Exact**, in rationals and integers, with PARI's p-adics at precision asserted, as in rank_fisher_theta4:

- the models, checked against I5's minimal models;
- rank_indep's descents;
- the level-2 curves `D_λ`, with local points at every place of S;
- the level-3 curves `C3`, checked locally soluble at R and every prime of S3 (rank_fisher.certify_C3, and
  t11_lib.recertify_xi for a twist);
- rank_fisher5's steps (1)–(8), with every conic point checked on its cone and l tangent along the ruling;
- the local terms. The term at v must not depend on the sampled point `P_v`, and that is checked. The negative test
  shows the check has teeth: an E′ class paired through an E-side `C3` gives a term that depends on the point.
- The controls: `Θ₃(y, y) = 0`; the pairing with the image of T is 0; symmetry (`Θ₃(14, 274) = Θ₃(274, 14)`); two
  lifts of the same class agree; `Θ₂` from either side agrees; and K2, K3, K5 of the `Θ₄` code.

**Written arguments, not machine-checked:**

- Outside S, the term is 1. `C3'` has good reduction, f0 mod p is a non-zero form on it, and a point off its at most
  8 zeros lifts for p ≥ 47.
- The `Θ₄` code's own written arguments (the evenness of div g0, the good-reduction step).
- The reading of Fisher's Theorem 2.1, Lemma 2.2 and Section 3, as in the repository. This note did not read the
  paper itself.
- The step from "the level-3 lift of a type-A class" to `Θ₃` uses the same side for `z`. That is the reading in
  rank_fisher.py's docstring (6) and rank_fisher_theta4_norm.py's obstruction ("Ob is 1 on S_2 … moves by
  Θ₁(s, ·)"). The point-independence check agrees with it on every curve here.

## Code (nothing in the repository was edited)

- **t11_lib.py** holds:
  - the models;
  - `level3_lift`: rank_fisher.analyse's loop, made general for any number of open classes;
  - `recertify_xi`: another member of the level-3 Selmer set, with checked local points;
  - `find_conic`: rank_fisher_theta4_find.find's body, taking R in place of a json index (a finder);
  - `find_conic_norm`: rank_fisher_theta4_norm.find_hint's body, likewise (a finder);
  - `theta3`: new;
  - `theta4_gen_module`.
- **rank_fisher_theta4_gen.py is used unchanged.** Module attributes are replaced at run time, and no file is
  edited:
  - `curve_index` and `base_curve` now return this curve's record in the rank_fisher_curves.json format;
  - `HINTS_GEN` points to t11_hints.json;
  - `DIMS` gets the dimensions of this curve;
  - **`choose_Xi_lin` is replaced** by t11_lib.choose_Xi_lin_images. It does the same F₂ linear algebra, but first
    closes the sampled classes of f0 at each place with the local image of the twisting side: rank_indep's certified
    image at bad primes, and the unit classes at good odd primes. The original closes them only to the coset they
    span. That was too small at the pencil prime 1332819889 of E1(8249,3720), where the sampler finds one class of
    a coset of two. This is a sampler: C4 is then checked locally soluble at every place, as before.
- **t11_theta.py** is the checked run, with hints from t11_hints.json; it writes t11_theta_<curve>.txt.
  **t11_negative.py** (`--negative`) runs 9 forged and true inputs and writes t11_theta_negative.txt.
  - A re-run prints `# matches`. This was checked for E1(10512,145)−, D1(12007,1799) and `--negative`.
  - times.txt holds the wall times. run_all.sh runs a list of curves with `--write`.
- **explore/** holds the finders and the exploration scripts. Run them from this folder, for example
  `python3 explore/try_find.py explore/l3_....json NAME "E'" DD U`.
  - Finders, which wrote t11_hints.json: try_l3.py (level-3 lifts), try_find.py and try_find2.py (conic points).
  - Exploration: try_t3.py, try_t3b.py, try_t4.py, try_ob.py, try_t2.py, try_xi.py, try_rexi.py, try_reg.py,
    try_sizes.py and try_probe.py.
  - The l3_*.json files are the level-3 records these scripts used. explore_*.py and explore_rank_all.json are the
    first survey.
  - Nothing rests on them: t11_theta.py recomputes every lift and re-checks every hint.

## Proposed CI lines (times on a loaded 2-core machine)

```
python3 -I experiments/conductor-t11/t11_theta.py "E1(10512,145)-"   # 70-140 s  (Theta'_3 = 0, two lifts)
python3 -I experiments/conductor-t11/t11_theta.py "Sig+(22057,2351)" # about 2 min (Theta_3 rank 2: 5 -> 3)
python3 -I experiments/conductor-t11/t11_theta.py "E1(4440,2231)"    # about 70 s
python3 -I experiments/conductor-t11/t11_theta.py "E1(8249,3720)"    # about 7 min (Theta_2 rank 1, Theta_4 = 0): nightly
python3 -I experiments/conductor-t11/t11_theta.py "E1(3540,2581)" "D1(12007,1799)" "E1(6480,4961)" "Sig+(12204,9853)" "E1(10512,145)+"   # about 2.5 min
python3 -I experiments/conductor-t11/t11_theta.py --negative        # about 30 s
```

Each line prints `# matches` against its committed output.

## What failed, exactly

- **Real quadratic fields with a large regulator.** A generator of an ideal of norm N there has about
  `R/(2 ln 10)` digits.
  - ℚ(√354384937) (h = 2, R ≈ 18900, GRH) gave conic coefficients of 5115 and 15443 digits after square reduction.
  - The lattice finder found no point in 20 trials. The norm finder stopped with a ZeroDivisionError in
    short_vectors (its 150-digit working precision).
  - That is D1(12007,1799), and E1(3540,2581) is worse: R ≥ 396790 for every field it could use.
- **Slow cases among small fields.** The finder over ℚ(√959) for the lift of +959 of Sig+(12204,9853) ran 20
  minutes without a point. The lift of −959 put ℚ(√959) on the other cone pair, and a point came in 13 s.
- **Unfinished.** The conic over ℚ(√−959) for the same lift, and the one over ℚ(√−23868431349) of E1(10512,145)+,
  were still running when the time ran out. The norm of the reduced coefficient was 134 digits; most of that time
  is probably PARI factoring norms.
- **E1(6480,4961)** was not attempted (time). Its real field ℚ(√34129302) has a regulator of about 490.
