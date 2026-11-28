# Fisher's higher descents on the five Magma-alone curves (note for `rank_fisher.py`)

Written by an AI model (Claude), 25 November 2026; not read by a person. Experiment, not counted. Fisher's paper
(arXiv:1509.03234v1) was read only through summaries made by a fetch tool; the PDF itself could not be downloaded.
So the statement numbers below are as those summaries gave them.

## The curves and the notation

There are five curves: `E+` of (760,39), (945,248), (72,65), (920,129), and `D1` of (631,71). Each is written
`E: y^2 = x(x^2 + a x + b)` with `T = (0,0)`. The isogenous curve is `E': Y^2 = X(X^2 + a'X + b')`, with `a' = -2a`
and `b' = a^2 - 4b`, and `phi: E -> E'` is the 2-isogeny. On all five curves `a^2 - 4b` is minus a square, so
`Q(E[2]) = Q(i)`. And `b` is not a square.

Fisher's notation is used (Section 2):
- `S_1 = S^(phi)(E/Q)` is rank_indep's "E' side": the classes `d' | b'`.
- `S'_1 = S^(phihat)(E'/Q)` is rank_indep's "E side": the classes `d | b`.
- `S_m` and `S'_m` are the images in `S_1` and `S'_1` of the higher Selmer groups `Sel(m)`. Here
  `Sel(2n) = S^(2^n)` and `Sel(2n+1) = S^(2^n phi)`, with E and E' alternating.

## Formulas used, and where they come from

1. **Fisher, Theorem 2.1.**
   - For `m` odd, `Theta_m` on `S_m x S_m` is alternating, with kernel `S_{m+1}`.
   - For `m` even, `Theta_m: S_m x S'_m` has left kernel `S_{m+1}` and right kernel `S'_{m+1}`.
   - The images of the rational points lie in every kernel.
2. **Fisher, Lemma 2.2.** The `2^n`-descent bound is
   `rank E(Q) <= dim S_{2n-1} + dim S'_{2n} - 2`, and the same holds with the primes swapped.
3. **Fisher, Section 3 (the Cassels–Tate pairing by pushouts).**
   - Let `x` be in `Sel(m)`, with covering curve `C`, and let `y` be in the first Selmer group of the other side when `m` is even, of the same side when `m` is odd (run 11, panel P1; `experiments/conductor-t11/` uses the odd case for `Θ₃`).
     Then `<x, y> = sum_v (f(P_v), y)_v`. Here `f` is in `Q(C)` with `div f = 2b`, where `b` is a rational divisor
     whose class corresponds to the rational 2-torsion point. The points `P_v` are in `C(Q_v)`.
   - Also `Theta_m(alpha_1, beta_1) = <alpha_m, beta_1>`, for any lift `alpha_m` (Section 4).
   - The Selmer set of the double covers `xi w^2 = f` is
     `S(C, f) = { xi : xi = f(P_v) mod Im(delta_v) for all v }`. It is empty, or a coset of a first Selmer group.
     It is non-empty exactly when `<x, eta> = 0` for every `eta` of the dual first Selmer group.
4. **Fisher, Lemma 7.1 (pushout forms).** If `l_1` and `l_2` are the tangent planes at points of two singular
   quadrics `Gamma_1` and `Gamma_2` of the pencil, and the pair `{Gamma_1, Gamma_2}` is Galois-stable, then `l_1 l_2`
   is a pushout. For the curves `D_lam: X^2 = lam U, Y^2 = lam V` of the second descent, both cones are rational, so
   `f = l_X l_Y`. This `f` is the one `rank_second.py` uses.
5. **The local images.** These come from rank_indep's first descent where every possible class is certified. At an
   odd prime of good reduction they are the unit classes. Where rank_indep leaves a class undecided (only at 2, on
   D1), the script uses `|Im delta_phi,p| |Im delta_phihat,p| = |E(Q_p)/2E(Q_p)| |E'(Q_p)[phihat]/phi(E(Q_p)[2])|`,
   which is 8 at p = 2 here. It combines this with lower bounds from the certified classes and from differences of
   values of `f`.

## What the script does, exactly

For each of the five curves:

1. **Levels 1 and 2** (rank_indep, rank_second with the v4 fix). The upper bounds are `S_2 <= <known, y0>` and
   `S'_2 <= <known, y1>`. Lifts `D_lam` of `y0` and of `y1` have a checked local point at every bad prime of `D_lam`
   and at R. So both dimensions are exact.
2. **`Theta_2(y0, y1) = +1`.** This is recomputed on every lift used. By (1), `S_3 = S_2` and `S'_3 = S'_2`.
3. **Level 3, constructively.** For each lift `D_lam`, the script takes one sampled value of `f = l_X l_Y` per place
   and the exact local image. It solves the linear system over F_2 for `S(D_lam, f)`. It checks that every further
   sampled value lies in the predicted coset. For one `xi` in the set, it finds a sampled point `P_v` with
   `xi f(P_v)` a non-zero square at every place of S3. Here S3 is R, 2, the bad primes of `D_lam`, E and E', the
   primes of `lam`, `delta` and `det`, and 3, 5, 7, 11, 13. The square tests are exact. Outside S3, the double
   cover has good reduction and is connected. A point exists by the Hasse bound (`p >= 17`, more than 8 points) and
   Hensel's lemma. So `C3 = {X^2 = lam U, Y^2 = lam V, xi w^2 = l_X l_Y}` is everywhere locally soluble, and the
   class lifts to `Sel(3)`.
   - **Cross-check.** The same Selmer set is predicted from Hilbert-symbol pairings with a basis of the dual first
     Selmer group. Both are computed, and they must agree.
   - **Control.** On E+ of (945,248) the lifts of the trivial class must give both outcomes. They give 12 empty and
     4 non-empty, and all agree.
4. **Level 4, by the theory.**
   - `Theta_3` is alternating on `S_3`, with kernel of codimension at most 1, so `Theta_3 = 0`.
   - The same holds for `Theta'_3`.
   - So `S_4 = S_3` and `S'_4 = S'_3`, and Lemma 2.2 with `n = 2` gives exactly `2 + 3 - 2 = 3`.
5. **The reduction of `Theta_4`.** Take the side whose first Selmer group equals its `S_2`: the E' side for
   (760,39), (72,65) and (920,129), the E side for (945,248) and for D1. Every level-3 lift of that side lifts to level 4: the obstruction pairs it with
   `<images of rational points, y>`, and `Theta_3(y, y) = 0`. With `<psi x, y> = <x, psihat y>`, this gives
   `Theta_4(y0, y1) = <[C3], [D]>`, the Cassels–Tate pairing of two explicit everywhere locally soluble curves.
   - `C3` is the level-3 lift of that side.
   - `D` is a level-2 lift `D_lam` of the other side.
   - Both curves are written to `rank_fisher_curves.json`.

## Result

| curve | exact bound | what it rests on | the 4-descent bound |
|---|---|---|---|
| E+ (760,39) | 3 | rank_indep (exact second descent) | exactly 3 (Thm 2.1, Lemma 2.2, levels 2–3 checked) |
| E+ (945,248) | 3 | same | exactly 3 |
| E+ (72,65) | 3 | same | exactly 3 |
| E+ (920,129) | 3 | same | exactly 3 |
| D1 (631,71) | 3 | same (via rank_second) | exactly 3 |

No bound below 3 is obtained, so none of the ten ratios leaves "Magma alone". Fisher's 4-descent cannot give less
than 3 on these curves. Magma's bound 1 needs `Theta_4(y0, y1) != 0` (Fisher's fifth descent, Section 7).

- If `Theta_4 != 0`, then `S_5` and `S'_5` are the images of the rational points. Lemma 2.2 with `n = 3` then gives
  `1 + 2 - 2 = 1`.
- If `Theta_4 = 0`, then `Theta_5` is again 0 by the same alternating argument, and `Theta_6` would be needed.

## Not done (not blocked)

`Theta_4` is not computed here, but it is not blocked. (An earlier version of this note said it was, and attributed
to Fisher a sentence that is not in his paper; that was wrong.) Fisher's Section 7 (arXiv:1509.03234) is exactly
this step, the fifth descent:

- The level-3 curve `C3` satisfies Fisher's hypothesis (19), `Im(rho) in <(12),(34)>`: the four singular quadrics
  of its pencil are the two over `Q(sqrt disc U)` and the two over `Q(sqrt disc V)`, and the pair attached to the
  rational 2-torsion point is Galois-stable.
- Lemma 7.1 gives a pushout form `l_1 l_2` on `C3`: `l_1` is the tangent plane at a point of one of the two
  conjugate cones, which needs a point of a conic over a quadratic field (non-split case, (20)).
- Theorem 7.2 gives the level-4 curve `C4` as a quadric intersection; Theorems 7.3 and 7.4 give a pushout form on
  `C4` from its companion `C4^dagger` (built from the other pair of cones); then
  `Theta_4(y0, y1) = sum_v (f(P_v), y)_v` by Section 3.
- Hypothesis (19) fails only one step later, at `C4` (Section 9). With one rational 2-torsion point this is where
  the method stops (`S_5`). It is presumably how Magma's `TwoPowerIsogenyDescentRankBound` reaches 1.

`Theta_4` is Fisher's fifth descent; the 8-descent is the sixth. `rank_fisher5.py` implements the first steps of
Section 7 on the level-3 curves (see its own output).
