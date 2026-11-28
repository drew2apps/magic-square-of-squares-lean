# Exploration 5: the pairing as Hilbert symbols in Q(sqrt d) (Redei form derived), and follow-ups (2026-11-27)

> **Read `NOTES.md` first.** This is the session's note, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.


Tier D. Written by an AI model (Claude), not read by a person. "Proved" below means: derived here step by step,
given Cassels' theorem, and not yet checked by anyone else. No counts change. Nothing was published.

## 1. Theorem B: classes with first component 1 on any curve with full 2-torsion

Setting. `E: y² = (x − e1)(x − e2)(x − e3)` over `ℚ`. `a = (1, d, d)` and `b = (1, β, β)` in Sel², `d` not a square,
`K = ℚ(√d)`. The 2-covering of `a` is `u1² = x − e1`, `d·u2² = x − e2`, `d·u3² = x − e3`. Put

    z2 = u1 + √d·u3   (norm c2 = e3 − e1),      z3 = u1 + √d·u2   (norm c3 = e2 − e1),      γ = z2·z3.

**Theorem B (given Cassels' theorem).** Suppose `β = N_{K/ℚ}(π)` for some `π` in `K` (that is, `x² − d·y² = β·z²` has a
rational point). Then

    ⟨a, b⟩ = ∏_w (γ_w, π)_{K_w}        (product over all places w of K, Hilbert symbols of K),

where for `w | v` the element `γ_w` is the value of `γ` at any point of the covering over `ℚ_v`.

**Proof.** Cassels: `⟨a, b⟩ = ∏_v (L2(P_v)·L3(P_v), β)_v`, where `L2`, `L3` are tangent forms, at rational points, of the
conics `u1² − d·u3² = c2` and `u1² − d·u2² = c3` (the conic that omits `u_j` goes with `b_j`; `b_1 = 1` drops out).
Write a rational point of the conic `k` as `q_k` in `K`, `N(q_k) = c_k`.

Step 1 (tangent form = norm of a difference). `L_k(P) = −N(z_k − q_k)/2` (identity on the conic; `tangent_norm.py`).
With `ε_k = z_k/q_k` (norm 1) and `ε_k = η̄_k/η_k` (Hilbert 90 in `K ⊗ ℚ_v`; `η_k` is unique up to `ℚ_v^×`),
`N(ε − 1) = −4d·y²/N(η)` for `η = x + y√d`. So `L_k(P_v) = 2·c_k·d·N(η_k)` up to squares of `ℚ_v`.

Step 2 (down to up). For `η` in `(K ⊗ ℚ_v)^×` and `β = π·π̄`:
`(N η, β)_v = ∏_{w|v} (η, β)_w = ∏_{w|v} (η, π)_w (η, π̄)_w = ∏_{w|v} (η, π)_w (η̄, π)_w = ∏_{w|v} (η·η̄, π)_w = ∏_{w|v} (ε, π)_w`,
using: `(N η, β)_v = ∏_{w|v} (η, β)_w` for `β` in `ℚ_v`; conjugation permutes the places above `v`; `η·η̄ = ε·η²`.

Step 3. `⟨a, b⟩ = ∏_v (4·c2·c3·d², β)_v · ∏_w (ε_2 ε_3, π)_w = ∏_w (γ_w, π)_w / ∏_w (q_2 q_3, π)_w = ∏_w (γ_w, π)_w`,
by the product formula over `ℚ` (first factor) and over `K` (the global element `q_2 q_3`). ∎

No reciprocity law is used, and the global conic points `q_k` have disappeared: the right side needs only local
points of the covering and one solution of `x² − d·y² = β·z²`. Without the hypothesis on `β` the same proof gives
`⟨a, b⟩ = ∏_v (N(η_{2,v})·N(η_{3,v}), β)_v`, which still needs the `q_k`.

**Reading as a Redei symbol.** `N(γ) = c2·c3`. At a place `w` not above 2 or ∞ where `π` is a unit,
`(γ_w, π)_w = (π/w)^{w(γ_w)}`: the residue symbol of `π` at the primes of `K` that divide `γ`. On the mixed curves
`c2·c3 = 4W`, `W = R+·R− = m⁸ + 14m⁴n⁴ + n⁸`: the main part is `∏ (π/𝔭)` over primes `𝔭` of `K` above the primes of odd
exponent in `W`, which is the symbol tested in exploration 3 (there on 14,160 and 16,638 cases). The remaining
factors sit at the places above 2, `d`, `β` and ∞. So the Redei form is now derived, not only observed, for classes with
first component 1 and `β` a norm.

**Test** (`redei_k.py`): right side computed with PARI's `nfhilbert` at all places above the primes of
`2·d·β·N(π)·∏(e_i − e_j)` and ∞, against the general engine: 784 of 784 pairs agree (507 on `E_mix`, 239 on `E_08`,
31 on `E_LU` at ten ratios, 7 on random curves); 4,200 further pairs were skipped because `β` is not a norm from
`ℚ(√d)` (most) or no local point was found by sampling. A first run showed 2 disagreements; the cause was a sign bug in
my test harness (products of negative classes), not the identity; `redei_dbg2.py` shows the two sides agree place
by place there.

**Not done.**
- Tuples on the mixed curves have classes `(2, κ·d_j·d_k)`, first component 2, not 1. Pairing such a class against
  `b = (1, β, β)` on the covering of `b` gives, by the same steps, a Theorem-B term (with `K = ℚ(√β)`) times
  `∏_v ((u2 − u3)·N(η_{3,v}), 2)_v`. I have not put this extra factor in closed form.
- The factors of Theorem B at the places above 2, `d`, `β`, ∞ are not yet written as explicit symbols for the family.
- Pairs where `β` is not a norm from `ℚ(√d)`.

## 2. The centre form is not in the stage-2 local test, and adding it removes most of what is left

Found by a sub-agent (another AI model instance); the self-test and the open-ratio run below are mine.

The repo's tuple test (`tuple_descent.magic`) uses four Gaussian forms: F0, F8, GL, GU. The centre
`F4 = 2(m² + n²)(u² + v²)` is the norm of `γ4 = (1 + i)(m + in)(u + iv)`. A solution needs F4 to be a square as well
(`L + U = 2·F4` does not make it one). Adding `γ4` as a fifth form, with a fifth class `d4` supported on the primes of
`2(m² + n²)`, and running the same code at the same bad primes (`f4_form.py`, `f4_all.py`):

| | four forms | five forms |
|---|---|---|
| ratios (of the 937) that keep a tuple | 937 | 211 |
| tuples | 6,904 | 1,370 extend to a 5-tuple |
| ratios with a tuple that also survives the six pairings | 5 | 1: (14903,9696) |
| sampled ratios beyond m = 15000 with a tuple (269) | 269 | 46 |
| pairing survivors in that sample | 6 tuples at 3 ratios | 0 |

Checks:
- Six degenerate ratios keep exactly their known solutions with five forms (`f4_form.txt`).
- Soundness (`f4_selftest.py`, mine): the repo's own `selftest`, applied to the five forms at seven ratios, finds 0
  failures on 37,974 points with the repo's `Local` and 0 on 27,873 points with `exact_local`.
- At the nine open ratios of the census (`f4_open.py`, mine) the five-form test leaves exactly 2 tuples at each
  (the number of known degenerate solutions there); with four forms it left 4 at (609,200), (809,409), (627,364), (991,263).

Not checked: that the 5-tuple code is complete in the sense the repo requires for a count (only the soundness
self-test was run), and whether the census already imposes F4 elsewhere. If it does not, this is a free strengthening
of stage 2: no new prime, no new curve. Tier D; no counts change.

Next on this line: pair curves with the centre, `y² = F_j·F4`, and the 5-tuple classes on them, for (14903,9696)
and its partner (24599,5207), the only ghost ratios left.

## 3. Ghost tuples: more data, two partial laws (sub-agent; not re-run by me)

- Prediction confirmed: the partner (24599,5207) of (14903,9696) has 8 residual tuples and 2 survivors, both of
  ghost shape.
- Beyond m = 15000: three samples, 321 ratios with residual tuples, 1,366 tuples; 6 survivors at 3 new ratios,
  (90384,58537), (126127,16889), (115513,18217); none of another shape. In all: 18 surviving tuples at 9 ratios,
  18 of ghost shape. The biased sample was stopped by hand and its number of ratios tried is only known as "at least 5,892".
- Law (a), necessary only: a ghost-shape tuple is residual only if every odd prime of odd exponent in `m² − n²`
  (shape S1) or in `mn` (shape S2) is ±1 mod 8: 333 of 333 and 166 of 166 ratios.
- Law (b): for `m, n` of opposite parity, `(1, 2mn, 2mn)` is in Sel²(`E_08`) exactly when no prime ≡ 5 mod 8 divides
  `m² − n²`: 1,994 of 1,994. For `m, n` both odd it is never in Sel²: 1,006 of 1,006.
- Law (c): no law found for when the ghost class is killed on `E_08`; four quartic-residue predicates failed.
- Ranks: at six ghost ratios the ghost class is a torsion class on one of `E_08`, `E_LU` and a non-torsion class in
  the kernel of the pairing on the other. At (3320,2303) and (5623,1017) PARI reports rank [1,1] on that curve.
  At the other four the rank is 0 or 2, not decided. No point of height up to 10⁵ on the ghost coverings.
- Not done: the class map to the genus-3 quotients E1, E2, E3; a 4-descent.

## 4. The (627,364) curve over Q(i): still not decided (sub-agent; not re-run by me)

- Model: `Y² = X(X − e2)(X − e3)`, `e2 = 273990441250`, `e3 = −102929464375 + 118967496648i`; torsion Z/4 × Z/2.
- dim Sel² = 4 over `ℚ(i)` (taken from an earlier run's `kdescent.py`, not re-derived).
- Searches: `X = e_j + b_j(U/W)²` for all 16 Selmer elements, `U, W` Gaussian integers with coordinates up to 100:
  only the two points of order 4. On the 2-isogenous curve `E/⟨(e3,0)⟩`, six classes, coordinates up to 280: nothing.
- Limits: the sieve was not tested on a planted non-torsion point, and no height bound for generators was derived.
  So this is "no small point found", not evidence for rank 0. Rank 2 is not proved; the route is not retired.
- Next: reduced 2-covering quartics and 4-coverings over `ℚ(i)`; Simon's `ell.gp` did not run under PARI 2.17.

## 5. Reviewer's eight fixes to the E_08/E_LU note

Done in `PAIRING_FORMULA_3_2026-11-27.md` (bundle `pairct3.tgz`, with `h2_check.py`). The good-prime step now needs no
bound on `p`; the reviewer's H2 counter-example is reproduced (32 of 64 pairs wrong at `E_LU` (725,497)).

## Files in explore5.tgz

`redei_k.py`, `redei_k.txt`, `redei_dbg2.py`, `tangent_norm.py`; `ghost/` (sub-agent: `GHOST_NOTE.md`, scripts, outputs,
plus my `f4_selftest.py`, `f4_open.py`); `kcurve/` (sub-agent: `KCURVE_NOTE.md`, scripts, logs).
