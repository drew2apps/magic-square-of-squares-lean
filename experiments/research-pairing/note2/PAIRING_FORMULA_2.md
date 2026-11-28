# The pairing on `E_08` and `E_LU`: every local term derived, and the part with first component 1 from factorisations alone

> **Read `NOTES.md` first.** This is the session's second note, committed unchanged apart from this line (and `shape.txt`, regenerated: the session's copy lacked one line the script prints); where they differ, `NOTES.md` (after review z13) holds.

27 November 2026 (UTC), second note. Written by an AI model (Claude, Fable 5.2); not read by a
person; not reviewed. Tier D; no count changes. It continues `PAIRING_FORMULA.md` and uses the
corrections of review z12 (coprime `m > n`; the convention of the three forms is forced).

## 0. Setting

`(A, B, C)` is a primitive Pythagorean triple with `C` odd and `A, C > 0`; `E: y² = x(x − A²)(x − C²)`,
`e_r = 0`, `e_m = A²`, `e_* = C²`. For coprime `m > n > 0`: `E_08` is the curve of
`(m² − n², 2mn, m² + n²)` and `E_LU` that of `(2mn, m² − n², m² + n²)`, each divided by 2 when
`m, n` are both odd (the curve does not change). A class is `(c_r, c_m, c_*)`, product a square.
`a = (1, δ, δ)`, `b = (b_r, b_m, b_*)` are in Sel². The 2-covering of `a` is

    u² = x,   δ·u_m² = u² − A²,   δ·u_*² = u² − C²,

and Cassels' pairing with the conic points `(A, 0)`, `(C, 0)`, `u_m = u_*` is

    ⟨a, b⟩ = ∏_v t_v,     t_v = (u − A, b_*)_v (u − C, b_m)_v (u_m − u_*, b_r)_v            (1)

at any point of the 2-covering over `ℚ_v` (Cassels 1998: for classes of Sel² each `t_v` does not
depend on the point). **Hypothesis H2:** if `A` is odd, 32 divides `C − A` or `C + A`. It holds for
`E_08` and `E_LU` whenever the ratio has a 2-adic point (`4 | mn` or `m ≡ ±n mod 8`), since
`C ∓ A` is `2n², 2m²`, resp. `(m ∓ n)²`, or half of these.

## 1. The formula

    ⟨a, b⟩ = ∏_{p | B, C ≡ −A (p)} (δ, b_*)_p  ·  ∏_{v | AC, and v = 2 if A is even} T_v ,            (2)

    T_v = (2√δ, b_r)_v                                       if δ is a square in ℚ_v,
    T_v = (−A, b_*)_v (−C, b_m)_v ((A − C)·√(−δ), b_r)_v      otherwise (then −δ is a square in ℚ_v).

For `p = 2` (when `A` is odd) "`C ≡ −A (p)`" means `C ≡ −A (mod 8)`. All other places, the real
one included, contribute 1. No local point and no conic is needed to evaluate (2).

## 2. Proofs of the local terms (given Cassels' theorem)

**Lemma 1 (local images where `e_m, e_*` meet).** Let `p | B` be odd, or `p = 2` with `A` odd and
H2. Every class in the image of `E(ℚ_p)` has `b_r` a square in `ℚ_p` and `b_m = b_*` in `ℚ_p^×/□`.
*Proof.* Let `x` be the abscissa of a point; `x(x − A²)(x − C²)` is a square. `(x − C²)/(x − A²) = 1 − B²/(x − A²)`.
If `v(x − A²) ≤ v(B²) − 1` (odd `p`), resp. `≤ v(B²) − 3` (`p = 2`), this ratio is `≡ 1` modulo `p`,
resp. modulo 8, so it is a square, `b_m = b_*`, and then `x` is a square. Otherwise `x ≡ A²` modulo
`p`, resp. modulo 8 (here H2 gives `v₂(B²) ≥ 6`), so `x` is a square and again `b_m = b_*`. The
points of order 2 have the images `(1, −1, −1)`, `(1, −1, −1)`, `(1, 1, 1)` there. ∎

**Lemma 2 (where `e_r` meets `e_m` or `e_*`).** At an odd `p | A`, a class in the local image has
`b_*` a unit that is `±1` modulo squares; at an odd `p | C` the same for `b_m`. At `p = 2` with `A`
even, `(1, δ, δ)` is locally soluble only if `δ` or `−δ` is a square in `ℚ_2`.
*Proof.* Odd `p | A`: if `x ≡ 0`, then `x − C² ≡ −C²`; if `x ≡ C²`, then `x, x − A²` are unit
squares and `x − C²` is a square; otherwise `x ≡ x − A²` and `x − C²` is a square. `p | C`: the
same with the roles exchanged. `p = 2`, `A ≡ 0 (4)`, `C` odd, a point `u`: if `v(u) ≤ −2` both
`1 − A²/u²`, `1 − C²/u²` are `≡ 1 (8)` and `δ` is a square; `v(u) = −1` gives the classes 1 and 5,
which differ; `u` odd gives `u² − A² ≡ 1 (8)`, so `δ` is a square; `u ≡ 2 (4)` gives `u² − C² ≡ 3 (8)`
and `u² − A² ∈ 4·{1, 5}`, which differ; `u ≡ 0 (4)` gives `u² − C² ≡ −1 (8)`, so `−δ` is a square. ∎

So for `a = (1, δ, δ) ∈ Sel²`: `δ` is `±` a square at every `v | AC` (and at 2 if `A` is even), which
is what (2) needs.

**(i) `δ` a square in `ℚ_v`** (any place). Take `u` with `|u|_v` large and `s = √δ`. Then
`u_m = (u/s)(1 − A²/2u² − …)`, `u_* = (u/s)(1 − C²/2u² − …)`, so `u_m − u_* = (B²/2su)(1 + ε)` and
`u − A`, `u − C` are `u` times squares. `t_v = (u, b_* b_m)(2su, b_r) = (2s, b_r)_v`.

**(ii) `−δ` a square in `ℚ_v`** (any place). Take `u` near 0 and `r = √(−δ)`: `u_m ≈ A/r`,
`u_* ≈ C/r`, `u − A ≈ −A`, `u − C ≈ −C`. `t_v = (−A, b_*)(−C, b_m)((A − C)r, b_r)_v`.
At the real place both give 1 (all three first entries are negative in (ii), and `b_r b_m b_*` is a square).

**(iii) `p | B` with `C ≡ A`** (odd `p`: mod `p`; `p = 2`: mod 32 by H2). By Lemma 1,
`t_p = ((u − A)(u − C), β)_p` with `β = b_*`. If `u − A` has valuation at most `v(C − A) − 1`
(resp. `− 3`), then `(u − C)/(u − A) ≡ 1` and the product is a square. Otherwise `u ≡ A`, `u + A` has
valuation 0 (resp. 1), and `(u + C)/(u + A) = 1 + (C − A)/(u + A) ≡ 1` modulo `p` (resp. 16), so
`(u + A)(u + C)` is a square; as
`(u² − A²)(u² − C²) = δ²u_m²u_*²` is a square, so is `(u − A)(u − C)`. So `t_p = 1`.

**(iv) `p | B` with `C ≡ −A`.** Again `t_p = ((u − A)(u − C), β)_p`. If `u + A` has small valuation,
`u − C = (u + A) − (C + A)` is `u + A` times a square, so `(u − A)(u − C) ∈ (u² − A²)·□ = δ·□`.
Otherwise `u + C = (u − A) + (A + C)` is `u − A` times a square, so `(u − A)(u + C)` is a square,
hence `(u − C)(u + A)` is one, and `(u − A)(u − C) ∈ (u² − C²)·□ = δ·□`. So `t_p = (δ, β)_p`.

(iii) and (iv) hold at every local point; (i) and (ii) use one point each and Cassels' theorem.
An odd `p ∤ ABC` gives 1 (all entries are units for a suitable point, or by (i), (ii)). ∎

## 3. The classes with first component 1, without any descent

    S′ = { β squarefree, a product of −1 and primes of B :  β is ± a square in ℚ_p for every p | AC, and in ℚ_2 if A is even }.

Then `{b ∈ Sel² : b_r = 1} = {(1, β, β) : β ∈ S′}`, and on it (from (2), with `b_r = 1`)

    Θ(δ, β) = ∏_{p | B, C ≡ −A (p)} (δ, β)_p · ∏_{v | AC or v = 2 | A :  δ not a square, −δ a square in ℚ_v} (AC, β)_v .            (3)

The inclusion "⊂" is Lemmas 1, 2. For "⊃": at `p | B` (and at 2 when `A` is odd) the local image
has 4 (resp. 8) elements, all of the form `(1, β, β)` by Lemma 1, and there are exactly 4 (resp. 8)
such classes; at `v | AC` (and at 2 when `A` is even) the classes `(1, ±1, ±1)` are images of points
of order 2; a good prime allows the unit classes; the real place allows both signs. `S′` is the
image in Sel² of the Selmer group of the 2-isogeny `E → E/⟨(0, 0)⟩`; `Θ` is the pairing on it
whose kernel is the next step of a descent by 2-isogenies (Fisher 2017). Both are now matrices of Legendre symbols of
the primes of `m, n, m ± n, m² + n²`, in the style of Monsky's matrix.

## 4. Tests (referee z12's tables `full0.jsonl`, `full1.jsonl`: every class against all of Sel²)

- `closed_formula.py`: formula (2) with no local point: 69,304 pairs on `E_08` (146 ratios) and
  35,096 on `E_LU` (133 ratios) agree with the committed values; no disagreement; never undefined.
  (316 rows of `E_LU` whose class is printed with `c_r ≠ 1` were skipped.)
- `two.py`, `two2.py` (the pilot's kill lines): when `A` is even `δ` is always `±` a square at 2 and
  (i), (ii) agree with the numerical term on 720 lines; when `A` is odd, `b_r` is a square and
  `b_m = b_*` at 2 on all 1,434 lines, `v₂(C ∓ A)` is 5 or at least 7, and (iii), (iv) agree on all.
- `isogeny_part.py`: `S′` from factorisations equals `{b ∈ Sel² : b_r = 1}` on all 146 + 133 curves
  (dimension 1 to 7); `Θ` equals the pairing on 29,484 + 11,928 pairs.
- `shape.py`: of the tuple classes in those tables that are killed on the curve, 172 of 202 on
  `E_08` and 86 of 134 on `E_LU` are killed by some `b` with `b_r = ±1`, that is by `Θ` alone.
  (Rows count a class and its negative.)

## 5. What this says about tuples with trivial pairing

On `E_08` a tuple enters through `δ = mn·d0·d8` only, on `E_LU` through `δ = (m² − n²)·dL·dU` only.
It passes the `Θ`-test on a curve exactly when its `δ` lies in the radical of `Θ` on `S′`. The
radical always contains `±1` (the class of the 2-torsion). When it is no more than that, a
surviving tuple has `d0·d8 = ±mn` (resp. `dL·dU = ±(m² − n²)`) up to squares: its class is a
torsion class, which no pairing can see, and which only a rank statement removes. This is the
shape of the pilot's exceptions: at (5623,1017) the two tuples have `d0·d8 = 4·635399` and the
squarefree part of `mn` is 635399. So the candidate answer to the question is: *a tuple that
survives the pairing on a pair curve of full square differences has a torsion class there, or the
radical of `Θ` is larger than `±1`*. Not proved; the full pairing (with `b_r ≠ 1`) and the mixed
curves cut further.

## 6. The four mixed curves and Fisher's 2022 formula (not solved)

In the triple the mixed curve is `y² = (x + 2C²)(x − 2AB)(x + 2AB)`, with differences
`2R± = 2(C² ± AB)` and `4AB`; only 232 of the 1,762 killed tuples of the pilot are killed on mixed
curves alone. Fisher's formula takes three binary quartics `g1, g2, g3` with the same invariants
and `z(g1)z(g2)z(g3)` a square in the cubic algebra; it saves the conics of the pairing itself,
not the one needed to write a class as a quartic or to find `g3`. Here the classes are given as
`(c_r, c_*)`, the conics are `2u² − κD·w² = 2R±` with `D = d_j·d_k`, and I found no point of them
that is a polynomial in `m, n` and the tuple (for `D = 1` one of the two has the point `(J, 2)`).
So this note has no formula for the mixed curves, and no proof that none exists.

## Files

`closed_formula.py`, `isogeny_part.py`, `shape.py` read `full0.jsonl`, `full1.jsonl` from the
current folder; `two.py`, `two2.py` read the pilot's `t10_*_[ab].txt` and import `parse.py`,
`ct08.py` of the first note. Outputs: the `.txt` files of the same names. Two paths in the scripts
(`/tmp/t10`) point to the folder of the first note and have to be adapted.

## 7. Where "degenerate" enters: the radical of `Θ` (`radical.py`, factorisations only)

1,500 random coprime non-degenerate ratios with `m ≤ 20000` and a 2-adic point, and 400 ratios of
square type; the dimension of the radical of `Θ` on `S′` (it always contains `−1`):

| | 1 | 2 | 3 | 4 | 5+ |
|---|---|---|---|---|---|
| `E_08`, non-degenerate | 482 | 636 | 323 | 57 | 2 |
| `E_LU`, non-degenerate | 494 | 686 | 270 | 46 | 4 |
| `E_08`, square type | 0 | 94 | 207 | 82 | 17 |
| `E_LU`, square type | 194 | 179 | 26 | 1 | 0 |

At a ratio of square type the radical on `E_08` is never `±1`. The reason is exact: the point of
order 4 has the class `(S, 2mn, 2mnS)`, and when `S` is a square this is `(1, 2mn, 2mn)`, a class of
`S′`. It is the class of the degenerate solutions at `t = ±1` (there `d0·d8 = 2`). When `S` is not a
square, `(1, 2mn, 2mn)` differs from a torsion class by `(S, 1, S)`, and formula (2) gives its
pairing with `β ∈ S′` as a product that contains `(2√β, S)_p` at each prime `p` of odd exponent in
`S`: the quartic residue symbol of `β` there. So on this curve the difference between the two
kinds of ratio is one class, and whether it is killed at a non-degenerate ratio is a sign of
quartic symbols at the primes of the odd squarefree part of `m² + n²`, not a uniform law. For
about a third of the non-degenerate ratios the radical is `±1` on a given curve, and then only a
tuple with `d0·d8 = ±mn` (resp. `dL·dU = ±(m² − n²)`) passes `Θ` there. These are counts on a
sample, not a density statement.
