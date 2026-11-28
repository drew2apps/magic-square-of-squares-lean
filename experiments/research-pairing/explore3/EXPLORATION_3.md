# Mixed curves, the modular picture, and where the difficulty sits (third note)

> **Read `NOTES.md` first.** This is the session's note, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.


27 November 2026 (UTC). Written by an AI model (Claude, Fable 5.2) with two sub-agents it started
(one computing, one acting as a panel); not read by a person; not reviewed. Tier D; no count changes.
Every count below was printed by a program in this folder.

## 1. The four mixed curves: a Rédei symbol with the octahedral invariant

With `(A, B, C) = (m² − n², 2mn, m² + n²)`, `R± = C² ± AB`, `W = R₊R₋ = m⁸ + 14m⁴n⁴ + n⁸`: the mixed
curve is `y² = (x − e1)(x − e2)(x − e3)` with `e2 − e1 = 2R₋`, `e3 − e1 = 2R₊`, `e3 − e2 = 4AB` exactly, on
all 804 curves of the pilot; in the triple it is `y² = (x + 2C²)(x − 2AB)(x + 2AB)`.

- **A general implementation of Cassels' pairing** (`cassels.py`: conics by `qfsolve`, local points
  by search, the real place by signs) gives −1 on all 6,776 mixed kill lines of the pilot, and on
  all 1,280 + 874 lines of `E_08`, `E_LU`. Controls on 300 lines: `⟨a, a⟩ = 1`, `⟨a, T_i⟩ = 1` for the
  three points of order 2, symmetry; four seeds give the same value and the same local factors on
  200 lines (`task1.log`, `extra_checks.log`).
- **Classes `(1, d, d)`** (first component at `e1`): their covering is
  `w² = d·t⁴ + 8C²·t²z² + (16A²B²/d)·z⁴`, whose conic becomes, with `s = dT + 4ABZ`, `D = dT − 4ABZ`,

      R₊·s² − R₋·D² = 2AB·d·w² ,

  so it has a point exactly when `d` is a norm from `ℚ(√W)`. The pairing of `(1, d, d)` with
  `(1, β, β)` computed from the tangent form of this conic equals the general pairing on
  2,995 + 2,784 pairs, no disagreement (`task23_*.log`). The set of such `d` has dimension 2 to 8
  on the 201 ratios; the same linear algebra reproduces the pilot's dimension of Sel² on 201 of 201.
- **Rédei form.** When `x² − d·y² = β·z²` has a rational point `(x0, y0, z0)`, the pairing equals

      ∏_{p | W, odd exponent}  ( (x0 + y0·√d) / p ) ,

  the Legendre symbol at `p` of `x0 + y0·s`, `s² ≡ d (mod p)`. This held on 14,160 of 14,160 pairs
  and on 16,638 of 16,638 in a second run; it did not depend on the sign of `s` or on the solution
  (14 solutions per pair). That is the Rédei symbol `[d, β, W]`; the conic over `ℚ(√W)` above is
  the same symbol read as `[W, d, β]`. For 2,248 of the first 2,995 pairs the conic
  `x² − dy² = βz²` has no point and only the second reading is defined.
- So the pairing on the mixed curves is not a product of Legendre symbols of the primes of the
  ratio; it is a Rédei symbol whose third entry is the octahedral form `W`. This is tested, not
  proved, and it covers the classes with first component 1; the killing classes `b` of the pilot
  have first component 1 on 2,760 lines, 2 on 488, another divisor of `2R₊R₋` on 3,528.
- **Where degenerate ratios show up here.** With `α1 = m(m + n)`, `β1 = n(m − n)`, `α2 = m(m − n)`,
  `β2 = n(m + n)`: `α1β1 = α2β2 = AB/2`, `R₊ = α1² + β1²`, `R₋ = α2² + β2²`, and
  `m² + n² = α1 − β1 = α2 + β2`. So `x² − α2·y² = β2·z²` and `x² − α1·y² = −β1·z²` have the point
  `(√(m² + n²), 1, 1)` exactly when `m² + n²` is a square. Not yet used.

## 2. The modular picture (checked by hand; elementary)

The `r`-line (`r = m/n`) is the modular curve `X(4)`: its six cusps are `r ∈ {0, ∞, ±1, ±i}`, the
zeros of `mn(m⁴ − n⁴)`; its automorphisms are the octahedral group; `E_08: y² = x(x − A²)(x − C²)` has
six fibres of type `I4` there and all its 4-torsion over `ℚ(i)(r)` (the halves of `(0, 0)` are
`x = ±AC`, `y = i·AC(C − A)`), so it is the elliptic modular surface of level 4, a K3 surface;
`E_LU` is its twist by −1. `E_mix` is another elliptic K3 over the same line, with twelve fibres
`I2` (at four cusps and at the eight zeros of `W`, the points with `j = 0`) and a section of
infinite order. The four corner forms vanish on the graphs of eight automorphisms of `X(4)`
(`t = ±i·r`, `t = ±i/r` and their conjugates by `x ↦ (x − 1)/(x + 1)`); the exchange of the two ratios
is `γ ↦ γ⁻¹`. Each pair surface `y² = F_j·F_k` is a K3 with two elliptic fibrations, one over each
ratio. I know no method that uses this for rational points.

## 3. The panel (one agent as statistician, surface geometer and referee; with web search)

- **A density theorem is plausible, "100 per cent" probably is not**: the number of locally allowed
  tuples and the size of the radicals are bounded in distribution, so a positive proportion of
  ratios should keep a tuple in every radical; the pilot already has two such ratios.
- Closest theorems are for twists of one curve (Smith; Kane) and do not apply: the family is not a
  twist family and has a rational point of order 4. Structurally nearest: Koymans–Pagano
  (arXiv:2412.01768). Needed and not found anywhere: joint equidistribution of Legendre symbols
  among primes of `m, n, m ± n, m² + n²` at one `(m, n)`, and of quartic symbols at the Gaussian
  primes of `m + ni` (sums of "spin" type).
- Rejected: Brauer–Manin on the surface; a Hilbert-modular reading (Chern numbers do not fit a
  quotient of a bidisc); Chabauty at degenerate ratios; Bombieri–Lang.
- Its honest answer on (R) for every ratio: no known route; through pair curves it needs the
  Mordell–Weil group modulo 2 of every member of a non-isotrivial family.

## 4. Where the difficulty sits: the parity of the two Frey curves (`rank0.py`, `rank0.txt`)

Rank 0 of `E_08` or of `E_LU` excludes a ratio (their eight torsion points lie over
`t ∈ {0, ±1, ∞}`). With the general code (Sel² and the pairing), on 150 random non-degenerate
ratios with `m ≤ 3000` and a 2-adic point: the bound is 0 on `E_08` for 37, on `E_LU` for 34, on at
least one for 62 (41 per cent). The parities of the two bounds are (even, even) 36, (even, odd) 47,
(odd, even) 41, (odd, odd) 26. So for about five ratios in six one of the two curves has even
parity and should have rank 0 or a large Ш; **the hard ratios are those where both bounds are odd**
(26 of 150 here; a quarter is the natural guess), where both curves have a point of infinite
order if Ш is finite. On 80 ratios of square type the numbers are similar (rank 0 on one curve for
38). A magic square of squares needs positive rank on `E_08`, `E_LU` at both of its ratios.

## 5. What I would do next

1. Write the bad-prime laws of tuples (I1) as linear conditions, next to `S′` and `Θ` of the second
   note and the Rédei symbol above. Then "which tuples survive" is linear algebra over `𝔽₂` with
   symbols as entries, and the failure rate `f(H)` the panel asks for can be computed to large
   height without any descent code.
2. Prove the Rédei form of section 1 (it should follow from Rédei reciprocity in Stevenhagen's
   form, arXiv:1806.06250) and find the classes with first component not 1.
3. Restrict attention to the both-odd ratios of section 4: describe them by root numbers (a
   count of primes `≡ 3 mod 4` of `mn(m² − n²)` and a term at 2), and test whether the local
   conditions of a solution correlate with that parity.
