# H6: the genus-1 and genus-2 quotients of the (ℤ/2)⁵-cover for (627,364) and (991,263)

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Machine: PARI 2.19.1 with cypari2, not CI's PARI 2.17.4. Real sympy 1.14.0. PARI's `ellrank` is used as a finder only.

## Bottom line

- **No Demjanenko–Manin or Chabauty route opens on these quotients.** *(computed: the genus-2 bounds are exact descents with certified local images; the genus-1 ranks are PARI's, consistent with PARI)*
  - The four genus-2 pieces F0·F4·L, F0·F4·U, F8·F4·L, F8·F4·U are one curve A up to isomorphism over ℚ, so J(cover) ~ A⁴ × B. Demjanenko–Manin through these four maps needs rank J(A)(ℚ) ≤ 3. The Richelot descent gives rank ≤ 8 for both ratios.
  - Every other genus-2 piece has rank bound 6 to 9. Chabauty on one piece needs rank ≤ 1.
  - Among the genus-1 pieces, the rank-4 class c3 has 4 maps, one short of what Demjanenko–Manin needs (rank < 4).
- **The two ratios mirror each other.** 991 = 627 + 364 and 263 = 627 − 364. Their pieces have the same Frobenius polynomials and the same rank bounds, in swapped positions. This is an internal consistency check. That the two covers are isomorphic is not proved here.

## The quotients (h6_split.py, h6_ranks.py)

- **The pieces:** for each subset S of {F0, F8, F4, L, U} with |S| ≥ 2, the curve y² = ∏_S 4F_j(t) has genus |S| − 1.
  - Each is split further by every involution of ℙ¹ that leaves the form exactly invariant (x → −x, x ↔ z). The quotients come from Kani–Rosen, which is quoted.
  - Every invariance is an exact polynomial identity, in pure Python.
- **Genus 1:** the Jacobian is built from the quartic's invariants I, J and must agree with PARI's `ellfromeqn`. PARI's `ellrank` is run on every curve of the 2-isogeny class, and the class gets the largest lower bound and the smallest upper bound.
- **Genus 2 with no involution:** the Frobenius polynomial is computed from point counts over F_p and F_{p²}. If it is irreducible at one good prime, the Jacobian is simple over ℚ. *(computed, exact)*

**(627,364)** (class labels as in h6_ranks_627_364.txt; [lo, hi] are PARI's bounds):

| rank | genus-1 classes (number of maps) |
|---|---|
| 1 (bounds meet) | c1 F0F8 (1), c5 LU (1), c14 F0F8LU (1) |
| [1,3] | c7 (2), c8 (2), c13 (1), c15 (1) |
| 2 | c2 (2), c6 (2), c10 (2), c11 (2), c16 (2) |
| 3 | c4 (2), c12 (2) |
| 4 | c3 F0L, F0U, F8L, F8U (4); c9 (2) |

The genus-2 pieces are all simple over ℚ: F0F4L/U and F8F4L/U at p = 31, the two F0F8F4L/U-wP at 23, F0F4LU-wP and F8F4LU-wP at 41, and F0F8F4LU P and wP at 17.

**(991,263):** the same shape, with the classes swapped as in the bottom line. For example, F0·F4 has rank 3 and F4·L rank 2 here, against 2 and 3 for (627,364). The F0F8L classes c7 and c8 are decided here (4 and 2), but F0LU's c9 and c10 are [1,3].

## Richelot descent on the genus-2 pieces (h6_richelot.py)

- **The curve:** each piece is y² = G1·G2·G3 with three rational factors of degree ≤ 2.
  - The descent goes through the Richelot isogeny to C′: Δy² = H1H2H3 (Bost–Mestre).
  - The twist of C′ is fixed by equal point counts over F_p and F_{p²}.
- **Local images:** found by sampling local points, near roots mod p too, with conjugate pairs over small quadratic fields. They are certified exactly by the product formula |W_p|·|W′_p| = 16 (p odd), 64 (p = 2), 4 (ℝ), which holds by Tate local duality, quoted.
- **The bound:** rank J(ℚ) ≤ dim Sel(δ) + dim Sel(δ′) − 2 − (the 2-torsion beyond T1, T2 that the images show). Every grouping of the factors is tried, and the best bound is kept.
- **Known points:** the differences of known global points are asserted to land in the Selmer groups.

| piece | maps in the cover | (627,364) | (991,263) |
|---|---|---|---|
| F0F4L, F0F4U, F8F4L, F8F4U | 4 (isomorphic) | ≤ 8 | ≤ 8 |
| F0F8F4L-wP, F0F8F4U-wP | 2 | ≤ 9 | ≤ 6 |
| F0F4LU-wP, F8F4LU-wP | 2 | ≤ 6 | ≤ 9 |
| F0F8F4LU-P, F0F8F4LU-wP | 1 each | ≤ 8 | ≤ 9 |

*(computed, exact given the quoted theory: the Richelot construction, the Kummer map by (G1, G2), local duality, unramified images at odd primes of good reduction)*

These are upper bounds only. The true ranks may be smaller, but closing the gap from 8 to 3 would need a much deeper descent than this.

## Why the four pieces are one curve

- σ: t → −t swaps L and U and fixes F0, F8, F4.
- τ: t → 1/t swaps F0 and F8 (after the factor t²) and fixes F4, L, U.
- So F0F4L ≅ F0F4U ≅ F8F4L ≅ F8F4U over ℚ, from these polynomial identities; their Frobenius polynomials agree in h6_ranks, as they must.
- The four quotients belong to distinct characters of the (ℤ/2)⁵ action, so they are four different factors of J(cover), each isogenous to J(A). That is the A⁴ in J(cover) ~ A⁴ × B. *(proved: the argument is standard, written out here, not machine-checked)*

## Checks

- **h6_richelot.py --controls (h6_richelot_controls.txt):**
  - y² = (x² + x + 1)(x² − x + 1)(x² + c) for 10 values of c: an even sextic, so its Jacobian is isogenous to E1 × E2;
  - the descent runs on a grouping with Δ ≠ 0, and its bound equals PARI's rank E1 + rank E2 (bounds meet) on 10 of 10. The bound need only be at least that sum, so equality is the stronger check.
- **h6_richelot.py --negative (h6_richelot_negative.txt):**
  - a wrong descent map is refused at p = 2 (local images too large);
  - dropping one generator at the bad prime 5 gives a product of 8, not 16, so it is not certified;
  - the wrong twist of the dual is refused by the point counts.
- **h6_ranks.py --negative:**
  - a forged Jacobian is refused;
  - a split curve is not shown simple;
  - a generic sextic is shown simple;
  - 37a1 and 11a1 get ranks 1 and 0.

## What failed

- **Hangs at large primes:** the first local-image sampler hung on nfinit of huge discriminants, and it could not certify p = 991, where L has a double root mod p. Fixed by sampling near roots mod p and by using small quadratic fields.
- **Negative test:** it first used p = 3, which is a good prime. Fixed by using 5.
- **times.txt:** the runs killed during these fixes are listed there with exit 143.
- **No candidate:** the result is negative.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** The genus-1 and genus-2 quotients of the (ℤ/2)⁵-cover for (627,364) and (991,263) (experiments/roadmap-h/h6): every genus-2 piece is simple over ℚ, and an exact Richelot descent bounds its rank by 6 to 9. The four pieces F0F4L, F0F4U, F8F4L, F8F4U are one curve, so Demjanenko–Manin would need rank ≤ 3, and the bound is 8. No Demjanenko–Manin or Chabauty candidate. Genus-1 ranks are PARI's bounds. Computed. Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-h/h6        # F = grep -v '^PARI stack size'
python3 h6_split.py 627 364 2>&1 | $F | tail -1 | grep -q '^# matches'      # about 1 s
python3 h6_ranks.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'   # about 1 s
python3 h6_richelot.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches' # about 7 s
python3 h6_richelot.py --controls 2>&1 | $F | tail -1 | grep -q '^# matches' # about 70 s
# not for CI: bash run_h6.sh m n (ranks, about 3 min) and bash run_h6_rich.sh m n (about 5 min); times.txt
```
