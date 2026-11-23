import Mathlib.Tactic
import Mathlib.RingTheory.Coprime.Basic
import Mathlib.Data.Real.Basic

/-!
# Algebra behind the argument with four points on `E0L`

`G1_FULL_SYSTEM.md`, "Rank 2 and 3 on `E0L`", and `scripts/g1_full/dm4.py` use three identities for
coordinates of sums and differences of points of the curve `y² = F0·GL`, and elementary bounds for
the heights of the resulting fractions. The identities themselves are statements about the group
law of the curve and are checked by the script, not here. This file contains the integer algebra
that the height bounds rest on. Notation: `t0 = p/q` in lowest terms, `S = m² + n²`,

* `Ψ0 = m²q² + n²p²`, `Ψ8 = n²q² + m²p²`, `Ψ = S(p² + q²) − 4mnpq`
  (`q²F0`, `q²F8`, `q²GL/2` at `t0`),
* `Π = 2mn(n²p² − m²q²) + (m⁴ − n⁴)pq`,
* `G = 2S(p² + q²)` and `r` with `r² = G² − 64m²n²p²q²`; `G′ = S(p² + q²)` and `r′` with
  `r′² = G′² − (m² − n²)²(p² − q²)²`.

Contents.

* `fourPoint_pm_fraction`, `fourPoint_inv_fraction`: the squares `((Y1 − Y2)/(2pq))²` and
  `((Y1 − Y3)/(p² − q²))²` of the identities are the fractions `32m²n²Ψ0/(G + r)` and
  `2(m² − n²)²Ψ/(G′ + r′)`.
* `fourPoint_pm_dvd`, `fourPoint_inv_dvd`: a common divisor of numerator and denominator of these
  fractions divides `64m⁴n⁴`, respectively `16(m² − n²)⁴`. (The denominators divide
  `64m²n²p²q²` and `(m² − n²)²(p² − q²)²`; the theorems are stated for any common divisor of the
  numerator and of that multiple.)
* `fourPoint_res_psi0_dvd`, `fourPoint_res_psi_dvd`: a common divisor of `Π` and `Ψ0` divides
  `m²n²(m⁸ + 14m⁴n⁴ + n⁸)`, and a common divisor of `Π` and `Ψ` divides
  `(m² − n²)²(m⁸ + 14m⁴n⁴ + n⁸)`. These two numbers are the resultants used for the height of
  `x(2Q − R0) = 2Π²/(Ψ0Ψ)`; the octic is the one of the singular fibres of Result 1.

* `fourPoint_rank2_bound`, `fourPoint_rank3_bound`: the step from the heights and height pairings
  of the four points to a point `z` of bounded height, for rank at most 2 and for rank 3. They are
  statements about real numbers: the entries of a singular Gram matrix. `N₁, …, N₄` are the heights
  `N(a_i)`, all in `[L, L + e]`; `x = ⟨a₁, a₄⟩`, `y = ⟨a₂, a₃⟩`; the other four pairings are at
  most `κ` in absolute value. That the Gram matrix of dependent vectors has a kernel vector is the
  hypothesis; `N₁ + N₄ − 2|x|` is the height of `a₁ − a₄` or of `a₁ + a₄`.

Not in Lean: the three identities on the curve, the comparison of naive and canonical heights,
the size estimates, the step from these divisibilities to the bound for the common factor of `2Π²`
and `Ψ0Ψ`, ranks, generators, and the two listings. The polynomial cofactors below were computed
with Singular (through Sage) and are checked here by `ring`.
-/

namespace G1

/-- If `g` divides `K·pᵏ` and `K·qᵏ` with `p`, `q` coprime, then `g` divides `K`. -/
theorem dvd_of_dvd_mul_pow_of_isCoprime {g K p q : ℤ} (k : ℕ) (hpq : IsCoprime p q)
    (hp : g ∣ K * p ^ k) (hq : g ∣ K * q ^ k) : g ∣ K := by
  obtain ⟨a, b, hab⟩ : IsCoprime (p ^ k) (q ^ k) := hpq.pow
  have h : K = a * (K * p ^ k) + b * (K * q ^ k) := by linear_combination (-K) * hab
  rw [h]
  exact dvd_add (Dvd.dvd.mul_left hp a) (Dvd.dvd.mul_left hq b)

/-- The identity `(I±)` as a fraction: `((Y1 − Y2)/(2pq))² = 32m²n²Ψ0/(G + r)`, cleared of
denominators. `Y1² = q⁴F0·GL`, `Y2² = q⁴F0·GU`, `Y1·Y2 = Ψ0·r`. -/
theorem fourPoint_pm_fraction {m n p q r Y1 Y2 : ℤ}
    (h1 : Y1 ^ 2 = (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) * (2 * (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
          - 8 * m * n * p * q))
    (h2 : Y2 ^ 2 = (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) * (2 * (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
          + 8 * m * n * p * q))
    (h3 : Y1 * Y2 = (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) * r)
    (hr : r ^ 2 = (2 * (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)) ^ 2
          - 64 * m ^ 2 * n ^ 2 * p ^ 2 * q ^ 2) :
    (Y1 - Y2) ^ 2 * (2 * (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
          + r) = 4 * p ^ 2 * q ^ 2 * (32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2)) := by
  linear_combination (2 * (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) + r) * h1 + (2 * (m ^ 2
        + n ^ 2) * (p ^ 2 + q ^ 2) + r) * h2 - 2 * (2 * (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
        + r) * h3 - 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) * hr

/-- The identity `(Iinv)` as a fraction: `((Y1 − Y3)/(p² − q²))² = 2(m² − n²)²Ψ/(G′ + r′)`, cleared
of denominators. `Y1² = q⁴F0·GL`, `Y3² = q⁴F8·GL`, `Y1·Y3 = Ψ·r′`. -/
theorem fourPoint_inv_fraction {m n p q r Y1 Y3 : ℤ}
    (h1 : Y1 ^ 2 = (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) * (2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
          - 4 * m * n * p * q)))
    (h2 : Y3 ^ 2 = (n ^ 2 * q ^ 2 + m ^ 2 * p ^ 2) * (2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
          - 4 * m * n * p * q)))
    (h3 : Y1 * Y3 = ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q) * r)
    (hr : r ^ 2 = ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)) ^ 2 - (m ^ 2 - n ^ 2) ^ 2 * (p ^ 2
          - q ^ 2) ^ 2) :
    (Y1 - Y3) ^ 2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) + r) = (p ^ 2 - q ^ 2) ^ 2 * (2 * (m ^ 2
          - n ^ 2) ^ 2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q)) := by
  linear_combination ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) + r) * h1 + ((m ^ 2 + n ^ 2) * (p ^ 2
        + q ^ 2) + r) * h2 - 2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) + r) * h3 - 2 * ((m ^ 2
        + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q) * hr

/-- A common divisor of `32m²n²Ψ0` and `64m²n²p²q²` divides `64m⁴n⁴`. -/
theorem fourPoint_pm_dvd {m n p q g : ℤ} (hpq : IsCoprime p q)
    (h1 : g ∣ 32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2))
    (h2 : g ∣ 64 * m ^ 2 * n ^ 2 * p ^ 2 * q ^ 2) : g ∣ 64 * m ^ 4 * n ^ 4 := by
  refine dvd_of_dvd_mul_pow_of_isCoprime 4 hpq ?_ ?_
  · have h := dvd_sub (Dvd.dvd.mul_left h1 (2 * m ^ 2 * p ^ 2)) (Dvd.dvd.mul_left h2 (m ^ 4))
    have e : 2 * m ^ 2 * p ^ 2 * (32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2))
        - m ^ 4 * (64 * m ^ 2 * n ^ 2 * p ^ 2 * q ^ 2) = 64 * m ^ 4 * n ^ 4 * p ^ 4 := by ring
    rwa [e] at h
  · have h := dvd_sub (Dvd.dvd.mul_left h1 (2 * n ^ 2 * q ^ 2)) (Dvd.dvd.mul_left h2 (n ^ 4))
    have e : 2 * n ^ 2 * q ^ 2 * (32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2))
        - n ^ 4 * (64 * m ^ 2 * n ^ 2 * p ^ 2 * q ^ 2) = 64 * m ^ 4 * n ^ 4 * q ^ 4 := by ring
    rwa [e] at h

/-- A common divisor of `2(m² − n²)²Ψ` and `(m² − n²)²(p² − q²)²` divides `16(m² − n²)⁴`. -/
theorem fourPoint_inv_dvd {m n p q g : ℤ} (hpq : IsCoprime p q)
    (h1 : g ∣ 2 * (m ^ 2 - n ^ 2) ^ 2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q))
    (h2 : g ∣ (m ^ 2 - n ^ 2) ^ 2 * (p ^ 2 - q ^ 2) ^ 2) : g ∣ 16 * (m ^ 2 - n ^ 2) ^ 4 := by
  refine dvd_of_dvd_mul_pow_of_isCoprime 5 hpq ?_ ?_
  · have h := dvd_add (Dvd.dvd.mul_left h1 (6 * m ^ 2 * p ^ 3 + 6 * n ^ 2 * p ^ 3
        + 16 * m * n * p ^ 2 * q - 2 * m ^ 2 * p * q ^ 2 - 2 * n ^ 2 * p * q ^ 2
        - 8 * m * n * q ^ 3)) (Dvd.dvd.mul_left h2 (4 * m ^ 4 * p - 56 * m ^ 2 * n ^ 2 * p
        + 4 * n ^ 4 * p + 16 * m ^ 3 * n * q + 16 * m * n ^ 3 * q))
    have e : (6 * m ^ 2 * p ^ 3 + 6 * n ^ 2 * p ^ 3 + 16 * m * n * p ^ 2 * q
          - 2 * m ^ 2 * p * q ^ 2 - 2 * n ^ 2 * p * q ^ 2 - 8 * m * n * q ^ 3) * (2 * (m ^ 2
          - n ^ 2) ^ 2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q))
        + (4 * m ^ 4 * p - 56 * m ^ 2 * n ^ 2 * p + 4 * n ^ 4 * p + 16 * m ^ 3 * n * q
              + 16 * m * n ^ 3 * q) * ((m ^ 2 - n ^ 2) ^ 2 * (p ^ 2 - q ^ 2) ^ 2)
        = 16 * (m ^ 2 - n ^ 2) ^ 4 * p ^ 5 := by ring
    rwa [e] at h
  · have h := dvd_add (Dvd.dvd.mul_left h1 (-8 * m * n * p ^ 3 - 2 * m ^ 2 * p ^ 2 * q
        - 2 * n ^ 2 * p ^ 2 * q + 16 * m * n * p * q ^ 2 + 6 * m ^ 2 * q ^ 3
        + 6 * n ^ 2 * q ^ 3)) (Dvd.dvd.mul_left h2 (16 * m ^ 3 * n * p + 16 * m * n ^ 3 * p
        + 4 * m ^ 4 * q - 56 * m ^ 2 * n ^ 2 * q + 4 * n ^ 4 * q))
    have e : (-8 * m * n * p ^ 3 - 2 * m ^ 2 * p ^ 2 * q - 2 * n ^ 2 * p ^ 2 * q
          + 16 * m * n * p * q ^ 2 + 6 * m ^ 2 * q ^ 3 + 6 * n ^ 2 * q ^ 3) * (2 * (m ^ 2
          - n ^ 2) ^ 2 * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q))
        + (16 * m ^ 3 * n * p + 16 * m * n ^ 3 * p + 4 * m ^ 4 * q - 56 * m ^ 2 * n ^ 2 * q
              + 4 * n ^ 4 * q) * ((m ^ 2 - n ^ 2) ^ 2 * (p ^ 2 - q ^ 2) ^ 2)
        = 16 * (m ^ 2 - n ^ 2) ^ 4 * q ^ 5 := by ring
    rwa [e] at h

/-- A common divisor of `Π(p, q)` and `Ψ0(p, q)` divides their resultant `m²n²(m⁸ + 14m⁴n⁴
      + n⁸)`. -/
theorem fourPoint_res_psi0_dvd {m n p q g : ℤ} (hpq : IsCoprime p q)
    (h1 : g ∣ 2 * m * n * (n ^ 2 * p ^ 2 - m ^ 2 * q ^ 2) + (m ^ 4
          - n ^ 4) * p * q) (h2 : g ∣ m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) : g ∣ m ^ 2 * n ^ 2 * (m ^ 8
          + 14 * m ^ 4 * n ^ 4 + n ^ 8) := by
  refine dvd_of_dvd_mul_pow_of_isCoprime 3 hpq ?_ ?_
  · have h := dvd_add (Dvd.dvd.mul_left h1 (4 * m ^ 5 * n ^ 3 * p - m ^ 8 * q
        + m ^ 4 * n ^ 4 * q)) (Dvd.dvd.mul_left h2 (m ^ 10 * p + 6 * m ^ 6 * n ^ 4 * p
        + m ^ 2 * n ^ 8 * p - 2 * m ^ 9 * n * q + 2 * m ^ 5 * n ^ 5 * q))
    have e : (4 * m ^ 5 * n ^ 3 * p - m ^ 8 * q + m ^ 4 * n ^ 4 * q) * (2 * m * n * (n ^ 2 * p ^ 2
          - m ^ 2 * q ^ 2) + (m ^ 4 - n ^ 4) * p * q) + (m ^ 10 * p + 6 * m ^ 6 * n ^ 4 * p
          + m ^ 2 * n ^ 8 * p - 2 * m ^ 9 * n * q + 2 * m ^ 5 * n ^ 5 * q) * (m ^ 2 * q ^ 2
          + n ^ 2 * p ^ 2)
        = m ^ 2 * n ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) * p ^ 3 := by ring
    rwa [e] at h
  · have h := dvd_add (Dvd.dvd.mul_left h1 (-m ^ 4 * n ^ 4 * p + n ^ 8 * p
        - 4 * m ^ 3 * n ^ 5 * q)) (Dvd.dvd.mul_left h2 (2 * m ^ 5 * n ^ 5 * p - 2 * m * n ^ 9 * p
        + m ^ 8 * n ^ 2 * q + 6 * m ^ 4 * n ^ 6 * q + n ^ 10 * q))
    have e : (-m ^ 4 * n ^ 4 * p + n ^ 8 * p
          - 4 * m ^ 3 * n ^ 5 * q) * (2 * m * n * (n ^ 2 * p ^ 2 - m ^ 2 * q ^ 2) + (m ^ 4
          - n ^ 4) * p * q) + (2 * m ^ 5 * n ^ 5 * p - 2 * m * n ^ 9 * p + m ^ 8 * n ^ 2 * q
          + 6 * m ^ 4 * n ^ 6 * q + n ^ 10 * q) * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2)
        = m ^ 2 * n ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) * q ^ 3 := by ring
    rwa [e] at h

/-- A common divisor of `Π(p, q)` and `Ψ(p, q)` divides their resultant
`(m² − n²)²(m⁸ + 14m⁴n⁴ + n⁸)`. -/
theorem fourPoint_res_psi_dvd {m n p q g : ℤ} (hpq : IsCoprime p q)
    (h1 : g ∣ 2 * m * n * (n ^ 2 * p ^ 2 - m ^ 2 * q ^ 2) + (m ^ 4
          - n ^ 4) * p * q) (h2 : g ∣ (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2)
          - 4 * m * n * p * q) : g ∣ (m ^ 2 - n ^ 2) ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4
          + n ^ 8) := by
  refine dvd_of_dvd_mul_pow_of_isCoprime 3 hpq ?_ ?_
  · have h := dvd_add (Dvd.dvd.mul_left h1 (6 * m ^ 7 * n * p - 22 * m ^ 5 * n ^ 3 * p
        + 2 * m ^ 3 * n ^ 5 * p - 2 * m * n ^ 7 * p - m ^ 8 * q + 6 * m ^ 6 * n ^ 2 * q
        + 8 * m ^ 4 * n ^ 4 * q + 2 * m ^ 2 * n ^ 6 * q
        + n ^ 8 * q)) (Dvd.dvd.mul_left h2 (m ^ 10 * p - 3 * m ^ 8 * n ^ 2 * p
        + 6 * m ^ 6 * n ^ 4 * p + 10 * m ^ 4 * n ^ 6 * p + m ^ 2 * n ^ 8 * p + n ^ 10 * p
        - 2 * m ^ 9 * n * q + 14 * m ^ 7 * n ^ 3 * q + 2 * m ^ 5 * n ^ 5 * q
        + 2 * m ^ 3 * n ^ 7 * q))
    have e : (6 * m ^ 7 * n * p - 22 * m ^ 5 * n ^ 3 * p + 2 * m ^ 3 * n ^ 5 * p
          - 2 * m * n ^ 7 * p - m ^ 8 * q + 6 * m ^ 6 * n ^ 2 * q + 8 * m ^ 4 * n ^ 4 * q
          + 2 * m ^ 2 * n ^ 6 * q + n ^ 8 * q) * (2 * m * n * (n ^ 2 * p ^ 2 - m ^ 2 * q ^ 2)
          + (m ^ 4 - n ^ 4) * p * q) + (m ^ 10 * p - 3 * m ^ 8 * n ^ 2 * p + 6 * m ^ 6 * n ^ 4 * p
          + 10 * m ^ 4 * n ^ 6 * p + m ^ 2 * n ^ 8 * p + n ^ 10 * p - 2 * m ^ 9 * n * q
          + 14 * m ^ 7 * n ^ 3 * q + 2 * m ^ 5 * n ^ 5 * q + 2 * m ^ 3 * n ^ 7 * q) * ((m ^ 2
          + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q)
        = (m ^ 2 - n ^ 2) ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) * p ^ 3 := by ring
    rwa [e] at h
  · have h := dvd_add (Dvd.dvd.mul_left h1 (-m ^ 8 * p - 2 * m ^ 6 * n ^ 2 * p
        - 8 * m ^ 4 * n ^ 4 * p - 6 * m ^ 2 * n ^ 6 * p + n ^ 8 * p + 2 * m ^ 7 * n * q
        - 2 * m ^ 5 * n ^ 3 * q + 22 * m ^ 3 * n ^ 5 * q
        - 6 * m * n ^ 7 * q)) (Dvd.dvd.mul_left h2 (2 * m ^ 7 * n ^ 3 * p + 2 * m ^ 5 * n ^ 5 * p
        + 14 * m ^ 3 * n ^ 7 * p - 2 * m * n ^ 9 * p + m ^ 10 * q + m ^ 8 * n ^ 2 * q
        + 10 * m ^ 6 * n ^ 4 * q + 6 * m ^ 4 * n ^ 6 * q - 3 * m ^ 2 * n ^ 8 * q + n ^ 10 * q))
    have e : (-m ^ 8 * p - 2 * m ^ 6 * n ^ 2 * p - 8 * m ^ 4 * n ^ 4 * p - 6 * m ^ 2 * n ^ 6 * p
          + n ^ 8 * p + 2 * m ^ 7 * n * q - 2 * m ^ 5 * n ^ 3 * q + 22 * m ^ 3 * n ^ 5 * q
          - 6 * m * n ^ 7 * q) * (2 * m * n * (n ^ 2 * p ^ 2 - m ^ 2 * q ^ 2) + (m ^ 4
          - n ^ 4) * p * q) + (2 * m ^ 7 * n ^ 3 * p + 2 * m ^ 5 * n ^ 5 * p
          + 14 * m ^ 3 * n ^ 7 * p - 2 * m * n ^ 9 * p + m ^ 10 * q + m ^ 8 * n ^ 2 * q
          + 10 * m ^ 6 * n ^ 4 * q + 6 * m ^ 4 * n ^ 6 * q - 3 * m ^ 2 * n ^ 8 * q
          + n ^ 10 * q) * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q)
        = (m ^ 2 - n ^ 2) ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8) * q ^ 3 := by ring
    rwa [e] at h


/-! ### The lattice step -/

/-- Rank at most 2. If the Gram matrix of `a₁, a₂, a₄` has a nonzero kernel vector, the pairings
`⟨a₁, a₂⟩`, `⟨a₂, a₄⟩` are at most `κ` in absolute value and the three norms lie in `[L, L + e]`,
then `(N₁ + N₄ − 2|x|)·L ≤ e²/2 + 4κ²`, where `x = ⟨a₁, a₄⟩`. -/
theorem fourPoint_rank2_bound {N1 N2 N4 g12 g24 x κ L e c1 c2 c4 : ℝ}
    (hc : c1 ≠ 0 ∨ c2 ≠ 0 ∨ c4 ≠ 0)
    (h1 : N1 * c1 + g12 * c2 + x * c4 = 0)
    (h2 : g12 * c1 + N2 * c2 + g24 * c4 = 0)
    (h3 : x * c1 + g24 * c2 + N4 * c4 = 0)
    (hg12 : |g12| ≤ κ) (hg24 : |g24| ≤ κ) (hL : 0 < L)
    (hN1 : L ≤ N1) (hN1' : N1 ≤ L + e) (hN2 : L ≤ N2) (hN4 : L ≤ N4) (hN4' : N4 ≤ L + e) :
    (N1 + N4 - 2 * |x|) * L ≤ e ^ 2 / 2 + 4 * κ ^ 2 := by
  have hκ : 0 ≤ κ := (abs_nonneg g12).trans hg12
  -- the determinant vanishes
  have hdet : N1 * N2 * N4 + 2 * g12 * g24 * x - N1 * g24 ^ 2 - N2 * x ^ 2 - N4 * g12 ^ 2 = 0 := by
    have d1 : (N1 * N2 * N4 + 2 * g12 * g24 * x - N1 * g24 ^ 2 - N2 * x ^ 2 - N4 * g12 ^ 2) * c1
        = 0 := by
      linear_combination (N2 * N4 - g24 ^ 2) * h1 - (g12 * N4 - x * g24) * h2
        + (g12 * g24 - x * N2) * h3
    have d2 : (N1 * N2 * N4 + 2 * g12 * g24 * x - N1 * g24 ^ 2 - N2 * x ^ 2 - N4 * g12 ^ 2) * c2
        = 0 := by
      linear_combination (-(g12 * N4 - g24 * x)) * h1 + (N1 * N4 - x ^ 2) * h2
        - (N1 * g24 - g12 * x) * h3
    have d4 : (N1 * N2 * N4 + 2 * g12 * g24 * x - N1 * g24 ^ 2 - N2 * x ^ 2 - N4 * g12 ^ 2) * c4
        = 0 := by
      linear_combination (g12 * g24 - N2 * x) * h1 - (N1 * g24 - x * g12) * h2
        + (N1 * N2 - g12 ^ 2) * h3
    rcases hc with hc | hc | hc
    · exact (mul_eq_zero.mp d1).resolve_right hc
    · exact (mul_eq_zero.mp d2).resolve_right hc
    · exact (mul_eq_zero.mp d4).resolve_right hc
  have g12sq : g12 ^ 2 ≤ κ ^ 2 := by
    rw [← sq_abs g12]; exact pow_le_pow_left₀ (abs_nonneg _) hg12 2
  have g24sq : g24 ^ 2 ≤ κ ^ 2 := by
    rw [← sq_abs g24]; exact pow_le_pow_left₀ (abs_nonneg _) hg24 2
  have cr : -(2 * g12 * g24 * x) ≤ 2 * κ ^ 2 * |x| := by
    have h : |g12 * g24 * x| ≤ κ * κ * |x| := by
      rw [abs_mul, abs_mul]
      exact mul_le_mul_of_nonneg_right (mul_le_mul hg12 hg24 (abs_nonneg _) hκ) (abs_nonneg _)
    have h' := neg_abs_le (g12 * g24 * x)
    nlinarith
  have hN1pos : 0 ≤ N1 := hL.le.trans hN1
  have hN4pos : 0 ≤ N4 := hL.le.trans hN4
  have hN2pos : 0 < N2 := hL.trans_le hN2
  have key : N2 * (N1 * N4 - x ^ 2) ≤ κ ^ 2 * (N1 + N4 + 2 * |x|) := by
    have a := mul_le_mul_of_nonneg_left g24sq hN1pos
    have b := mul_le_mul_of_nonneg_left g12sq hN4pos
    nlinarith
  have hT : 2 * L ≤ N1 + N4 + 2 * |x| := by
    have := abs_nonneg x
    linarith
  have hTpos : 0 < N1 + N4 + 2 * |x| := by linarith
  have mT : (N1 + N4 - 2 * |x|) * (N1 + N4 + 2 * |x|)
      = (N1 - N4) ^ 2 + 4 * (N1 * N4 - x ^ 2) := by
    have h := sq_abs x
    calc (N1 + N4 - 2 * |x|) * (N1 + N4 + 2 * |x|) = (N1 + N4) ^ 2 - 4 * |x| ^ 2 := by ring
      _ = (N1 - N4) ^ 2 + 4 * (N1 * N4 - x ^ 2) := by rw [h]; ring
  have hd : (N1 - N4) ^ 2 ≤ e ^ 2 := by
    have h1' : N1 - N4 ≤ e := by linarith
    have h2' : -e ≤ N1 - N4 := by linarith
    nlinarith
  have hA : (N1 + N4 - 2 * |x|) * (N1 + N4 + 2 * |x|) * N2
      ≤ e ^ 2 * N2 + 4 * κ ^ 2 * (N1 + N4 + 2 * |x|) := by
    rw [mT]
    nlinarith [mul_le_mul_of_nonneg_right hd hN2pos.le]
  have hpos : 0 < (N1 + N4 + 2 * |x|) * N2 := mul_pos hTpos hN2pos
  refine le_of_mul_le_mul_right ?_ hpos
  have hB := mul_le_mul_of_nonneg_right hA hL.le
  have hC : e ^ 2 * N2 * (2 * L) ≤ e ^ 2 * N2 * (N1 + N4 + 2 * |x|) :=
    mul_le_mul_of_nonneg_left hT (mul_nonneg (sq_nonneg e) hN2pos.le)
  have hD : 4 * κ ^ 2 * (N1 + N4 + 2 * |x|) * L ≤ 4 * κ ^ 2 * (N1 + N4 + 2 * |x|) * N2 :=
    mul_le_mul_of_nonneg_left hN2 (mul_nonneg (by positivity) hTpos.le)
  nlinarith

/-- A positive definite binary form: if `0 < p` and `x² < p·r`, then `p·a² + 2x·a·b + r·b²` is
non-negative and vanishes only at `a = b = 0`. -/
theorem binary_form_pos {p r x a b : ℝ} (hp : 0 < p) (h : x ^ 2 < p * r) :
    0 ≤ p * a ^ 2 + 2 * x * a * b + r * b ^ 2 ∧
      (p * a ^ 2 + 2 * x * a * b + r * b ^ 2 = 0 → a = 0 ∧ b = 0) := by
  have e : p * (p * a ^ 2 + 2 * x * a * b + r * b ^ 2)
      = (p * a + x * b) ^ 2 + (p * r - x ^ 2) * b ^ 2 := by ring
  have hpos : 0 < p * r - x ^ 2 := by linarith
  have h0 : 0 ≤ p * (p * a ^ 2 + 2 * x * a * b + r * b ^ 2) := by
    rw [e]; exact add_nonneg (sq_nonneg _) (mul_nonneg hpos.le (sq_nonneg _))
  refine ⟨nonneg_of_mul_nonneg_right h0 hp, fun hz => ?_⟩
  rw [hz, mul_zero] at e
  have hb2 : (p * r - x ^ 2) * b ^ 2 = 0 := by
    have := sq_nonneg (p * a + x * b)
    have := mul_nonneg hpos.le (sq_nonneg b)
    linarith
  have hb : b = 0 := by
    rcases mul_eq_zero.mp hb2 with h' | h'
    · exact absurd h' hpos.ne'
    · exact (pow_eq_zero_iff two_ne_zero).mp h'
  have ha2 : (p * a + x * b) ^ 2 = 0 := by linarith
  have ha : p * a = 0 := by
    have h' := (pow_eq_zero_iff two_ne_zero).mp ha2
    rw [hb] at h'
    linarith
  exact ⟨(mul_eq_zero.mp ha).resolve_left hp.ne', hb⟩

/-- The estimate for one block: if `(N − 2κ)(N' − 2κ) ≤ z²` with `N, N'` in `[L, L + e]` and
`2κ < L`, then `(N + N' − 2|z|)·4(L − κ) ≤ e² + 16κ(L + e) − 16κ²`. -/
theorem fourPoint_block_bound {N N' z κ L e : ℝ} (hκ : 0 ≤ κ) (hκL : 2 * κ < L)
    (hN : L ≤ N) (hN' : N ≤ L + e) (hM : L ≤ N') (hM' : N' ≤ L + e)
    (hz : (N - 2 * κ) * (N' - 2 * κ) ≤ z ^ 2) :
    (N + N' - 2 * |z|) * (4 * (L - κ)) ≤ e ^ 2 + 16 * κ * (L + e) - 16 * κ ^ 2 := by
  have hz0 := abs_nonneg z
  have hzL : L - 2 * κ ≤ |z| := by
    by_contra hlt
    have hlt' := not_le.mp hlt
    have h1' : (L - 2 * κ) * (L - 2 * κ) ≤ (N - 2 * κ) * (N' - 2 * κ) :=
      mul_le_mul (by linarith) (by linarith) (by linarith) (by linarith)
    have h2' : |z| * |z| < (L - 2 * κ) * (L - 2 * κ) := mul_self_lt_mul_self hz0 hlt'
    have h3' := sq_abs z
    nlinarith
  have mT : (N + N' - 2 * |z|) * (N + N' + 2 * |z|)
      = (N - N') ^ 2 + 4 * (N * N' - z ^ 2) := by
    have h := sq_abs z
    calc (N + N' - 2 * |z|) * (N + N' + 2 * |z|) = (N + N') ^ 2 - 4 * |z| ^ 2 := by ring
      _ = (N - N') ^ 2 + 4 * (N * N' - z ^ 2) := by rw [h]; ring
  have hd : (N - N') ^ 2 ≤ e ^ 2 := by
    have h1' : N - N' ≤ e := by linarith
    have h2' : -e ≤ N - N' := by linarith
    nlinarith
  have hR : (N + N' - 2 * |z|) * (N + N' + 2 * |z|)
      ≤ e ^ 2 + 16 * κ * (L + e) - 16 * κ ^ 2 := by
    rw [mT]
    have hs' : κ * (N + N') ≤ κ * (2 * (L + e)) :=
      mul_le_mul_of_nonneg_left (by linarith) hκ
    nlinarith
  have hT : 4 * (L - κ) ≤ N + N' + 2 * |z| := by linarith
  have hTpos : 0 < 4 * (L - κ) := by linarith
  by_cases hm : 0 ≤ N + N' - 2 * |z|
  · exact (mul_le_mul_of_nonneg_left hT hm).trans hR
  · have h1' : (N + N' - 2 * |z|) * (4 * (L - κ)) ≤ 0 :=
      mul_nonpos_of_nonpos_of_nonneg (not_le.mp hm).le hTpos.le
    have h2' : 0 ≤ κ * (L + e - κ) := mul_nonneg hκ (by linarith)
    nlinarith [sq_nonneg e]

/-- `(|a| + |b|)² ≤ 2(a² + b²)`. -/
theorem abs_add_abs_sq_le (a b : ℝ) : (|a| + |b|) ^ 2 ≤ 2 * (a ^ 2 + b ^ 2) := by
  nlinarith [sq_nonneg (|a| - |b|), sq_abs a, sq_abs b]

/-- A bilinear form with coefficients at most `κ` in absolute value. -/
theorem bilinear_abs_le {g12 g13 g42 g43 κ u1 u4 v2 v3 : ℝ}
    (hg12 : |g12| ≤ κ) (hg13 : |g13| ≤ κ) (hg42 : |g42| ≤ κ) (hg43 : |g43| ≤ κ) :
    |g12 * u1 * v2 + g13 * u1 * v3 + g42 * u4 * v2 + g43 * u4 * v3|
      ≤ κ * ((|u1| + |u4|) * (|v2| + |v3|)) := by
  have t1 : |g12 * u1 * v2| ≤ κ * (|u1| * |v2|) := by
    rw [abs_mul, abs_mul, mul_assoc]
    exact mul_le_mul_of_nonneg_right hg12 (by positivity)
  have t2 : |g13 * u1 * v3| ≤ κ * (|u1| * |v3|) := by
    rw [abs_mul, abs_mul, mul_assoc]
    exact mul_le_mul_of_nonneg_right hg13 (by positivity)
  have t3 : |g42 * u4 * v2| ≤ κ * (|u4| * |v2|) := by
    rw [abs_mul, abs_mul, mul_assoc]
    exact mul_le_mul_of_nonneg_right hg42 (by positivity)
  have t4 : |g43 * u4 * v3| ≤ κ * (|u4| * |v3|) := by
    rw [abs_mul, abs_mul, mul_assoc]
    exact mul_le_mul_of_nonneg_right hg43 (by positivity)
  have a1 := abs_add_le (g12 * u1 * v2 + g13 * u1 * v3 + g42 * u4 * v2) (g43 * u4 * v3)
  have a2 := abs_add_le (g12 * u1 * v2 + g13 * u1 * v3) (g42 * u4 * v2)
  have a3 := abs_add_le (g12 * u1 * v2) (g13 * u1 * v3)
  have e' : κ * ((|u1| + |u4|) * (|v2| + |v3|))
      = κ * (|u1| * |v2|) + κ * (|u1| * |v3|) + κ * (|u4| * |v2|) + κ * (|u4| * |v3|) := by
    ring
  rw [e']
  linarith

/-- If `s = 2κ·Pu + dA = 2κ·Pv + dD` with all terms non-negative and `s² ≤ 4κ²·Pu·Pv`, then
`dA = 0` or `dD = 0`. -/
theorem excess_eq_zero {s κ Pu Pv dA dD : ℝ} (hκ : 0 ≤ κ) (hPu : 0 ≤ Pu) (hPv : 0 ≤ Pv)
    (hdA : 0 ≤ dA) (hdD : 0 ≤ dD) (s1 : s = 2 * κ * Pu + dA) (s2 : s = 2 * κ * Pv + dD)
    (ssq : s ^ 2 ≤ 4 * κ ^ 2 * (Pu * Pv)) : dA = 0 ∨ dD = 0 := by
  have h2κ : 0 ≤ 2 * κ := by linarith
  have e1 : s * s = (2 * κ * Pu + dA) * (2 * κ * Pv + dD) := by rw [← s1, ← s2]
  have t1 : 0 ≤ 2 * κ * Pu * dD := mul_nonneg (mul_nonneg h2κ hPu) hdD
  have t2 : 0 ≤ 2 * κ * Pv * dA := mul_nonneg (mul_nonneg h2κ hPv) hdA
  have hprod0 : dA * dD ≤ 0 := by nlinarith
  by_contra hne
  obtain ⟨hne1, hne2⟩ := not_or.mp hne
  have hpos' : 0 < dA * dD :=
    mul_pos (lt_of_le_of_ne hdA (Ne.symm hne1)) (lt_of_le_of_ne hdD (Ne.symm hne2))
  linarith

/-- The two diagonal blocks of a singular Gram matrix cannot both stay positive definite after
subtracting `2κ`, when the entries of the off-diagonal block are at most `κ` in absolute value. -/
theorem fourPoint_blocks_not_both_definite
    {N1 N2 N3 N4 g12 g13 g42 g43 x y κ u1 u4 v2 v3 : ℝ}
    (hw : u1 ≠ 0 ∨ u4 ≠ 0 ∨ v2 ≠ 0 ∨ v3 ≠ 0)
    (r1 : N1 * u1 + x * u4 + g12 * v2 + g13 * v3 = 0)
    (r4 : x * u1 + N4 * u4 + g42 * v2 + g43 * v3 = 0)
    (r2 : g12 * u1 + g42 * u4 + N2 * v2 + y * v3 = 0)
    (r3 : g13 * u1 + g43 * u4 + y * v2 + N3 * v3 = 0)
    (hg12 : |g12| ≤ κ) (hg13 : |g13| ≤ κ) (hg42 : |g42| ≤ κ) (hg43 : |g43| ≤ κ)
    (hp1 : 0 < N1 - 2 * κ) (hp2 : 0 < N2 - 2 * κ)
    (hx : x ^ 2 < (N1 - 2 * κ) * (N4 - 2 * κ)) (hy : y ^ 2 < (N2 - 2 * κ) * (N3 - 2 * κ)) :
    False := by
  have hκ : 0 ≤ κ := (abs_nonneg g12).trans hg12
  have h2κ : 0 ≤ 2 * κ := by linarith
  obtain ⟨dA0, dAz⟩ := binary_form_pos (a := u1) (b := u4) hp1 hx
  obtain ⟨dD0, dDz⟩ := binary_form_pos (a := v2) (b := v3) hp2 hy
  -- `s = uᵀAu = vᵀDv = −uᵀBv`
  have sA : N1 * u1 ^ 2 + 2 * x * u1 * u4 + N4 * u4 ^ 2
      = -(g12 * u1 * v2 + g13 * u1 * v3 + g42 * u4 * v2 + g43 * u4 * v3) := by
    linear_combination u1 * r1 + u4 * r4
  have sD : N2 * v2 ^ 2 + 2 * y * v2 * v3 + N3 * v3 ^ 2
      = -(g12 * u1 * v2 + g13 * u1 * v3 + g42 * u4 * v2 + g43 * u4 * v3) := by
    linear_combination v2 * r2 + v3 * r3
  have hB := bilinear_abs_le (u1 := u1) (u4 := u4) (v2 := v2) (v3 := v3) hg12 hg13 hg42 hg43
  have hPu0 : 0 ≤ u1 ^ 2 + u4 ^ 2 := by positivity
  have hPv0 : 0 ≤ v2 ^ 2 + v3 ^ 2 := by positivity
  have s1 : N1 * u1 ^ 2 + 2 * x * u1 * u4 + N4 * u4 ^ 2 = 2 * κ * (u1 ^ 2 + u4 ^ 2)
      + ((N1 - 2 * κ) * u1 ^ 2 + 2 * x * u1 * u4 + (N4 - 2 * κ) * u4 ^ 2) := by ring
  have s2 : N1 * u1 ^ 2 + 2 * x * u1 * u4 + N4 * u4 ^ 2 = 2 * κ * (v2 ^ 2 + v3 ^ 2)
      + ((N2 - 2 * κ) * v2 ^ 2 + 2 * y * v2 * v3 + (N3 - 2 * κ) * v3 ^ 2) := by
    rw [sA, ← sD]; ring
  have ssq : (N1 * u1 ^ 2 + 2 * x * u1 * u4 + N4 * u4 ^ 2) ^ 2
      ≤ 4 * κ ^ 2 * ((u1 ^ 2 + u4 ^ 2) * (v2 ^ 2 + v3 ^ 2)) := by
    have habs : |N1 * u1 ^ 2 + 2 * x * u1 * u4 + N4 * u4 ^ 2|
        ≤ κ * ((|u1| + |u4|) * (|v2| + |v3|)) := by
      rw [sA, abs_neg]; exact hB
    have hsq := pow_le_pow_left₀ (abs_nonneg _) habs 2
    rw [sq_abs] at hsq
    have hprod : (|u1| + |u4|) ^ 2 * (|v2| + |v3|) ^ 2
        ≤ 2 * (u1 ^ 2 + u4 ^ 2) * (2 * (v2 ^ 2 + v3 ^ 2)) :=
      mul_le_mul (abs_add_abs_sq_le u1 u4) (abs_add_abs_sq_le v2 v3) (by positivity)
        (mul_nonneg (by norm_num) hPu0)
    calc (N1 * u1 ^ 2 + 2 * x * u1 * u4 + N4 * u4 ^ 2) ^ 2
        ≤ (κ * ((|u1| + |u4|) * (|v2| + |v3|))) ^ 2 := hsq
      _ = κ ^ 2 * ((|u1| + |u4|) ^ 2 * (|v2| + |v3|) ^ 2) := by ring
      _ ≤ κ ^ 2 * (2 * (u1 ^ 2 + u4 ^ 2) * (2 * (v2 ^ 2 + v3 ^ 2))) :=
          mul_le_mul_of_nonneg_left hprod (sq_nonneg κ)
      _ = 4 * κ ^ 2 * ((u1 ^ 2 + u4 ^ 2) * (v2 ^ 2 + v3 ^ 2)) := by ring
  have hall : (u1 = 0 ∧ u4 = 0) ∧ (v2 = 0 ∧ v3 = 0) := by
    rcases excess_eq_zero hκ hPu0 hPv0 dA0 dD0 s1 s2 ssq with hz | hz
    · have hu := dAz hz
      have hdD0 : (N2 - 2 * κ) * v2 ^ 2 + 2 * y * v2 * v3 + (N3 - 2 * κ) * v3 ^ 2 = 0 := by
        have hk : 2 * κ * (u1 ^ 2 + u4 ^ 2) = 0 := by rw [hu.1, hu.2]; ring
        have t : 0 ≤ 2 * κ * (v2 ^ 2 + v3 ^ 2) := mul_nonneg h2κ hPv0
        refine le_antisymm ?_ dD0
        linarith
      exact ⟨hu, dDz hdD0⟩
    · have hv := dDz hz
      have hdA0 : (N1 - 2 * κ) * u1 ^ 2 + 2 * x * u1 * u4 + (N4 - 2 * κ) * u4 ^ 2 = 0 := by
        have hk : 2 * κ * (v2 ^ 2 + v3 ^ 2) = 0 := by rw [hv.1, hv.2]; ring
        have t : 0 ≤ 2 * κ * (u1 ^ 2 + u4 ^ 2) := mul_nonneg h2κ hPu0
        refine le_antisymm ?_ dA0
        linarith
      exact ⟨dAz hdA0, hv⟩
  rcases hw with h | h | h | h
  · exact h hall.1.1
  · exact h hall.1.2
  · exact h hall.2.1
  · exact h hall.2.2

/-- Rank 3. If the Gram matrix of `a₁, a₄, a₂, a₃` has a nonzero kernel vector, the four pairings
between `{a₁, a₄}` and `{a₂, a₃}` are at most `κ` in absolute value, the four norms lie in
`[L, L + e]` and `2κ < L`, then `(N₁ + N₄ − 2|x|)·4(L − κ) ≤ e² + 16κ(L + e) − 16κ²` or the same
holds for `N₂, N₃, y`; here `x = ⟨a₁, a₄⟩`, `y = ⟨a₂, a₃⟩`. -/
theorem fourPoint_rank3_bound {N1 N2 N3 N4 g12 g13 g42 g43 x y κ L e u1 u4 v2 v3 : ℝ}
    (hw : u1 ≠ 0 ∨ u4 ≠ 0 ∨ v2 ≠ 0 ∨ v3 ≠ 0)
    (r1 : N1 * u1 + x * u4 + g12 * v2 + g13 * v3 = 0)
    (r4 : x * u1 + N4 * u4 + g42 * v2 + g43 * v3 = 0)
    (r2 : g12 * u1 + g42 * u4 + N2 * v2 + y * v3 = 0)
    (r3 : g13 * u1 + g43 * u4 + y * v2 + N3 * v3 = 0)
    (hg12 : |g12| ≤ κ) (hg13 : |g13| ≤ κ) (hg42 : |g42| ≤ κ) (hg43 : |g43| ≤ κ)
    (hκL : 2 * κ < L)
    (hN1 : L ≤ N1) (hN1' : N1 ≤ L + e) (hN2 : L ≤ N2) (hN2' : N2 ≤ L + e)
    (hN3 : L ≤ N3) (hN3' : N3 ≤ L + e) (hN4 : L ≤ N4) (hN4' : N4 ≤ L + e) :
    (N1 + N4 - 2 * |x|) * (4 * (L - κ)) ≤ e ^ 2 + 16 * κ * (L + e) - 16 * κ ^ 2 ∨
      (N2 + N3 - 2 * |y|) * (4 * (L - κ)) ≤ e ^ 2 + 16 * κ * (L + e) - 16 * κ ^ 2 := by
  have hκ : 0 ≤ κ := (abs_nonneg g12).trans hg12
  by_cases hx : (N1 - 2 * κ) * (N4 - 2 * κ) ≤ x ^ 2
  · exact Or.inl (fourPoint_block_bound hκ hκL hN1 hN1' hN4 hN4' hx)
  by_cases hy : (N2 - 2 * κ) * (N3 - 2 * κ) ≤ y ^ 2
  · exact Or.inr (fourPoint_block_bound hκ hκL hN2 hN2' hN3 hN3' hy)
  exact (fourPoint_blocks_not_both_definite hw r1 r4 r2 r3 hg12 hg13 hg42 hg43
    (by linarith) (by linarith) (not_le.mp hx) (not_le.mp hy)).elim

/-! ## The common factor in Lemma C and the lower bounds of Lemmas A and B

`FOUR_POINTS.md`, Section 3. The fractions of the three identities are not in lowest terms; these
theorems bound what can cancel, and give the lower bounds for the heights that follow. -/

/-- A common divisor of `2a²` and `bc` divides `2·gcd(a, b)²·gcd(a, c)²` (natural numbers). -/
theorem gcd_two_sq_mul_dvd (a b c : ℕ) :
    Nat.gcd (2 * a ^ 2) (b * c) ∣ 2 * Nat.gcd a b ^ 2 * Nat.gcd a c ^ 2 := by
  have sq : ∀ d : ℕ, Nat.gcd (a ^ 2) d ∣ Nat.gcd a d ^ 2 := by
    intro d
    have h : Nat.gcd d (a * a) ∣ Nat.gcd d a * Nat.gcd d a := Nat.gcd_mul_right_dvd_mul_gcd d a a
    rw [Nat.gcd_comm (a ^ 2) d, Nat.gcd_comm a d, pow_two, pow_two]
    exact h
  have h1 : Nat.gcd (2 * a ^ 2) (b * c) ∣ Nat.gcd (b * c) 2 * Nat.gcd (b * c) (a ^ 2) := by
    rw [Nat.gcd_comm (2 * a ^ 2) (b * c)]
    exact Nat.gcd_mul_right_dvd_mul_gcd (b * c) 2 (a ^ 2)
  have h2 : Nat.gcd (b * c) 2 ∣ 2 := Nat.gcd_dvd_right _ _
  have h3 : Nat.gcd (b * c) (a ^ 2) ∣ Nat.gcd a b ^ 2 * Nat.gcd a c ^ 2 := by
    rw [Nat.gcd_comm (b * c) (a ^ 2)]
    exact (Nat.gcd_mul_right_dvd_mul_gcd (a ^ 2) b c).trans (Nat.mul_dvd_mul (sq b) (sq c))
  have h4 : Nat.gcd (b * c) 2 * Nat.gcd (b * c) (a ^ 2)
      ∣ 2 * (Nat.gcd a b ^ 2 * Nat.gcd a c ^ 2) := Nat.mul_dvd_mul h2 h3
  rw [mul_assoc]
  exact h1.trans h4

/-- Integers: if every common divisor of `P` and `A` divides `R0` and every common divisor of `P`
and `B` divides `R1`, then a common divisor of `2P²` and `AB` divides `2·R0²·R1²`. -/
theorem dvd_two_sq_mul_of_common {g P A B R0 R1 : ℤ} (h1 : g ∣ 2 * P ^ 2) (h2 : g ∣ A * B)
    (hA : ∀ d : ℤ, d ∣ P → d ∣ A → d ∣ R0) (hB : ∀ d : ℤ, d ∣ P → d ∣ B → d ∣ R1) :
    g ∣ 2 * R0 ^ 2 * R1 ^ 2 := by
  have hg : g ∣ ((Int.gcd (2 * P ^ 2) (A * B) : ℕ) : ℤ) := Int.dvd_coe_gcd h1 h2
  have e1 : (2 * P ^ 2).natAbs = 2 * P.natAbs ^ 2 := by
    rw [Int.natAbs_mul, Int.natAbs_pow]; rfl
  have e2 : (A * B).natAbs = A.natAbs * B.natAbs := Int.natAbs_mul A B
  have hnat : Int.gcd (2 * P ^ 2) (A * B) ∣ 2 * Int.gcd P A ^ 2 * Int.gcd P B ^ 2 := by
    show Nat.gcd (2 * P ^ 2).natAbs (A * B).natAbs
      ∣ 2 * Nat.gcd P.natAbs A.natAbs ^ 2 * Nat.gcd P.natAbs B.natAbs ^ 2
    rw [e1, e2]
    exact gcd_two_sq_mul_dvd _ _ _
  have hA' : ((Int.gcd P A : ℕ) : ℤ) ∣ R0 :=
    hA _ (Int.gcd_dvd_left (a := P) (b := A)) (Int.gcd_dvd_right (a := P) (b := A))
  have hB' : ((Int.gcd P B : ℕ) : ℤ) ∣ R1 :=
    hB _ (Int.gcd_dvd_left (a := P) (b := B)) (Int.gcd_dvd_right (a := P) (b := B))
  have hcast : ((Int.gcd (2 * P ^ 2) (A * B) : ℕ) : ℤ)
      ∣ 2 * ((Int.gcd P A : ℕ) : ℤ) ^ 2 * ((Int.gcd P B : ℕ) : ℤ) ^ 2 := by
    have h := Int.natCast_dvd_natCast.mpr hnat
    push_cast at h
    exact h
  exact hg.trans (hcast.trans (mul_dvd_mul (mul_dvd_mul (dvd_refl 2) (pow_dvd_pow_of_dvd hA' 2))
    (pow_dvd_pow_of_dvd hB' 2)))

/-- **Lemma C, the common factor.** For coprime `p`, `q`, a common divisor of `2Π²` and `Ψ0·Ψ`
divides `2·R0²·R1²` with `R0 = m²n²(m⁸ + 14m⁴n⁴ + n⁸)`, `R1 = (m² − n²)²(m⁸ + 14m⁴n⁴ + n⁸)`. -/
theorem fourPoint_lemmaC_dvd {m n p q g : ℤ} (hpq : IsCoprime p q)
    (h1 : g ∣ 2 * (2 * m * n * (n ^ 2 * p ^ 2 - m ^ 2 * q ^ 2) + (m ^ 4 - n ^ 4) * p * q) ^ 2)
    (h2 : g ∣ (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2)
        * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q)) :
    g ∣ 2 * (m ^ 2 * n ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8)) ^ 2
        * ((m ^ 2 - n ^ 2) ^ 2 * (m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8)) ^ 2 :=
  dvd_two_sq_mul_of_common h1 h2 (fun _ a b => fourPoint_res_psi0_dvd hpq a b)
    (fun _ a b => fourPoint_res_psi_dvd hpq a b)

/-- **Lemma A, the lower bound.** If `g > 0` divides `64m⁴n⁴` (as every common factor does, by
`fourPoint_pm_dvd`) and `N` is the numerator `32m²n²Ψ0` divided by `g`, then
`max(p², q²) ≤ 2m²·N`. Here `m > n > 0`. -/
theorem fourPoint_pm_lower {m n p q g N : ℤ} (hn : 0 < n) (hnm : n < m) (hg : 0 < g)
    (hgd : g ∣ 64 * m ^ 4 * n ^ 4)
    (hN : g * N = 32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2)) :
    p ^ 2 ≤ 2 * m ^ 2 * N ∧ q ^ 2 ≤ 2 * m ^ 2 * N := by
  have hm : 0 < m := lt_trans hn hnm
  have hK : 0 < 64 * m ^ 4 * n ^ 4 := by positivity
  have hgle : g ≤ 64 * m ^ 4 * n ^ 4 := Int.le_of_dvd hK hgd
  have hR : 0 ≤ 32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2 + n ^ 2 * p ^ 2) := by positivity
  have hN0 : 0 ≤ N := by
    rcases le_or_gt 0 N with h | h
    · exact h
    · exfalso
      have := mul_pos hg (neg_pos.mpr h)
      linarith
  have h64 : g * N ≤ 64 * m ^ 4 * n ^ 4 * N := mul_le_mul_of_nonneg_right hgle hN0
  have hc : 0 < 32 * m ^ 2 * n ^ 4 := by positivity
  have hmn : 0 ≤ m ^ 2 - n ^ 2 := by nlinarith
  constructor
  · have h0 : 0 ≤ 32 * m ^ 2 * n ^ 2 * (m ^ 2 * q ^ 2) := by positivity
    have h : 32 * m ^ 2 * n ^ 4 * p ^ 2 ≤ 32 * m ^ 2 * n ^ 4 * (2 * m ^ 2 * N) := by linarith
    exact le_of_mul_le_mul_left h hc
  · have h0 : 0 ≤ 32 * m ^ 2 * n ^ 2 * (n ^ 2 * p ^ 2) := by positivity
    have h1 : 0 ≤ 32 * m ^ 2 * n ^ 2 * q ^ 2 * (m ^ 2 - n ^ 2) :=
      mul_nonneg (by positivity) hmn
    have h : 32 * m ^ 2 * n ^ 4 * q ^ 2 ≤ 32 * m ^ 2 * n ^ 4 * (2 * m ^ 2 * N) := by linarith
    exact le_of_mul_le_mul_left h hc

/-- **Lemma B, the lower bound.** If `g > 0` divides `16(m² − n²)⁴` (as every common factor does,
by `fourPoint_inv_dvd`) and `N` is the numerator `2(m² − n²)²Ψ` divided by `g`, then
`max(p², q²) ≤ 8(m + n)²·N`. Here `m > n > 0`. -/
theorem fourPoint_inv_lower {m n p q g N : ℤ} (hn : 0 < n) (hnm : n < m) (hg : 0 < g)
    (hgd : g ∣ 16 * (m ^ 2 - n ^ 2) ^ 4)
    (hN : g * N = 2 * (m ^ 2 - n ^ 2) ^ 2
        * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q)) :
    p ^ 2 ≤ 8 * (m + n) ^ 2 * N ∧ q ^ 2 ≤ 8 * (m + n) ^ 2 * N := by
  have hm : 0 < m := lt_trans hn hnm
  have hd : 0 < m - n := sub_pos.mpr hnm
  have hs : 0 < m ^ 2 - n ^ 2 := by nlinarith
  have hK : 0 < 16 * (m ^ 2 - n ^ 2) ^ 4 := mul_pos (by norm_num) (pow_pos hs 4)
  have hgle : g ≤ 16 * (m ^ 2 - n ^ 2) ^ 4 := Int.le_of_dvd hK hgd
  have hsq : 0 ≤ (m ^ 2 - n ^ 2) ^ 2 := sq_nonneg _
  have hpq2 : 0 ≤ 2 * m * n * (p - q) ^ 2 :=
    mul_nonneg (mul_nonneg (mul_nonneg (by norm_num) hm.le) hn.le) (sq_nonneg _)
  have hR : 0 ≤ 2 * (m ^ 2 - n ^ 2) ^ 2
      * ((m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q) := by
    have e : (m ^ 2 + n ^ 2) * (p ^ 2 + q ^ 2) - 4 * m * n * p * q
        = (m - n) ^ 2 * (p ^ 2 + q ^ 2) + 2 * m * n * (p - q) ^ 2 := by ring
    rw [e]
    exact mul_nonneg (mul_nonneg (by norm_num) hsq)
      (add_nonneg (mul_nonneg (sq_nonneg _) (add_nonneg (sq_nonneg _) (sq_nonneg _))) hpq2)
  have hN0 : 0 ≤ N := by
    rcases le_or_gt 0 N with h | h
    · exact h
    · exfalso
      have := mul_pos hg (neg_pos.mpr h)
      linarith
  have h16 : g * N ≤ 16 * (m ^ 2 - n ^ 2) ^ 4 * N := mul_le_mul_of_nonneg_right hgle hN0
  have hc : 0 < 2 * (m ^ 2 - n ^ 2) ^ 2 * (m - n) ^ 2 :=
    mul_pos (mul_pos (by norm_num) (pow_pos hs 2)) (pow_pos hd 2)
  have hA : 0 ≤ 2 * (m ^ 2 - n ^ 2) ^ 2 * (2 * m * n * (p - q) ^ 2) :=
    mul_nonneg (mul_nonneg (by norm_num) hsq) hpq2
  constructor
  · have h0 : 0 ≤ 2 * (m ^ 2 - n ^ 2) ^ 2 * ((m - n) ^ 2 * q ^ 2) :=
      mul_nonneg (mul_nonneg (by norm_num) hsq) (mul_nonneg (sq_nonneg _) (sq_nonneg _))
    have h : 2 * (m ^ 2 - n ^ 2) ^ 2 * (m - n) ^ 2 * p ^ 2
        ≤ 2 * (m ^ 2 - n ^ 2) ^ 2 * (m - n) ^ 2 * (8 * (m + n) ^ 2 * N) := by linarith
    exact le_of_mul_le_mul_left h hc
  · have h0 : 0 ≤ 2 * (m ^ 2 - n ^ 2) ^ 2 * ((m - n) ^ 2 * p ^ 2) :=
      mul_nonneg (mul_nonneg (by norm_num) hsq) (mul_nonneg (sq_nonneg _) (sq_nonneg _))
    have h : 2 * (m ^ 2 - n ^ 2) ^ 2 * (m - n) ^ 2 * q ^ 2
        ≤ 2 * (m ^ 2 - n ^ 2) ^ 2 * (m - n) ^ 2 * (8 * (m + n) ^ 2 * N) := by linarith
    exact le_of_mul_le_mul_left h hc

end G1
