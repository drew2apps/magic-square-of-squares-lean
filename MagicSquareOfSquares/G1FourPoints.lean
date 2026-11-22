import Mathlib.Tactic
import Mathlib.RingTheory.Coprime.Basic

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

end G1
