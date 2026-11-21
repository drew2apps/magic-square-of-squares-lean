import MagicSquareOfSquares.Defs
import MagicSquareOfSquares.G1FullCertificates

/-!
# From a magic square of squares to the G1 system

Take a 3×3 magic square of squares with roots

```
a b c
d e f
g h i
```

The four edge cells satisfy `b² + h² = d² + f²`. If `(b + d) : (h + f) = m : n`, then
`b = mA + nB`, `d = mA − nB`, `h = nA − mB`, `f = nA + mB` for rationals `A`, `B`, and the
corners and the centre are the five forms of the G1 family at `t = A/B = (h + f)/(b − d)`.

## Main results

* `IsMSS.fullSquares`: such a square gives `G1.FullSquares m n ((h + f)/(b − d))`.
* `IsMSS.ratio_ne_of_le_200`: **in a magic square of squares, `(b + d) : (h + f)` is not
  `m : n` for any coprime `1 ≤ n < m ≤ 200` outside the 60 exceptions.** The signs of the
  roots are arbitrary, so this covers every choice of sign.
-/

namespace MagicSquareOfSquares

open G1

private theorem isSquare_of_mul_sq {x D R : ℚ} (hD : D ≠ 0) (h : x * D ^ 2 = R ^ 2) :
    IsSquare x :=
  ⟨R / D, by
    rw [div_mul_div_comm, eq_div_iff (mul_ne_zero hD hD)]
    linear_combination h⟩

private theorem fullSquares_of {m n : ℕ} {τ A B a c e g i : ℚ} (hB : B ≠ 0) (hτ : τ * B = A)
    (ha : a ^ 2 = (n : ℚ) ^ 2 * A ^ 2 + (m : ℚ) ^ 2 * B ^ 2)
    (he : e ^ 2 = ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (A ^ 2 + B ^ 2) / 2)
    (hi : i ^ 2 = (m : ℚ) ^ 2 * A ^ 2 + (n : ℚ) ^ 2 * B ^ 2)
    (hc : c ^ 2 = ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (A ^ 2 + B ^ 2) / 2 - 2 * m * n * A * B)
    (hg : g ^ 2 = ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (A ^ 2 + B ^ 2) / 2 + 2 * m * n * A * B) :
    FullSquares m n τ := by
  refine ⟨⟨?_, ?_, ?_⟩, ?_, ?_⟩
  · exact isSquare_of_mul_sq hB (R := a)
      (by linear_combination ((n : ℚ) ^ 2 * (τ * B + A)) * hτ - ha)
  · exact isSquare_of_mul_sq hB (R := e)
      (by linear_combination (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) / 2 * (τ * B + A)) * hτ - he)
  · exact isSquare_of_mul_sq hB (R := i)
      (by linear_combination ((m : ℚ) ^ 2 * (τ * B + A)) * hτ - hi)
  · exact isSquare_of_mul_sq hB (R := c)
      (by
        linear_combination (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) / 2 * (τ * B + A) - 2 * m * n * B) * hτ -
          hc)
  · exact isSquare_of_mul_sq hB (R := g)
      (by
        linear_combination (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) / 2 * (τ * B + A) + 2 * m * n * B) * hτ -
          hg)

/-- A magic square of squares whose edge roots satisfy `(b + d)·n = (h + f)·m` gives a solution
of the five-form system of the G1 family with parameters `m`, `n`. -/
theorem IsMSS.fullSquares {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hn : 0 < n) (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) (hbd : b ≠ d) :
    FullSquares m n (((h : ℚ) + f) / (b - d)) := by
  obtain ⟨r1, r2, r3, c1, c2, c3, dg, ad⟩ := hs
  have e1 : b ^ 2 + h ^ 2 = d ^ 2 + f ^ 2 := by linarith
  have ea : 2 * a ^ 2 = f ^ 2 + h ^ 2 := by linarith
  have ec : 2 * c ^ 2 = d ^ 2 + h ^ 2 := by linarith
  have eg : 2 * g ^ 2 = b ^ 2 + f ^ 2 := by linarith
  have ei : 2 * i ^ 2 = b ^ 2 + d ^ 2 := by linarith
  have ee : 2 * e ^ 2 = b ^ 2 + h ^ 2 := by linarith
  have e1q : (b : ℚ) ^ 2 + (h : ℚ) ^ 2 = (d : ℚ) ^ 2 + (f : ℚ) ^ 2 := by exact_mod_cast e1
  have eaq : 2 * (a : ℚ) ^ 2 = (f : ℚ) ^ 2 + (h : ℚ) ^ 2 := by exact_mod_cast ea
  have ecq : 2 * (c : ℚ) ^ 2 = (d : ℚ) ^ 2 + (h : ℚ) ^ 2 := by exact_mod_cast ec
  have egq : 2 * (g : ℚ) ^ 2 = (b : ℚ) ^ 2 + (f : ℚ) ^ 2 := by exact_mod_cast eg
  have eiq : 2 * (i : ℚ) ^ 2 = (b : ℚ) ^ 2 + (d : ℚ) ^ 2 := by exact_mod_cast ei
  have eeq : 2 * (e : ℚ) ^ 2 = (b : ℚ) ^ 2 + (h : ℚ) ^ 2 := by exact_mod_cast ee
  have hmnq : ((b : ℚ) + d) * n = ((h : ℚ) + f) * m := by exact_mod_cast hmn
  have hneq : (h : ℚ) + f ≠ 0 := by exact_mod_cast hne
  have hbdq : (b : ℚ) - d ≠ 0 := sub_ne_zero.mpr (by exact_mod_cast hbd)
  have hnq : (n : ℚ) ≠ 0 := by exact_mod_cast hn.ne'
  have h2n : (2 : ℚ) * n ≠ 0 := mul_ne_zero two_ne_zero hnq
  have key : (m : ℚ) * ((b : ℚ) - d) = n * ((f : ℚ) - h) := by
    apply mul_left_cancel₀ hneq
    linear_combination (-((b : ℚ) - d)) * hmnq + (n : ℚ) * e1q
  obtain ⟨A, hA⟩ : ∃ A : ℚ, 2 * n * A = (h : ℚ) + f := ⟨((h : ℚ) + f) / (2 * n), by field_simp⟩
  obtain ⟨B, hB⟩ : ∃ B : ℚ, 2 * n * B = (b : ℚ) - d := ⟨((b : ℚ) - d) / (2 * n), by field_simp⟩
  have hB0 : B ≠ 0 := by
    intro h0
    rw [h0, mul_zero] at hB
    exact hbdq hB.symm
  have hb : (b : ℚ) = m * A + n * B := by
    apply mul_left_cancel₀ h2n
    linear_combination hmnq - (m : ℚ) * hA - (n : ℚ) * hB
  have hd : (d : ℚ) = m * A - n * B := by
    apply mul_left_cancel₀ h2n
    linear_combination hmnq - (m : ℚ) * hA + (n : ℚ) * hB
  have hh : (h : ℚ) = n * A - m * B := by
    apply mul_left_cancel₀ h2n
    linear_combination key - (n : ℚ) * hA + (m : ℚ) * hB
  have hf : (f : ℚ) = n * A + m * B := by
    apply mul_left_cancel₀ h2n
    linear_combination -key - (n : ℚ) * hA - (m : ℚ) * hB
  have hτ : ((h : ℚ) + f) / (b - d) * B = A := by
    rw [div_mul_eq_mul_div, div_eq_iff hbdq]
    linear_combination A * hB - B * hA
  refine fullSquares_of (a := a) (c := c) (e := e) (g := g) (i := i) hB0 hτ ?_ ?_ ?_ ?_ ?_
  · rw [hf, hh] at eaq
    linear_combination (1 / 2 : ℚ) * eaq
  · rw [hb, hh] at eeq
    linear_combination (1 / 2 : ℚ) * eeq
  · rw [hb, hd] at eiq
    linear_combination (1 / 2 : ℚ) * eiq
  · rw [hd, hh] at ecq
    linear_combination (1 / 2 : ℚ) * ecq
  · rw [hb, hf] at egq
    linear_combination (1 / 2 : ℚ) * egq

/-- **In a 3×3 magic square of squares with roots `a b c / d e f / g h i`, the ratio
`(b + d) : (h + f)` of the edge roots is not `m : n` for any coprime `1 ≤ n < m ≤ 200` outside
the 60 pairs of `G1.exceptions200`.** No distinctness is assumed beyond `b ≠ d`, and the signs
of the roots are arbitrary. -/
theorem IsMSS.ratio_ne_of_le_200 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hn0 : 0 < n) (hnm : n < m) (hm : m ≤ 200) (hcop : Nat.Coprime m n)
    (hex : (m, n) ∉ exceptions200) (hne : h + f ≠ 0) (hbd : b ≠ d) :
    (b + d) * n ≠ (h + f) * m :=
  fun hmn => no_solution_le_200 hn0 hnm hm hcop hex _ (hs.fullSquares hn0 hmn hne hbd)

end MagicSquareOfSquares
