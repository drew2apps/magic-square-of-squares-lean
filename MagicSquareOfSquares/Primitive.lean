import MagicSquareOfSquares.Basic

/-!
# Reduction to primitive squares

Dividing the nine roots of a magic square of distinct squares by their greatest common
divisor gives a primitive magic square of distinct squares. So if the open problem has a
solution, it has a primitive one, and the necessary conditions for primitive squares apply.

## Main results

* `IsMSS.of_mul`: a common factor can be removed from the nine roots.
* `IsMSS.exists_primitive`: a magic square of distinct squares yields a primitive one.
-/

namespace MagicSquareOfSquares

variable {a b c d e f g h i T : ℤ}

/-- A common nonzero factor `k` can be removed from the nine roots. -/
theorem IsMSS.of_mul {k : ℤ} (hk : k ≠ 0)
    (H : IsMSS (k * a) (k * b) (k * c) (k * d) (k * e) (k * f) (k * g) (k * h) (k * i)
      (k ^ 2 * T)) :
    IsMSS a b c d e f g h i T := by
  have HM : IsMagic ((k * a) ^ 2) ((k * b) ^ 2) ((k * c) ^ 2) ((k * d) ^ 2) ((k * e) ^ 2)
      ((k * f) ^ 2) ((k * g) ^ 2) ((k * h) ^ 2) ((k * i) ^ 2) (k ^ 2 * T) := H
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := HM
  have hk2 : k ^ 2 ≠ 0 := pow_ne_zero 2 hk
  exact IsMagic.mk (mul_left_cancel₀ hk2 (by linear_combination h1))
    (mul_left_cancel₀ hk2 (by linear_combination h2))
    (mul_left_cancel₀ hk2 (by linear_combination h3))
    (mul_left_cancel₀ hk2 (by linear_combination h4))
    (mul_left_cancel₀ hk2 (by linear_combination h5))
    (mul_left_cancel₀ hk2 (by linear_combination h6))
    (mul_left_cancel₀ hk2 (by linear_combination h7))
    (mul_left_cancel₀ hk2 (by linear_combination h8))

/-- From a magic square of distinct squares one gets a primitive magic square of distinct
squares, by dividing the nine roots by their greatest common divisor. -/
theorem IsMSS.exists_primitive (H : IsMSS a b c d e f g h i T)
    (hnd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    ∃ a' b' c' d' e' f' g' h' i' T' : ℤ, IsMSS a' b' c' d' e' f' g' h' i' T' ∧
      [a' ^ 2, b' ^ 2, c' ^ 2, d' ^ 2, e' ^ 2, f' ^ 2, g' ^ 2, h' ^ 2, i' ^ 2].Nodup ∧
      Primitive a' b' c' d' e' f' g' h' i' := by
  obtain ⟨σ, hσ, hgcd⟩ := Finset.extract_gcd (s := (Finset.univ : Finset (Fin 9)))
    ![a, b, c, d, e, f, g, h, i] Finset.univ_nonempty
  generalize (Finset.univ : Finset (Fin 9)).gcd ![a, b, c, d, e, f, g, h, i] = G at hσ
  have ha : a = G * σ 0 := hσ 0 (Finset.mem_univ _)
  have hb : b = G * σ 1 := hσ 1 (Finset.mem_univ _)
  have hc : c = G * σ 2 := hσ 2 (Finset.mem_univ _)
  have hd : d = G * σ 3 := hσ 3 (Finset.mem_univ _)
  have he : e = G * σ 4 := hσ 4 (Finset.mem_univ _)
  have hf : f = G * σ 5 := hσ 5 (Finset.mem_univ _)
  have hg : g = G * σ 6 := hσ 6 (Finset.mem_univ _)
  have hh : h = G * σ 7 := hσ 7 (Finset.mem_univ _)
  have hi : i = G * σ 8 := hσ 8 (Finset.mem_univ _)
  have hab : a ^ 2 ≠ b ^ 2 := by
    intro hEq
    have hnot := (List.nodup_cons.mp hnd).1
    rw [hEq] at hnot
    exact hnot (by simp)
  have hG : G ≠ 0 := by
    rintro rfl
    apply hab
    rw [ha, hb]
    simp
  have hT : G ^ 2 * (3 * σ 4 ^ 2) = T := by
    rw [H.sum_eq, he]
    ring
  have H' : IsMSS (G * σ 0) (G * σ 1) (G * σ 2) (G * σ 3) (G * σ 4) (G * σ 5) (G * σ 6)
      (G * σ 7) (G * σ 8) (G ^ 2 * (3 * σ 4 ^ 2)) := by
    rw [← ha, ← hb, ← hc, ← hd, ← he, ← hf, ← hg, ← hh, ← hi, hT]
    exact H
  refine ⟨σ 0, σ 1, σ 2, σ 3, σ 4, σ 5, σ 6, σ 7, σ 8, 3 * σ 4 ^ 2, IsMSS.of_mul hG H', ?_, ?_⟩
  · have hl : [σ 0 ^ 2, σ 1 ^ 2, σ 2 ^ 2, σ 3 ^ 2, σ 4 ^ 2, σ 5 ^ 2, σ 6 ^ 2, σ 7 ^ 2,
        σ 8 ^ 2].map (fun x => G ^ 2 * x) =
        [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2] := by
      simp only [List.map_cons, List.map_nil]
      rw [ha, hb, hc, hd, he, hf, hg, hh, hi]
      simp only [mul_pow]
    exact List.Nodup.of_map _ (hl ▸ hnd)
  · intro p hp ⟨h0, h1, h2, h3, h4, h5, h6, h7, h8⟩
    have hall : ∀ k ∈ (Finset.univ : Finset (Fin 9)), (p : ℤ) ∣ σ k := by
      intro k _
      fin_cases k <;> assumption
    have h1' : (p : ℤ) ∣ 1 := hgcd ▸ Finset.dvd_gcd hall
    have hp1 : p ∣ 1 := Int.natCast_dvd_natCast.mp (by exact_mod_cast h1')
    exact hp.one_lt.ne' (Nat.dvd_one.mp hp1)

end MagicSquareOfSquares
