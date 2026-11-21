import MagicSquareOfSquares.G1

/-!
# The G1 family: primes that are 5 modulo 8

If a prime `p ≡ 5 (mod 8)` divides `m² + n²` to an odd power and does not divide `n`, then no
rational `t` makes the three diagonal cells of the G1 family squares. This is the general
form of the obstruction that the audit sees at the moduli `25`, `169`, `841`, ….

The argument: clear denominators, `t = U/V` in lowest terms. The middle cell gives
`2(m² + n²)(U² + V²) = b²`, so `p` divides `U² + V²` as well. Modulo `p` the first cell then
reads `a² = -2n²V²` with `n`, `V` invertible, so `-2` would be a square modulo `p`, which
fails for `p ≡ 5 (mod 8)`.

The statement is the one used in the Kepler71 notes; this file is an independent Lean proof.
-/

namespace MagicSquareOfSquares.G1

/-- A prime `p ≡ 5 (mod 8)` that divides `m² + n²` to an odd power and does not divide `n`
rules out every rational `t`. -/
theorem no_point_of_five_mod_eight {m n p : ℕ} (hp : p.Prime) (hp8 : p % 8 = 5)
    (hodd : Odd (padicValNat p (m * m + n * n))) (hn : ¬p ∣ n) (t : ℚ) :
    ¬DiagonalSquares m n t := by
  intro hd
  have : Fact p.Prime := ⟨hp⟩
  obtain ⟨U, V, hV, hcop, ⟨a, ha⟩, ⟨b, hb⟩, -⟩ := exists_coprime_squares hd
  obtain ⟨u₀, hu₀⟩ : ∃ u₀ : ℕ, u₀ = U.natAbs := ⟨_, rfl⟩
  have hUU : U * U = (u₀ : ℤ) * u₀ := by rw [hu₀, Int.natCast_natAbs, abs_mul_abs_self]
  have ha' : (m : ℤ) * m * V * V + n * n * ((u₀ : ℤ) * u₀) = a * a := by
    rw [← hUU, ← ha]
    ring
  have hb' : 2 * ((m : ℤ) * m + n * n) * ((u₀ : ℤ) * u₀ + V * V) = b * b := by
    rw [← hUU, ← hb]
  have hnat : 2 * (m * m + n * n) * (u₀ * u₀ + V * V) = b.natAbs * b.natAbs := by
    have h2 : ((2 * (m * m + n * n) * (u₀ * u₀ + V * V) : ℕ) : ℤ) =
        ((b.natAbs * b.natAbs : ℕ) : ℤ) := by
      push_cast
      rw [abs_mul_abs_self, ← hb']
    exact_mod_cast h2
  have hN : m * m + n * n ≠ 0 := by
    intro h0
    rw [h0] at hodd
    simp at hodd
  have hW : u₀ * u₀ + V * V ≠ 0 := by
    have := Nat.mul_pos hV hV
    omega
  have hL : 2 * (m * m + n * n) ≠ 0 := Nat.mul_ne_zero two_ne_zero hN
  have hβ : b.natAbs ≠ 0 := by
    intro h0
    exact Nat.mul_ne_zero hL hW (hnat.trans (by simp [h0]))
  have hp2 : ¬p ∣ 2 := by
    intro h
    have h2 : p = 2 := le_antisymm (Nat.le_of_dvd (by norm_num) h) hp.two_le
    omega
  have hval : padicValNat p (m * m + n * n) + padicValNat p (u₀ * u₀ + V * V) =
      2 * padicValNat p b.natAbs := by
    have h1 := congrArg (padicValNat p) hnat
    rw [padicValNat.mul hL hW, padicValNat.mul two_ne_zero hN, padicValNat.mul hβ hβ,
      padicValNat.eq_zero_of_not_dvd hp2] at h1
    omega
  obtain ⟨k, hk⟩ := hodd
  have hpW : p ∣ u₀ * u₀ + V * V := dvd_of_one_le_padicValNat (by omega)
  have hpN : p ∣ m * m + n * n := dvd_of_one_le_padicValNat (by omega)
  have hVp : ¬p ∣ V := by
    intro hpV
    have h1 : p ∣ u₀ * u₀ := (Nat.dvd_add_left (Dvd.dvd.mul_right hpV V)).mp hpW
    have h2 : p ∣ u₀ := (hp.dvd_mul.mp h1).elim id id
    have hg := Nat.dvd_gcd h2 hpV
    rw [hu₀, hcop.gcd_eq_one] at hg
    exact hp.one_lt.ne' (Nat.dvd_one.mp hg)
  have hm0 : (m : ZMod p) * m + n * n = 0 := by
    have h := (ZMod.natCast_eq_zero_iff (m * m + n * n) p).mpr hpN
    push_cast at h
    exact h
  have hU0 : (u₀ : ZMod p) * u₀ + V * V = 0 := by
    have h := (ZMod.natCast_eq_zero_iff (u₀ * u₀ + V * V) p).mpr hpW
    push_cast at h
    exact h
  have hAc : (m : ZMod p) * m * V * V + n * n * (u₀ * u₀) = a * a := by
    have h := congrArg (Int.cast : ℤ → ZMod p) ha'
    push_cast at h
    exact h
  have hsq : (a : ZMod p) * a = -2 * (n * n) * (V * V) := by
    linear_combination (-1 : ZMod p) * hAc + (V * V : ZMod p) * hm0 + (n * n : ZMod p) * hU0
  have hn0 : (n : ZMod p) ≠ 0 := fun h => hn ((ZMod.natCast_eq_zero_iff n p).mp h)
  have hV0 : (V : ZMod p) ≠ 0 := fun h => hVp ((ZMod.natCast_eq_zero_iff V p).mp h)
  have hnV : (n : ZMod p) * V ≠ 0 := mul_ne_zero hn0 hV0
  have hsquare : IsSquare (-2 : ZMod p) := by
    refine ⟨a * ((n : ZMod p) * V)⁻¹, ?_⟩
    have h : (a : ZMod p) * ((n : ZMod p) * V)⁻¹ * (a * ((n : ZMod p) * V)⁻¹) = -2 := by
      calc (a : ZMod p) * ((n : ZMod p) * V)⁻¹ * (a * ((n : ZMod p) * V)⁻¹)
          = a * a * (((n : ZMod p) * V)⁻¹ * ((n : ZMod p) * V)⁻¹) := by ring
        _ = -2 * ((n : ZMod p) * V * ((n : ZMod p) * V)⁻¹) *
              ((n : ZMod p) * V * ((n : ZMod p) * V)⁻¹) := by
            rw [hsq]
            ring
        _ = -2 := by
            rw [mul_inv_cancel₀ hnV]
            ring
    exact h.symm
  have h := (ZMod.exists_sq_eq_neg_two_iff (p := p) (by omega)).mp hsquare
  omega

/-- The same for coprime `m`, `n`, where `p ∤ n` is automatic. -/
theorem no_point_of_five_mod_eight_of_coprime {m n p : ℕ} (hp : p.Prime) (hp8 : p % 8 = 5)
    (hodd : Odd (padicValNat p (m * m + n * n))) (hmn : Nat.Coprime m n) (t : ℚ) :
    ¬DiagonalSquares m n t := by
  refine no_point_of_five_mod_eight hp hp8 hodd ?_ t
  intro hpn
  have hpN : p ∣ m * m + n * n :=
    dvd_of_one_le_padicValNat (by obtain ⟨k, hk⟩ := hodd; omega)
  have h1 : p ∣ m * m := (Nat.dvd_add_left (Dvd.dvd.mul_right hpn n)).mp hpN
  have h2 : p ∣ m := (hp.dvd_mul.mp h1).elim id id
  have hg := Nat.dvd_gcd h2 hpn
  rw [hmn.gcd_eq_one] at hg
  exact hp.one_lt.ne' (Nat.dvd_one.mp hg)

/-- The case seen most often: `p ≡ 5 (mod 8)` divides `m² + n²` exactly once. -/
theorem no_point_of_exactly_divides {m n p : ℕ} (hp : p.Prime) (hp8 : p % 8 = 5)
    (h1 : p ∣ m * m + n * n) (h2 : ¬p * p ∣ m * m + n * n) (hn : ¬p ∣ n) (t : ℚ) :
    ¬DiagonalSquares m n t := by
  have : Fact p.Prime := ⟨hp⟩
  obtain ⟨q, hq⟩ := h1
  have hq0 : ¬p ∣ q := fun ⟨r, hr⟩ => h2 ⟨r, by rw [hq, hr]; ring⟩
  have hqne : q ≠ 0 := fun h => hq0 (h ▸ dvd_zero p)
  refine no_point_of_five_mod_eight hp hp8 ?_ hn t
  rw [hq, padicValNat.mul hp.ne_zero hqne, padicValNat_self,
    padicValNat.eq_zero_of_not_dvd hq0]
  exact ⟨0, by norm_num⟩

end MagicSquareOfSquares.G1
