import MagicSquareOfSquares.Basic

/-!
# Magic squares of squares: prime divisors of the entries

The layout is `a b c`, `d e f`, `g h i`: the centre is `e`, the corners are `a c g i` and the
(non-corner) edges are `b d f h`.

## Main results

* `IsMSS.not_dvd_center`: in a primitive square no prime `p ≡ 3 (mod 4)` divides the centre.
* `IsMSS.center_prime_factors`: in a primitive square every prime factor of the centre is
  `≡ 1 (mod 4)`.
* `IsMSS.dvd_edges_of_dvd_corner`: a prime `p ≡ 3 (mod 4)` dividing a corner divides the two
  edges that do not touch it.
* `IsMSS.dvd_center_of_dvd_corner`: a prime `p ≡ 3, 5 (mod 8)` dividing a corner divides the
  centre and the opposite corner.
* `IsMSS.not_dvd_edge`: in a primitive square no prime `p ≡ 3, 5 (mod 8)` divides an edge.
* `IsMSS.not_dvd_of_three_mod_eight`: in a primitive square no prime `p ≡ 3 (mod 8)` divides
  any entry.
-/

namespace MagicSquareOfSquares

variable {a b c d e f g h i T : ℤ}

/-- A prime `p ≡ 3 (mod 4)` dividing a sum of two squares divides both roots. -/
theorem dvd_of_dvd_sq_add_sq {p : ℕ} [Fact p.Prime] (hp : p % 4 = 3) {y z : ℤ}
    (hd : (p : ℤ) ∣ y ^ 2 + z ^ 2) : (p : ℤ) ∣ y ∧ (p : ℤ) ∣ z := by
  have h0 : (y : ZMod p) ^ 2 + (z : ZMod p) ^ 2 = 0 := by
    have := (ZMod.intCast_zmod_eq_zero_iff_dvd (y ^ 2 + z ^ 2) p).mpr hd
    push_cast at this
    exact this
  have hz : (z : ZMod p) = 0 := by
    by_contra hz
    exact ZMod.mod_four_ne_three_of_sq_eq_neg_sq' hz (eq_neg_of_add_eq_zero_left h0) hp
  have hy : (y : ZMod p) = 0 := by
    rw [hz] at h0
    simpa using h0
  exact ⟨(ZMod.intCast_zmod_eq_zero_iff_dvd y p).mp hy,
    (ZMod.intCast_zmod_eq_zero_iff_dvd z p).mp hz⟩

/-- A prime `p ≡ 3, 5 (mod 8)` dividing `x² - 2y²` divides both `x` and `y`,
because `2` is not a square modulo `p`. -/
theorem dvd_of_dvd_sq_sub_two_mul_sq {p : ℕ} [Fact p.Prime] (hp : p % 8 = 3 ∨ p % 8 = 5)
    {x y : ℤ} (hd : (p : ℤ) ∣ x ^ 2 - 2 * y ^ 2) : (p : ℤ) ∣ x ∧ (p : ℤ) ∣ y := by
  have hp2 : p ≠ 2 := by omega
  have h0 : (x : ZMod p) ^ 2 = 2 * (y : ZMod p) ^ 2 := by
    have := (ZMod.intCast_zmod_eq_zero_iff_dvd (x ^ 2 - 2 * y ^ 2) p).mpr hd
    push_cast at this
    exact sub_eq_zero.mp this
  have hy : (y : ZMod p) = 0 := by
    by_contra hy
    have hsq : IsSquare (2 : ZMod p) := by
      refine ⟨(x : ZMod p) * (y : ZMod p)⁻¹, ?_⟩
      have hmul : (x : ZMod p) * (y : ZMod p)⁻¹ * ((x : ZMod p) * (y : ZMod p)⁻¹) =
          (x : ZMod p) ^ 2 * ((y : ZMod p) ^ 2)⁻¹ := by ring
      rw [hmul, h0, mul_assoc, mul_inv_cancel₀ (pow_ne_zero 2 hy), mul_one]
    rcases (ZMod.exists_sq_eq_two_iff hp2).mp hsq with h8 | h8 <;> omega
  have hx : (x : ZMod p) = 0 := by
    rw [hy] at h0
    simpa using h0
  exact ⟨(ZMod.intCast_zmod_eq_zero_iff_dvd x p).mp hx,
    (ZMod.intCast_zmod_eq_zero_iff_dvd y p).mp hy⟩

private theorem dvd_sq {p : ℕ} {x : ℤ} (hx : (p : ℤ) ∣ x) : (p : ℤ) ∣ x ^ 2 :=
  Dvd.dvd.pow hx two_ne_zero

namespace IsMSS

/-- A prime `p ≡ 3 (mod 4)` dividing the centre divides all nine roots. -/
theorem dvd_all_of_dvd_center (H : IsMSS a b c d e f g h i T) {p : ℕ} [Fact p.Prime]
    (hp : p % 4 = 3) (he : (p : ℤ) ∣ e) :
    (p : ℤ) ∣ a ∧ (p : ℤ) ∣ b ∧ (p : ℤ) ∣ c ∧ (p : ℤ) ∣ d ∧ (p : ℤ) ∣ e ∧
    (p : ℤ) ∣ f ∧ (p : ℤ) ∣ g ∧ (p : ℤ) ∣ h ∧ (p : ℤ) ∣ i := by
  obtain ⟨h1, h2, h3, h4⟩ := H.opposite
  have key : ∀ {y z : ℤ}, y ^ 2 + z ^ 2 = 2 * e ^ 2 → (p : ℤ) ∣ y ∧ (p : ℤ) ∣ z := by
    intro y z hyz
    apply dvd_of_dvd_sq_add_sq hp
    rw [hyz]
    exact Dvd.dvd.mul_left (dvd_sq he) 2
  exact ⟨(key h1).1, (key h3).1, (key h2).1, (key h4).1, he, (key h4).2, (key h2).2,
    (key h3).2, (key h1).2⟩

/-- In a primitive magic square of squares no prime `p ≡ 3 (mod 4)` divides the centre. -/
theorem not_dvd_center (H : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i)
    {p : ℕ} [Fact p.Prime] (hp : p % 4 = 3) : ¬(p : ℤ) ∣ e :=
  fun he => hprim p Fact.out (H.dvd_all_of_dvd_center hp he)

/-- In a primitive magic square of squares every prime factor of the centre is
`≡ 1 (mod 4)`. -/
theorem center_prime_factors (H : IsMSS a b c d e f g h i T)
    (hprim : Primitive a b c d e f g h i) {p : ℕ} (hpp : p.Prime) (he : (p : ℤ) ∣ e) :
    p % 4 = 1 := by
  have : Fact p.Prime := ⟨hpp⟩
  have h3 : p % 4 ≠ 3 := fun h3 => H.not_dvd_center hprim h3 he
  have h2 : p ≠ 2 := by
    rintro rfl
    have hodd : Odd e := (H.odd hprim).2.2.2.2.1
    have heven : Even e := even_iff_two_dvd.mpr (by exact_mod_cast he)
    exact (Int.not_even_iff_odd.mpr hodd) heven
  have := hpp.eq_two_or_odd
  omega

/-- A prime `p ≡ 3 (mod 4)` dividing a corner divides the two edges that do not touch
that corner. -/
theorem dvd_edges_of_dvd_corner (H : IsMSS a b c d e f g h i T) {p : ℕ} [Fact p.Prime]
    (hp : p % 4 = 3) :
    ((p : ℤ) ∣ a → (p : ℤ) ∣ f ∧ (p : ℤ) ∣ h) ∧ ((p : ℤ) ∣ c → (p : ℤ) ∣ d ∧ (p : ℤ) ∣ h) ∧
    ((p : ℤ) ∣ g → (p : ℤ) ∣ b ∧ (p : ℤ) ∣ f) ∧ ((p : ℤ) ∣ i → (p : ℤ) ∣ b ∧ (p : ℤ) ∣ d) := by
  obtain ⟨c1, c2, c3, c4⟩ := H.corner
  have key : ∀ {x y z : ℤ}, 2 * x ^ 2 = y ^ 2 + z ^ 2 → (p : ℤ) ∣ x →
      (p : ℤ) ∣ y ∧ (p : ℤ) ∣ z := by
    intro x y z hxyz hx
    apply dvd_of_dvd_sq_add_sq hp
    rw [← hxyz]
    exact Dvd.dvd.mul_left (dvd_sq hx) 2
  exact ⟨key c1, key c2, key c3, key c4⟩

/-- A prime `p ≡ 3, 5 (mod 8)` dividing a corner divides the centre and the opposite
corner. -/
theorem dvd_center_of_dvd_corner (H : IsMSS a b c d e f g h i T) {p : ℕ} [Fact p.Prime]
    (hp : p % 8 = 3 ∨ p % 8 = 5) :
    ((p : ℤ) ∣ a → (p : ℤ) ∣ e ∧ (p : ℤ) ∣ i) ∧ ((p : ℤ) ∣ c → (p : ℤ) ∣ e ∧ (p : ℤ) ∣ g) ∧
    ((p : ℤ) ∣ g → (p : ℤ) ∣ e ∧ (p : ℤ) ∣ c) ∧ ((p : ℤ) ∣ i → (p : ℤ) ∣ e ∧ (p : ℤ) ∣ a) := by
  obtain ⟨h1, h2, h3, h4⟩ := H.opposite
  have key : ∀ {x y : ℤ}, x ^ 2 + y ^ 2 = 2 * e ^ 2 → (p : ℤ) ∣ x →
      (p : ℤ) ∣ e ∧ (p : ℤ) ∣ y := by
    intro x y hxy hx
    have hd : (p : ℤ) ∣ y ^ 2 - 2 * e ^ 2 := by
      have hrw : y ^ 2 - 2 * e ^ 2 = -(x ^ 2) := by linarith
      rw [hrw]
      exact (dvd_neg).mpr (dvd_sq hx)
    exact (dvd_of_dvd_sq_sub_two_mul_sq hp hd).symm
  exact ⟨key h1, key h2, key (by linarith), key (by linarith)⟩

/-- A prime `p ≡ 3, 5 (mod 8)` dividing the top edge `b` divides all nine roots. -/
theorem dvd_all_of_dvd_edge (H : IsMSS a b c d e f g h i T) {p : ℕ} [Fact p.Prime]
    (hp : p % 8 = 3 ∨ p % 8 = 5) (hb : (p : ℤ) ∣ b) :
    (p : ℤ) ∣ a ∧ (p : ℤ) ∣ b ∧ (p : ℤ) ∣ c ∧ (p : ℤ) ∣ d ∧ (p : ℤ) ∣ e ∧
    (p : ℤ) ∣ f ∧ (p : ℤ) ∣ g ∧ (p : ℤ) ∣ h ∧ (p : ℤ) ∣ i := by
  obtain ⟨o1, o2, o3, o4⟩ := H.opposite
  obtain ⟨c1, c2, c3, c4⟩ := H.corner
  have hb2 : (p : ℤ) ∣ -(b ^ 2) := (dvd_neg).mpr (dvd_sq hb)
  have hfg : (p : ℤ) ∣ f ∧ (p : ℤ) ∣ g := by
    apply dvd_of_dvd_sq_sub_two_mul_sq hp
    have hrw : f ^ 2 - 2 * g ^ 2 = -(b ^ 2) := by linarith
    rw [hrw]
    exact hb2
  have hdi : (p : ℤ) ∣ d ∧ (p : ℤ) ∣ i := by
    apply dvd_of_dvd_sq_sub_two_mul_sq hp
    have hrw : d ^ 2 - 2 * i ^ 2 = -(b ^ 2) := by linarith
    rw [hrw]
    exact hb2
  have hhe : (p : ℤ) ∣ h ∧ (p : ℤ) ∣ e := by
    apply dvd_of_dvd_sq_sub_two_mul_sq hp
    have hrw : h ^ 2 - 2 * e ^ 2 = -(b ^ 2) := by linarith
    rw [hrw]
    exact hb2
  have ha : (p : ℤ) ∣ a := by
    apply Int.Prime.dvd_pow' (Fact.out : p.Prime) (k := 2)
    have hrw : a ^ 2 = 2 * e ^ 2 - i ^ 2 := by linarith
    rw [hrw]
    exact dvd_sub (Dvd.dvd.mul_left (dvd_sq hhe.2) 2) (dvd_sq hdi.2)
  have hc : (p : ℤ) ∣ c := by
    apply Int.Prime.dvd_pow' (Fact.out : p.Prime) (k := 2)
    have hrw : c ^ 2 = 2 * e ^ 2 - g ^ 2 := by linarith
    rw [hrw]
    exact dvd_sub (Dvd.dvd.mul_left (dvd_sq hhe.2) 2) (dvd_sq hfg.2)
  exact ⟨ha, hb, hc, hdi.1, hhe.2, hfg.1, hfg.2, hhe.1, hdi.2⟩

/-- In a primitive magic square of squares no prime `p ≡ 3, 5 (mod 8)` divides an edge. -/
theorem not_dvd_edge (H : IsMSS a b c d e f g h i T) (hprim : Primitive a b c d e f g h i)
    {p : ℕ} [Fact p.Prime] (hp : p % 8 = 3 ∨ p % 8 = 5) :
    ¬(p : ℤ) ∣ b ∧ ¬(p : ℤ) ∣ d ∧ ¬(p : ℤ) ∣ f ∧ ¬(p : ℤ) ∣ h := by
  refine ⟨fun hb => ?_, fun hd => ?_, fun hf => ?_, fun hh => ?_⟩
  · exact hprim p Fact.out (H.dvd_all_of_dvd_edge hp hb)
  · obtain ⟨xg, xd, xa, xh, xe, xb, xi, xf, xc⟩ := H.rotate.dvd_all_of_dvd_edge hp hd
    exact hprim p Fact.out ⟨xa, xb, xc, xd, xe, xf, xg, xh, xi⟩
  · obtain ⟨xc, xf, xi, xb, xe, xh, xa, xd, xg⟩ :=
      H.rotate.rotate.rotate.dvd_all_of_dvd_edge hp hf
    exact hprim p Fact.out ⟨xa, xb, xc, xd, xe, xf, xg, xh, xi⟩
  · obtain ⟨xi, xh, xg, xf, xe, xd, xc, xb, xa⟩ := H.rotate.rotate.dvd_all_of_dvd_edge hp hh
    exact hprim p Fact.out ⟨xa, xb, xc, xd, xe, xf, xg, xh, xi⟩

/-- In a primitive magic square of squares no prime `p ≡ 3 (mod 8)` divides any entry. -/
theorem not_dvd_of_three_mod_eight (H : IsMSS a b c d e f g h i T)
    (hprim : Primitive a b c d e f g h i) {p : ℕ} [Fact p.Prime] (hp : p % 8 = 3) :
    ¬(p : ℤ) ∣ a ∧ ¬(p : ℤ) ∣ b ∧ ¬(p : ℤ) ∣ c ∧ ¬(p : ℤ) ∣ d ∧ ¬(p : ℤ) ∣ e ∧
    ¬(p : ℤ) ∣ f ∧ ¬(p : ℤ) ∣ g ∧ ¬(p : ℤ) ∣ h ∧ ¬(p : ℤ) ∣ i := by
  have hp4 : p % 4 = 3 := by omega
  have he : ¬(p : ℤ) ∣ e := H.not_dvd_center hprim hp4
  obtain ⟨hb, hd, hf, hh⟩ := H.not_dvd_edge hprim (Or.inl hp)
  obtain ⟨ka, kc, kg, ki⟩ := H.dvd_center_of_dvd_corner (p := p) (Or.inl hp)
  exact ⟨fun ha => he (ka ha).1, hb, fun hc => he (kc hc).1, hd, he, hf,
    fun hg => he (kg hg).1, hh, fun hi => he (ki hi).1⟩

end IsMSS

end MagicSquareOfSquares
