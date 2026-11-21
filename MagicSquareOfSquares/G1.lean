import Mathlib

/-!
# Local obstructions for the diagonal of the G1 family

In the two-parameter family of 3×3 magic squares called G1 in
[Kepler71/magic-square-of-squares](https://github.com/Kepler71/magic-square-of-squares), six
cells are squares or products of the parameters, and the three diagonal cells are

`F0 = m² + n²t²`, `F4 = (m² + n²)(t² + 1)/2`, `F8 = n² + m²t²`.

A magic square of nine squares in the family with parameters `m, n` needs a rational `t` for
which `F0`, `F4`, `F8` are all squares of rationals.

Writing `t = U/V` in lowest terms, this needs the integers
`m²V² + n²U²`, `2(m² + n²)(U² + V²)`, `n²V² + m²U²` to be squares. This file proves that a
finite check modulo `M` rules that out.

## Main results

* `sqMask M`: the squares modulo `M`, as a bit mask.
* `obstructs m n M S`: a computable check. With `S = sqMask M`, it says that no pair `(u, v)`
  modulo `M` with `gcd(u, v, M) = 1` makes the three integers squares modulo `M`.
* `no_point_of_obstructs`: if the check succeeds then no rational `t` makes `F0`, `F4`, `F8`
  all squares.
-/

namespace MagicSquareOfSquares.G1

/-- `sqMask M`: the natural number whose bit `i` is set when `i` is a square modulo `M`. -/
def sqMask (M : ℕ) : ℕ :=
  (List.range M).foldl (fun acc r => acc ||| 1 <<< (r * r % M)) 0

/-- `obstructs m n M S`: with `S` the mask of squares modulo `M`, no pair `(u, v)` modulo `M`
with `gcd(u, v, M) = 1` makes `m²v² + n²u²`, `2(m² + n²)(u² + v²)` and `n²v² + m²u²` all
squares modulo `M`. -/
def obstructs (m n M S : ℕ) : Bool :=
  (List.range M).all fun u => (List.range M).all fun v =>
    !(Nat.gcd (Nat.gcd u v) M == 1) ||
      !(S.testBit ((m * m * v * v + n * n * u * u) % M) &&
        S.testBit ((2 * (m * m + n * n) * (u * u + v * v)) % M) &&
        S.testBit ((n * n * v * v + m * m * u * u) % M))

private theorem testBit_foldl {M : ℕ} (l : List ℕ) (acc i : ℕ)
    (h : acc.testBit i = true ∨ ∃ r ∈ l, r * r % M = i) :
    (l.foldl (fun acc r => acc ||| 1 <<< (r * r % M)) acc).testBit i = true := by
  induction l generalizing acc with
  | nil => simpa using h
  | cons r l ih =>
    rw [List.foldl_cons]
    apply ih
    rcases h with h | ⟨ρ, hρ, hi⟩
    · left
      rw [Nat.testBit_or, h, Bool.true_or]
    · rcases List.mem_cons.mp hρ with hr | hρ'
      · left
        subst hi
        subst hr
        rw [Nat.testBit_or, Nat.one_shiftLeft, Nat.testBit_two_pow_self, Bool.or_true]
      · right
        exact ⟨ρ, hρ', hi⟩

private theorem exists_nat_rep {M : ℕ} (hM : 0 < M) (x : ℤ) :
    ∃ k : ℕ, k < M ∧ (k : ℤ) ≡ x [ZMOD M] := by
  have hM' : (0 : ℤ) < M := by exact_mod_cast hM
  refine ⟨(x % M).toNat, ?_, ?_⟩
  · have h1 := Int.emod_lt_of_pos x hM'
    have h2 := Int.emod_nonneg x hM'.ne'
    omega
  · rw [Int.toNat_of_nonneg (Int.emod_nonneg x hM'.ne')]
    exact Int.mod_modEq x M

private theorem testBit_sqMask_of_int {M : ℕ} (hM : 0 < M) {x : ℕ} {r : ℤ}
    (h : r * r ≡ (x : ℤ) [ZMOD M]) : (sqMask M).testBit (x % M) = true := by
  obtain ⟨ρ, hρlt, hρ⟩ := exists_nat_rep hM r
  have h1 : ((ρ * ρ : ℕ) : ℤ) ≡ (x : ℤ) [ZMOD M] := by
    push_cast
    exact (hρ.mul hρ).trans h
  have h2 : ρ * ρ ≡ x [MOD M] := Int.natCast_modEq_iff.mp h1
  exact testBit_foldl _ _ _ (Or.inr ⟨ρ, List.mem_range.mpr hρlt, h2⟩)

/-- If the finite check `obstructs m n M (sqMask M)` succeeds, then no rational `t` makes the three
diagonal cells `m² + n²t²`, `(m² + n²)(t² + 1)/2`, `n² + m²t²` all squares. -/
theorem no_point_of_obstructs {m n M S : ℕ} (hM : 0 < M) (hS : S = sqMask M)
    (h : obstructs m n M S = true) (t : ℚ) :
    ¬(IsSquare ((m : ℚ) ^ 2 + (n : ℚ) ^ 2 * t ^ 2) ∧
      IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2) ∧
      IsSquare ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t ^ 2)) := by
  subst hS
  rintro ⟨h0, h4, h8⟩
  obtain ⟨U, hUdef⟩ : ∃ U : ℤ, U = t.num := ⟨_, rfl⟩
  obtain ⟨V, hVdef⟩ : ∃ V : ℕ, V = t.den := ⟨_, rfl⟩
  have hVpos : 0 < V := by rw [hVdef]; exact t.den_pos
  have hcop : Nat.Coprime U.natAbs V := by rw [hUdef, hVdef]; exact t.reduced
  have hV0 : (V : ℚ) ≠ 0 := by exact_mod_cast hVpos.ne'
  have hUt : (U : ℚ) = t * V := by
    rw [hUdef, hVdef]
    exact (Rat.mul_den_eq_num t).symm
  -- the three integers are squares
  have hA : IsSquare ((m : ℤ) * m * V * V + n * n * U * U) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : (((m : ℤ) * m * V * V + n * n * U * U : ℤ) : ℚ) =
        ((m : ℚ) ^ 2 + (n : ℚ) ^ 2 * t ^ 2) * ((V : ℚ) * V) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h0.mul ⟨(V : ℚ), rfl⟩
  have hB : IsSquare (2 * ((m : ℤ) * m + n * n) * (U * U + V * V)) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : ((2 * ((m : ℤ) * m + n * n) * (U * U + V * V) : ℤ) : ℚ) =
        (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2) * ((2 * (V : ℚ)) * (2 * V)) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h4.mul ⟨2 * (V : ℚ), rfl⟩
  have hC : IsSquare ((n : ℤ) * n * V * V + m * m * U * U) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : (((n : ℤ) * n * V * V + m * m * U * U : ℤ) : ℚ) =
        ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t ^ 2) * ((V : ℚ) * V) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h8.mul ⟨(V : ℚ), rfl⟩
  obtain ⟨a, ha⟩ := hA
  obtain ⟨b, hb⟩ := hB
  obtain ⟨c, hc⟩ := hC
  -- residues of U and V
  obtain ⟨u, hult, hu⟩ := exists_nat_rep hM U
  obtain ⟨v, hvlt, hv⟩ := exists_nat_rep hM (V : ℤ)
  have hbody := List.all_eq_true.mp (List.all_eq_true.mp h u (List.mem_range.mpr hult)) v
    (List.mem_range.mpr hvlt)
  have hg : Nat.gcd (Nat.gcd u v) M = 1 := by
    have hdM : Nat.gcd (Nat.gcd u v) M ∣ M := Nat.gcd_dvd_right _ _
    have hdu : Nat.gcd (Nat.gcd u v) M ∣ u := (Nat.gcd_dvd_left _ _).trans (Nat.gcd_dvd_left _ _)
    have hdv : Nat.gcd (Nat.gcd u v) M ∣ v := (Nat.gcd_dvd_left _ _).trans (Nat.gcd_dvd_right _ _)
    have key : ∀ {x : ℤ} {k : ℕ}, (k : ℤ) ≡ x [ZMOD M] → Nat.gcd (Nat.gcd u v) M ∣ k →
        ((Nat.gcd (Nat.gcd u v) M : ℕ) : ℤ) ∣ x := by
      intro x k hk hdk
      have h1 : (M : ℤ) ∣ x - k := Int.modEq_iff_dvd.mp hk
      have h2 := (Int.natCast_dvd_natCast.mpr hdM).trans h1
      have h3 := dvd_add h2 (Int.natCast_dvd_natCast.mpr hdk)
      simpa using h3
    have hdU : Nat.gcd (Nat.gcd u v) M ∣ U.natAbs := Int.natCast_dvd.mp (key hu hdu)
    have hdV : Nat.gcd (Nat.gcd u v) M ∣ V := Int.natCast_dvd_natCast.mp (key hv hdv)
    have hd1 := Nat.dvd_gcd hdU hdV
    rw [hcop.gcd_eq_one] at hd1
    exact Nat.dvd_one.mp hd1
  have s0 : (sqMask M).testBit ((m * m * v * v + n * n * u * u) % M) = true := by
    apply testBit_sqMask_of_int hM (r := a)
    rw [← ha]
    push_cast
    exact (((Int.ModEq.refl _).mul hv.symm).mul hv.symm).add
      (((Int.ModEq.refl _).mul hu.symm).mul hu.symm)
  have s4 : (sqMask M).testBit ((2 * (m * m + n * n) * (u * u + v * v)) % M) = true := by
    apply testBit_sqMask_of_int hM (r := b)
    rw [← hb]
    push_cast
    exact (Int.ModEq.refl _).mul ((hu.symm.mul hu.symm).add (hv.symm.mul hv.symm))
  have s8 : (sqMask M).testBit ((n * n * v * v + m * m * u * u) % M) = true := by
    apply testBit_sqMask_of_int hM (r := c)
    rw [← hc]
    push_cast
    exact (((Int.ModEq.refl _).mul hv.symm).mul hv.symm).add
      (((Int.ModEq.refl _).mul hu.symm).mul hu.symm)
  simp [hg, s0, s4, s8] at hbody

end MagicSquareOfSquares.G1
