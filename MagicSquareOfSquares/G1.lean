import Mathlib

/-!
# Local obstructions for the diagonal of the G1 family

In the two-parameter family of 3×3 magic squares called G1 in
[Kepler71/magic-square-of-squares](https://github.com/Kepler71/magic-square-of-squares), six
cells are squares or products of the parameters, and the three diagonal cells are

`F0 = m² + n²t²`, `F4 = (m² + n²)(t² + 1)/2`, `F8 = n² + m²t²`.

A magic square of nine squares in the family with parameters `m, n` needs a rational `t` for
which `F0`, `F4`, `F8` are all squares of rationals (`DiagonalSquares m n t`).

Writing `t = U/V` in lowest terms, this needs the integers
`m²V² + n²U²`, `2(m² + n²)(U² + V²)`, `n²V² + m²U²` to be squares. This file proves that a
finite check modulo `M` rules that out.

## Main results

* `DiagonalSquares m n t`: the three diagonal cells are squares.
* `exists_coprime_squares`: the reduction to three integer squares.
* `sqMask M`: the squares modulo `M`, as a bit mask.
* `obstructs m n M S`: a computable check. With `S = sqMask M`, it says that no pair `(u, v)`
  modulo `M` with `gcd(u, v, M) = 1` makes the three integers squares modulo `M`.
* `no_point_of_obstructs`: if the check succeeds then `DiagonalSquares m n t` fails for every
  rational `t`.
* `no_point_of_obstructs_mod`: the same when the check succeeds for `m % M` and `n % M`.
* `obstructsProj`, `no_point_of_obstructsProj_mod`: a cheaper check for prime powers.
-/

namespace MagicSquareOfSquares.G1

/-- The three diagonal cells `m² + n²t²`, `(m² + n²)(t² + 1)/2`, `n² + m²t²` of the G1 family
are squares of rationals. -/
def DiagonalSquares (m n : ℕ) (t : ℚ) : Prop :=
  IsSquare ((m : ℚ) ^ 2 + (n : ℚ) ^ 2 * t ^ 2) ∧
    IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2) ∧
    IsSquare ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t ^ 2)

/-- If the three diagonal cells are squares for `t = U/V` in lowest terms, then the integers
`m²V² + n²U²`, `2(m² + n²)(U² + V²)` and `n²V² + m²U²` are squares. -/
theorem exists_coprime_squares {m n : ℕ} {t : ℚ} (h : DiagonalSquares m n t) :
    ∃ (U : ℤ) (V : ℕ), 0 < V ∧ Nat.Coprime U.natAbs V ∧
      IsSquare ((m : ℤ) * m * V * V + n * n * U * U) ∧
      IsSquare (2 * ((m : ℤ) * m + n * n) * (U * U + V * V)) ∧
      IsSquare ((n : ℤ) * n * V * V + m * m * U * U) := by
  obtain ⟨h0, h4, h8⟩ := h
  obtain ⟨U, hUdef⟩ : ∃ U : ℤ, U = t.num := ⟨_, rfl⟩
  obtain ⟨V, hVdef⟩ : ∃ V : ℕ, V = t.den := ⟨_, rfl⟩
  have hVpos : 0 < V := by rw [hVdef]; exact t.den_pos
  have hcop : Nat.Coprime U.natAbs V := by rw [hUdef, hVdef]; exact t.reduced
  have hUt : (U : ℚ) = t * V := by
    rw [hUdef, hVdef]
    exact (Rat.mul_den_eq_num t).symm
  refine ⟨U, V, hVpos, hcop, ?_, ?_, ?_⟩
  · rw [← Rat.isSquare_intCast_iff]
    have e : (((m : ℤ) * m * V * V + n * n * U * U : ℤ) : ℚ) =
        ((m : ℚ) ^ 2 + (n : ℚ) ^ 2 * t ^ 2) * ((V : ℚ) * V) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h0.mul ⟨(V : ℚ), rfl⟩
  · rw [← Rat.isSquare_intCast_iff]
    have e : ((2 * ((m : ℤ) * m + n * n) * (U * U + V * V) : ℤ) : ℚ) =
        (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2) * ((2 * (V : ℚ)) * (2 * V)) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h4.mul ⟨2 * (V : ℚ), rfl⟩
  · rw [← Rat.isSquare_intCast_iff]
    have e : (((n : ℤ) * n * V * V + m * m * U * U : ℤ) : ℚ) =
        ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t ^ 2) * ((V : ℚ) * V) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h8.mul ⟨(V : ℚ), rfl⟩

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

theorem testBit_sqMask_of_int {M : ℕ} (hM : 0 < M) {x : ℕ} {r : ℤ}
    (h : r * r ≡ (x : ℤ) [ZMOD M]) : (sqMask M).testBit (x % M) = true := by
  obtain ⟨ρ, hρlt, hρ⟩ := exists_nat_rep hM r
  have h1 : ((ρ * ρ : ℕ) : ℤ) ≡ (x : ℤ) [ZMOD M] := by
    push_cast
    exact (hρ.mul hρ).trans h
  have h2 : ρ * ρ ≡ x [MOD M] := Int.natCast_modEq_iff.mp h1
  exact testBit_foldl _ _ _ (Or.inr ⟨ρ, List.mem_range.mpr hρlt, h2⟩)

private theorem core {m n m' n' M S : ℕ} (hM : 0 < M) (hS : S = sqMask M)
    (hm : (m' : ℤ) ≡ (m : ℤ) [ZMOD M]) (hn : (n' : ℤ) ≡ (n : ℤ) [ZMOD M])
    (h : obstructs m' n' M S = true) (t : ℚ) : ¬DiagonalSquares m n t := by
  subst hS
  intro hd
  obtain ⟨U, V, -, hcop, ⟨a, ha⟩, ⟨b, hb⟩, ⟨c, hc⟩⟩ := exists_coprime_squares hd
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
  have s0 : (sqMask M).testBit ((m' * m' * v * v + n' * n' * u * u) % M) = true := by
    apply testBit_sqMask_of_int hM (r := a)
    rw [← ha]
    push_cast
    exact (((hm.symm.mul hm.symm).mul hv.symm).mul hv.symm).add
      (((hn.symm.mul hn.symm).mul hu.symm).mul hu.symm)
  have s4 : (sqMask M).testBit ((2 * (m' * m' + n' * n') * (u * u + v * v)) % M) = true := by
    apply testBit_sqMask_of_int hM (r := b)
    rw [← hb]
    push_cast
    exact ((Int.ModEq.refl 2).mul ((hm.symm.mul hm.symm).add (hn.symm.mul hn.symm))).mul
      ((hu.symm.mul hu.symm).add (hv.symm.mul hv.symm))
  have s8 : (sqMask M).testBit ((n' * n' * v * v + m' * m' * u * u) % M) = true := by
    apply testBit_sqMask_of_int hM (r := c)
    rw [← hc]
    push_cast
    exact (((hn.symm.mul hn.symm).mul hv.symm).mul hv.symm).add
      (((hm.symm.mul hm.symm).mul hu.symm).mul hu.symm)
  simp [hg, s0, s4, s8] at hbody

/-- If the finite check `obstructs m n M (sqMask M)` succeeds, then no rational `t` makes the
three diagonal cells of the G1 family squares. -/
theorem no_point_of_obstructs {m n M S : ℕ} (hM : 0 < M) (hS : S = sqMask M)
    (h : obstructs m n M S = true) (t : ℚ) : ¬DiagonalSquares m n t :=
  core hM hS (Int.ModEq.refl _) (Int.ModEq.refl _) h t

/-- The same conclusion when the check succeeds for the residues of `m` and `n` modulo `M`.
So the exclusion depends only on `m` and `n` modulo `M`, and holds for parameters of any
size. -/
theorem no_point_of_obstructs_mod {m n M S : ℕ} (hM : 0 < M) (hS : S = sqMask M)
    (h : obstructs (m % M) (n % M) M S = true) (t : ℚ) : ¬DiagonalSquares m n t :=
  core hM hS (by push_cast; exact Int.mod_modEq _ _) (by push_cast; exact Int.mod_modEq _ _) h t

/-! ## A cheaper check for prime powers

For `M` a power of a prime `p`, every pair `(u, v)` with `gcd(u, v, M) = 1` is a unit multiple
of `(1, v')` or of `(u', 1)` with `p ∣ u'`. Scaling by a unit multiplies the three integers by
a square, so it is enough to test those pairs. -/

/-- `sqTriple m n M S u v`: with `S` the mask of squares modulo `M`, the three integers
`m²v² + n²u²`, `2(m² + n²)(u² + v²)`, `n²v² + m²u²` are squares modulo `M`. -/
def sqTriple (m n M S u v : ℕ) : Bool :=
  S.testBit ((m * m * v * v + n * n * u * u) % M) &&
    S.testBit ((2 * (m * m + n * n) * (u * u + v * v)) % M) &&
    S.testBit ((n * n * v * v + m * m * u * u) % M)

/-- `obstructsProj m n M p S`: no pair `(1, v)`, and no pair `(u, 1)` with `p ∣ u`, makes the
three integers squares modulo `M`. -/
def obstructsProj (m n M p S : ℕ) : Bool :=
  ((List.range M).all fun v => !sqTriple m n M S 1 v) &&
    ((List.range M).all fun u => !(u % p == 0) || !sqTriple m n M S u 1)

private theorem sqTriple_of_int {m n m' n' M u v : ℕ} (hM : 0 < M)
    (hm : (m' : ℤ) ≡ (m : ℤ) [ZMOD M]) (hn : (n' : ℤ) ≡ (n : ℤ) [ZMOD M])
    {X Y a b c : ℤ} (hu : (u : ℤ) ≡ X [ZMOD M]) (hv : (v : ℤ) ≡ Y [ZMOD M])
    (ha : (m : ℤ) * m * Y * Y + n * n * X * X = a * a)
    (hb : 2 * ((m : ℤ) * m + n * n) * (X * X + Y * Y) = b * b)
    (hc : (n : ℤ) * n * Y * Y + m * m * X * X = c * c) :
    sqTriple m' n' M (sqMask M) u v = true := by
  have s0 : (sqMask M).testBit ((m' * m' * v * v + n' * n' * u * u) % M) = true := by
    apply testBit_sqMask_of_int hM (r := a)
    rw [← ha]
    push_cast
    exact (((hm.symm.mul hm.symm).mul hv.symm).mul hv.symm).add
      (((hn.symm.mul hn.symm).mul hu.symm).mul hu.symm)
  have s4 : (sqMask M).testBit ((2 * (m' * m' + n' * n') * (u * u + v * v)) % M) = true := by
    apply testBit_sqMask_of_int hM (r := b)
    rw [← hb]
    push_cast
    exact ((Int.ModEq.refl 2).mul ((hm.symm.mul hm.symm).add (hn.symm.mul hn.symm))).mul
      ((hu.symm.mul hu.symm).add (hv.symm.mul hv.symm))
  have s8 : (sqMask M).testBit ((n' * n' * v * v + m' * m' * u * u) % M) = true := by
    apply testBit_sqMask_of_int hM (r := c)
    rw [← hc]
    push_cast
    exact (((hn.symm.mul hn.symm).mul hv.symm).mul hv.symm).add
      (((hm.symm.mul hm.symm).mul hu.symm).mul hu.symm)
  simp [sqTriple, s0, s4, s8]

private theorem core_proj {m n m' n' M p S : ℕ} (hp : p.Prime) (hM1 : 1 < M) (hpM : p ∣ M)
    (hcopM : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x M) (hS : S = sqMask M)
    (hm : (m' : ℤ) ≡ (m : ℤ) [ZMOD M]) (hn : (n' : ℤ) ≡ (n : ℤ) [ZMOD M])
    (h : obstructsProj m' n' M p S = true) (t : ℚ) : ¬DiagonalSquares m n t := by
  subst hS
  intro hd
  have hM : 0 < M := by omega
  obtain ⟨U, V, -, hcop, ⟨a, ha⟩, ⟨b, hb⟩, ⟨c, hc⟩⟩ := exists_coprime_squares hd
  obtain ⟨u₀, hu₀⟩ : ∃ u₀ : ℕ, u₀ = U.natAbs := ⟨_, rfl⟩
  have hUU : U * U = (u₀ : ℤ) * u₀ := by rw [hu₀, Int.natCast_natAbs, abs_mul_abs_self]
  have ha' : (m : ℤ) * m * V * V + n * n * ((u₀ : ℤ) * u₀) = a * a := by
    rw [← hUU, ← ha]
    ring
  have hb' : 2 * ((m : ℤ) * m + n * n) * ((u₀ : ℤ) * u₀ + V * V) = b * b := by
    rw [← hUU, ← hb]
  have hc' : (n : ℤ) * n * V * V + m * m * ((u₀ : ℤ) * u₀) = c * c := by
    rw [← hUU, ← hc]
    ring
  simp only [obstructsProj, Bool.and_eq_true] at h
  obtain ⟨h1, h2⟩ := h
  by_cases hpu : p ∣ u₀
  · have hpV : ¬p ∣ V := by
      intro hpV
      have hg := Nat.dvd_gcd hpu hpV
      rw [hu₀, hcop.gcd_eq_one] at hg
      exact hp.one_lt.ne' (Nat.dvd_one.mp hg)
    obtain ⟨w, -, hw⟩ := Nat.exists_mul_mod_eq_one_of_coprime (hcopM V hpV) hM1
    have hone : (V : ℤ) * w ≡ 1 [ZMOD M] := by
      have hnat : V * w ≡ 1 [MOD M] := by
        show V * w % M = 1 % M
        rw [hw, Nat.mod_eq_of_lt hM1]
      exact_mod_cast Int.natCast_modEq_iff.mpr hnat
    have hu' : (((u₀ * w) % M : ℕ) : ℤ) ≡ (u₀ : ℤ) * w [ZMOD M] := by
      push_cast
      exact Int.mod_modEq _ _
    have hv' : ((1 : ℕ) : ℤ) ≡ (V : ℤ) * w [ZMOD M] := by
      push_cast
      exact hone.symm
    have hdiv : (u₀ * w) % M % p = 0 :=
      Nat.mod_eq_zero_of_dvd ((Nat.dvd_mod_iff hpM).mpr (Dvd.dvd.mul_right hpu w))
    have hs : sqTriple m' n' M (sqMask M) ((u₀ * w) % M) 1 = true :=
      sqTriple_of_int hM hm hn (a := a * w) (b := b * w) (c := c * w) hu' hv'
        (by linear_combination ((w : ℤ) * w) * ha') (by linear_combination ((w : ℤ) * w) * hb')
        (by linear_combination ((w : ℤ) * w) * hc')
    have hbody := List.all_eq_true.mp h2 ((u₀ * w) % M) (List.mem_range.mpr (Nat.mod_lt _ hM))
    simp [hdiv, hs] at hbody
  · obtain ⟨w, -, hw⟩ := Nat.exists_mul_mod_eq_one_of_coprime (hcopM u₀ hpu) hM1
    have hone : (u₀ : ℤ) * w ≡ 1 [ZMOD M] := by
      have hnat : u₀ * w ≡ 1 [MOD M] := by
        show u₀ * w % M = 1 % M
        rw [hw, Nat.mod_eq_of_lt hM1]
      exact_mod_cast Int.natCast_modEq_iff.mpr hnat
    have hv' : (((V * w) % M : ℕ) : ℤ) ≡ (V : ℤ) * w [ZMOD M] := by
      push_cast
      exact Int.mod_modEq _ _
    have hu' : ((1 : ℕ) : ℤ) ≡ (u₀ : ℤ) * w [ZMOD M] := by
      push_cast
      exact hone.symm
    have hs : sqTriple m' n' M (sqMask M) 1 ((V * w) % M) = true :=
      sqTriple_of_int hM hm hn (a := a * w) (b := b * w) (c := c * w) hu' hv'
        (by linear_combination ((w : ℤ) * w) * ha') (by linear_combination ((w : ℤ) * w) * hb')
        (by linear_combination ((w : ℤ) * w) * hc')
    have hbody := List.all_eq_true.mp h1 ((V * w) % M) (List.mem_range.mpr (Nat.mod_lt _ hM))
    simp [hs] at hbody

/-- For a prime `p` and `M` a power of `p` (expressed by `hcopM`): if the cheaper check
succeeds for the residues of `m` and `n` modulo `M`, then no rational `t` makes the three
diagonal cells of the G1 family squares. -/
theorem no_point_of_obstructsProj_mod {m n M p S : ℕ} (hp : p.Prime) (hM1 : 1 < M)
    (hpM : p ∣ M) (hcopM : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x M) (hS : S = sqMask M)
    (h : obstructsProj (m % M) (n % M) M p S = true) (t : ℚ) : ¬DiagonalSquares m n t :=
  core_proj hp hM1 hpM hcopM hS (by push_cast; exact Int.mod_modEq _ _)
    (by push_cast; exact Int.mod_modEq _ _) h t

/-- A number not divisible by the prime `p` is coprime to every power of `p`. -/
theorem coprime_pow_of_not_dvd {p : ℕ} (hp : p.Prime) (k : ℕ) :
    ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x (p ^ k) :=
  fun _ hx => Nat.Coprime.pow_right k (hp.coprime_iff_not_dvd.mpr hx).symm

end MagicSquareOfSquares.G1
