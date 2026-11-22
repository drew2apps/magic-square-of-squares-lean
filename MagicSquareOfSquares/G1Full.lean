import MagicSquareOfSquares.G1

/-!
# The G1 family with all nine cells

`G1.lean` treats the three diagonal cells. Here all five free cells are required to be squares:

```
F0 = m² + n²t²,  F4 = (m² + n²)(t² + 1)/2,  F8 = n² + m²t²,  L = F4 − 2mnt,  U = F4 + 2mnt.
```

The other four cells of the square, `(mt ± n)²` and `(nt ± m)²`, are squares identically.
The five-form test is the one used in the Kepler71 notes
(`PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md`).

## Main definitions and results

* `FullSquares m n t`: the five cells are squares of rationals.
* `obstructsFull m n M p S`: a finite check modulo `M`, a power of the prime `p`.
* `no_solution_of_obstructsFull_mod`: if the check succeeds for `m % M`, `n % M`, no rational
  `t` satisfies `FullSquares m n t`.
* `no_solution_of_ratio`: for a prime `p`, if the check succeeds for `(r, 1)`, the same holds
  for every `(m, n)` with `m ≡ r·n (mod p)` and `p ∤ n`.
-/

namespace MagicSquareOfSquares.G1

/-- The five free cells of the G1 family are squares of rationals. -/
def FullSquares (m n : ℕ) (t : ℚ) : Prop :=
  DiagonalSquares m n t ∧
    IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2 - 2 * m * n * t) ∧
    IsSquare (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2 + 2 * m * n * t)

theorem FullSquares.diagonal {m n : ℕ} {t : ℚ} (h : FullSquares m n t) : DiagonalSquares m n t :=
  h.1

/-- Clearing denominators: coprime natural numbers `u`, `V` with five integer squares. -/
theorem exists_coprime_squares_full {m n : ℕ} {t : ℚ} (h : FullSquares m n t) :
    ∃ (u V : ℕ), 0 < V ∧ Nat.Coprime u V ∧
      IsSquare ((m : ℤ) * m * V * V + n * n * (u * u)) ∧
      IsSquare (2 * ((m : ℤ) * m + n * n) * (u * u + V * V)) ∧
      IsSquare ((n : ℤ) * n * V * V + m * m * (u * u)) ∧
      IsSquare (2 * ((m : ℤ) * m + n * n) * (u * u + V * V) - 8 * m * n * u * V) ∧
      IsSquare (2 * ((m : ℤ) * m + n * n) * (u * u + V * V) + 8 * m * n * u * V) := by
  obtain ⟨⟨h0, h4, h8⟩, hL, hU⟩ := h
  obtain ⟨U, hUdef⟩ : ∃ U : ℤ, U = t.num := ⟨_, rfl⟩
  obtain ⟨V, hVdef⟩ : ∃ V : ℕ, V = t.den := ⟨_, rfl⟩
  have hVpos : 0 < V := by rw [hVdef]; exact t.den_pos
  have hcop : Nat.Coprime U.natAbs V := by rw [hUdef, hVdef]; exact t.reduced
  have hUt : (U : ℚ) = t * V := by
    rw [hUdef, hVdef]
    exact (Rat.mul_den_eq_num t).symm
  have hA : IsSquare ((m : ℤ) * m * V * V + n * n * (U * U)) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : (((m : ℤ) * m * V * V + n * n * (U * U) : ℤ) : ℚ) =
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
  have hC : IsSquare ((n : ℤ) * n * V * V + m * m * (U * U)) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : (((n : ℤ) * n * V * V + m * m * (U * U) : ℤ) : ℚ) =
        ((n : ℚ) ^ 2 + (m : ℚ) ^ 2 * t ^ 2) * ((V : ℚ) * V) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact h8.mul ⟨(V : ℚ), rfl⟩
  have hD : IsSquare (2 * ((m : ℤ) * m + n * n) * (U * U + V * V) - 8 * m * n * U * V) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : ((2 * ((m : ℤ) * m + n * n) * (U * U + V * V) - 8 * m * n * U * V : ℤ) : ℚ) =
        (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2 - 2 * m * n * t) *
          ((2 * (V : ℚ)) * (2 * V)) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact hL.mul ⟨2 * (V : ℚ), rfl⟩
  have hE : IsSquare (2 * ((m : ℤ) * m + n * n) * (U * U + V * V) + 8 * m * n * U * V) := by
    rw [← Rat.isSquare_intCast_iff]
    have e : ((2 * ((m : ℤ) * m + n * n) * (U * U + V * V) + 8 * m * n * U * V : ℤ) : ℚ) =
        (((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (t ^ 2 + 1) / 2 + 2 * m * n * t) *
          ((2 * (V : ℚ)) * (2 * V)) := by
      push_cast
      rw [hUt]
      ring
    rw [e]
    exact hU.mul ⟨2 * (V : ℚ), rfl⟩
  obtain ⟨u, hu⟩ : ∃ u : ℕ, u = U.natAbs := ⟨_, rfl⟩
  have hUU : U * U = (u : ℤ) * u := by rw [hu, Int.natCast_natAbs, abs_mul_abs_self]
  have hsign : U = (u : ℤ) ∨ U = -(u : ℤ) := by rw [hu]; exact Int.natAbs_eq U
  rw [hUU] at hA hB hC
  refine ⟨u, V, hVpos, by rw [hu]; exact hcop, hA, hB, hC, ?_, ?_⟩
  · rcases hsign with hs | hs
    · rw [hs] at hD
      exact hD
    · have e : 2 * ((m : ℤ) * m + n * n) * ((u : ℤ) * u + V * V) - 8 * m * n * u * V =
          2 * ((m : ℤ) * m + n * n) * (U * U + V * V) + 8 * m * n * U * V := by
        rw [hs]
        ring
      rw [e]
      exact hE
  · rcases hsign with hs | hs
    · rw [hs] at hE
      exact hE
    · have e : 2 * ((m : ℤ) * m + n * n) * ((u : ℤ) * u + V * V) + 8 * m * n * u * V =
          2 * ((m : ℤ) * m + n * n) * (U * U + V * V) - 8 * m * n * U * V := by
        rw [hs]
        ring
      rw [e]
      exact hD

/-! ## The finite check -/

/-- `sqFive m n M S u v`: with `S` the mask of squares modulo `M`, the five integers are squares
modulo `M`. The fourth one, `2(m² + n²)(u² + v²) − 8mnuv`, is computed without subtraction. -/
def sqFive (m n M S u v : ℕ) : Bool :=
  sqTriple m n M S u v &&
    S.testBit ((2 * (m * m + n * n) * (u * u + v * v) + (M - 8 * m * n * u * v % M)) % M) &&
    S.testBit ((2 * (m * m + n * n) * (u * u + v * v) + 8 * m * n * u * v) % M)

/-- `obstructsFull m n M p S`: no pair `(1, v)`, and no pair `(u, 1)` with `p ∣ u`, makes the
five integers squares modulo `M`. -/
def obstructsFull (m n M p S : ℕ) : Bool :=
  ((List.range M).all fun v => !sqFive m n M S 1 v) &&
    ((List.range M).all fun u => !(u % p == 0) || !sqFive m n M S u 1)

private theorem modEq_rhs {M : ℕ} {x y z : ℤ} (h : x ≡ y [ZMOD M]) (e : y = z) :
    x ≡ z [ZMOD M] := e ▸ h

/-- A value of a quadratic form in `(m, n)` that is a square stays a square modulo `M` when
`(m, n)` is replaced by `(m', n')` with `(m, n) ≡ c·(m', n')` and `c` invertible. -/
private theorem quad_scale {M : ℕ} {P Q R m n m' n' c c' a0 : ℤ}
    (h : P * (m * m) + Q * (m * n) + R * (n * n) = a0 * a0)
    (hm : m ≡ c * m' [ZMOD M]) (hn : n ≡ c * n' [ZMOD M]) (hc : c * c' ≡ 1 [ZMOD M]) :
    (a0 * c') * (a0 * c') ≡ P * (m' * m') + Q * (m' * n') + R * (n' * n') [ZMOD M] := by
  have e1 : (a0 * c') * (a0 * c') =
      (P * (m * m) + Q * (m * n) + R * (n * n)) * (c' * c') := by
    rw [h]
    ring
  have h1 : P * (m * m) + Q * (m * n) + R * (n * n) ≡
      P * ((c * m') * (c * m')) + Q * ((c * m') * (c * n')) + R * ((c * n') * (c * n'))
      [ZMOD M] :=
    (((Int.ModEq.refl P).mul (hm.mul hm)).add ((Int.ModEq.refl Q).mul (hm.mul hn))).add
      ((Int.ModEq.refl R).mul (hn.mul hn))
  have e2 : (P * ((c * m') * (c * m')) + Q * ((c * m') * (c * n')) +
      R * ((c * n') * (c * n'))) * (c' * c') =
      (c * c') * (c * c') * (P * (m' * m') + Q * (m' * n') + R * (n' * n')) := by ring
  have h2 : (c * c') * (c * c') * (P * (m' * m') + Q * (m' * n') + R * (n' * n')) ≡
      1 * 1 * (P * (m' * m') + Q * (m' * n') + R * (n' * n')) [ZMOD M] :=
    (hc.mul hc).mul_right _
  rw [e1]
  exact modEq_rhs ((modEq_rhs (h1.mul_right _) e2).trans h2) (by ring)

theorem sqFive_of_solution {m n m' n' M u v u₀ V : ℕ} (hM : 0 < M)
    {c c' w a b c3 d e : ℤ}
    (hm : (m : ℤ) ≡ c * m' [ZMOD M]) (hn : (n : ℤ) ≡ c * n' [ZMOD M])
    (hc : c * c' ≡ 1 [ZMOD M])
    (hu : (u : ℤ) ≡ (u₀ : ℤ) * w [ZMOD M]) (hv : (v : ℤ) ≡ (V : ℤ) * w [ZMOD M])
    (ha : (m : ℤ) * m * V * V + n * n * (u₀ * u₀) = a * a)
    (hb : 2 * ((m : ℤ) * m + n * n) * (u₀ * u₀ + V * V) = b * b)
    (hc3 : (n : ℤ) * n * V * V + m * m * (u₀ * u₀) = c3 * c3)
    (hd : 2 * ((m : ℤ) * m + n * n) * (u₀ * u₀ + V * V) - 8 * m * n * u₀ * V = d * d)
    (he : 2 * ((m : ℤ) * m + n * n) * (u₀ * u₀ + V * V) + 8 * m * n * u₀ * V = e * e) :
    sqFive m' n' M (sqMask M) u v = true := by
  obtain ⟨X, hX⟩ : ∃ X : ℤ, X = (u₀ : ℤ) * w := ⟨_, rfl⟩
  obtain ⟨Y, hY⟩ : ∃ Y : ℤ, Y = (V : ℤ) * w := ⟨_, rfl⟩
  rw [← hX] at hu
  rw [← hY] at hv
  have qa : (a * w * c') * (a * w * c') ≡ (m' : ℤ) * m' * Y * Y + n' * n' * X * X [ZMOD M] :=
    modEq_rhs (quad_scale (P := Y * Y) (Q := 0) (R := X * X) (a0 := a * w)
      (by rw [hX, hY]; linear_combination (w * w) * ha) hm hn hc) (by ring)
  have qb : (b * w * c') * (b * w * c') ≡
      2 * ((m' : ℤ) * m' + n' * n') * (X * X + Y * Y) [ZMOD M] :=
    modEq_rhs (quad_scale (P := 2 * (X * X + Y * Y)) (Q := 0) (R := 2 * (X * X + Y * Y))
      (a0 := b * w) (by rw [hX, hY]; linear_combination (w * w) * hb) hm hn hc) (by ring)
  have qc : (c3 * w * c') * (c3 * w * c') ≡ (n' : ℤ) * n' * Y * Y + m' * m' * X * X [ZMOD M] :=
    modEq_rhs (quad_scale (P := X * X) (Q := 0) (R := Y * Y) (a0 := c3 * w)
      (by rw [hX, hY]; linear_combination (w * w) * hc3) hm hn hc) (by ring)
  have qd : (d * w * c') * (d * w * c') ≡
      2 * ((m' : ℤ) * m' + n' * n') * (X * X + Y * Y) + -(8 * m' * n' * X * Y) [ZMOD M] :=
    modEq_rhs (quad_scale (P := 2 * (X * X + Y * Y)) (Q := -(8 * X * Y))
      (R := 2 * (X * X + Y * Y)) (a0 := d * w)
      (by rw [hX, hY]; linear_combination (w * w) * hd) hm hn hc) (by ring)
  have qe : (e * w * c') * (e * w * c') ≡
      2 * ((m' : ℤ) * m' + n' * n') * (X * X + Y * Y) + 8 * m' * n' * X * Y [ZMOD M] :=
    modEq_rhs (quad_scale (P := 2 * (X * X + Y * Y)) (Q := 8 * X * Y)
      (R := 2 * (X * X + Y * Y)) (a0 := e * w)
      (by rw [hX, hY]; linear_combination (w * w) * he) hm hn hc) (by ring)
  have s0 : (sqMask M).testBit ((m' * m' * v * v + n' * n' * u * u) % M) = true :=
    testBit_sqMask_of_int hM (qa.trans (by
      push_cast
      exact ((((Int.ModEq.refl _).mul hv).mul hv).add
        (((Int.ModEq.refl _).mul hu).mul hu)).symm))
  have s4 : (sqMask M).testBit ((2 * (m' * m' + n' * n') * (u * u + v * v)) % M) = true :=
    testBit_sqMask_of_int hM (qb.trans (by
      push_cast
      exact ((Int.ModEq.refl _).mul ((hu.mul hu).add (hv.mul hv))).symm))
  have s8 : (sqMask M).testBit ((n' * n' * v * v + m' * m' * u * u) % M) = true :=
    testBit_sqMask_of_int hM (qc.trans (by
      push_cast
      exact ((((Int.ModEq.refl _).mul hv).mul hv).add
        (((Int.ModEq.refl _).mul hu).mul hu)).symm))
  have sU : (sqMask M).testBit
      ((2 * (m' * m' + n' * n') * (u * u + v * v) + 8 * m' * n' * u * v) % M) = true :=
    testBit_sqMask_of_int hM (qe.trans (by
      push_cast
      exact (((Int.ModEq.refl _).mul ((hu.mul hu).add (hv.mul hv))).add
        (((Int.ModEq.refl _).mul hu).mul hv)).symm))
  have sL : (sqMask M).testBit
      ((2 * (m' * m' + n' * n') * (u * u + v * v) + (M - 8 * m' * n' * u * v % M)) % M) =
      true := by
    have hsub : ((M - 8 * m' * n' * u * v % M : ℕ) : ℤ) =
        (M : ℤ) - ((8 * m' * n' * u * v % M : ℕ) : ℤ) :=
      Nat.cast_sub (Nat.mod_lt _ hM).le
    have hmod : ((8 * m' * n' * u * v % M : ℕ) : ℤ) ≡ 8 * m' * n' * X * Y [ZMOD M] := by
      push_cast
      exact (Int.mod_modEq _ _).trans (((Int.ModEq.refl _).mul hu).mul hv)
    have h1 : ((M - 8 * m' * n' * u * v % M : ℕ) : ℤ) ≡ -(8 * m' * n' * X * Y) [ZMOD M] := by
      rw [hsub, ← zero_sub]
      exact (Int.modEq_zero_iff_dvd.mpr (dvd_refl (M : ℤ))).sub hmod
    have h2 : ((2 * (m' * m' + n' * n') * (u * u + v * v) : ℕ) : ℤ) ≡
        2 * ((m' : ℤ) * m' + n' * n') * (X * X + Y * Y) [ZMOD M] := by
      push_cast
      exact (Int.ModEq.refl _).mul ((hu.mul hu).add (hv.mul hv))
    refine testBit_sqMask_of_int hM (qd.trans ?_)
    rw [Nat.cast_add]
    exact (h2.add h1).symm
  unfold sqFive sqTriple
  rw [s0, s4, s8, sL, sU]
  rfl

private theorem core_full {m n m' n' M p S : ℕ} (hp : p.Prime) (hM1 : 1 < M) (hpM : p ∣ M)
    (hcopM : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x M) (hS : S = sqMask M) {c c' : ℤ}
    (hm : (m : ℤ) ≡ c * m' [ZMOD M]) (hn : (n : ℤ) ≡ c * n' [ZMOD M])
    (hc : c * c' ≡ 1 [ZMOD M])
    (h : obstructsFull m' n' M p S = true) (t : ℚ) : ¬FullSquares m n t := by
  subst hS
  intro hf
  have hM : 0 < M := by omega
  obtain ⟨u₀, V, -, hcop, ⟨a, ha⟩, ⟨b, hb⟩, ⟨c3, hc3⟩, ⟨d, hd⟩, ⟨e, he⟩⟩ :=
    exists_coprime_squares_full hf
  simp only [obstructsFull, Bool.and_eq_true] at h
  obtain ⟨h1, h2⟩ := h
  by_cases hpu : p ∣ u₀
  · have hpV : ¬p ∣ V := by
      intro hpV
      have hg := Nat.dvd_gcd hpu hpV
      rw [hcop.gcd_eq_one] at hg
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
    have hs := sqFive_of_solution hM hm hn hc hu' hv' ha hb hc3 hd he
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
    have hs := sqFive_of_solution hM hm hn hc hu' hv' ha hb hc3 hd he
    have hbody := List.all_eq_true.mp h1 ((V * w) % M) (List.mem_range.mpr (Nat.mod_lt _ hM))
    simp [hs] at hbody

/-- For `M` a power of the prime `p`: if the check succeeds for the remainders of `m` and `n`,
no rational `t` makes the five cells squares. -/
theorem no_solution_of_obstructsFull_mod {m n M p S : ℕ} (hp : p.Prime) (hM1 : 1 < M)
    (hpM : p ∣ M) (hcopM : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x M) (hS : S = sqMask M)
    (h : obstructsFull (m % M) (n % M) M p S = true) (t : ℚ) : ¬FullSquares m n t := by
  have hm : (m : ℤ) ≡ 1 * ((m % M : ℕ) : ℤ) [ZMOD M] := by
    rw [one_mul]
    push_cast
    exact (Int.mod_modEq _ _).symm
  have hn : (n : ℤ) ≡ 1 * ((n % M : ℕ) : ℤ) [ZMOD M] := by
    rw [one_mul]
    push_cast
    exact (Int.mod_modEq _ _).symm
  exact core_full hp hM1 hpM hcopM hS hm hn (c' := 1) (by rw [one_mul]) h t

/-- For a prime `p`: if the check succeeds for the ratio `(r : 1)`, then no rational `t` makes
the five cells squares for any `(m, n)` with `m ≡ r·n (mod p)` and `p ∤ n`. -/
theorem no_solution_of_ratio {m n p r S : ℕ} (hp : p.Prime) (hS : S = sqMask p)
    (h : obstructsFull r 1 p p S = true) (hmn : m % p = r * n % p) (hn : n % p ≠ 0) (t : ℚ) :
    ¬FullSquares m n t := by
  have hnp : ¬p ∣ n := fun hd => hn (Nat.mod_eq_zero_of_dvd hd)
  have hcopp : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x p :=
    fun _ hx => (hp.coprime_iff_not_dvd.mpr hx).symm
  obtain ⟨w, -, hw⟩ := Nat.exists_mul_mod_eq_one_of_coprime (hcopp n hnp) hp.one_lt
  have hc : (n : ℤ) * w ≡ 1 [ZMOD p] := by
    have hnat : n * w ≡ 1 [MOD p] := by
      show n * w % p = 1 % p
      rw [hw, Nat.mod_eq_of_lt hp.one_lt]
    exact_mod_cast Int.natCast_modEq_iff.mpr hnat
  have hm : (m : ℤ) ≡ (n : ℤ) * (r : ℤ) [ZMOD p] := by
    have h1 : m ≡ r * n [MOD p] := hmn
    have h2 := Int.natCast_modEq_iff.mpr h1
    push_cast at h2
    have e : (r : ℤ) * n = n * r := mul_comm _ _
    rw [e] at h2
    exact h2
  have hn' : (n : ℤ) ≡ (n : ℤ) * ((1 : ℕ) : ℤ) [ZMOD p] := by
    rw [Nat.cast_one, mul_one]
  exact core_full hp hp.one_lt (dvd_refl p) hcopp hS hm hn' hc h t

end MagicSquareOfSquares.G1
