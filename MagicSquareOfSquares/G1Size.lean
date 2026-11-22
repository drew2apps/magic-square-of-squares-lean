import MagicSquareOfSquares.G1SizeCheck

/-!
# A lower bound for the edge roots

The five-form system of the G1 family is symmetric: exchanging `(m, n)` with `(v, u)` leaves the
five integers unchanged. So in a solution both ratios must avoid the excluded classes. For ratios
of height at most 1000 that leaves the 425 exceptional pairs for each, and `G1SizeCheck.lean`
shows that no two exceptional pairs fit together.

## Main result

* `IsMSS.edge_root_gt_500`: in a 3×3 magic square of squares with roots `a b c / d e f / g h i`
  and `b² ≠ d²`, `h² ≠ f²`, `b² ≠ f²`, one of the edge roots `b`, `d`, `f`, `h` is larger than
  500 in absolute value.

Published computer searches exclude far larger squares. This bound is only the one that follows
from the theorems of this repository.
-/

namespace MagicSquareOfSquares.G1

/-- The five integers of the G1 family, cleared of denominators, are squares. -/
def IntSq (m n u v : ℤ) : Prop :=
  IsSquare (m * m * v * v + n * n * (u * u)) ∧
    IsSquare (2 * (m * m + n * n) * (u * u + v * v)) ∧
    IsSquare (n * n * v * v + m * m * (u * u)) ∧
    IsSquare (2 * (m * m + n * n) * (u * u + v * v) - 8 * m * n * u * v) ∧
    IsSquare (2 * (m * m + n * n) * (u * u + v * v) + 8 * m * n * u * v)

private theorem sqeq {x y : ℤ} (h : IsSquare x) (e : y = x) : IsSquare y := e ▸ h

/-- The symmetry between the two ratios. -/
theorem IntSq.swap {m n u v : ℤ} (h : IntSq m n u v) : IntSq v u n m := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := h
  exact ⟨sqeq h1 (by ring), sqeq h2 (by ring), sqeq h3 (by ring), sqeq h4 (by ring),
    sqeq h5 (by ring)⟩

theorem IntSq.flip {m n u v : ℤ} (h : IntSq m n u v) : IntSq n m v u := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := h
  exact ⟨sqeq h1 (by ring), sqeq h2 (by ring), sqeq h3 (by ring), sqeq h4 (by ring),
    sqeq h5 (by ring)⟩

theorem IntSq.neg_right {m n u v : ℤ} (h : IntSq m n u v) : IntSq m n u (-v) := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := h
  exact ⟨sqeq h1 (by ring), sqeq h2 (by ring), sqeq h3 (by ring), sqeq h5 (by ring),
    sqeq h4 (by ring)⟩

theorem IntSq.neg_left {m n u v : ℤ} (h : IntSq m n u v) : IntSq m n (-u) v := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := h
  exact ⟨sqeq h1 (by ring), sqeq h2 (by ring), sqeq h3 (by ring), sqeq h5 (by ring),
    sqeq h4 (by ring)⟩

theorem IntSq.abs_right {m n u v : ℤ} (h : IntSq m n u v) : IntSq m n u |v| := by
  rcases abs_choice v with e | e
  · rw [e]
    exact h
  · rw [e]
    exact h.neg_right

theorem IntSq.abs_left {m n u v : ℤ} (h : IntSq m n u v) : IntSq m n |u| v := by
  rcases abs_choice u with e | e
  · rw [e]
    exact h
  · rw [e]
    exact h.neg_left

private theorem isSquare_of_sq_mul {g k : ℤ} (hg : g ≠ 0) (h : IsSquare (g * g * k)) :
    IsSquare k := by
  obtain ⟨r, hr⟩ := h
  have hd : g ∣ r := by
    have h2 : g.natAbs ^ 2 ∣ r.natAbs ^ 2 := ⟨k.natAbs, by
      have h3 := congrArg Int.natAbs hr
      simp only [Int.natAbs_mul] at h3
      rw [sq, sq, ← h3]⟩
    exact Int.natAbs_dvd_natAbs.mp ((Nat.pow_dvd_pow_iff two_ne_zero).mp h2)
  obtain ⟨q, rfl⟩ := hd
  exact ⟨q, mul_left_cancel₀ (mul_ne_zero hg hg) (by rw [hr]; ring)⟩

/-- A common factor of `m` and `n` can be removed. -/
theorem IntSq.of_mul_left {g m n u v : ℤ} (hg : g ≠ 0) (h : IntSq (g * m) (g * n) u v) :
    IntSq m n u v := by
  obtain ⟨h1, h2, h3, h4, h5⟩ := h
  exact ⟨isSquare_of_sq_mul hg (sqeq h1 (by ring)), isSquare_of_sq_mul hg (sqeq h2 (by ring)),
    isSquare_of_sq_mul hg (sqeq h3 (by ring)), isSquare_of_sq_mul hg (sqeq h4 (by ring)),
    isSquare_of_sq_mul hg (sqeq h5 (by ring))⟩

private theorem isSquare_of_mul_sq' {x D R : ℚ} (hD : D ≠ 0) (h : x * D ^ 2 = R ^ 2) :
    IsSquare x :=
  ⟨R / D, by
    rw [div_mul_div_comm, eq_div_iff (mul_ne_zero hD hD)]
    linear_combination h⟩

/-- An integer solution gives a rational one. -/
theorem IntSq.fullSquares {m n : ℕ} {u v : ℤ} (hv : v ≠ 0) (h : IntSq m n u v) :
    FullSquares m n ((u : ℚ) / v) := by
  obtain ⟨⟨r1, h1⟩, ⟨r2, h2⟩, ⟨r3, h3⟩, ⟨r4, h4⟩, ⟨r5, h5⟩⟩ := h
  have hvq : (v : ℚ) ≠ 0 := by exact_mod_cast hv
  have h2v : (2 : ℚ) * v ≠ 0 := mul_ne_zero two_ne_zero hvq
  obtain ⟨τ, hτdef⟩ : ∃ τ : ℚ, τ = (u : ℚ) / v := ⟨_, rfl⟩
  have hτ : τ * v = u := by
    rw [hτdef, div_mul_eq_mul_div, div_eq_iff hvq]
  rw [← hτdef]
  have q1 : (m : ℚ) * m * v * v + n * n * (u * u) = r1 * r1 := by exact_mod_cast h1
  have q2 : 2 * ((m : ℚ) * m + n * n) * (u * u + v * v) = r2 * r2 := by exact_mod_cast h2
  have q3 : (n : ℚ) * n * v * v + m * m * (u * u) = r3 * r3 := by exact_mod_cast h3
  have q4 : 2 * ((m : ℚ) * m + n * n) * (u * u + v * v) - 8 * m * n * u * v = r4 * r4 := by
    exact_mod_cast h4
  have q5 : 2 * ((m : ℚ) * m + n * n) * (u * u + v * v) + 8 * m * n * u * v = r5 * r5 := by
    exact_mod_cast h5
  refine ⟨⟨?_, ?_, ?_⟩, ?_, ?_⟩
  · exact isSquare_of_mul_sq' hvq (R := r1)
      (by linear_combination ((n : ℚ) ^ 2 * (τ * v + u)) * hτ + q1)
  · exact isSquare_of_mul_sq' h2v (R := r2)
      (by linear_combination (2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (τ * v + u)) * hτ + q2)
  · exact isSquare_of_mul_sq' hvq (R := r3)
      (by linear_combination ((m : ℚ) ^ 2 * (τ * v + u)) * hτ + q3)
  · exact isSquare_of_mul_sq' h2v (R := r4)
      (by
        linear_combination (2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (τ * v + u) - 8 * m * n * v) * hτ +
          q4)
  · exact isSquare_of_mul_sq' h2v (R := r5)
      (by
        linear_combination (2 * ((m : ℚ) ^ 2 + (n : ℚ) ^ 2) * (τ * v + u) + 8 * m * n * v) * hτ +
          q5)

/-- An integer solution passes the check modulo every `M`. -/
theorem sqFive_of_intSq {m n u v M : ℕ} (hM : 0 < M) (h : IntSq m n u v) :
    sqFive m n M (sqMask M) u v = true := by
  obtain ⟨⟨r1, h1⟩, ⟨r2, h2⟩, ⟨r3, h3⟩, ⟨r4, h4⟩, ⟨r5, h5⟩⟩ := h
  exact sqFive_of_solution hM (c := 1) (c' := 1) (w := 1) (u₀ := u) (V := v)
    (by rw [one_mul]) (by rw [one_mul]) (by rw [one_mul]) (by rw [mul_one]) (by rw [mul_one])
    h1 h2 h3 h4 h5

theorem not_intSq_of_badQuad {p q u v : ℕ} (h : badQuad p q u v = true) : ¬IntSq p q u v := by
  intro hI
  unfold badQuad at h
  obtain ⟨⟨P, S⟩, he, h2⟩ := List.any_eq_true.mp h
  have hv := List.all_eq_true.mp sizePrimes_valid _ he
  simp only [decide_eq_true_eq] at hv
  obtain ⟨hP, hS⟩ := hv
  have hs := sqFive_of_intSq (M := P) hP hI
  rw [← hS] at hs
  simp [hs] at h2

/-- Two exceptional pairs do not fit together. -/
private theorem stage2 {M N x y : ℕ} (hMN : (M, N) ∈ exceptions1000) (hx : 0 < x) (hy : 0 < y)
    (hx' : x ≤ 1000) (hy' : y ≤ 1000) (hxy : x ≠ y) (h : IntSq M N x y) : False := by
  have hpos := List.all_eq_true.mp exceptions1000_pos _ hMN
  simp only [decide_eq_true_eq] at hpos
  obtain ⟨hN0, hNM⟩ := hpos
  have hMne : (M : ℤ) ≠ 0 := by
    have : 0 < M := by omega
    exact_mod_cast this.ne'
  have hNne : (N : ℤ) ≠ 0 := by exact_mod_cast hN0.ne'
  have hs := h.swap
  obtain ⟨y0, x0, hcop, hy0, hx0⟩ := Nat.exists_coprime y x
  obtain ⟨G, hG⟩ : ∃ G, G = Nat.gcd y x := ⟨_, rfl⟩
  rw [← hG] at hy0 hx0
  have hGpos : 0 < G := by
    rw [hG]
    exact Nat.gcd_pos_of_pos_left _ hy
  have ey : (y : ℤ) = G * y0 := by
    have h3 := congrArg (Nat.cast : ℕ → ℤ) hy0
    push_cast at h3
    rw [h3]
    ring
  have ex : (x : ℤ) = G * x0 := by
    have h3 := congrArg (Nat.cast : ℕ → ℤ) hx0
    push_cast at h3
    rw [h3]
    ring
  have hs0 : IntSq (y0 : ℤ) x0 N M := by
    refine IntSq.of_mul_left (g := (G : ℤ)) (by exact_mod_cast hGpos.ne') ?_
    rw [← ey, ← ex]
    exact hs
  have hy0le : y0 ≤ 1000 :=
    le_trans (Nat.le_mul_of_pos_right y0 hGpos) (by rw [← hy0]; exact hy')
  have hx0le : x0 ≤ 1000 :=
    le_trans (Nat.le_mul_of_pos_right x0 hGpos) (by rw [← hx0]; exact hx')
  have hy0pos : 0 < y0 := Nat.pos_of_ne_zero (fun h0 => by
    rw [h0, zero_mul] at hy0
    omega)
  have hx0pos : 0 < x0 := Nat.pos_of_ne_zero (fun h0 => by
    rw [h0, zero_mul] at hx0
    omega)
  rcases Nat.lt_trichotomy x0 y0 with hlt | heq | hgt
  · by_cases hex : (y0, x0) ∈ exceptions1000
    · have hq := quadOk_of_mem hMN hex
      simp only [quadOk, Bool.and_eq_true] at hq
      exact not_intSq_of_badQuad hq.1 hs0
    · exact no_solution_le_1000 hx0pos hlt hy0le hcop hex _ (hs0.fullSquares hMne)
  · apply hxy
    rw [hx0, hy0, heq]
  · have hf := hs0.flip
    by_cases hex : (x0, y0) ∈ exceptions1000
    · have hq := quadOk_of_mem hMN hex
      simp only [quadOk, Bool.and_eq_true] at hq
      exact not_intSq_of_badQuad hq.2 hf
    · exact no_solution_le_1000 hy0pos hgt hx0le hcop.symm hex _ (hf.fullSquares hNne)

end MagicSquareOfSquares.G1

namespace MagicSquareOfSquares

open G1

/-- The edge roots of a magic square of squares give an integer solution of the five-form
system, with `(m, n) = (|b| + |d|, |h| + |f|)` and `(u, v) = (|h| + |f|, |b| − |d|)`. -/
theorem IsMSS.intSq {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T) :
    IntSq (|b| + |d|) (|h| + |f|) (|h| + |f|) (|b| - |d|) := by
  obtain ⟨r1, r2, r3, c1, c2, c3, dg, ad⟩ := hs
  obtain ⟨B, hB⟩ : ∃ B : ℤ, B = |b| := ⟨_, rfl⟩
  obtain ⟨D, hD⟩ : ∃ D : ℤ, D = |d| := ⟨_, rfl⟩
  obtain ⟨H, hH⟩ : ∃ H : ℤ, H = |h| := ⟨_, rfl⟩
  obtain ⟨F, hF⟩ : ∃ F : ℤ, F = |f| := ⟨_, rfl⟩
  have hB2 : B ^ 2 = b ^ 2 := by rw [hB, sq_abs]
  have hD2 : D ^ 2 = d ^ 2 := by rw [hD, sq_abs]
  have hH2 : H ^ 2 = h ^ 2 := by rw [hH, sq_abs]
  have hF2 : F ^ 2 = f ^ 2 := by rw [hF, sq_abs]
  rw [← hB, ← hD, ← hH, ← hF]
  have e1 : B ^ 2 + H ^ 2 = D ^ 2 + F ^ 2 := by rw [hB2, hH2, hD2, hF2]; linarith
  have ea : 2 * a ^ 2 = F ^ 2 + H ^ 2 := by rw [hF2, hH2]; linarith
  have ec : 2 * c ^ 2 = D ^ 2 + H ^ 2 := by rw [hD2, hH2]; linarith
  have eg : 2 * g ^ 2 = B ^ 2 + F ^ 2 := by rw [hB2, hF2]; linarith
  have ei : 2 * i ^ 2 = B ^ 2 + D ^ 2 := by rw [hB2, hD2]; linarith
  have ee : 2 * e ^ 2 = B ^ 2 + H ^ 2 := by rw [hB2, hH2]; linarith
  refine ⟨⟨2 * (H + F) * a, ?_⟩, ⟨4 * (H + F) * e, ?_⟩, ⟨2 * (H + F) * i, ?_⟩,
    ⟨4 * (H + F) * c, ?_⟩, ⟨4 * (H + F) * g, ?_⟩⟩
  · linear_combination ((B + D) * (B - D) + (F - H) * (H + F)) * e1 - 2 * (H + F) ^ 2 * ea
  · linear_combination (2 * ((B + D) * (B - D) + (F - H) * (H + F)) - 4 * (H + F) ^ 2) * e1 -
      8 * (H + F) ^ 2 * ee
  · linear_combination (-2 * (H + F) ^ 2) * ei
  · linear_combination (2 * ((B + D) * (B - D) + (F - H) * (H + F)) - 4 * (H + F) ^ 2) * e1 -
      8 * (H + F) ^ 2 * ec
  · linear_combination (2 * ((B + D) * (B - D) + (F - H) * (H + F)) + 4 * (H + F) ^ 2) * e1 -
      8 * (H + F) ^ 2 * eg

/-- The relation between the four edge roots. -/
theorem IsMSS.abs_edges {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T) :
    |b| ^ 2 + |h| ^ 2 = |d| ^ 2 + |f| ^ 2 := by
  obtain ⟨r1, r2, r3, c1, c2, c3, dg, ad⟩ := hs
  rw [sq_abs, sq_abs, sq_abs, sq_abs]
  linarith

private theorem degenerate {B D H F : ℤ} (e1 : B ^ 2 + H ^ 2 = D ^ 2 + F ^ 2)
    (hds : (H + F) ^ 2 = (B - D) ^ 2) (hs : H + F ≠ 0) : B ^ 2 = F ^ 2 := by
  have h1 : (H + F) ^ 2 * ((B + D) ^ 2 - (F - H) ^ 2) = 0 := by
    linear_combination ((B + D) ^ 2) * hds +
      (2 * (B - D) * (B + D) - (B ^ 2 + H ^ 2 - D ^ 2 - F ^ 2)) * e1
  have h2 : (B + D) ^ 2 - (F - H) ^ 2 = 0 := by
    rcases mul_eq_zero.mp h1 with h3 | h3
    · exact absurd (pow_eq_zero_iff two_ne_zero |>.mp h3) hs
    · exact h3
  linarith

/-- **In a 3×3 magic square of squares with roots `a b c / d e f / g h i` whose edge cells
satisfy `b² ≠ d²`, `h² ≠ f²` and `b² ≠ f²`, one of the edge roots is larger than 500 in absolute
value.** -/
theorem IsMSS.edge_root_gt_500 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hbd : b ^ 2 ≠ d ^ 2) (hhf : h ^ 2 ≠ f ^ 2) (hbf : b ^ 2 ≠ f ^ 2) :
    500 < |b| ∨ 500 < |d| ∨ 500 < |f| ∨ 500 < |h| := by
  by_contra hcon
  simp only [not_or, not_lt] at hcon
  obtain ⟨cb, cd, cf, ch⟩ := hcon
  obtain ⟨B, hB⟩ : ∃ B : ℕ, (B : ℤ) = |b| := ⟨b.natAbs, Int.natCast_natAbs b⟩
  obtain ⟨D, hD⟩ : ∃ D : ℕ, (D : ℤ) = |d| := ⟨d.natAbs, Int.natCast_natAbs d⟩
  obtain ⟨H, hH⟩ : ∃ H : ℕ, (H : ℤ) = |h| := ⟨h.natAbs, Int.natCast_natAbs h⟩
  obtain ⟨F, hF⟩ : ∃ F : ℕ, (F : ℤ) = |f| := ⟨f.natAbs, Int.natCast_natAbs f⟩
  have hI := hs.intSq
  have e1 := hs.abs_edges
  rw [← hB, ← hD, ← hH, ← hF] at hI e1
  rw [← hB] at cb
  rw [← hD] at cd
  rw [← hF] at cf
  rw [← hH] at ch
  have hBle : B ≤ 500 := by exact_mod_cast cb
  have hDle : D ≤ 500 := by exact_mod_cast cd
  have hFle : F ≤ 500 := by exact_mod_cast cf
  have hHle : H ≤ 500 := by exact_mod_cast ch
  -- distinctness in terms of the absolute values
  have nBD : (B : ℤ) ^ 2 ≠ (D : ℤ) ^ 2 := by rw [hB, hD, sq_abs, sq_abs]; exact hbd
  have nHF : (H : ℤ) ^ 2 ≠ (F : ℤ) ^ 2 := by rw [hH, hF, sq_abs, sq_abs]; exact hhf
  have nBF : (B : ℤ) ^ 2 ≠ (F : ℤ) ^ 2 := by rw [hB, hF, sq_abs, sq_abs]; exact hbf
  have hBD : B ≠ D := fun h0 => nBD (by rw [h0])
  have hHF : H ≠ F := fun h0 => nHF (by rw [h0])
  have hQpos : 0 < H + F := by omega
  have hPpos : 0 < B + D := by omega
  have hsne : (H : ℤ) + F ≠ 0 := by exact_mod_cast hQpos.ne'
  have hδne : (B : ℤ) - D ≠ 0 := sub_ne_zero.mpr (by exact_mod_cast hBD)
  obtain ⟨y, hy⟩ : ∃ y : ℕ, (y : ℤ) = |(B : ℤ) - D| :=
    ⟨((B : ℤ) - D).natAbs, Int.natCast_natAbs _⟩
  have hypos : 0 < y := by
    have : (0 : ℤ) < y := by rw [hy]; exact abs_pos.mpr hδne
    exact_mod_cast this
  have hyle : y ≤ 1000 := by
    have : (y : ℤ) ≤ 1000 := by
      rw [hy, abs_le]
      constructor <;> omega
    exact_mod_cast this
  have hxy : H + F ≠ y := by
    intro heq
    have h3 : ((H : ℤ) + F) ^ 2 = ((B : ℤ) - D) ^ 2 := by
      have h4 := congrArg (Nat.cast : ℕ → ℤ) heq
      push_cast at h4
      rw [h4, hy, sq_abs]
    exact nBF (degenerate e1 h3 hsne)
  -- reduce the first ratio
  obtain ⟨m0, n0, hcop, hP, hQ⟩ := Nat.exists_coprime (B + D) (H + F)
  obtain ⟨G, hG⟩ : ∃ G, G = Nat.gcd (B + D) (H + F) := ⟨_, rfl⟩
  rw [← hG] at hP hQ
  have hGpos : 0 < G := by
    rw [hG]
    exact Nat.gcd_pos_of_pos_left _ hPpos
  have eP : (B : ℤ) + D = G * m0 := by
    have h3 := congrArg (Nat.cast : ℕ → ℤ) hP
    push_cast at h3
    rw [h3]
    ring
  have eQ : (H : ℤ) + F = G * n0 := by
    have h3 := congrArg (Nat.cast : ℕ → ℤ) hQ
    push_cast at h3
    rw [h3]
    ring
  obtain ⟨Qz, hQz⟩ : ∃ Qz : ℤ, Qz = (H : ℤ) + F := ⟨_, rfl⟩
  have hQzne : Qz ≠ 0 := by rw [hQz]; exact hsne
  have hQzcast : ((H + F : ℕ) : ℤ) = Qz := by rw [hQz]; push_cast; rfl
  have hI0 : IntSq (m0 : ℤ) n0 Qz ((B : ℤ) - D) := by
    refine IntSq.of_mul_left (g := (G : ℤ)) (by exact_mod_cast hGpos.ne') ?_
    rw [← eP, ← eQ, hQz]
    exact hI
  have hm0le : m0 ≤ 1000 :=
    le_trans (Nat.le_mul_of_pos_right m0 hGpos) (by rw [← hP]; omega)
  have hn0le : n0 ≤ 1000 :=
    le_trans (Nat.le_mul_of_pos_right n0 hGpos) (by rw [← hQ]; omega)
  have hm0pos : 0 < m0 := Nat.pos_of_ne_zero (fun h0 => by
    rw [h0, zero_mul] at hP
    omega)
  have hn0pos : 0 < n0 := Nat.pos_of_ne_zero (fun h0 => by
    rw [h0, zero_mul] at hQ
    omega)
  rcases Nat.lt_trichotomy n0 m0 with hlt | heq | hgt
  · by_cases hex : (m0, n0) ∈ exceptions1000
    · refine stage2 hex hQpos hypos (by omega) hyle hxy ?_
      rw [hy, hQzcast]
      exact hI0.abs_right
    · exact no_solution_le_1000 hn0pos hlt hm0le hcop hex _ (hI0.fullSquares hδne)
  · have hPQ : (B : ℤ) + D = H + F := by rw [eP, eQ, heq]
    have h3 : ((B : ℤ) - F) * (B + F + D + H) = 0 := by
      linear_combination e1 + ((D : ℤ) + H) * hPQ
    rcases mul_eq_zero.mp h3 with h4 | h4
    · exact nBF (by rw [sub_eq_zero.mp h4])
    · have : (0 : ℤ) < B + F + D + H := by exact_mod_cast (by omega : 0 < B + F + D + H)
      omega
  · have hf := hI0.flip
    by_cases hex : (n0, m0) ∈ exceptions1000
    · refine stage2 hex hypos hQpos hyle (by omega) (fun h0 => hxy h0.symm) ?_
      rw [hy, hQzcast]
      exact hf.abs_left
    · exact no_solution_le_1000 hm0pos hgt hn0le hcop.symm hex _ (hf.fullSquares hQzne)

end MagicSquareOfSquares
