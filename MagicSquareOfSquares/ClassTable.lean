import MagicSquareOfSquares.ClassFive

/-!
# The table of Theorem C as lower bounds, at every odd prime of `mn(m² − n²)`

* `Passes p k d m n`: there is a point `(u, v)`, primitive modulo `p`, at which the four norms have
  square roots modulo `p^(k+1)` not divisible by `p³` and `d` times each of the eight products is a
  Gaussian square modulo `p^(k+1)` (`ClassAt`, in `ClassLocal.lean`).
* `Row p m d`: the condition of the table of Theorem C of `COMPLETE.md` on a class `d` at an odd
  prime `p` dividing `m`: `d` prime to `p` (and then, if `p = 5`, either `d ≡ ±2 (mod 5)` or 25
  divides `m`), or `d = p·d₁` with `p ∤ d₁` and `p ≡ ±1 (mod 8)`.
* `passes_of_row`: **every class of the table passes, modulo every power of `p`**, for `p | m`.
* `Passes.swap`, `Passes.rot`: the symmetries `(m, n) ↦ (n, m)` and
  `(m, n, u, v) ↦ (m + n, m − n, u + v, u − v)`; `passes_of_row_n`, `passes_of_row_add`,
  `passes_of_row_sub`: the same statement for `p` dividing `n`, `m + n`, `m − n`.
-/

namespace MagicSquareOfSquares.CoverDescent

/-- `d` passes for the ratio `m : n` modulo `p^(k+1)`: there is a point `(u, v)`, primitive modulo
`p`, at which the four norms have square roots modulo `p^(k+1)` that are not divisible by `p³`, and
`d` times each of the eight products is a Gaussian square modulo `p^(k+1)` (`ClassAt`).

The first version of this definition (commit `fdaf0b5`) had only the last condition. That was
weaker than "`d` is the class of a point of `X_p`": for `p ≡ 3 (mod 8)` with `p ∥ m`, `d = p`
satisfied it at `(0, 1)` for every `k`, although `p` is not a class at such a prime
(`COMPLETE.md`, section 5). With the roots it does not: `not_passes_three` below. -/
def Passes (p k : ℕ) (d m n : ℤ) : Prop :=
  ∃ u v : ℤ, (¬(p : ℤ) ∣ u ∨ ¬(p : ℤ) ∣ v) ∧
    ClassAt ((p : ℤ) ^ (k + 1)) ((p : ℤ) ^ 3) d m n u v

/-- The classes of the table of Theorem C at an odd prime `p` dividing `m`. -/
def Row (p : ℕ) (m d : ℤ) : Prop :=
  (¬(p : ℤ) ∣ d ∧ (p ≠ 5 ∨ d % 5 = 2 ∨ d % 5 = 3 ∨ (25 : ℤ) ∣ m)) ∨
    ∃ d1 : ℤ, d = (p : ℤ) * d1 ∧ ¬(p : ℤ) ∣ d1 ∧ (p % 8 = 1 ∨ p % 8 = 7)

/-- **Every class of the table passes at an odd prime `p` dividing `m`, modulo every power of
`p`.** -/
theorem passes_of_row {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n d : ℤ} (hm : (p : ℤ) ∣ m)
    (hn : ¬(p : ℤ) ∣ n) (hrow : Row p m d) (k : ℕ) : Passes p k d m n := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h1 : ¬(p : ℤ) ∣ 1 := fun h => hp'.not_isUnit (isUnit_of_dvd_one h)
  rcases hrow with ⟨hd, h5⟩ | ⟨d1, rfl, hd1, h8⟩
  · by_cases hp5 : p = 5
    · subst hp5
      have hn5 : ¬(5 : ℤ) ∣ n := by exact_mod_cast hn
      have hd5 : ¬(5 : ℤ) ∣ d := by exact_mod_cast hd
      have hm5 : (5 : ℤ) ∣ m := by exact_mod_cast hm
      by_cases h23 : d % 5 = 2 ∨ d % 5 = 3
      · exact ⟨1, 1, Or.inl h1, by exact_mod_cast passes_at_five hm5 hn5 h23 k⟩
      · have h14 : d % 5 = 1 ∨ d % 5 = 4 := by omega
        have h25 : (25 : ℤ) ∣ m := by
          rcases h5 with h | h | h | h
          · exact absurd rfl h
          · exact absurd (Or.inl h) h23
          · exact absurd (Or.inr h) h23
          · exact h
        obtain ⟨m1, rfl⟩ := h25
        by_cases hm1 : (5 : ℤ) ∣ m1
        · obtain ⟨m3, rfl⟩ := hm1
          refine ⟨1, 7, Or.inl h1, ?_⟩
          have e : (25 : ℤ) * (5 * m3) = 125 * m3 := by ring
          rw [e]
          exact_mod_cast class_one_passes_at_five_of_125 (m3 := m3) hn5 h14 k
        · have hq : ¬(5 : ℤ) ∣ m1 * n := by
            have hp5' : Prime (5 : ℤ) := by exact_mod_cast hp'
            exact not_dvd_mul_of_prime hp5' hm1 hn5
          generalize hq' : m1 * n = q at hq
          have hcases : q % 5 = 1 ∨ q % 5 = 4 ∨ q % 5 = 2 ∨ q % 5 = 3 := by omega
          rcases hcases with h | h | h | h
          · refine ⟨1, 57, Or.inl h1, ?_⟩
            have hτ : (1 : ℤ) ^ 2 ≡ m1 * n [ZMOD 5] := by
              rw [hq']
              unfold Int.ModEq
              omega
            exact_mod_cast class_one_passes_at_five_of_25 hn5 hτ (by norm_num) h14 k
          · refine ⟨1, 57, Or.inl h1, ?_⟩
            have hτ : (2 : ℤ) ^ 2 ≡ m1 * n [ZMOD 5] := by
              rw [hq']
              unfold Int.ModEq
              omega
            exact_mod_cast class_one_passes_at_five_of_25 hn5 hτ (by norm_num) h14 k
          · have hτ : (1 : ℤ) ^ 2 ≡ 3 * (m1 * n) [ZMOD 5] := by
              rw [hq']
              unfold Int.ModEq
              omega
            obtain ⟨v, hv⟩ := class_one_passes_at_five_of_25_nonres hn5 hτ (by norm_num) h14 k
            exact ⟨1, v, Or.inl h1, by exact_mod_cast hv⟩
          · have hτ : (2 : ℤ) ^ 2 ≡ 3 * (m1 * n) [ZMOD 5] := by
              rw [hq']
              unfold Int.ModEq
              omega
            obtain ⟨v, hv⟩ := class_one_passes_at_five_of_25_nonres hn5 hτ (by norm_num) h14 k
            exact ⟨1, v, Or.inl h1, by exact_mod_cast hv⟩
    · have h4 : p % 4 = 1 ∨ p % 4 = 3 := by
        have h2 := hp.eq_one_or_self_of_dvd 2
        have h3 : p % 2 = 1 := by
          by_contra hc
          have h6 : 2 ∣ p := by omega
          rcases h2 h6 with h7 | h7 <;> omega
        omega
      rcases h4 with h4 | h4
      · obtain ⟨v, hv⟩ := unit_passes_of_one_mod_four_ne_five hp h4 hp5 hm hn hd k
        exact ⟨1, v, Or.inl h1, hv⟩
      · exact ⟨1, 1, Or.inl h1, unit_passes_of_three_mod_four hp h4 hm hn hd k⟩
  · obtain ⟨m1, rfl⟩ := hm
    by_cases hm1 : (p : ℤ) ∣ m1
    · obtain ⟨m2, rfl⟩ := hm1
      refine ⟨(p : ℤ) * (2 * d1), 1, Or.inr h1, ?_⟩
      have e : (p : ℤ) * ((p : ℤ) * m2) = (p : ℤ) ^ 2 * m2 := by ring
      rw [e]
      exact prime_mul_unit_passes_mod_eight hp h8 hn hd1 k
    · rcases h8 with h8 | h8
      · obtain ⟨u1, hu1⟩ := prime_mul_unit_passes_one_mod_eight hp h8 hn hm1 hd1 k
        exact ⟨(p : ℤ) * u1, 1, Or.inr h1, hu1⟩
      · exact ⟨(p : ℤ) * 0, 1, Or.inr h1,
          prime_mul_unit_passes_seven_mod_eight hp h8 hn hm1 hd1 k⟩

/-! ### The symmetries -/

theorem GSqM.negfst {M d : ℤ} {p : ℤ × ℤ} (h : GSqM M d p) : GSqM M d (-p.1, p.2) := by
  obtain ⟨a, b, e, f⟩ := h
  refine ⟨b, a, ?_, ?_⟩
  · show d * -p.1 ≡ b ^ 2 - a ^ 2 [ZMOD M]
    calc d * -p.1 = -(d * p.1) := by ring
      _ ≡ -(a ^ 2 - b ^ 2) [ZMOD M] := e.neg
      _ = b ^ 2 - a ^ 2 := by ring
  · show d * p.2 ≡ 2 * b * a [ZMOD M]
    calc d * p.2 ≡ 2 * a * b [ZMOD M] := f
      _ = 2 * b * a := by ring

theorem GSqM.negboth {M d : ℤ} {p : ℤ × ℤ} (h : GSqM M d p) : GSqM M d (-p.1, -p.2) := by
  obtain ⟨a, b, e, f⟩ := h
  refine ⟨b, -a, ?_, ?_⟩
  · show d * -p.1 ≡ b ^ 2 - (-a) ^ 2 [ZMOD M]
    calc d * -p.1 = -(d * p.1) := by ring
      _ ≡ -(a ^ 2 - b ^ 2) [ZMOD M] := e.neg
      _ = b ^ 2 - (-a) ^ 2 := by ring
  · show d * -p.2 ≡ 2 * b * -a [ZMOD M]
    calc d * -p.2 = -(d * p.2) := by ring
      _ ≡ -(2 * a * b) [ZMOD M] := f.neg
      _ = 2 * b * -a := by ring

theorem GSqM.scale16 {M d : ℤ} {p : ℤ × ℤ} (h : GSqM M d p) : GSqM M d (16 * p.1, 16 * p.2) := by
  obtain ⟨a, b, e, f⟩ := h
  refine ⟨4 * a, 4 * b, ?_, ?_⟩
  · show d * (16 * p.1) ≡ (4 * a) ^ 2 - (4 * b) ^ 2 [ZMOD M]
    calc d * (16 * p.1) = 16 * (d * p.1) := by ring
      _ ≡ 16 * (a ^ 2 - b ^ 2) [ZMOD M] := e.mul_left _
      _ = (4 * a) ^ 2 - (4 * b) ^ 2 := by ring
  · show d * (16 * p.2) ≡ 2 * (4 * a) * (4 * b) [ZMOD M]
    calc d * (16 * p.2) = 16 * (d * p.2) := by ring
      _ ≡ 16 * (2 * a * b) [ZMOD M] := f.mul_left _
      _ = 2 * (4 * a) * (4 * b) := by ring

theorem GSqM.scale16neg {M d : ℤ} {p : ℤ × ℤ} (h : GSqM M d p) :
    GSqM M d (16 * p.1, -(16 * p.2)) := by
  obtain ⟨a, b, e, f⟩ := h
  refine ⟨4 * a, -(4 * b), ?_, ?_⟩
  · show d * (16 * p.1) ≡ (4 * a) ^ 2 - (-(4 * b)) ^ 2 [ZMOD M]
    calc d * (16 * p.1) = 16 * (d * p.1) := by ring
      _ ≡ 16 * (a ^ 2 - b ^ 2) [ZMOD M] := e.mul_left _
      _ = (4 * a) ^ 2 - (-(4 * b)) ^ 2 := by ring
  · show d * -(16 * p.2) ≡ 2 * (4 * a) * -(4 * b) [ZMOD M]
    calc d * -(16 * p.2) = -16 * (d * p.2) := by ring
      _ ≡ -16 * (2 * a * b) [ZMOD M] := f.mul_left _
      _ = 2 * (4 * a) * -(4 * b) := by ring

/-- Exchanging `m` and `n` permutes the eight products up to sign and conjugation. -/
theorem gsqM_swap {M d m n u v : ℤ} (h : ∀ f : Bool × Bool × Bool, GSqM M d (Q m n u v f)) :
    ∀ f : Bool × Bool × Bool, GSqM M d (Q n m u v f) := by
  have key1 : ∀ b c : Bool, Q n m u v (false, b, c) =
      (-(Q m n u v (false, !b, !c)).1, -(Q m n u v (false, !b, !c)).2) := by
    intro b c
    cases b <;> cases c <;>
      simp only [Q, gm, cj, cjIf, P0, P8, PL, PU, Bool.not_false, Bool.not_true, if_true,
        if_false, Bool.false_eq_true, Prod.mk.injEq] <;>
      constructor <;> ring
  have key2 : ∀ b c : Bool, Q n m u v (true, b, c) =
      (-(Q m n u v (true, b, c)).1, (Q m n u v (true, b, c)).2) := by
    intro b c
    cases b <;> cases c <;>
      simp only [Q, gm, cj, cjIf, P0, P8, PL, PU, if_true, if_false, Bool.false_eq_true,
        Prod.mk.injEq] <;>
      constructor <;> ring
  rintro ⟨a, b, c⟩
  cases a
  · rw [key1]
    exact (h _).negboth
  · rw [key2]
    exact (h _).negfst

/-- The substitution `(m, n, u, v) ↦ (m + n, m − n, u + v, u − v)` permutes the eight products up
to the factor 16 and conjugation. -/
theorem gsqM_rot {M d m n u v : ℤ} (h : ∀ f : Bool × Bool × Bool, GSqM M d (Q m n u v f)) :
    ∀ f : Bool × Bool × Bool, GSqM M d (Q (m + n) (m - n) (u + v) (u - v) f) := by
  have key1 : ∀ a c : Bool, Q (m + n) (m - n) (u + v) (u - v) (a, false, c) =
      (16 * (Q m n u v (c, false, a)).1, 16 * (Q m n u v (c, false, a)).2) := by
    intro a c
    cases a <;> cases c <;>
      simp only [Q, gm, cj, cjIf, P0, P8, PL, PU, if_true, if_false, Bool.false_eq_true,
        Prod.mk.injEq] <;>
      constructor <;> ring
  have key2 : ∀ a c : Bool, Q (m + n) (m - n) (u + v) (u - v) (a, true, c) =
      (16 * (Q m n u v (!c, true, !a)).1, -(16 * (Q m n u v (!c, true, !a)).2)) := by
    intro a c
    cases a <;> cases c <;>
      simp only [Q, gm, cj, cjIf, P0, P8, PL, PU, Bool.not_false, Bool.not_true, if_true,
        if_false, Bool.false_eq_true, Prod.mk.injEq] <;>
      constructor <;> ring
  rintro ⟨a, b, c⟩
  cases b
  · rw [key1]
    exact (h _).scale16
  · rw [key2]
    exact (h _).scale16neg

theorem ClassAt.swap {M q d m n u v : ℤ} (h : ClassAt M q d m n u v) :
    ClassAt M q d n m u v := by
  obtain ⟨w0, w8, wL, wU, ⟨h0, h8, hL, hU⟩, ⟨n0, n8, nL, nU⟩, g⟩ := h
  have eL : nrm (PL n m u v) = nrm (PL m n u v) := by simp only [nrm, PL]; ring
  have eU : nrm (PU n m u v) = nrm (PU m n u v) := by simp only [nrm, PU]; ring
  refine ⟨w8, w0, wL, wU, ⟨h8, h0, ?_, ?_⟩, ⟨n8, n0, nL, nU⟩, gsqM_swap g⟩
  · rw [eL]
    exact hL
  · rw [eU]
    exact hU

/-- The substitution `(m, n, u, v) ↦ (m + n, m − n, u + v, u − v)` sends the four norms to
`N_L, N_U, 16·N_0, 16·N_8`. -/
theorem ClassAt.rot {M q d m n u v : ℤ} (hq : IsCoprime (4 : ℤ) q)
    (h : ClassAt M q d m n u v) : ClassAt M q d (m + n) (m - n) (u + v) (u - v) := by
  obtain ⟨w0, w8, wL, wU, ⟨h0, h8, hL, hU⟩, ⟨n0, n8, nL, nU⟩, g⟩ := h
  have eL : nrm (PL (m + n) (m - n) (u + v) (u - v)) = 16 * nrm (P0 m n u v) := by
    simp only [nrm, PL, P0]; ring
  have eU : nrm (PU (m + n) (m - n) (u + v) (u - v)) = 16 * nrm (P8 m n u v) := by
    simp only [nrm, PU, P8]; ring
  refine ⟨wL, wU, 4 * w0, 4 * w8, ⟨hL, hU, ?_, ?_⟩,
    ⟨nL, nU, fun hd => n0 (hq.symm.dvd_of_dvd_mul_left hd),
      fun hd => n8 (hq.symm.dvd_of_dvd_mul_left hd)⟩, gsqM_rot g⟩
  · rw [eL]
    calc 16 * nrm (P0 m n u v) ≡ 16 * w0 ^ 2 [ZMOD M] := h0.mul_left 16
      _ = (4 * w0) ^ 2 := by ring
  · rw [eU]
    calc 16 * nrm (P8 m n u v) ≡ 16 * w8 ^ 2 [ZMOD M] := h8.mul_left 16
      _ = (4 * w8) ^ 2 := by ring

/-- From the ratio `2m : 2n` back to `m : n`, when 4 has an inverse `σ` modulo `M` and `2σ` is
prime to `q`: the norms are divided by 4 and the products by 16. -/
theorem ClassAt.unscale {M q d m n u v σ : ℤ} (hσ : σ * 4 ≡ 1 [ZMOD M])
    (hq : IsCoprime (2 * σ) q) (h : ClassAt M q d (2 * m) (2 * n) u v) :
    ClassAt M q d m n u v := by
  obtain ⟨w0, w8, wL, wU, ⟨h0, h8, hL, hU⟩, ⟨n0, n8, nL, nU⟩, g⟩ := h
  have key : ∀ (N N' w : ℤ), N' = 4 * N → N' ≡ w ^ 2 [ZMOD M] →
      N ≡ (2 * σ * w) ^ 2 [ZMOD M] := by
    intro N N' w e hw
    rw [e] at hw
    calc N = N * 1 ^ 2 := by ring
      _ ≡ N * (σ * 4) ^ 2 [ZMOD M] := ((hσ.pow 2).symm).mul_left _
      _ = 4 * σ ^ 2 * (4 * N) := by ring
      _ ≡ 4 * σ ^ 2 * w ^ 2 [ZMOD M] := hw.mul_left _
      _ = (2 * σ * w) ^ 2 := by ring
  have nd : ∀ w : ℤ, ¬q ∣ w → ¬q ∣ 2 * σ * w :=
    fun w hw hd => hw (hq.symm.dvd_of_dvd_mul_left hd)
  refine ⟨2 * σ * w0, 2 * σ * w8, 2 * σ * wL, 2 * σ * wU,
    ⟨key _ _ w0 (by simp only [nrm, P0]; ring) h0, key _ _ w8 (by simp only [nrm, P8]; ring) h8,
      key _ _ wL (by simp only [nrm, PL]; ring) hL, key _ _ wU (by simp only [nrm, PU]; ring) hU⟩,
    ⟨nd w0 n0, nd w8 n8, nd wL nL, nd wU nU⟩, fun f => ?_⟩
  have h3 := g f
  have e : Q (2 * m) (2 * n) u v f = (16 * (Q m n u v f).1, 16 * (Q m n u v f).2) := by
    have h4 := Q_scale 2 1 m n u v f
    rw [one_mul, one_mul] at h4
    rw [h4]
    norm_num
  rw [e] at h3
  obtain ⟨a, b, ea, eb⟩ := h3
  refine ⟨σ * a, σ * b, ?_, ?_⟩
  · calc d * (Q m n u v f).1 = d * (Q m n u v f).1 * 1 ^ 2 := by ring
      _ ≡ d * (Q m n u v f).1 * (σ * 4) ^ 2 [ZMOD M] := ((hσ.pow 2).symm).mul_left _
      _ = σ ^ 2 * (d * (16 * (Q m n u v f).1)) := by ring
      _ ≡ σ ^ 2 * (a ^ 2 - b ^ 2) [ZMOD M] := ea.mul_left _
      _ = (σ * a) ^ 2 - (σ * b) ^ 2 := by ring
  · calc d * (Q m n u v f).2 = d * (Q m n u v f).2 * 1 ^ 2 := by ring
      _ ≡ d * (Q m n u v f).2 * (σ * 4) ^ 2 [ZMOD M] := ((hσ.pow 2).symm).mul_left _
      _ = σ ^ 2 * (d * (16 * (Q m n u v f).2)) := by ring
      _ ≡ σ ^ 2 * (2 * a * b) [ZMOD M] := eb.mul_left _
      _ = 2 * (σ * a) * (σ * b) := by ring

theorem Passes.swap {p k : ℕ} {d m n : ℤ} (h : Passes p k d m n) : Passes p k d n m := by
  obtain ⟨u, v, hprim, g⟩ := h
  exact ⟨u, v, hprim, g.swap⟩

/-- From the ratio `(m + n) : (m − n)` back to `m : n`, for an odd prime: the substitution applied
twice gives `(2m, 2n, 2u, 2v)`; here it is applied once more and the factor 2 of `m, n` removed. -/
theorem Passes.rot {p k : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {d m n : ℤ}
    (h : Passes p k d (m + n) (m - n)) : Passes p k d m n := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  obtain ⟨u, v, hprim, g⟩ := h
  have h16 : ¬(p : ℤ) ∣ 4 := by
    have e4 : (4 : ℤ) = 2 * 2 := by norm_num
    rw [e4]
    exact not_dvd_mul_of_prime hp' h2 h2
  have hq4 : IsCoprime (4 : ℤ) ((p : ℤ) ^ 3) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr h16).symm.pow_right
  have g2 := g.rot hq4
  have e1 : m + n + (m - n) = 2 * m := by ring
  have e2 : m + n - (m - n) = 2 * n := by ring
  rw [e1, e2] at g2
  obtain ⟨σ, μ, hσ⟩ : IsCoprime (4 : ℤ) ((p : ℤ) ^ (k + 1)) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr h16).symm.pow_right
  have hσ1 : σ * 4 ≡ 1 [ZMOD (p : ℤ) ^ (k + 1)] :=
    Int.modEq_iff_dvd.mpr ⟨μ, by linear_combination -hσ⟩
  have hnd : ¬(p : ℤ) ∣ 2 * σ := by
    intro hd
    have e : σ * 4 = 2 * (2 * σ) := by ring
    have h1 : (p : ℤ) ∣ σ * 4 + μ * (p : ℤ) ^ (k + 1) := by
      rw [e]
      exact dvd_add (Dvd.dvd.mul_left hd 2)
        (Dvd.dvd.mul_left (dvd_pow_self _ (Nat.succ_ne_zero k)) μ)
    rw [hσ] at h1
    exact hp'.not_isUnit (isUnit_of_dvd_one h1)
  have hq2 : IsCoprime (2 * σ) ((p : ℤ) ^ 3) :=
    (hp'.irreducible.coprime_iff_not_dvd.mpr hnd).symm.pow_right
  refine ⟨u + v, u - v, ?_, g2.unscale hσ1 hq2⟩
  by_contra hc
  rw [not_or, not_not, not_not] at hc
  obtain ⟨c1, c2⟩ := hc
  have hu : (p : ℤ) ∣ 2 * u := by
    have := dvd_add c1 c2
    have e : u + v + (u - v) = 2 * u := by ring
    rwa [e] at this
  have hv : (p : ℤ) ∣ 2 * v := by
    have := dvd_sub c1 c2
    have e : u + v - (u - v) = 2 * v := by ring
    rwa [e] at this
  have hu' : (p : ℤ) ∣ u := (hp'.dvd_or_dvd hu).resolve_left h2
  have hv' : (p : ℤ) ∣ v := (hp'.dvd_or_dvd hv).resolve_left h2
  rcases hprim with h | h
  · exact h hu'
  · exact h hv'

/-- The table for `p` dividing `n`. -/
theorem passes_of_row_n {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n d : ℤ} (hn : (p : ℤ) ∣ n)
    (hm : ¬(p : ℤ) ∣ m) (hrow : Row p n d) (k : ℕ) : Passes p k d m n :=
  (passes_of_row hp hp2 hn hm hrow k).swap

/-- The table for `p` dividing `m + n`. -/
theorem passes_of_row_add {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n d : ℤ}
    (h1 : (p : ℤ) ∣ m + n) (h2 : ¬(p : ℤ) ∣ m - n) (hrow : Row p (m + n) d) (k : ℕ) :
    Passes p k d m n :=
  Passes.rot hp hp2 (passes_of_row hp hp2 h1 h2 hrow k)

/-- The table for `p` dividing `m − n`. -/
theorem passes_of_row_sub {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n d : ℤ}
    (h1 : (p : ℤ) ∣ m - n) (h2 : ¬(p : ℤ) ∣ m + n) (hrow : Row p (m - n) d) (k : ℕ) :
    Passes p k d m n :=
  Passes.rot hp hp2 (passes_of_row hp hp2 h1 h2 hrow k).swap

/-- **The table of Theorem C as lower bounds, in one statement**: for coprime `m`, `n` and an odd
prime `p` dividing one of `m`, `n`, `m + n`, `m − n`, every class of the row of that factor passes
modulo every power of `p`. -/
theorem passes_of_table {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {m n d : ℤ} (hc : IsCoprime m n)
    (h : ((p : ℤ) ∣ m ∧ Row p m d) ∨ ((p : ℤ) ∣ n ∧ Row p n d) ∨
      ((p : ℤ) ∣ m + n ∧ Row p (m + n) d) ∨ ((p : ℤ) ∣ m - n ∧ Row p (m - n) d)) (k : ℕ) :
    Passes p k d m n := by
  have hp' : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have h2 := not_dvd_two_of_odd hp hp2
  have hboth : ¬((p : ℤ) ∣ m ∧ (p : ℤ) ∣ n) := by
    rintro ⟨a, b⟩
    obtain ⟨x, y, e⟩ := hc
    have h1 : (p : ℤ) ∣ x * m + y * n := dvd_add (Dvd.dvd.mul_left a x) (Dvd.dvd.mul_left b y)
    rw [e] at h1
    exact hp'.not_isUnit (isUnit_of_dvd_one h1)
  have hpm : ¬((p : ℤ) ∣ m + n ∧ (p : ℤ) ∣ m - n) := by
    rintro ⟨a, b⟩
    have hm : (p : ℤ) ∣ 2 * m := by
      have := dvd_add a b
      have e : m + n + (m - n) = 2 * m := by ring
      rwa [e] at this
    have hn : (p : ℤ) ∣ 2 * n := by
      have := dvd_sub a b
      have e : m + n - (m - n) = 2 * n := by ring
      rwa [e] at this
    exact hboth ⟨(hp'.dvd_or_dvd hm).resolve_left h2, (hp'.dvd_or_dvd hn).resolve_left h2⟩
  rcases h with ⟨a, r⟩ | ⟨a, r⟩ | ⟨a, r⟩ | ⟨a, r⟩
  · exact passes_of_row hp hp2 a (fun b => hboth ⟨a, b⟩) r k
  · exact passes_of_row_n hp hp2 a (fun b => hboth ⟨b, a⟩) r k
  · exact passes_of_row_add hp hp2 a (fun b => hpm ⟨a, b⟩) r k
  · exact passes_of_row_sub hp hp2 a (fun b => hpm ⟨b, a⟩) r k

/-! ### The predicate is not the weak one -/

/-- **`d = 3` does not pass at `p = 3` for the ratio `3 : 1`**, already modulo 9. With the eight
products alone it would: at `(0, 1)` each product is `3` times a Gaussian unit of square norm, and
3 times that is a Gaussian square modulo every power of 3 (`class_gaps.py`, part 1). The class 3
is not in the table at a prime `≡ 3 (mod 8)`, so this is what the table asks. -/
theorem not_passes_three : ¬Passes 3 1 3 3 1 := by
  rintro ⟨u, v, hprim, w0, w8, wL, wU, ⟨-, -, hL, -⟩, -, g⟩
  have e9 : (((3 : ℕ) : ℤ)) ^ (1 + 1) = 9 := by norm_num
  rw [e9] at hL g
  have hprim' : ¬(3 : ℤ) ∣ u ∨ ¬(3 : ℤ) ∣ v := by exact_mod_cast hprim
  have h39 : (3 : ℤ) ∣ 9 := by norm_num
  have h30 : (3 : ZMod 3) = 0 := by decide
  -- `u` and `v` are units modulo 3, because `N_L` is a square
  have hA : ∀ x y w : ZMod 3, (x ≠ 0 ∨ y ≠ 0) → (4 * (x - y)) ^ 2 + (2 * (x + y)) ^ 2 = w ^ 2 →
      x ≠ 0 ∧ y ≠ 0 := by decide
  have hL3 : nrm (PL 3 1 u v) ≡ wL ^ 2 [ZMOD ((3 : ℕ) : ℤ)] := by exact_mod_cast hL.of_dvd h39
  have hLz := (ZMod.intCast_eq_intCast_iff _ _ 3).mpr hL3
  simp only [nrm, PL] at hLz
  push_cast at hLz
  have hx : ((u : ℤ) : ZMod 3) ≠ 0 ∨ ((v : ℤ) : ZMod 3) ≠ 0 := by
    rcases hprim' with h | h
    · exact Or.inl (fun hz => h (by exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd u 3).mp hz))
    · exact Or.inr (fun hz => h (by exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd v 3).mp hz))
  obtain ⟨hu, hv⟩ := hA _ _ ((wL : ℤ) : ZMod 3) hx (by linear_combination hLz)
  -- `3·Q` a Gaussian square modulo 9 forces `3 | Q`
  obtain ⟨a, b, ea, eb⟩ := g (false, false, false)
  have hB : ∀ x y : ZMod 3, x ^ 2 - y ^ 2 = 0 → 2 * x * y = 0 → x = 0 ∧ y = 0 := by decide
  have ea3 : 3 * (Q 3 1 u v (false, false, false)).1 ≡ a ^ 2 - b ^ 2 [ZMOD ((3 : ℕ) : ℤ)] := by
    exact_mod_cast ea.of_dvd h39
  have eb3 : 3 * (Q 3 1 u v (false, false, false)).2 ≡ 2 * a * b [ZMOD ((3 : ℕ) : ℤ)] := by
    exact_mod_cast eb.of_dvd h39
  have eaz := (ZMod.intCast_eq_intCast_iff _ _ 3).mpr ea3
  have ebz := (ZMod.intCast_eq_intCast_iff _ _ 3).mpr eb3
  push_cast at eaz ebz
  rw [h30, zero_mul] at eaz ebz
  obtain ⟨ha, hb⟩ := hB _ _ eaz.symm ebz.symm
  have ha3 : (3 : ℤ) ∣ a := by exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd a 3).mp ha
  have hb3 : (3 : ℤ) ∣ b := by exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd b 3).mp hb
  obtain ⟨a', rfl⟩ := ha3
  obtain ⟨b', rfl⟩ := hb3
  have e1 : (3 * a') ^ 2 - (3 * b') ^ 2 = 9 * (a' ^ 2 - b' ^ 2) := by ring
  have e2 : 2 * (3 * a') * (3 * b') = 9 * (2 * a' * b') := by ring
  rw [e1] at ea
  rw [e2] at eb
  have d1 : (9 : ℤ) ∣ 3 * (Q 3 1 u v (false, false, false)).1 :=
    Int.modEq_zero_iff_dvd.mp (ea.trans (Int.modEq_zero_iff_dvd.mpr (dvd_mul_right 9 _)))
  have d2 : (9 : ℤ) ∣ 3 * (Q 3 1 u v (false, false, false)).2 :=
    Int.modEq_zero_iff_dvd.mp (eb.trans (Int.modEq_zero_iff_dvd.mpr (dvd_mul_right 9 _)))
  have t1 : (3 : ℤ) ∣ (Q 3 1 u v (false, false, false)).1 := by omega
  have t2 : (3 : ℤ) ∣ (Q 3 1 u v (false, false, false)).2 := by omega
  -- so 3 divides the product of the four norms, which are units
  have hN : (3 : ℤ) ∣ nrm (Q 3 1 u v (false, false, false)) :=
    dvd_add (Dvd.dvd.pow t1 two_ne_zero) (Dvd.dvd.pow t2 two_ne_zero)
  rw [nrm_Q] at hN
  have hNz := (ZMod.intCast_zmod_eq_zero_iff_dvd _ 3).mpr hN
  simp only [nrm, P0, P8, PL, PU] at hNz
  push_cast at hNz
  have hC : ∀ x y : ZMod 3, x ≠ 0 → y ≠ 0 →
      ((3 * y) ^ 2 + x ^ 2) * (y ^ 2 + (3 * x) ^ 2) * ((4 * (x - y)) ^ 2 + (2 * (x + y)) ^ 2) *
        ((2 * (x - y)) ^ 2 + (4 * (x + y)) ^ 2) ≠ 0 := by decide
  exact hC _ _ hu hv (by linear_combination hNz)

end MagicSquareOfSquares.CoverDescent
