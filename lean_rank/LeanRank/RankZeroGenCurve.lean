import LeanRank.RankZeroGen

/-!
# `rank E(ℚ) = 0` from a certificate, for `E : y² = x(x − e₂)(x − e₃)`

Written by an AI model; not read by a person.  Experiment, not counted.

`finrank_eq_zero_of`: if a certificate `c` passes `Cert.check` (and its tables of squares are
right), the primes `S₁`, `S₂` cover `e₂e₃` and `e₂(e₂ − e₃)`, `c.good` has at most four classes
including `(1, 1)`, and the three points of order 2 have their classes in `c.good` (explicit
witnesses), then `Module.finrank ℤ E(ℚ) = 0`.  The proof is that of `LeanRank.RankZero`: Stoll's
`card_range_μ` (`#(im μ) = 2 ^ rank · #E[2]`), `card_ker_nsmul_two` (`#E[2] = 4`), the injective
map `Ψ` to `(ℚˣ/(ℚˣ)²)³` and `descent_gen` with `Cert.good_of_check`.
-/

open WeierstrassCurve WeierstrassCurve.Affine Polynomial
open LeanRank.RankZero (int_sq_of_rat_sq cls cls_eq cls_one)

namespace LeanRank.RankZeroGen

variable (e₂ e₃ : ℤ)

/-- `y² = x(x − e₂)(x − e₃)`. -/
def curve : Affine ℚ := ⟨0, -((e₂ : ℚ) + e₃), 0, (e₂ : ℚ) * e₃, 0⟩

instance : (curve e₂ e₃).IsCharNeTwoNF := ⟨rfl, rfl⟩

variable {e₂ e₃}

theorem curve_Δ : (curve e₂ e₃).Δ = 16 * ((e₂ : ℚ) * e₃ * (e₂ - e₃)) ^ 2 := by
  simp only [curve, WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆, WeierstrassCurve.b₈]
  ring

theorem isElliptic (h2 : e₂ ≠ 0) (h3 : e₃ ≠ 0) (h23 : e₂ ≠ e₃) : (curve e₂ e₃).IsElliptic := by
  rw [isElliptic_iff, isUnit_iff_ne_zero, curve_Δ]
  have a2 : (e₂ : ℚ) ≠ 0 := by exact_mod_cast h2
  have a3 : (e₃ : ℚ) ≠ 0 := by exact_mod_cast h3
  have a23 : (e₂ : ℚ) - e₃ ≠ 0 := sub_ne_zero.mpr (by exact_mod_cast h23)
  exact mul_ne_zero (by norm_num) (pow_ne_zero 2 (mul_ne_zero (mul_ne_zero a2 a3) a23))

theorem f_eq : (curve e₂ e₃).f = X * (X - C (e₂ : ℚ)) * (X - C (e₃ : ℚ)) := by
  simp only [Affine.f, curve, map_neg, map_add, map_mul, map_zero]
  ring

theorem eval_f (x : ℚ) : (curve e₂ e₃).f.eval x = x * (x - e₂) * (x - e₃) := by
  rw [f_eq]; simp

theorem equation_iff' (x y : ℚ) :
    (curve e₂ e₃).Equation x y ↔ y ^ 2 = x * (x - e₂) * (x - e₃) := by
  rw [equation_iff_eval_f_eq_sq, eval_f, eq_comm]

theorem r0 : (curve e₂ e₃).f.eval 0 = 0 := by rw [eval_f]; ring
theorem r2 : (curve e₂ e₃).f.eval (e₂ : ℚ) = 0 := by rw [eval_f]; ring
theorem r3 : (curve e₂ e₃).f.eval (e₃ : ℚ) = 0 := by rw [eval_f]; ring

/-- Evaluation at a root of `f`. -/
noncomputable def ev (r : ℚ) (hr : (curve e₂ e₃).f.eval r = 0) : (curve e₂ e₃).A →+* ℚ :=
  AdjoinRoot.lift (RingHom.id ℚ) r hr

theorem ev_mk (r : ℚ) (hr : (curve e₂ e₃).f.eval r = 0) (p : ℚ[X]) :
    ev r hr (AdjoinRoot.mk (curve e₂ e₃).f p) = p.eval r :=
  AdjoinRoot.lift_mk hr p

/-- Evaluation at the three roots, on square classes. -/
noncomputable def Ψ : (curve e₂ e₃).M →* Units.modPow ℚ 2 × Units.modPow ℚ 2 × Units.modPow ℚ 2 :=
  (Units.modPow.map (ev 0 r0).toMonoidHom 2).prod
    ((Units.modPow.map (ev (e₂ : ℚ) r2).toMonoidHom 2).prod
      (Units.modPow.map (ev (e₃ : ℚ) r3).toMonoidHom 2))

theorem dvd_of_roots (h2 : (e₂ : ℚ) ≠ 0) (h3 : (e₃ : ℚ) ≠ 0) (h23 : (e₂ : ℚ) ≠ e₃) {p : ℚ[X]}
    (h0 : p.eval 0 = 0) (h2' : p.eval (e₂ : ℚ) = 0) (h3' : p.eval (e₃ : ℚ) = 0) :
    (curve e₂ e₃).f ∣ p := by
  rw [f_eq]
  obtain ⟨p₁, rfl⟩ : X - C 0 ∣ p := dvd_iff_isRoot.mpr h0
  obtain ⟨p₂, rfl⟩ : X - C (e₂ : ℚ) ∣ p₁ := by
    apply dvd_iff_isRoot.mpr
    simp only [eval_mul, eval_sub, eval_X, eval_C, sub_zero] at h2'
    exact (mul_eq_zero.mp h2').resolve_left h2
  obtain ⟨p₃, rfl⟩ : X - C (e₃ : ℚ) ∣ p₂ := by
    apply dvd_iff_isRoot.mpr
    simp only [eval_mul, eval_sub, eval_X, eval_C, sub_zero] at h3'
    rcases mul_eq_zero.mp h3' with h | h
    · exact absurd h h3
    · exact (mul_eq_zero.mp h).resolve_left (sub_ne_zero.mpr h23.symm)
  exact ⟨p₃, by simp only [map_zero, sub_zero]; ring⟩

theorem Ψ_injective (h2 : e₂ ≠ 0) (h3 : e₃ ≠ 0) (h23 : e₂ ≠ e₃) :
    Function.Injective (Ψ (e₂ := e₂) (e₃ := e₃)) := by
  have a2 : (e₂ : ℚ) ≠ 0 := by exact_mod_cast h2
  have a3 : (e₃ : ℚ) ≠ 0 := by exact_mod_cast h3
  have a23 : (e₂ : ℚ) - e₃ ≠ 0 := sub_ne_zero.mpr (by exact_mod_cast h23)
  have a32 : (e₃ : ℚ) - e₂ ≠ 0 := sub_ne_zero.mpr (by exact_mod_cast h23.symm)
  refine (injective_iff_map_eq_one _).mpr fun m hm => ?_
  induction m using QuotientGroup.induction_on with
  | H u =>
    simp only [Ψ, MonoidHom.prod_apply, Units.modPow.map_mk, Prod.mk_eq_one] at hm
    obtain ⟨h0, h12, h16⟩ := hm
    rw [Units.modPow.mk_eq_one_iff_isSquare] at h0 h12 h16 ⊢
    simp only [Units.coe_map, RingHom.toMonoidHom_eq_coe, MonoidHom.coe_ofClass] at h0 h12 h16
    obtain ⟨g, hg⟩ := AdjoinRoot.mk_surjective (u : (curve e₂ e₃).A)
    rw [← hg, ev_mk] at h0 h12 h16
    obtain ⟨s₀, hs₀⟩ := h0
    obtain ⟨s₁, hs₁⟩ := h12
    obtain ⟨s₂, hs₂⟩ := h16
    set L : ℚ[X] := C (s₀ / (e₂ * e₃)) * ((X - C (e₂ : ℚ)) * (X - C (e₃ : ℚ))) +
      C (s₁ / (e₂ * (e₂ - e₃))) * (X * (X - C (e₃ : ℚ))) +
      C (s₂ / (e₃ * (e₃ - e₂))) * (X * (X - C (e₂ : ℚ))) with hL
    have hL0 : L.eval 0 = s₀ := by
      simp only [hL, eval_add, eval_mul, eval_C, eval_X, eval_sub]
      field_simp
      try ring
    have hL12 : L.eval (e₂ : ℚ) = s₁ := by
      simp only [hL, eval_add, eval_mul, eval_C, eval_X, eval_sub]
      field_simp
      try ring
    have hL16 : L.eval (e₃ : ℚ) = s₂ := by
      simp only [hL, eval_add, eval_mul, eval_C, eval_X, eval_sub]
      field_simp
      try ring
    refine ⟨AdjoinRoot.mk (curve e₂ e₃).f L, ?_⟩
    rw [← hg, ← map_mul, AdjoinRoot.mk_eq_mk]
    apply dvd_of_roots a2 a3 (by exact_mod_cast h23) <;>
      simp [eval_mul, hL0, hL12, hL16, hs₀, hs₁, hs₂]

/-- The triple of values `(v₁, v₂, v₃)` has the classes `(d₁, d₂, d₁d₂)` of a good class. -/
def RealizesTriple (good : List (ℤ × ℤ)) (v₁ v₂ v₃ : ℚ) : Prop :=
  ∃ d ∈ good, ∃ w₁ w₂ w₃ : ℚ, w₁ ≠ 0 ∧ w₂ ≠ 0 ∧ w₃ ≠ 0 ∧ v₁ = (d.1 : ℚ) * w₁ ^ 2 ∧
    v₂ = (d.2 : ℚ) * w₂ ^ 2 ∧ v₃ = ((d.1 * d.2 : ℤ) : ℚ) * w₃ ^ 2

open Classical in
/-- The target set. -/
noncomputable def Tset (good : List (ℤ × ℤ)) :
    Finset (Units.modPow ℚ 2 × Units.modPow ℚ 2 × Units.modPow ℚ 2) :=
  (good.map fun d => (cls (d.1 : ℚ), cls (d.2 : ℚ), cls ((d.1 * d.2 : ℤ) : ℚ))).toFinset

open Classical in
theorem card_Tset (good : List (ℤ × ℤ)) : (Tset good).card ≤ good.length := by
  unfold Tset
  exact (List.toFinset_card_le _).trans (List.length_map _).le

theorem cls_of {u : ℚˣ} {d w : ℚ} (hw : w ≠ 0) (h : (u : ℚ) = d * w ^ 2) :
    (QuotientGroup.mk u : Units.modPow ℚ 2) = cls d := by
  have hd : d ≠ 0 := by
    rintro rfl; rw [zero_mul] at h; exact u.ne_zero h
  exact cls_eq hd hw h

open Classical in
theorem mem_Tset_of {good : List (ℤ × ℤ)} {a : (curve e₂ e₃).A} (ha : IsUnit a) {v₁ v₂ v₃ : ℚ}
    (h₁ : ev 0 r0 a = v₁) (h₂ : ev (e₂ : ℚ) r2 a = v₂) (h₃ : ev (e₃ : ℚ) r3 a = v₃)
    (hR : RealizesTriple good v₁ v₂ v₃) : Ψ (ha.unit : (curve e₂ e₃).M) ∈ Tset good := by
  obtain ⟨d, hd, w₁, w₂, w₃, hw₁, hw₂, hw₃, e1, e2, e3⟩ := hR
  have hΨ : Ψ (ha.unit : (curve e₂ e₃).M) =
      (QuotientGroup.mk (Units.map (ev 0 r0).toMonoidHom ha.unit),
        QuotientGroup.mk (Units.map (ev (e₂ : ℚ) r2).toMonoidHom ha.unit),
        QuotientGroup.mk (Units.map (ev (e₃ : ℚ) r3).toMonoidHom ha.unit)) := rfl
  have c1 : (QuotientGroup.mk (Units.map (ev 0 r0).toMonoidHom ha.unit) :
      Units.modPow ℚ 2) = cls (d.1 : ℚ) := cls_of hw₁ (by simpa [h₁] using e1)
  have c2 : (QuotientGroup.mk (Units.map (ev (e₂ : ℚ) r2).toMonoidHom ha.unit) :
      Units.modPow ℚ 2) = cls (d.2 : ℚ) := cls_of hw₂ (by simpa [h₂] using e2)
  have c3 : (QuotientGroup.mk (Units.map (ev (e₃ : ℚ) r3).toMonoidHom ha.unit) :
      Units.modPow ℚ 2) = cls ((d.1 * d.2 : ℤ) : ℚ) := cls_of hw₃ (by simpa [h₃] using e3)
  rw [hΨ, c1, c2, c3]
  unfold Tset
  rw [List.mem_toFinset]
  exact List.mem_map.mpr ⟨d, hd, rfl⟩

/-- **`rank E(ℚ) = 0` from a certificate.** -/
theorem finrank_eq_zero_of (c : Cert) (h2 : c.e₂ ≠ 0) (h3 : c.e₃ ≠ 0) (h23 : c.e₂ ≠ c.e₃)
    (hS₁ : ∀ p ∈ c.S₁, p.Prime) (hS₂ : ∀ p ∈ c.S₂, p.Prime)
    (hd₁ : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ c.e₂ * c.e₃ → p ∈ c.S₁)
    (hd₂ : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ c.e₂ * (c.e₂ - c.e₃) → p ∈ c.S₂)
    (hc : c.check = true) (hT : c.TablesOK) (hlen : c.good.length ≤ 4)
    (h11 : ((1 : ℤ), (1 : ℤ)) ∈ c.good)
    (t0 : RealizesTriple c.good ((c.e₂ : ℚ) * c.e₃) (-(c.e₂ : ℚ)) (-(c.e₃ : ℚ)))
    (t2 : RealizesTriple c.good (c.e₂ : ℚ) ((c.e₂ : ℚ) * (c.e₂ - c.e₃)) ((c.e₂ : ℚ) - c.e₃))
    (t3 : RealizesTriple c.good (c.e₃ : ℚ) ((c.e₃ : ℚ) - c.e₂) ((c.e₃ : ℚ) * (c.e₃ - c.e₂))) :
    Module.finrank ℤ (curve c.e₂ c.e₃).Point = 0 := by
  classical
  haveI := isElliptic h2 h3 h23
  haveI : IsAddFG (curve c.e₂ c.e₃).Point := fg_point_of_numberField
  have a2 : (c.e₂ : ℚ) ≠ 0 := by exact_mod_cast h2
  have a3 : (c.e₃ : ℚ) ≠ 0 := by exact_mod_cast h3
  have a23 : (c.e₂ : ℚ) ≠ c.e₃ := by exact_mod_cast h23
  -- every point lands in `Tset c.good`
  have hmemT : ∀ P : (curve c.e₂ c.e₃).Point, Ψ (μ (Multiplicative.ofAdd P)) ∈ Tset c.good := by
    intro P
    rw [μ_apply]
    rcases P with _ | ⟨x, y, h⟩
    · have : Ψ (μ₀ (W := curve c.e₂ c.e₃) 0) = 1 := by rw [μ₀_zero, map_one]
      rw [show (0 : (curve c.e₂ c.e₃).Point) = .zero from rfl] at this
      rw [this]
      unfold Tset
      rw [List.mem_toFinset]
      refine List.mem_map.mpr ⟨(1, 1), h11, ?_⟩
      simp [cls_one]
    · rw [μ₀_some]
      by_cases hx : (curve c.e₂ c.e₃).f.eval x = 0
      · rw [μX_of_eval_f_eq_zero hx]
        rw [eval_f] at hx
        have hx' : x = 0 ∨ x = c.e₂ ∨ x = c.e₃ := by
          rcases mul_eq_zero.mp hx with h1 | h1
          · rcases mul_eq_zero.mp h1 with h2 | h2
            · exact Or.inl h2
            · exact Or.inr (Or.inl (by linarith))
          · exact Or.inr (Or.inr (by linarith))
        rcases hx' with rfl | rfl | rfl
        · exact mem_Tset_of _ (by rw [ev_mk]; simp [curve]; try ring)
            (by rw [ev_mk]; simp [curve]; try ring)
            (by rw [ev_mk]; simp [curve]; try ring) t0
        · exact mem_Tset_of _ (by rw [ev_mk]; simp [curve]; try ring)
            (by rw [ev_mk]; simp [curve]; try ring)
            (by rw [ev_mk]; simp [curve]; try ring) t2
        · exact mem_Tset_of _ (by rw [ev_mk]; simp [curve]; try ring)
            (by rw [ev_mk]; simp [curve]; try ring)
            (by rw [ev_mk]; simp [curve]; try ring) t3
      · rw [μX_of_eval_f_ne_zero hx]
        have hxe : y ^ 2 = x * (x - c.e₂) * (x - c.e₃) := (equation_iff' x y).mp h.1
        rw [eval_f] at hx
        have hx0 : x ≠ 0 := fun h0 => hx (by rw [h0]; ring)
        have hx2 : x - c.e₂ ≠ 0 := fun h0 => hx (by rw [h0]; ring)
        have hx3 : x - c.e₃ ≠ 0 := fun h0 => hx (by rw [h0]; ring)
        obtain ⟨n, e, a, b, t, d₁, d₂, hd₁m, hd₂m, he, ha, hb, ht, hxn, hna, hgb, hgt, hne⟩ :=
          descent_gen c.S₁ c.S₂ hS₁ hS₂ hd₁ hd₂ hxe hx0 hx2 hx3
        have hgood := c.good_of_check hc hT hd₁m hd₂m he ha hb ht hna hgb hgt hne
        have heQ : (e : ℚ) ≠ 0 := by exact_mod_cast he
        have hd1 : (d₁ : ℚ) ≠ 0 := by
          rintro h0; apply hx0
          have : (n : ℚ) = 0 := by rw [show (n : ℚ) = d₁ * a ^ 2 by exact_mod_cast hna, h0]; ring
          rw [← hxn] at this
          exact (mul_eq_zero.mp this).resolve_right (pow_ne_zero 2 heQ)
        have hd2 : (d₂ : ℚ) ≠ 0 := by
          rintro h0; apply hx2
          have : (n : ℚ) - c.e₂ * e ^ 2 = 0 := by
            rw [show (n : ℚ) - c.e₂ * e ^ 2 = d₂ * b ^ 2 by exact_mod_cast hgb, h0]; ring
          have : (x - c.e₂) * (e : ℚ) ^ 2 = 0 := by linear_combination this + hxn
          exact (mul_eq_zero.mp this).resolve_right (pow_ne_zero 2 heQ)
        have hd12 : ((d₁ * d₂ : ℤ) : ℚ) ≠ 0 := by push_cast; exact mul_ne_zero hd1 hd2
        refine mem_Tset_of _ (v₁ := x) (v₂ := x - c.e₂) (v₃ := x - c.e₃)
          (by rw [ev_mk]; simp)
          (by rw [ev_mk]; simp)
          (by rw [ev_mk]; simp)
          ⟨(d₁, d₂), hgood, a / e, b / e, t / ((d₁ * d₂ : ℤ) * e),
            div_ne_zero (by exact_mod_cast ha) heQ, div_ne_zero (by exact_mod_cast hb) heQ,
            div_ne_zero (by exact_mod_cast ht) (mul_ne_zero hd12 heQ), ?_, ?_, ?_⟩
        · have hnQ : (n : ℚ) = d₁ * a ^ 2 := by exact_mod_cast hna
          field_simp
          linear_combination hxn + hnQ
        · have hgQ : (n : ℚ) - c.e₂ * e ^ 2 = d₂ * b ^ 2 := by exact_mod_cast hgb
          field_simp
          linear_combination hxn + hgQ
        · have htQ : (d₁ : ℚ) * d₂ * (n - c.e₃ * e ^ 2) = t ^ 2 := by exact_mod_cast hgt
          push_cast
          field_simp
          linear_combination (d₁ : ℚ) * d₂ * hxn + htQ
  -- counting
  have hcard : Nat.card (μ (W := curve c.e₂ c.e₃)).range ≤ 4 := by
    have hmem : ∀ m ∈ (μ (W := curve c.e₂ c.e₃)).range,
        Ψ (e₂ := c.e₂) (e₃ := c.e₃) m ∈ Tset c.good := by
      rintro _ ⟨P, rfl⟩
      exact hmemT (Multiplicative.toAdd P)
    have hinj : Function.Injective
        (fun m : (μ (W := curve c.e₂ c.e₃)).range =>
          (⟨Ψ (e₂ := c.e₂) (e₃ := c.e₃) (m : (curve c.e₂ c.e₃).M), hmem m m.2⟩ :
            Tset c.good)) := by
      intro m m' hmm
      exact Subtype.ext (Ψ_injective h2 h3 h23 (congrArg Subtype.val hmm))
    calc Nat.card (μ (W := curve c.e₂ c.e₃)).range ≤ Nat.card (Tset c.good) :=
          Nat.card_le_card_of_injective _ hinj
      _ = (Tset c.good).card := by rw [Nat.card_eq_fintype_card, Fintype.card_coe]
      _ ≤ 4 := (card_Tset _).trans hlen
  have hroots : Nat.card {x : ℚ // (curve c.e₂ c.e₃).f.eval x = 0} = 3 := by
    have : {x : ℚ // (curve c.e₂ c.e₃).f.eval x = 0} ≃
        ({0, (c.e₂ : ℚ), (c.e₃ : ℚ)} : Finset ℚ) := by
      refine Equiv.subtypeEquivRight fun x => ?_
      rw [eval_f]
      simp only [Finset.mem_insert, Finset.mem_singleton]
      constructor
      · intro hx
        rcases mul_eq_zero.mp hx with h1 | h1
        · rcases mul_eq_zero.mp h1 with h2 | h2
          · exact Or.inl h2
          · exact Or.inr (Or.inl (by linarith))
        · exact Or.inr (Or.inr (by linarith))
      · rintro (rfl | rfl | rfl) <;> ring
    rw [Nat.card_congr this, Nat.card_eq_fintype_card, Fintype.card_coe]
    rw [Finset.card_insert_of_notMem, Finset.card_pair a23]
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or]
    exact ⟨Ne.symm a2, Ne.symm a3⟩
  have h := card_range_μ (W := curve c.e₂ c.e₃)
  rw [card_ker_nsmul_two, hroots] at h
  rw [h] at hcard
  by_contra hr
  have := Nat.le_self_pow hr 2
  omega

/-- Primes dividing `N` lie in `S` when `|N|` is a product of primes of `S`. -/
theorem cover_of_factor {N : ℤ} {F S : List ℕ} (hF : ∀ q ∈ F, q.Prime) (hFS : ∀ q ∈ F, q ∈ S)
    (hN : N.natAbs = F.prod) : ∀ p : ℕ, p.Prime → (p : ℤ) ∣ N → p ∈ S := by
  intro p hp hpN
  have hpN' : p ∣ F.prod := by rw [← hN]; exact Int.natCast_dvd.mp hpN
  clear hN hpN
  induction F with
  | nil => exact absurd (Nat.dvd_one.mp hpN') hp.one_lt.ne'
  | cons q F ih =>
    rw [List.prod_cons] at hpN'
    rcases (Nat.Prime.dvd_mul hp).mp hpN' with h | h
    · have := (Nat.prime_dvd_prime_iff_eq hp (hF q (by simp))).mp h
      rw [this]; exact hFS q (by simp)
    · exact ih (fun r hr => hF r (by simp [hr])) (fun r hr => hFS r (by simp [hr])) h

theorem int_prime (p : ℕ) (hp : p.Prime) : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp

end LeanRank.RankZeroGen

#print axioms LeanRank.RankZeroGen.Ψ_injective
#print axioms LeanRank.RankZeroGen.finrank_eq_zero_of
#print axioms LeanRank.RankZeroGen.cover_of_factor
