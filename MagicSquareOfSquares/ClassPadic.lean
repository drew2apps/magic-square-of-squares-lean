import MagicSquareOfSquares.ClassHensel

/-!
# The step to `ℤ_p` (the last step of Proposition C of `COMPLETE.md`, section 5)

`ClassHensel.lean` ends with: every class of the table of Theorem C at an odd prime of
`mn(m² − n²)` is the class modulo every power of `p` at one primitive integer point
(`exists_point_of_table`). Here that becomes a statement about `ℤ_p`.

* `padic_sq_of_forall`, `padic_gsq_of_forall`: an integer that is a square modulo every power of
  `p` is a square in `ℤ_p`, and the same for a Gaussian integer `X + iY` and the form
  `(a² − b², 2ab)`. The proofs use that `ℤ_p` is compact: a sequence of roots modulo `p^k` has a
  convergent subsequence, and its limit is a root. (`COMPLETE.md` does this step by Hensel's lemma;
  compactness needs no case distinction on `−1` being a square.)
* `PadicClassAt p d m n u v`: at the integer point `(u, v)`, the four norms are nonzero squares in
  `ℤ_p`, and `d` times each of the eight products is `(a² − b², 2ab)` with `a, b ∈ ℤ_p`. By
  Lemma 1 of `COMPLETE.md` read over `ℤ_p` (by hand, not in Lean) this says that `(u, v)` is a
  point of `X_p` with class `d`, for a point with integer coordinates.
* `padicClassAt_of_classAt`: `ClassAt` modulo every power of `p` gives `PadicClassAt`.
* `exists_padic_point_of_table`: **for every class `d` of the table of Theorem C at an odd prime
  of `mn(m² − n²)` there is an integer point, primitive at `p`, with `PadicClassAt p d`.**
* `not_padicClassAt_three`: a negative control. At the prime 3, for the ratio `3 : 1`, the point
  `(0, 1)` is not a point of `X_3` (`N_L = 20` is not a square in `ℤ_3`), for any `d`. This is the
  point at which the first version of `Passes` accepted the class 3.
* `classAt_sq_of_padicClassAt`: the other direction at the modulus `p²` (reduction by
  `PadicInt.toZModPow`), and `not_padicClassAt_three_class`, a negative control for the class: for
  `3 : 1` no integer point that is primitive at 3 has `PadicClassAt 3 3`.

Not here: Lemma 1 over `ℤ_p`, the upper bounds for points of `X_p` (Lean has them for integer
solutions), points of `X_p` that do not have integer coordinates, and the primes 2 and those of
`(m² + n²)W`.
-/

open Filter Topology

namespace MagicSquareOfSquares.CoverDescent

variable {p : ℕ} [hp : Fact p.Prime]

/-- A sequence of integers whose `k`-th term is divisible by `p^(φ k)`, `φ` strictly increasing,
tends to 0 in `ℤ_p`. -/
theorem padic_tendsto_zero_of_dvd {f : ℕ → ℤ} {φ : ℕ → ℕ} (hφ : StrictMono φ)
    (h : ∀ k, (p : ℤ) ^ φ k ∣ f k) :
    Tendsto (fun k => ((f k : ℤ) : ℤ_[p])) atTop (𝓝 0) := by
  have hp1 : (1 : ℝ) < p := by exact_mod_cast hp.out.one_lt
  have hgeo : Tendsto (fun k : ℕ => ((p : ℝ)⁻¹) ^ k) atTop (𝓝 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one (by positivity) (inv_lt_one_of_one_lt₀ hp1)
  have hl : Tendsto (fun k => (p : ℝ) ^ (-(φ k : ℤ))) atTop (𝓝 0) := by
    refine (hgeo.comp hφ.tendsto_atTop).congr (fun k => ?_)
    simp [zpow_neg, inv_pow]
  refine squeeze_zero_norm (fun k => ?_) hl
  exact PadicInt.norm_int_le_pow_iff_dvd.mpr (by exact_mod_cast h k)

/-- **An integer that is a square modulo every power of `p` is a square in `ℤ_p`.** -/
theorem padic_sq_of_forall {x : ℤ} (h : ∀ k : ℕ, ∃ y : ℤ, x ≡ y ^ 2 [ZMOD (p : ℤ) ^ k]) :
    ∃ w : ℤ_[p], (x : ℤ_[p]) = w ^ 2 := by
  choose y hy using h
  obtain ⟨a, -, φ, hφ, ha⟩ := isCompact_univ.tendsto_subseq
    (x := fun k => ((y k : ℤ) : ℤ_[p])) (fun _ => Set.mem_univ _)
  refine ⟨a, ?_⟩
  have h1 : Tendsto (fun k => (((x - y (φ k) ^ 2 : ℤ)) : ℤ_[p])) atTop
      (𝓝 ((x : ℤ_[p]) - a ^ 2)) := by
    refine (tendsto_const_nhds.sub (ha.pow 2)).congr (fun k => ?_)
    simp
  have h2 := padic_tendsto_zero_of_dvd (p := p) (f := fun k => x - y (φ k) ^ 2) hφ
    (fun k => Int.modEq_iff_dvd.mp (hy (φ k)).symm)
  exact sub_eq_zero.mp (tendsto_nhds_unique h1 h2)

/-- **A Gaussian integer `X + iY` that is a square modulo every power of `p` is a square in
`ℤ_p[i]`**, in coordinates: `X = a² − b²`, `Y = 2ab` with `a, b ∈ ℤ_p`. -/
theorem padic_gsq_of_forall {X Y : ℤ}
    (h : ∀ k : ℕ, ∃ a b : ℤ, X ≡ a ^ 2 - b ^ 2 [ZMOD (p : ℤ) ^ k] ∧
      Y ≡ 2 * a * b [ZMOD (p : ℤ) ^ k]) :
    ∃ a b : ℤ_[p], (X : ℤ_[p]) = a ^ 2 - b ^ 2 ∧ (Y : ℤ_[p]) = 2 * a * b := by
  choose a b hab using h
  obtain ⟨c, -, φ, hφ, hc⟩ := isCompact_univ.tendsto_subseq
    (x := fun k => (((a k : ℤ) : ℤ_[p]), ((b k : ℤ) : ℤ_[p]))) (fun _ => Set.mem_univ _)
  have ha : Tendsto (fun k => ((a (φ k) : ℤ) : ℤ_[p])) atTop (𝓝 c.1) :=
    (continuous_fst.tendsto c).comp hc
  have hb : Tendsto (fun k => ((b (φ k) : ℤ) : ℤ_[p])) atTop (𝓝 c.2) :=
    (continuous_snd.tendsto c).comp hc
  refine ⟨c.1, c.2, ?_, ?_⟩
  · have h1 : Tendsto (fun k => (((X - (a (φ k) ^ 2 - b (φ k) ^ 2) : ℤ)) : ℤ_[p])) atTop
        (𝓝 ((X : ℤ_[p]) - (c.1 ^ 2 - c.2 ^ 2))) := by
      refine (tendsto_const_nhds.sub ((ha.pow 2).sub (hb.pow 2))).congr (fun k => ?_)
      simp
    have h2 := padic_tendsto_zero_of_dvd (p := p)
      (f := fun k => X - (a (φ k) ^ 2 - b (φ k) ^ 2)) hφ
      (fun k => Int.modEq_iff_dvd.mp (hab (φ k)).1.symm)
    exact sub_eq_zero.mp (tendsto_nhds_unique h1 h2)
  · have h1 : Tendsto (fun k => (((Y - 2 * a (φ k) * b (φ k) : ℤ)) : ℤ_[p])) atTop
        (𝓝 ((Y : ℤ_[p]) - 2 * c.1 * c.2)) := by
      refine (tendsto_const_nhds.sub ((tendsto_const_nhds.mul ha).mul hb)).congr (fun k => ?_)
      simp
    have h2 := padic_tendsto_zero_of_dvd (p := p)
      (f := fun k => Y - 2 * a (φ k) * b (φ k)) hφ
      (fun k => Int.modEq_iff_dvd.mp (hab (φ k)).2.symm)
    exact sub_eq_zero.mp (tendsto_nhds_unique h1 h2)

omit hp in
/-- An integer that is, modulo `p^6`, the square of an integer not divisible by `p³` is not 0. -/
theorem ne_zero_of_sq_modEq {x y : ℤ} (h : x ≡ y ^ 2 [ZMOD (p : ℤ) ^ 6]) (hy : ¬(p : ℤ) ^ 3 ∣ y) :
    x ≠ 0 := by
  rintro rfl
  have h1 : (p : ℤ) ^ 6 ∣ y ^ 2 := by
    have := Int.modEq_iff_dvd.mp h
    simpa using this
  have e : (p : ℤ) ^ 6 = ((p : ℤ) ^ 3) ^ 2 := by ring
  rw [e] at h1
  exact hy ((Int.pow_dvd_pow_iff two_ne_zero).mp h1)

/-- **`d` is the class of the integer point `(u, v)` in `ℤ_p`**: the four norms are nonzero squares
in `ℤ_p`, and `d` times each of the eight products is a Gaussian square over `ℤ_p`, that is, of the
form `(a² − b², 2ab)` with `a, b ∈ ℤ_p`. -/
def PadicClassAt (p : ℕ) [Fact p.Prime] (d m n u v : ℤ) : Prop :=
  (∃ w0 w8 wL wU : ℤ_[p], (w0 ≠ 0 ∧ w8 ≠ 0 ∧ wL ≠ 0 ∧ wU ≠ 0) ∧
    ((nrm (P0 m n u v) : ℤ) : ℤ_[p]) = w0 ^ 2 ∧ ((nrm (P8 m n u v) : ℤ) : ℤ_[p]) = w8 ^ 2 ∧
    ((nrm (PL m n u v) : ℤ) : ℤ_[p]) = wL ^ 2 ∧ ((nrm (PU m n u v) : ℤ) : ℤ_[p]) = wU ^ 2) ∧
  ∀ f : Bool × Bool × Bool, ∃ a b : ℤ_[p],
    ((d * (Q m n u v f).1 : ℤ) : ℤ_[p]) = a ^ 2 - b ^ 2 ∧
    ((d * (Q m n u v f).2 : ℤ) : ℤ_[p]) = 2 * a * b

/-- One norm: a square modulo every power of `p` with roots not divisible by `p³` is a nonzero
square in `ℤ_p`. -/
theorem padic_norm_sq {x : ℤ}
    (h : ∀ k : ℕ, ∃ y : ℤ, x ≡ y ^ 2 [ZMOD (p : ℤ) ^ k] ∧ ¬(p : ℤ) ^ 3 ∣ y) :
    ∃ w : ℤ_[p], w ≠ 0 ∧ (x : ℤ_[p]) = w ^ 2 := by
  obtain ⟨y6, h6, n6⟩ := h 6
  have hx : x ≠ 0 := ne_zero_of_sq_modEq h6 n6
  obtain ⟨w, hw⟩ := padic_sq_of_forall (p := p) (x := x) (fun k => ⟨(h k).choose, (h k).choose_spec.1⟩)
  refine ⟨w, ?_, hw⟩
  rintro rfl
  have : ((x : ℤ) : ℤ_[p]) = 0 := by simpa using hw
  exact hx (by exact_mod_cast this)

/-- **From every power of `p` to `ℤ_p`**, at one integer point. -/
theorem padicClassAt_of_classAt {d m n u v : ℤ}
    (h : ∀ k : ℕ, ClassAt ((p : ℤ) ^ k) ((p : ℤ) ^ 3) d m n u v) : PadicClassAt p d m n u v := by
  obtain ⟨w0, n0, e0⟩ := padic_norm_sq (p := p) (x := nrm (P0 m n u v)) (fun k => by
    obtain ⟨a0, _, _, _, ⟨h0, _, _, _⟩, ⟨q0, _, _, _⟩, _⟩ := h k
    exact ⟨a0, h0, q0⟩)
  obtain ⟨w8, n8, e8⟩ := padic_norm_sq (p := p) (x := nrm (P8 m n u v)) (fun k => by
    obtain ⟨_, a8, _, _, ⟨_, h8, _, _⟩, ⟨_, q8, _, _⟩, _⟩ := h k
    exact ⟨a8, h8, q8⟩)
  obtain ⟨wL, nL, eL⟩ := padic_norm_sq (p := p) (x := nrm (PL m n u v)) (fun k => by
    obtain ⟨_, _, aL, _, ⟨_, _, hL, _⟩, ⟨_, _, qL, _⟩, _⟩ := h k
    exact ⟨aL, hL, qL⟩)
  obtain ⟨wU, nU, eU⟩ := padic_norm_sq (p := p) (x := nrm (PU m n u v)) (fun k => by
    obtain ⟨_, _, _, aU, ⟨_, _, _, hU⟩, ⟨_, _, _, qU⟩, _⟩ := h k
    exact ⟨aU, hU, qU⟩)
  refine ⟨⟨w0, w8, wL, wU, ⟨n0, n8, nL, nU⟩, e0, e8, eL, eU⟩, fun f => ?_⟩
  exact padic_gsq_of_forall (p := p) (fun k => (h k).gsq f)

/-- **The lower bounds of Theorem C at an odd prime of `mn(m² − n²)`, as a statement about
`ℤ_p`**: for coprime `m, n` and an odd prime `p` dividing one of `m`, `n`, `m + n`, `m − n`, every
class `d` of the row of that factor has `PadicClassAt p d` at an integer point that is primitive
at `p` (so, by Lemma 1 over `ℤ_p`, which is by hand, `d` is the class of a point of `X_p`). -/
theorem exists_padic_point_of_table (hp2 : p ≠ 2) {m n d : ℤ} (hc : IsCoprime m n)
    (h : ((p : ℤ) ∣ m ∧ Row p m d) ∨ ((p : ℤ) ∣ n ∧ Row p n d) ∨
      ((p : ℤ) ∣ m + n ∧ Row p (m + n) d) ∨ ((p : ℤ) ∣ m - n ∧ Row p (m - n) d)) :
    ∃ u v : ℤ, (¬(p : ℤ) ∣ u ∨ ¬(p : ℤ) ∣ v) ∧ PadicClassAt p d m n u v := by
  obtain ⟨u, v, huv, hcl⟩ := exists_point_of_table hp.out hp2 hc h
  exact ⟨u, v, huv, padicClassAt_of_classAt hcl⟩

/-- **A negative control**: at the prime 3 and the ratio `3 : 1`, the point `(0, 1)` is not a
point of `X_3`, whatever `d` is, because `N_L(0, 1) = 20` is not a square in `ℤ_3`. The eight
products alone accept the class 3 there (`COMPLETE.md`, section 5). -/
theorem not_padicClassAt_three (d : ℤ) : ¬PadicClassAt 3 d 3 1 0 1 := by
  rintro ⟨⟨w0, w8, wL, wU, -, -, -, hL, -⟩, -⟩
  have e : ((nrm (PL 3 1 0 1) : ℤ) : ℤ_[3]) = 20 := by norm_num [nrm, PL]
  rw [e] at hL
  have h := congrArg (PadicInt.toZMod (p := 3)) hL
  rw [map_pow, map_ofNat] at h
  have key : ∀ z : ZMod 3, (20 : ZMod 3) ≠ z ^ 2 := by decide
  exact key _ h

/-! ### Back from `ℤ_p` to a modulus

The other direction, at the modulus `p²` only: a root in `ℤ_p` reduces to a root modulo `p^k`
(`PadicInt.toZModPow`). At `p²` a root can be chosen not divisible by `p³`, so `PadicClassAt` gives
`ClassAt (p²) (p³)`, and a class that fails modulo `p²` at every primitive point is not the class of
any primitive integer point over `ℤ_p`. -/

/-- An integer that is a square in `ℤ_p` is a square modulo every power of `p`. -/
theorem modEq_sq_of_padic (k : ℕ) {x : ℤ} {w : ℤ_[p]} (h : (x : ℤ_[p]) = w ^ 2) :
    ∃ z : ℤ, x ≡ z ^ 2 [ZMOD (p : ℤ) ^ k] := by
  obtain ⟨z, hz⟩ := ZMod.intCast_surjective (PadicInt.toZModPow k w)
  refine ⟨z, ?_⟩
  have h' := congrArg (PadicInt.toZModPow k) h
  rw [map_pow, map_intCast, ← hz] at h'
  have e : ((x : ℤ) : ZMod (p ^ k)) = ((z ^ 2 : ℤ) : ZMod (p ^ k)) := by
    push_cast; exact h'
  have := (ZMod.intCast_eq_intCast_iff x (z ^ 2) (p ^ k)).mp e
  simpa using this

/-- A pair of integers that is a Gaussian square over `ℤ_p` is one modulo every power of `p`. -/
theorem modEq_gsq_of_padic (k : ℕ) {X Y : ℤ} {a b : ℤ_[p]}
    (h1 : (X : ℤ_[p]) = a ^ 2 - b ^ 2) (h2 : (Y : ℤ_[p]) = 2 * a * b) :
    ∃ a' b' : ℤ, X ≡ a' ^ 2 - b' ^ 2 [ZMOD (p : ℤ) ^ k] ∧ Y ≡ 2 * a' * b' [ZMOD (p : ℤ) ^ k] := by
  obtain ⟨a', ha⟩ := ZMod.intCast_surjective (PadicInt.toZModPow k a)
  obtain ⟨b', hb⟩ := ZMod.intCast_surjective (PadicInt.toZModPow k b)
  refine ⟨a', b', ?_, ?_⟩
  · have h' := congrArg (PadicInt.toZModPow k) h1
    rw [map_sub, map_pow, map_pow, map_intCast, ← ha, ← hb] at h'
    have e : ((X : ℤ) : ZMod (p ^ k)) = ((a' ^ 2 - b' ^ 2 : ℤ) : ZMod (p ^ k)) := by
      push_cast; exact h'
    have := (ZMod.intCast_eq_intCast_iff X (a' ^ 2 - b' ^ 2) (p ^ k)).mp e
    simpa using this
  · have h' := congrArg (PadicInt.toZModPow k) h2
    rw [map_mul, map_mul, map_ofNat, map_intCast, ← ha, ← hb] at h'
    have e : ((Y : ℤ) : ZMod (p ^ k)) = ((2 * a' * b' : ℤ) : ZMod (p ^ k)) := by
      push_cast; exact h'
    have := (ZMod.intCast_eq_intCast_iff Y (2 * a' * b') (p ^ k)).mp e
    simpa using this

/-- Modulo `p²` a square root can be chosen not divisible by `p³`. -/
theorem exists_sq_root_not_dvd {x z : ℤ} (h : x ≡ z ^ 2 [ZMOD (p : ℤ) ^ 2]) :
    ∃ z' : ℤ, x ≡ z' ^ 2 [ZMOD (p : ℤ) ^ 2] ∧ ¬(p : ℤ) ^ 3 ∣ z' := by
  by_cases hz : (p : ℤ) ^ 3 ∣ z
  · refine ⟨z + (p : ℤ) ^ 2, h.trans (Int.ModEq.pow 2 ?_), fun hd => ?_⟩
    · exact Int.modEq_iff_dvd.mpr ⟨1, by ring⟩
    · have h32 : (p : ℤ) ^ 3 ∣ (p : ℤ) ^ 2 := (Int.dvd_add_right hz).mp hd
      have hp2 : (2 : ℤ) ≤ p := by exact_mod_cast hp.out.two_le
      have hpos : (0 : ℤ) < (p : ℤ) ^ 2 := by positivity
      have hle := Int.le_of_dvd hpos h32
      have e : (p : ℤ) ^ 3 = (p : ℤ) ^ 2 * p := by ring
      nlinarith [mul_le_mul_of_nonneg_left hp2 hpos.le]
  · exact ⟨z, h, hz⟩

/-- `PadicClassAt` gives `ClassAt` modulo `p²`. -/
theorem classAt_sq_of_padicClassAt {d m n u v : ℤ} (h : PadicClassAt p d m n u v) :
    ClassAt ((p : ℤ) ^ 2) ((p : ℤ) ^ 3) d m n u v := by
  obtain ⟨⟨w0, w8, wL, wU, -, h0, h8, hL, hU⟩, hg⟩ := h
  obtain ⟨z0, e0⟩ := modEq_sq_of_padic 2 h0
  obtain ⟨z0, e0, n0⟩ := exists_sq_root_not_dvd e0
  obtain ⟨z8, e8⟩ := modEq_sq_of_padic 2 h8
  obtain ⟨z8, e8, n8⟩ := exists_sq_root_not_dvd e8
  obtain ⟨zL, eL⟩ := modEq_sq_of_padic 2 hL
  obtain ⟨zL, eL, nL⟩ := exists_sq_root_not_dvd eL
  obtain ⟨zU, eU⟩ := modEq_sq_of_padic 2 hU
  obtain ⟨zU, eU, nU⟩ := exists_sq_root_not_dvd eU
  refine ⟨z0, z8, zL, zU, ⟨e0, e8, eL, eU⟩, ⟨n0, n8, nL, nU⟩, fun f => ?_⟩
  obtain ⟨a, b, h1, h2⟩ := hg f
  exact modEq_gsq_of_padic 2 h1 h2

/-- Negative control for the class. For the ratio `3 : 1` the class 3 is outside the table at 3
(`not_passes_three`), and it is not the class over `ℤ_3` of any primitive integer point. -/
theorem not_padicClassAt_three_class (u v : ℤ) (huv : ¬(3 : ℤ) ∣ u ∨ ¬(3 : ℤ) ∣ v) :
    ¬PadicClassAt 3 3 3 1 u v := fun h =>
  not_passes_three ⟨u, v, by exact_mod_cast huv, classAt_sq_of_padicClassAt h⟩

end MagicSquareOfSquares.CoverDescent
