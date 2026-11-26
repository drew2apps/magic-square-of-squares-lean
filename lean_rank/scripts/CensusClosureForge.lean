import LeanRank.CensusClosure

/-!
# Forged claims against `LeanRank.CensusClosure` (each must FAIL)

Written by an AI model; not read by a person.  Not part of the library.  Run with
`lake env lean scripts/CensusClosureForge.lean` from `lean_rank/`.

**Expected: exactly 13 errors, one per forged claim F1–F13, and no other error.**  (Each forged
claim is one `example`; the error is the one named in its comment.)

What is refused, precisely: a wrong point count over `𝔽₁₁` (F1, F2); a wrong reduction model (F3);
a prime of bad reduction (F4); the prime `2`, where reduction need not be injective on torsion
(F5); a point that is not on the curve (F6); a list of 16 points with a repeat (F7); a list of the
affine points with one missing (F8); a bound `#E(ℚ) ≤ 8` for a curve with 16 points (F9); a
non-square certificate for a square (F10); the magic-square theorem applied to the degenerate
square at the ratio `4 : 3` (`τ = 1`), where `b² = f²` (F11); a `τ`-statement without `τ = −1`
(F12); a "partner" `18 : 7` of `12 : 5` (the partner is `17 : 7`) (F13).  Certificates are not unique, so this does not say that every wrong certificate is
refused; it says these are.
-/

open WeierstrassCurve WeierstrassCurve.Affine
open LeanRank.RankZeroGen LeanRank.CensusClosure
open LeanRank.RankZero (intPrime)

-- F1: the reduction of `(4, 3)` modulo 11 does not have 15 points (`decide` evaluates to false).
example : Fintype.card Ep_4_3.toAffine.Point = 15 := by decide +kernel

-- F2: the reduction of `(7, 1)` modulo 11 does not have 16 points.
example : Fintype.card Ep_7_1.toAffine.Point = 16 := by decide +kernel

-- F3: a wrong reduced model (`a₄ = 2` instead of `1`).
example : (((W₀ 256 81).toAffine ⁄ (ZMod 11)) : WeierstrassCurve _).toAffine =
    (⟨0, 7, 0, 2, 0⟩ : WeierstrassCurve (ZMod 11)).toAffine := by
  ext <;> decide +kernel

-- F4: `5` is a prime of bad reduction of `W₀ 256 81` (`5 ∣ 256 − 81 = 175`).
example : (redCurve (intPrime 5) (W₀ 256 81)).IsElliptic :=
  redCurve_isElliptic 5 256 81 (by norm_num)

-- F5: at `p = 2` the condition `p ∉ (p)^(p−1)` fails (`norm_num` cannot prove `¬ 2 ∣ 2`).
example : ¬ ((2 : ℤ) ^ (2 - 1) ∣ (2 : ℤ)) := by norm_num

-- F6: `(24, 841)` is not on `Y² = X(X + 256)(X + 81)` (`840² = 24·280·105`).
example : (curve (-256) (-81)).Nonsingular 24 841 := ns (by norm_num)

-- F7: sixteen points with `(864, 30240)` twice and `(864, −30240)` missing are not distinct.
noncomputable def forged_pts : Fin 16 → (curve (-256) (-81)).Point :=
  ![0, .some 0 0 (ns (by norm_num)), .some (-81) 0 (ns (by norm_num)),
    .some (-256) 0 (ns (by norm_num)), .some 24 840 (ns (by norm_num)),
    .some 24 (-840) (ns (by norm_num)), .some (-96) 480 (ns (by norm_num)),
    .some (-96) (-480) (ns (by norm_num)), .some (-144) 1008 (ns (by norm_num)),
    .some (-144) (-1008) (ns (by norm_num)), .some 144 3600 (ns (by norm_num)),
    .some 144 (-3600) (ns (by norm_num)), .some (-216) 1080 (ns (by norm_num)),
    .some (-216) (-1080) (ns (by norm_num)), .some 864 30240 (ns (by norm_num)),
    .some 864 30240 (ns (by norm_num))]

example : Function.Injective (key ∘ forged_pts) := by decide +kernel

-- F8: the affine points of `(4, 3)` without `(864, −30240)`: the case of that point is left.
example {x y : ℚ} (h : y ^ 2 = x * (x + 256) * (x + 81)) :
    (x, y) ∈ ([(0, 0), ((-81), 0), ((-256), 0), (24, 840), (24, (-840)), ((-96), 480),
      ((-96), (-480)), ((-144), 1008), ((-144), (-1008)), (144, 3600), (144, (-3600)),
      ((-216), 1080), ((-216), (-1080)), (864, 30240)] : List (ℚ × ℚ)) := by
  obtain ⟨i, hi⟩ := points_4_3 (.some x y (ns (by push_cast; linear_combination h)))
  have hk := congrArg key hi
  fin_cases i <;> simp [pts_4_3, key] at hk <;> obtain ⟨rfl, rfl⟩ := hk <;> simp

-- F9: `#E(ℚ) ≤ 8` for `(4, 3)` through the same lemma: the count `8` over `𝔽₁₁` is refused.
example : Nat.card (curve (-256) (-81)).Point ≤ 8 :=
  (card_le 11 256 81 (by norm_num) LeanRank.RankZeroCensus.finrank_4_3 Ep_4_3 base_4_3
    (N := 8) (by decide +kernel)).2

-- F10: `576 = 144·2²` is `144` times a square; no bracket `k² < 4 < (k + 1)²` exists.
example {τ : ℚ} (h : 144 * τ ^ 2 = 576) : False :=
  tau_nonsq (by norm_num) h (N := 4) (k := 1) (Q := 1) (by norm_num) (by norm_num) (by norm_num)

-- F11: the degenerate square `25 49 1 / 1 25 49 / 49 1 25` (roots `5 7 1 / 1 5 7 / 7 −1 5`) has
-- `(b + d) : (h + f) = 8 : 6 = 4 : 3`; the hypothesis `b² ≠ f²` fails (`49 = 49`).
example : ((7 : ℤ) + 1) * 3 ≠ ((-1 : ℤ) + 7) * 4 :=
  ratio_ne_4_3 (a := 5) (b := 7) (c := 1) (d := 1) (e := 5) (f := 7) (g := 7) (h := -1) (i := 5)
    (T := 75) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by norm_num)
    (by norm_num) (by norm_num) (by norm_num) (by norm_num)

-- F12: `τ = −1` cannot be dropped from `tau_4_3` (`F0 = F8 = 25` there).
example {τ r0 r8 : ℚ} (h0 : r0 ^ 2 = 16 + 9 * τ ^ 2) (h8 : r8 ^ 2 = 9 + 16 * τ ^ 2) :
    τ = 0 ∨ τ = 1 := by
  rcases tau_4_3 h0 h8 with h | h | h
  · exact Or.inl h
  · exact Or.inr h
  · exact absurd h (by norm_num)

-- F13: `18 : 7` is not the partner of `12 : 5` (`(18 + 7)·5 = 125 ≠ 132 = (18 − 7)·12`).
example {a b c d e f g h i T : ℤ}
    (r₁ : a ^ 2 + b ^ 2 + c ^ 2 = T) (r₂ : d ^ 2 + e ^ 2 + f ^ 2 = T)
    (r₃ : g ^ 2 + h ^ 2 + i ^ 2 = T) (c₁ : a ^ 2 + d ^ 2 + g ^ 2 = T)
    (c₂ : b ^ 2 + e ^ 2 + h ^ 2 = T) (c₃ : c ^ 2 + f ^ 2 + i ^ 2 = T)
    (dg : a ^ 2 + e ^ 2 + i ^ 2 = T) (ad : c ^ 2 + e ^ 2 + g ^ 2 = T)
    (hbd : b ^ 2 ≠ d ^ 2) (hbf : b ^ 2 ≠ f ^ 2) : (b + d) * (7 : ℕ) ≠ (h + f) * (18 : ℕ) :=
  ratio_ne_partner (m := 18) (n := 7) (M := 12) (N := 5) (by norm_num) (by norm_num)
    (fun a b c d e f g h i T r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf => by
      have := ratio_ne_12_5 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf
      exact_mod_cast this)
    r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf
