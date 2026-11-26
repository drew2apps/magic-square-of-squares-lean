import LeanRank.TorsionE1

/-!
# Forged claims for `TorsionE1.lean` (negative tests)

Written by an AI model (Claude, made by Anthropic) on 26 November 2026 (UTC); not read by a person.

**Expected: exactly 7 errors**, one per forged claim (F1–F7), and no error in the accepted controls (C1–C3).

* F1 a wrong residue map `ω ↦ 2` into `ZMod 7` (`2² = 4 ≠ 2`);
* F2 a wrong reduced curve (`a₄ = 4` instead of `3`);
* F3 a wrong point count `#Ẽ(𝔽₇) = 11`;
* F4 the list of `W7_points` with `(6, 6)` dropped;
* F5 a wrong doubling `2A = C` (the true one is `2A = −C`);
* F6 the ramification certificate with `ZMod 21` in place of `ZMod 49` (`7 · 7 ≠ 0` there);
* F7 a forged point `(4, 4√2)` on `E1K`.
* C1 `#Ẽ(𝔽₇) = 10`; C2 dropping `−D` from the torsion list makes `torsion_E1K` false (the list is sharp);
  C3 the torsion statement itself, restated.
-/

namespace MagicSquareOfSquares.Family.TorsionE1.Forge

open WeierstrassCurve WeierstrassCurve.Affine

-- F1: wrong residue map
noncomputable def badPsi : O →ₐ[ℤ] ZMod 7 := QuadraticAlgebra.lift ⟨2, by decide⟩

-- F2: wrong reduced curve
example : ((W0.toAffine ⁄ (ZMod 7)) : WeierstrassCurve _).toAffine =
    (⟨0, 5, 0, 4, 0⟩ : WeierstrassCurve (ZMod 7)).toAffine := by
  change W0.map (algebraMap O (ZMod 7)) = _
  ext <;> simp [W0, algebraMap_omega_7] <;> decide

-- F3: wrong point count
example : Fintype.card W7.toAffine.Point = 11 := by decide

-- F4: a point of `W7(𝔽₇)` dropped
example : ∀ x y : ZMod 7, W7.toAffine.Equation x y →
    (x = 0 ∧ y = 0) ∨ (x = 4 ∧ y = 3) ∨ (x = 4 ∧ y = 4) ∨ (x = 3 ∧ y = 5) ∨ (x = 3 ∧ y = 2) ∨
    (x = 1 ∧ y = 3) ∨ (x = 1 ∧ y = 4) ∨ (x = 6 ∧ y = 1) := by decide

-- F5: wrong doubling `2A = C`
example : pA + pA = pC := double_some (L := 2 + 4 * s2) nsA nsC (q_ne_zero (Or.inr (by norm_num)))
  (by simp only [q_eq]; push_cast; linear_combination (-72) * s2_sq)
  (by simp only [q_eq]; push_cast; linear_combination 16 * s2_sq)
  (by simp only [q_eq]; push_cast; linear_combination (-32) * s2_sq)

-- F6: the ramification certificate fails modulo 21
example : ∀ a b : ZMod 21, ZMod.castHom (show 7 ∣ 21 by norm_num) (ZMod 7) a = 0 →
    ZMod.castHom (show 7 ∣ 21 by norm_num) (ZMod 7) b = 0 → a * b = 0 := by decide

-- F7: forged point `(4, 4√2)`
example : E1K.Nonsingular (q 4 0) (q 0 4) :=
  E1K_ns (by simp only [q_eq]; push_cast; linear_combination 96 * s2_sq)

-- C1 (accepted): `#Ẽ(𝔽₇) = 10`
example : Fintype.card W7.toAffine.Point = 10 := card_W7

-- C2 (accepted): the list is sharp: without `−D` the torsion statement is false
theorem drop_Dm_false : ¬ ∀ P : E1K.Point, IsOfFinAddOrder P →
    P = 0 ∨ P = pT ∨ P = pA ∨ P = pAm ∨ P = pB ∨ P = pBm ∨ P = pC ∨ P = pCm ∨ P = pD := by
  intro hall
  have hx (P : E1K.Point) (a b c e : ℤ) (h : E1K.Nonsingular (q a b) (q c e)) (hP : P = .some _ _ h)
      (hne : -12 ≠ a) : pDm ≠ P := by
    rw [hP]; intro h'
    exact hne (q_inj_re ((Point.some.injEq _ _ _ _ _ _).mp h').1)
  rcases hall pDm tors_Dm with h | h | h | h | h | h | h | h | h
  · exact Point.some_ne_zero _ h
  · exact hx pT 0 0 0 0 nsT rfl (by norm_num) h
  · exact hx pA 4 0 0 8 nsA rfl (by norm_num) h
  · exact hx pAm 4 0 0 (-8) nsAm rfl (by norm_num) h
  · exact hx pB (-4) 0 16 8 nsB rfl (by norm_num) h
  · exact hx pBm (-4) 0 (-16) (-8) nsBm rfl (by norm_num) h
  · exact hx pC 12 8 80 56 nsC rfl (by norm_num) h
  · exact hx pCm 12 8 (-80) (-56) nsCm rfl (by norm_num) h
  · -- same abscissa as `−D`: compare ordinates (`32 ≠ −32`)
    have hy := ((Point.some.injEq _ _ _ _ _ _).mp h).2
    exact absurd (q_inj_re hy) (by norm_num)

-- C3 (accepted): the headline, restated
example (P : E1K.Point) (hP : IsOfFinAddOrder P) :
    P = 0 ∨ P = pT ∨ P = pA ∨ P = pAm ∨ P = pB ∨ P = pBm ∨ P = pC ∨ P = pCm ∨ P = pD ∨ P = pDm :=
  torsion_E1K P hP

end MagicSquareOfSquares.Family.TorsionE1.Forge
