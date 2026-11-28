import MagicSquareOfSquares.CoverDescent
import MagicSquareOfSquares.G1Full

/-!
# The exchange symmetry: the edge ratio `m : n` and the parameter `u : v` can be exchanged

The four corner pairs `P0, P8, PL, PU` of the G1 family and the centre form are unchanged, up to the
order of the two coordinates of a pair, by `(m, n, u, v) ↦ (u, v, m, n)`. So a solution of the
five-form system at the ratio `m : n` with parameter `t = u/v` is a solution at the ratio `u : v`
with parameter `m/n`.

* `P0_exchange`, `P8_exchange`, `PL_exchange`, `PU_exchange`: the pairs (`P0` and `PL` have their
  coordinates exchanged, `x + iy ↦ y + ix`; `P8` and `PU` are unchanged).
* `nrm_P0_exchange` and the three others: the norms are equal.
* `nrm_PL_add_nrm_PU`, `nrm_P0_add_nrm_P8`: the centre `2(m² + n²)(u² + v²)` is half of `L + U` and
  twice `F0 + F8`; it is symmetric by inspection.
* `fiveSquaresAt_exchange`: the five integer forms are squares together at `(m, n; u, v)` exactly
  when they are at `(u, v; m, n)`.
* `G1.fullSquares_exchange`: the same for the rational statement `G1.FullSquares` of the library:
  `FullSquares m n (u/v) ↔ FullSquares u v (m/n)` for positive `n, v`.

Written by an AI model (Claude); not read by a person. Ring identities only.
-/

namespace MagicSquareOfSquares.CoverDescent

theorem P0_exchange (m n u v : ℤ) : P0 u v m n = ((P0 m n u v).2, (P0 m n u v).1) := by
  ext <;> simp [P0, mul_comm]

theorem P8_exchange (m n u v : ℤ) : P8 u v m n = P8 m n u v := by
  ext <;> simp [P8, mul_comm]

theorem PL_exchange (m n u v : ℤ) : PL u v m n = ((PL m n u v).2, (PL m n u v).1) := by
  ext <;> simp [PL, mul_comm]

theorem PU_exchange (m n u v : ℤ) : PU u v m n = PU m n u v := by
  ext <;> simp [PU, mul_comm]

theorem nrm_P0_exchange (m n u v : ℤ) : nrm (P0 u v m n) = nrm (P0 m n u v) := by
  simp only [nrm, P0]; ring

theorem nrm_P8_exchange (m n u v : ℤ) : nrm (P8 u v m n) = nrm (P8 m n u v) := by
  simp only [nrm, P8]; ring

theorem nrm_PL_exchange (m n u v : ℤ) : nrm (PL u v m n) = nrm (PL m n u v) := by
  simp only [nrm, PL]; ring

theorem nrm_PU_exchange (m n u v : ℤ) : nrm (PU u v m n) = nrm (PU m n u v) := by
  simp only [nrm, PU]; ring

/-- `L + U` is twice the centre form `2(m² + n²)(u² + v²)`. -/
theorem nrm_PL_add_nrm_PU (m n u v : ℤ) :
    nrm (PL m n u v) + nrm (PU m n u v) = 4 * ((m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2)) := by
  simp only [nrm, PL, PU]; ring

/-- `F0 + F8` is half the centre form. -/
theorem nrm_P0_add_nrm_P8 (m n u v : ℤ) :
    nrm (P0 m n u v) + nrm (P8 m n u v) = (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) := by
  simp only [nrm, P0, P8]; ring

/-- The four corner norms and the centre form are squares of integers. -/
def FiveSquaresAt (m n u v : ℤ) : Prop :=
  IsSquare (nrm (P0 m n u v)) ∧ IsSquare (nrm (P8 m n u v)) ∧ IsSquare (nrm (PL m n u v)) ∧
    IsSquare (nrm (PU m n u v)) ∧ IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2))

/-- **The exchange symmetry**, for the integer forms. -/
theorem fiveSquaresAt_exchange (m n u v : ℤ) : FiveSquaresAt u v m n ↔ FiveSquaresAt m n u v := by
  unfold FiveSquaresAt
  rw [nrm_P0_exchange, nrm_P8_exchange, nrm_PL_exchange, nrm_PU_exchange,
    show 2 * (u ^ 2 + v ^ 2) * (m ^ 2 + n ^ 2) = 2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) by ring]

end MagicSquareOfSquares.CoverDescent

namespace MagicSquareOfSquares.G1

theorem isSquare_mul_sq_iff {a c : ℚ} (hc : c ≠ 0) : IsSquare (a * c ^ 2) ↔ IsSquare a := by
  constructor
  · rintro ⟨r, hr⟩
    refine ⟨r / c, ?_⟩
    field_simp
    linarith [hr]
  · rintro ⟨r, hr⟩
    exact ⟨r * c, by rw [hr]; ring⟩

theorem isSquare_congr_of_mul_sq {a b c d : ℚ} (hc : c ≠ 0) (hd : d ≠ 0)
    (h : a * c ^ 2 = b * d ^ 2) : IsSquare a ↔ IsSquare b := by
  rw [← isSquare_mul_sq_iff hc, h, isSquare_mul_sq_iff hd]

/-- **The exchange symmetry**, for the statement `FullSquares` of the library: the five cells of
the G1 family at the ratio `m : n` and the parameter `u/v` are squares exactly when those at the
ratio `u : v` and the parameter `m/n` are. -/
theorem fullSquares_exchange {m n u v : ℕ} (hn : n ≠ 0) (hv : v ≠ 0) :
    FullSquares m n ((u : ℚ) / v) ↔ FullSquares u v ((m : ℚ) / n) := by
  have hn' : (n : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hn
  have hv' : (v : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hv
  unfold FullSquares DiagonalSquares
  refine and_congr (and_congr ?_ (and_congr ?_ ?_)) (and_congr ?_ ?_) <;>
    refine isSquare_congr_of_mul_sq hv' hn' ?_ <;>
      first | (field_simp; ring) | field_simp

end MagicSquareOfSquares.G1
