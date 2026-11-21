import MagicSquareOfSquares.Defs

/-!
# Known examples and near misses

These are sanity checks of the definitions against the published examples.

* `bremner_isMagic`, `bremner_not_isSquare`: Bremner's square is magic and has seven square
  entries out of nine.
* `sallows`: in Sallows' square of squares the rows, columns and one diagonal sum to `147²`,
  and the other diagonal does not.
* `parker`: in the Parker square the rows, columns and one diagonal sum to `3051`, the other
  diagonal does not, and two entries repeat.
* `exists_magic_square_squares_four`: Euler's 4×4 magic square of distinct squares (1770).
-/

namespace MagicSquareOfSquares

/-- Bremner's square
```
373²    289²   565²
360721  425²   23²
205²    527²   222121
```
is magic with sum `541875`. -/
theorem bremner_isMagic :
    IsMagic (373 ^ 2) (289 ^ 2) (565 ^ 2) 360721 (425 ^ 2) (23 ^ 2) (205 ^ 2) (527 ^ 2) 222121
      541875 := by
  constructor <;> norm_num

/-- The two remaining entries of Bremner's square are not squares. -/
theorem bremner_not_isSquare : ¬IsSquare (360721 : ℕ) ∧ ¬IsSquare (222121 : ℕ) := by
  constructor
  · rintro ⟨r, hr⟩
    have h1 : r < 601 := by nlinarith
    have h2 : 600 < r := by nlinarith
    omega
  · rintro ⟨r, hr⟩
    have h1 : r < 472 := by nlinarith
    have h2 : 471 < r := by nlinarith
    omega

/-- Sallows' square
```
127²  46²   58²
2²    113²  94²
74²   82²   97²
```
has rows, columns and one diagonal summing to `147²`. The other diagonal does not. -/
theorem sallows :
    127 ^ 2 + 46 ^ 2 + 58 ^ 2 = 147 ^ 2 ∧ 2 ^ 2 + 113 ^ 2 + 94 ^ 2 = 147 ^ 2 ∧
    74 ^ 2 + 82 ^ 2 + 97 ^ 2 = 147 ^ 2 ∧ 127 ^ 2 + 2 ^ 2 + 74 ^ 2 = 147 ^ 2 ∧
    46 ^ 2 + 113 ^ 2 + 82 ^ 2 = 147 ^ 2 ∧ 58 ^ 2 + 94 ^ 2 + 97 ^ 2 = 147 ^ 2 ∧
    58 ^ 2 + 113 ^ 2 + 74 ^ 2 = 147 ^ 2 ∧ 127 ^ 2 + 113 ^ 2 + 97 ^ 2 ≠ 147 ^ 2 := by
  norm_num

/-- The Parker square
```
29²  1²   47²
41²  37²  1²
23²  41²  29²
```
has rows, columns and one diagonal summing to `3051`. The other diagonal sums to `4107`. -/
theorem parker :
    29 ^ 2 + 1 ^ 2 + 47 ^ 2 = 3051 ∧ 41 ^ 2 + 37 ^ 2 + 1 ^ 2 = 3051 ∧
    23 ^ 2 + 41 ^ 2 + 29 ^ 2 = 3051 ∧ 29 ^ 2 + 41 ^ 2 + 23 ^ 2 = 3051 ∧
    1 ^ 2 + 37 ^ 2 + 41 ^ 2 = 3051 ∧ 47 ^ 2 + 1 ^ 2 + 29 ^ 2 = 3051 ∧
    29 ^ 2 + 37 ^ 2 + 29 ^ 2 = 3051 ∧ 47 ^ 2 + 37 ^ 2 + 23 ^ 2 = 4107 := by
  norm_num

/-- The roots of Euler's 4×4 magic square of squares. -/
def euler : Fin 4 → Fin 4 → ℕ :=
  ![![68, 29, 41, 37], ![17, 31, 79, 32], ![59, 28, 23, 61], ![11, 77, 8, 49]]

/-- A 4×4 magic square of distinct positive squares exists (Euler, 1770). The statement has
the same shape as the open 3×3 statement in `formal-conjectures`. -/
theorem exists_magic_square_squares_four :
    ∃ m : Fin 4 → Fin 4 → ℕ, ∃ t : ℕ,
      m.Injective2 ∧
      (∀ i j, 0 < m i j ∧ IsSquare (m i j)) ∧
      (∀ i, ∑ j, m i j = t) ∧
      (∀ j, ∑ i, m i j = t) ∧
      ∑ i, m i i = t ∧
      ∑ i, m i i.rev = t := by
  refine ⟨fun i j => euler i j ^ 2, 8515, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · unfold Function.Injective2
    decide
  · intro i j
    have hpos : ∀ i j, 0 < euler i j := by decide
    exact ⟨pow_pos (hpos i j) 2, euler i j, sq _⟩
  · decide
  · decide
  · decide
  · decide

end MagicSquareOfSquares
