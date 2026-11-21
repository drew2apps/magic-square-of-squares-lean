import Mathlib

/-!
# 3×3 magic squares: definitions and linear structure

A 3×3 array of integers

```
A B C
D E F
G H I
```

is *magic* with sum `T` if its three rows, three columns and two diagonals all sum to `T`.
This file defines `IsMagic`, the squared version `IsMSS`, and `Primitive`, and proves the
facts that hold for every 3×3 magic square (no squares needed).

## Main results

* `IsMagic.sum_eq`: the magic sum is three times the centre.
* `IsMagic.opposite`: two cells opposite each other through the centre sum to twice the centre.
* `IsMagic.corner`: twice a corner is the sum of the two edge cells that do not touch it.
* `IsMagic.iff_param`: the parametrisation of all 3×3 magic squares by the centre and two steps.
* `IsMagic.transpose`, `IsMagic.rotate`: symmetries.
-/

namespace MagicSquareOfSquares

/-- `IsMagic A B C D E F G H I T`: the 3×3 integer array with rows `A B C`, `D E F`, `G H I`
has all three rows, all three columns and both diagonals summing to `T`. -/
structure IsMagic (A B C D E F G H I T : ℤ) : Prop where
  row₁ : A + B + C = T
  row₂ : D + E + F = T
  row₃ : G + H + I = T
  col₁ : A + D + G = T
  col₂ : B + E + H = T
  col₃ : C + F + I = T
  diag : A + E + I = T
  anti : C + E + G = T

/-- `IsMSS a b c d e f g h i T`: the squares of the nine integers, arranged in rows
`a² b² c²`, `d² e² f²`, `g² h² i²`, form a magic square with sum `T`.
The nine integers are called the *roots*. No distinctness is assumed. -/
def IsMSS (a b c d e f g h i T : ℤ) : Prop :=
  IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T

/-- The nine roots are *primitive* if no prime divides all of them. This is the same as the
greatest common divisor of the nine entries being `1`. -/
def Primitive (a b c d e f g h i : ℤ) : Prop :=
  ∀ p : ℕ, p.Prime → ¬((p : ℤ) ∣ a ∧ (p : ℤ) ∣ b ∧ (p : ℤ) ∣ c ∧ (p : ℤ) ∣ d ∧ (p : ℤ) ∣ e ∧
    (p : ℤ) ∣ f ∧ (p : ℤ) ∣ g ∧ (p : ℤ) ∣ h ∧ (p : ℤ) ∣ i)

namespace IsMagic

variable {A B C D E F G H I T : ℤ}

/-- The magic sum is three times the centre. -/
theorem sum_eq (h : IsMagic A B C D E F G H I T) : T = 3 * E := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := h
  linarith

/-- Two cells opposite each other through the centre sum to twice the centre. -/
theorem opposite (h : IsMagic A B C D E F G H I T) :
    A + I = 2 * E ∧ C + G = 2 * E ∧ B + H = 2 * E ∧ D + F = 2 * E := by
  have hT := h.sum_eq
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := h
  refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith

/-- Twice a corner is the sum of the two edge cells that do not touch it. -/
theorem corner (h : IsMagic A B C D E F G H I T) :
    2 * A = F + H ∧ 2 * C = D + H ∧ 2 * G = B + F ∧ 2 * I = B + D := by
  have hT := h.sum_eq
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := h
  refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith

/-- Every 3×3 magic square is determined by its centre `E` and two steps `u`, `v`. -/
theorem iff_param :
    IsMagic A B C D E F G H I T ↔
      ∃ u v : ℤ, A = E + u ∧ B = E - u - v ∧ C = E + v ∧ D = E - u + v ∧ F = E + u - v ∧
        G = E - v ∧ H = E + u + v ∧ I = E - u ∧ T = 3 * E := by
  constructor
  · intro h
    have hT := h.sum_eq
    obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := h
    refine ⟨A - E, C - E, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩ <;> linarith
  · rintro ⟨u, v, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl, rfl⟩
    constructor <;> ring

/-- Reflecting a magic square in its main diagonal gives a magic square. -/
theorem transpose (h : IsMagic A B C D E F G H I T) : IsMagic A D G B E H C F I T := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := h
  constructor <;> linarith

/-- Rotating a magic square by a quarter turn gives a magic square. -/
theorem rotate (h : IsMagic A B C D E F G H I T) : IsMagic G D A H E B I F C T := by
  obtain ⟨h1, h2, h3, h4, h5, h6, h7, h8⟩ := h
  constructor <;> linarith

end IsMagic

end MagicSquareOfSquares
