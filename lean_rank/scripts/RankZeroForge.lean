import LeanRank.RankZero

/-!
# Negative controls for `LeanRank.RankZero` (not part of the library)

Written by an AI model; not read by a person.  Run with `lake env lean scripts/RankZeroForge.lean`;
every `example` must compile.  What is refused, precisely:

1. the 2-adic checker does not exclude the four classes of the torsion points (so it is not a
   checker that excludes everything);
2. a forged certificate list containing a good class is refused (`all` returns `false`);
3. the depth matters: modulo `32` the class `(1, −3)` is not excluded, modulo `64` it is;
4. a forged table of squares modulo `64` (with `17` removed) fails the table check that
   `sq64_mem` rests on;
5. without the primitivity condition `¬ (2 ∣ n ∧ 2 ∣ e)` the test excludes nothing (the residues
   `a ≡ e ≡ 0` pass every class), so the condition is what carries the 2-adic information.
-/

open LeanRank.RankZero

-- 1. the four good classes pass
example : excluded 1 1 = false := by decide +kernel
example : excluded 3 (-3) = false := by decide +kernel
example : excluded 2 (-1) = false := by decide +kernel
example : excluded 6 3 = false := by decide +kernel

-- 2. a forged certificate list with a good class added is refused
example : ((1, 1) :: badClasses).all (fun c => excluded c.1 c.2) = false := by decide +kernel
example : ((6, 3) :: badClasses).all (fun c => excluded c.1 c.2) = false := by decide +kernel

-- 3. modulo 32 the class (1, -3) survives; modulo 64 it does not
def sq32 : List ℤ := [0, 1, 4, 9, 16, 17, 25]

def pass32 (d₁ d₂ s q : ℤ) : Bool :=
  ((d₁ * s) % 2 == 1 || q % 2 == 1) &&
    sq32.any (fun r => (d₁ * s - 12 * q - d₂ * r) % 32 == 0) &&
    sq32.any (fun r => (d₁ * d₂ * (d₁ * s - 16 * q) - r) % 32 == 0)

def excluded32 (d₁ d₂ : ℤ) : Bool := !(sq32.any fun s => sq32.any fun q => pass32 d₁ d₂ s q)

example : ∀ k : ℕ, k < 32 → ((k : ℤ) * k) % 32 ∈ sq32 := by decide +kernel
example : excluded32 1 (-3) = false := by decide +kernel
example : excluded 1 (-3) = true := by decide +kernel

-- 4. a forged table of squares (17 removed) fails the table check
example : ¬ (∀ k : ℕ, k < 64 →
    ((k : ℤ) * k) % 64 ∈ ([0, 1, 4, 9, 16, 25, 33, 36, 41, 49, 57] : List ℤ)) := by
  decide +kernel

-- 5. without the primitivity condition nothing is excluded
def passNoPrim (d₁ d₂ s q : ℤ) : Bool :=
  sq64.any (fun r => (d₁ * s - 12 * q - d₂ * r) % 64 == 0) &&
    sq64.any (fun r => (d₁ * d₂ * (d₁ * s - 16 * q) - r) % 64 == 0)

example : badClasses.all (fun c => sq64.any fun s => sq64.any fun q => passNoPrim c.1 c.2 s q)
    = true := by decide +kernel
