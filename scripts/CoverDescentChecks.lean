import MagicSquareOfSquares.CoverDescentTable

/-!
# Checks on the table of `CoverDescentTable.lean`

Not part of the library: nothing here is counted among its theorems.  CI runs this file with
`lake env lean scripts/CoverDescentChecks.lean`; it fails if any check below fails.

1. The count "97 of the 425": the table has 97 pairs, without repetition, all of them among the
   425 pairs of `G1.exceptions1000`; so the 328 pairs of `CoverDescent.exceptions` are exactly the
   other ones.
2. Which named ratios are in the table and which are not.
3. Forged certificates are refused by the checker `entryOk`, and a forged list of masks is refused
   by `MasksOk`.  A certificate is not unique: a class may be moved to another line when that
   modulus excludes it as well; what is refused is a line that claims an exclusion which does not
   hold, a class that no line excludes, a wrong ratio, a missing prime, an unknown modulus.
-/

open MagicSquareOfSquares MagicSquareOfSquares.CoverDescent

/-! ## 1. The count -/

example : ratios.length = 97 := by decide +kernel
example : ratios.Nodup := by decide +kernel
example : ∀ x ∈ ratios, x ∈ G1.exceptions1000 := by decide +kernel
example : G1.exceptions1000.length = 425 := by decide +kernel
example : G1.exceptions1000.Nodup := by decide +kernel
example : exceptions.length + ratios.length = G1.exceptions1000.length := by decide +kernel

/-! ## 2. Named ratios -/

example : (820, 799) ∈ ratios ∧ (889, 320) ∈ ratios := by decide +kernel
example : (627, 364) ∉ ratios ∧ (991, 263) ∉ ratios := by decide +kernel
/-- No pair of the table has `m² + n²` a square or twice a square (the two types with a
degenerate solution, for which no certificate can exist). -/
example : ∀ x ∈ ratios, ¬(Nat.sqrt (x.1 ^ 2 + x.2 ^ 2) ^ 2 = x.1 ^ 2 + x.2 ^ 2) ∧
    ¬(Nat.sqrt (2 * (x.1 ^ 2 + x.2 ^ 2)) ^ 2 = 2 * (x.1 ^ 2 + x.2 ^ 2)) := by decide +kernel

/-! ## 3. Forged certificates -/

/-- A genuine entry, copied from `chunk3`. -/
def genuine : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ) :=
  (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xffffffffffffffff), (3, 2, 0xffffffff000000000000000000000000), (13, 1, 0xc33c3cc30000000000000000), (17, 1, 0x244242240000000000000000), (5, 2, 0x188101100000000000000000), (41, 1, 0x80000000000000000000), (53, 1, 0x80000000000000000)])

example : genuine ∈ chunk3 := by decide +kernel
example : entryOk masks genuine = true := by decide +kernel

/-- One class moved from the line modulo 13 to the line modulo 8, which does not exclude it. -/
example : entryOk masks (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xffffffffffffffff ||| (1 <<< 64)), (3, 2, 0xffffffff000000000000000000000000), (13, 1, 0xc33c3cc30000000000000000 ^^^ (1 <<< 64)), (17, 1, 0x244242240000000000000000), (5, 2, 0x188101100000000000000000), (41, 1, 0x80000000000000000000), (53, 1, 0x80000000000000000)]) = false := by
  decide +kernel
/-- One class dropped from its line: it is then excluded by no line. -/
example : entryOk masks (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xffffffffffffffff), (3, 2, 0xffffffff000000000000000000000000), (13, 1, 0xc33c3cc30000000000000000), (17, 1, 0x244242240000000000000000), (5, 2, 0x188101100000000000000000), (41, 1, 0x80000000000000000000), (53, 1, 0)]) = false := by
  decide +kernel
/-- The empty certificate. -/
example : entryOk masks (555, 196, [2, 3, 5, 7, 37, 359, 751], []) = false := by decide +kernel
/-- The same certificate for another ratio. -/
example : entryOk masks (555, 197, genuine.2.2.1, genuine.2.2.2) = false := by decide +kernel
example : entryOk masks (196, 555, genuine.2.2.1, genuine.2.2.2) = false := by decide +kernel
/-- A prime of `2mn(m − n)(m + n)` left out of the list. -/
example : entryOk masks (555, 196, [2, 3, 5, 7, 37, 359], genuine.2.2.2) = false := by decide +kernel
/-- A composite number in the list of primes. -/
example : entryOk masks (555, 196, [2, 3, 5, 7, 37, 359, 751, 4], genuine.2.2.2) = false := by
  decide +kernel
/-- A modulus for which no mask was proved (7), an exponent 0, a wrong exponent for a mask. -/
example : entryOk masks (555, 196, genuine.2.2.1, (7, 1, 0) :: genuine.2.2.2) = false := by
  decide +kernel
example : entryOk masks (555, 196, genuine.2.2.1, (2, 0, 0) :: genuine.2.2.2) = false := by
  decide +kernel
/-- A ratio with a solution, `(4, 3)` (`4² + 3² = 5²`): claiming that every class is excluded
modulo 8 is refused, as it must be. -/
example : entryOk masks (4, 3, [2, 3, 7], [(2, 3, 0xff)]) = false := by decide +kernel
/-- A forged list of masks ("no residue is a Gaussian square modulo 8") would make the checker
accept that claim; the theorem `no_fullSquares_of_entry` asks for `MasksOk`, which fails for it. -/
example : entryOk [(8, S8, 0)] (4, 3, [2, 3, 7], [(2, 3, 0xff)]) = true := by decide +kernel
example : ¬MasksOk [(8, S8, 0)] := by
  intro h
  have h2 := (h (8, S8, 0) (by simp)).2
  rw [← mask8.2] at h2
  revert h2
  decide +kernel
