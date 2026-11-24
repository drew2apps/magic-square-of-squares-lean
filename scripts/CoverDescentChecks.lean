import MagicSquareOfSquares.CoverDescentTable

/-!
# Checks on the table of `CoverDescentTable.lean`

Not part of the library: nothing here is counted among its theorems.  CI runs this file with
`lake env lean scripts/CoverDescentChecks.lean`; it fails if any check below fails.

1. The count "119 of the 425": the two tables have 119 pairs (115 with one class, 4 with four
   classes), without repetition, all of them among the 425 pairs of `G1.exceptions1000`; so the
   306 pairs of `CoverDescent.exceptions` are exactly the other ones, and the 119 are exactly the
   pairs of the 425 for which `m² + n²` is neither a square nor twice a square.
2. Which named ratios are in the table and which are not.
3. Forged certificates are refused by the checker `entryOk8` (the one the table uses: candidates
   are the products of 2 and of the primes `≡ ±1 (mod 8)`, by `exists_class_prod8`), and a forged
   list of masks is refused by `MasksOk`.  A certificate is not unique: a class may be moved to another line when that
   modulus excludes it as well; what is refused is a line that claims an exclusion which does not
   hold, a class that no line excludes, a wrong ratio, a missing prime, an unknown modulus.
4. A line of the second kind, `(p, 0, c)` (Euler's criterion modulo a prime `p ≡ 1 (mod 4)`), is
   refused when the prime does not exclude the class, when `p` is composite, 2 or `≡ 3 (mod 4)`.
5. Certificates with four classes (`entryOk4`): a line dropped, a line that cannot be read, another
   ratio, a missing prime, and a ratio that has a solution are refused.
6. The class from its structure alone (`structOk`): another class, a prime less, a prime that does
   not divide `(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)`, a composite number, a missing prime are refused.
7. The earlier checker `entryOk` (all products of the primes as candidates) is still in the
   library; a certificate written for one checker is refused by the other.
-/

open MagicSquareOfSquares MagicSquareOfSquares.CoverDescent

/-! ## 1. The count -/

example : ratios.length = 119 := by decide +kernel
example : table.length = 115 ∧ table4.length = 4 := by decide +kernel
example : ratios.Nodup := by decide +kernel
example : ∀ x ∈ ratios, x ∈ G1.exceptions1000 := by decide +kernel
example : G1.exceptions1000.length = 425 := by decide +kernel
example : G1.exceptions1000.Nodup := by decide +kernel
example : exceptions.length + ratios.length = G1.exceptions1000.length := by decide +kernel

/-! ## 2. Named ratios -/

example : (820, 799) ∈ ratios ∧ (889, 320) ∈ ratios := by decide +kernel
example : (627, 364) ∉ ratios ∧ (991, 263) ∉ ratios := by decide +kernel
/-- The six ratios that have a certificate only with the candidates of `exists_class_prod8`. -/
example : (427, 320) ∈ ratios ∧ (595, 452) ∈ ratios ∧ (747, 107) ∈ ratios ∧ (884, 485) ∈ ratios ∧
    (913, 892) ∈ ratios ∧ (984, 655) ∈ ratios := by decide +kernel
/-- Twelve more have a certificate with a line by Euler's criterion at 137, 157, 233 or 337. -/
example : (388, 287) ∈ ratios ∧ (457, 100) ∈ ratios ∧ (557, 357) ∈ ratios ∧ (564, 31) ∈ ratios ∧
    (576, 257) ∈ ratios ∧ (595, 533) ∈ ratios ∧ (644, 535) ∈ ratios ∧ (675, 101) ∈ ratios ∧
    (784, 391) ∈ ratios ∧ (788, 91) ∈ ratios ∧ (833, 319) ∈ ratios ∧ (879, 697) ∈ ratios := by
  decide +kernel
/-- The four ratios for which one class is not enough have a certificate with four classes. -/
example : table4.map (fun e => (e.1, e.2.1)) = [(308, 185), (340, 53), (393, 287), (493, 123)] := by
  decide +kernel
/-- The 119 are exactly the pairs of the 425 with `m² + n²` neither a square nor twice a square. -/
example : ∀ x ∈ G1.exceptions1000, (x ∈ ratios ↔
    ¬(Nat.sqrt (x.1 ^ 2 + x.2 ^ 2) ^ 2 = x.1 ^ 2 + x.2 ^ 2) ∧
    ¬(Nat.sqrt (2 * (x.1 ^ 2 + x.2 ^ 2)) ^ 2 = 2 * (x.1 ^ 2 + x.2 ^ 2))) := by decide +kernel
/-- No pair of the tables has `m² + n²` a square or twice a square (the two types with a
degenerate solution, for which no certificate can exist). -/
example : ∀ x ∈ ratios, ¬(Nat.sqrt (x.1 ^ 2 + x.2 ^ 2) ^ 2 = x.1 ^ 2 + x.2 ^ 2) ∧
    ¬(Nat.sqrt (2 * (x.1 ^ 2 + x.2 ^ 2)) ^ 2 = 2 * (x.1 ^ 2 + x.2 ^ 2)) := by decide +kernel

/-! ## 3. Forged certificates -/

/-- A genuine entry, copied from `chunk3`. The primes kept as candidates are `2, 7, 359, 751`. -/
def genuine : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ) :=
  (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xff), (13, 1, 0xc300), (17, 1, 0x2400), (5, 2, 0x1000), (53, 1, 0x800)])

example : genuine ∈ chunk3 := by decide +kernel
example : entryOk8 masks genuine = true := by decide +kernel
example : genuine.2.2.1.filter keep8 = [2, 7, 359, 751] := by decide +kernel

/-- One class moved from the line modulo 13 to the line modulo 8, which does not exclude it. -/
example : entryOk8 masks (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xff ||| (1 <<< 8)), (13, 1, 0xc300 ^^^ (1 <<< 8)), (17, 1, 0x2400), (5, 2, 0x1000), (53, 1, 0x800)]) = false := by
  decide +kernel
/-- One class dropped from its line: it is then excluded by no line. -/
example : entryOk8 masks (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xff), (13, 1, 0xc300), (17, 1, 0x2400), (5, 2, 0x1000), (53, 1, 0)]) = false := by
  decide +kernel
/-- The empty certificate. -/
example : entryOk8 masks (555, 196, [2, 3, 5, 7, 37, 359, 751], []) = false := by decide +kernel
/-- The same certificate for another ratio. -/
example : entryOk8 masks (555, 197, genuine.2.2.1, genuine.2.2.2) = false := by decide +kernel
example : entryOk8 masks (196, 555, genuine.2.2.1, genuine.2.2.2) = false := by decide +kernel
/-- A prime of `2mn(m − n)(m + n)` left out of the list: one that is a candidate (359), and one
that is not (37). The list must be complete in both cases. -/
example : entryOk8 masks (555, 196, [2, 3, 5, 7, 37, 751], genuine.2.2.2) = false := by
  decide +kernel
example : entryOk8 masks (555, 196, [2, 3, 5, 7, 359, 751], genuine.2.2.2) = false := by
  decide +kernel
/-- A composite number in the list of primes. -/
example : entryOk8 masks (555, 196, [2, 3, 5, 7, 37, 359, 751, 4], genuine.2.2.2) = false := by
  decide +kernel
/-- A modulus for which no mask was proved (7), an exponent 0, a wrong exponent for a mask. -/
example : entryOk8 masks (555, 196, genuine.2.2.1, (7, 1, 0) :: genuine.2.2.2) = false := by
  decide +kernel
example : entryOk8 masks (555, 196, genuine.2.2.1, (2, 0, 0) :: genuine.2.2.2) = false := by
  decide +kernel
/-- A ratio with a solution, `(4, 3)` (`4² + 3² = 5²`): claiming that every class is excluded
modulo 8 is refused, as it must be. The candidates are `1, 7, 2, 14`. -/
example : entryOk8 masks (4, 3, [2, 3, 7], [(2, 3, 0xf)]) = false := by decide +kernel
/-- A forged list of masks ("no residue is a Gaussian square modulo 8") would make the checker
accept that claim; the theorem `no_fullSquares_of_entry` asks for `MasksOk`, which fails for it. -/
example : entryOk8 [(8, S8, 0)] (4, 3, [2, 3, 7], [(2, 3, 0xf)]) = true := by decide +kernel
example : ¬MasksOk [(8, S8, 0)] := by
  intro h
  have h2 := (h (8, S8, 0) (by simp)).2
  rw [← mask8.2] at h2
  revert h2
  decide +kernel

/-! ## 4. Lines by Euler's criterion -/

/-- A genuine entry with such a line, copied from `chunk8`: the last class is excluded modulo 233. -/
def genuineE : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ) :=
  (784, 391, [2, 3, 5, 7, 17, 23, 47, 131], [(2, 3, 0xffff), (13, 1, 0xaa550000), (29, 1, 0x5a00000), (37, 1, 0x40020000), (41, 1, 0x80000), (233, 0, 0x10000000)])

example : genuineE ∈ chunk8 := by decide +kernel
example : entryOk8 masks genuineE = true := by decide +kernel
/-- The same with another prime, which does not exclude that class (229), with a composite
number that has a square root of `−1` (221 = 13·17), with the prime 2, with a prime `≡ 3 (mod 4)`
(239), and with the class dropped from the line. -/
example : entryOk8 masks (784, 391, [2, 3, 5, 7, 17, 23, 47, 131], [(2, 3, 0xffff), (13, 1, 0xaa550000), (29, 1, 0x5a00000), (37, 1, 0x40020000), (41, 1, 0x80000), (229, 0, 0x10000000)]) = false := by decide +kernel
example : entryOk8 masks (784, 391, [2, 3, 5, 7, 17, 23, 47, 131], [(2, 3, 0xffff), (13, 1, 0xaa550000), (29, 1, 0x5a00000), (37, 1, 0x40020000), (41, 1, 0x80000), (221, 0, 0x10000000)]) = false := by decide +kernel
example : entryOk8 masks (784, 391, [2, 3, 5, 7, 17, 23, 47, 131], [(2, 3, 0xffff), (13, 1, 0xaa550000), (29, 1, 0x5a00000), (37, 1, 0x40020000), (41, 1, 0x80000), (2, 0, 0x10000000)]) = false := by decide +kernel
example : entryOk8 masks (784, 391, [2, 3, 5, 7, 17, 23, 47, 131], [(2, 3, 0xffff), (13, 1, 0xaa550000), (29, 1, 0x5a00000), (37, 1, 0x40020000), (41, 1, 0x80000), (239, 0, 0x10000000)]) = false := by decide +kernel
example : entryOk8 masks (784, 391, [2, 3, 5, 7, 17, 23, 47, 131], [(2, 3, 0xffff), (13, 1, 0xaa550000), (29, 1, 0x5a00000), (37, 1, 0x40020000), (41, 1, 0x80000), (233, 0, 0)]) = false := by decide +kernel

/-! ## 5. Certificates with four classes -/

/-- A genuine entry of `table4`. -/
def genuine4 : ℕ × ℕ × List ℕ × List (ℕ × ℕ) :=
  (308, 185, [2, 3, 5, 7, 11, 17, 29, 37, 41], [(2, 3), (5, 1), (13, 1), (29, 1)])

example : genuine4 ∈ table4 := by decide +kernel
example : entryOk4 masks genuine4 = true := by decide +kernel
/-- With one line less tuples are left; with a line that cannot be read (a modulus without a mask,
an exponent 0, a composite number) the same; the certificate of another ratio; a prime left out. -/
example : entryOk4 masks (308, 185, genuine4.2.2.1, [(2, 3), (5, 1), (13, 1)]) = false := by
  decide +kernel
example : entryOk4 masks (308, 185, genuine4.2.2.1, [(2, 3), (5, 1), (13, 1), (7, 1)]) = false := by
  decide +kernel
example : entryOk4 masks (308, 185, genuine4.2.2.1, [(2, 3), (5, 1), (13, 1), (29, 0)]) = false := by
  decide +kernel
example : entryOk4 masks (308, 187, genuine4.2.2.1, genuine4.2.2.2) = false := by decide +kernel
example : entryOk4 masks (308, 185, [2, 3, 5, 7, 11, 17, 29, 37], genuine4.2.2.2) = false := by
  decide +kernel
/-- A ratio with a solution, `(4, 3)`: with all the lines used in `table4` a tuple is left. -/
example : entryOk4 masks (4, 3, [2, 3, 7], [(2, 3), (5, 1), (13, 1), (17, 1), (29, 1)]) = false := by
  decide +kernel

/-! ## 6. The class from its structure alone -/

example : structTable.map (fun e => (e.1, e.2.1, e.2.2.2.1)) =
    [(627, 364, 521266), (991, 263, 521266)] := by decide +kernel
example : structOk (627, 364, [2, 3, 7, 11, 13, 19, 263, 991], 521266, [5, 29, 37]) = true := by
  decide +kernel
/-- Another class; one prime less (then two candidates are left); a prime that divides neither
`m² + n²` nor the octic (41); a composite number; a prime of `2mn(m² − n²)` left out. -/
example : structOk (627, 364, [2, 3, 7, 11, 13, 19, 263, 991], 2, [5, 29, 37]) = false := by
  decide +kernel
example : structOk (627, 364, [2, 3, 7, 11, 13, 19, 263, 991], 521266, [5, 29]) = false := by
  decide +kernel
example : structOk (627, 364, [2, 3, 7, 11, 13, 19, 263, 991], 521266, [5, 29, 37, 41]) = false := by
  decide +kernel
example : structOk (627, 364, [2, 3, 7, 11, 13, 19, 263, 991], 521266, [5, 29, 37, 25]) = false := by
  decide +kernel
example : structOk (627, 364, [2, 3, 7, 11, 13, 19, 263], 521266, [5, 29, 37]) = false := by
  decide +kernel

/-! ## 7. The earlier checker -/

/-- A certificate over all products of the primes (an entry of the table of the previous
version, which `entryOk` accepted with the masks of that table). -/
def genuineAll : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ) :=
  (555, 196, [2, 3, 5, 7, 37, 359, 751], [(2, 3, 0xffffffffffffffff), (3, 2, 0xffffffff000000000000000000000000), (13, 1, 0xc33c3cc30000000000000000), (17, 1, 0x244242240000000000000000), (5, 2, 0x188101100000000000000000), (41, 1, 0x80000000000000000000), (53, 1, 0x80000000000000000)])

/-- Read with the other set of candidates, each certificate is refused. (`genuineAll` also uses
the modulus 9, for which the present table has no mask.) -/
example : entryOk8 masks genuineAll = false := by decide +kernel
example : entryOk masks genuine = false := by decide +kernel
