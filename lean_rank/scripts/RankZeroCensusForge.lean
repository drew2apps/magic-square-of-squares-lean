import LeanRank.RankZeroCensus

/-!
# Negative controls for the generic certificate checker (not part of the library)

Written by an AI model; not read by a person.  Run with
`lake env lean scripts/RankZeroCensusForge.lean`; every `example` must compile.  Refused:

1. the certificate of `(4, 3)` with the good class `(6, 70)` removed (it is the class of a point of
   order 4, so no test can exclude it);
2. the same certificate without the table modulo `3`, or modulo `5`, or modulo `16`, or with no
   table at all (the table modulo `4` is redundant: removing it alone is accepted);
3. a certificate whose prime list `S₂` omits `7`: `Cert.check` accepts it, but the coverage
   hypothesis of `finrank_eq_zero_of` cannot be discharged, since `7 ∣ e₂(e₂ − e₃) = 256·175`;
4. a certificate for `(15, 8)` (rank 1 by PARI) with four claimed classes.
-/

open LeanRank.RankZeroGen LeanRank.RankZeroCensus

set_option maxHeartbeats 400000

-- 1. a good class removed
example : ({ c_4_3 with good := [(1, 1), (-1, 7), (-6, 10)] } : Cert).check = false := by
  decide +kernel

-- 2. tables removed
example : ({ c_4_3 with tables := [(2, 4, sq4), (5, 5, sq5), (2, 16, sq16)] } : Cert).check
    = false := by decide +kernel
example : ({ c_4_3 with tables := [(2, 4, sq4), (3, 3, sq3), (2, 16, sq16)] } : Cert).check
    = false := by decide +kernel
example : ({ c_4_3 with tables := [(2, 4, sq4), (3, 3, sq3), (5, 5, sq5)] } : Cert).check
    = false := by decide +kernel
example : ({ c_4_3 with tables := [] } : Cert).check = false := by decide +kernel
example : ({ c_4_3 with tables := [(3, 3, sq3), (5, 5, sq5), (2, 16, sq16)] } : Cert).check
    = true := by decide +kernel

-- 3. a prime missing from `S₂`: the check passes, the coverage fails
example : ({ c_4_3 with S₂ := [2, 5] } : Cert).check = true := by decide +kernel
example : ¬ (∀ q ∈ ([2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 7] : List ℕ), q ∈ ([2, 5] : List ℕ)) := by
  decide +kernel
example : ((-256 : ℤ) * (-256 - -81)).natAbs = ([2, 2, 2, 2, 2, 2, 2, 2, 5, 5, 7] : List ℕ).prod := by
  decide +kernel

-- 4. (15, 8): Y² = X(X + 50625)(X + 4096), rank 1
def c_15_8_forged : Cert where
  e₂ := -50625
  e₃ := -4096
  S₁ := [2, 3, 5]
  S₂ := [3, 5, 7, 17, 23]
  good := [(1, 1), (-1, 7 * 17 * 23), (-1, 1), (1, 7 * 17 * 23)]
  tables := [(2, 4, sq4), (3, 3, sq3), (5, 5, sq5), (2, 16, sq16)]

example : c_15_8_forged.check = false := by decide +kernel
