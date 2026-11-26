import MagicSquareOfSquares.Defs
import LeanRank.CensusClosure

/-!
# The six census closures of `CensusClosure.lean`, stated with the main library's `IsMSS`

Written by an AI model; not read by a person.

`CensusClosure.lean` states its exclusions with the eight line equations written out (`r₁ … ad`), because
`lean_rank/` did not compile `MagicSquareOfSquares/Defs.lean`. `Defs.lean` imports only Mathlib and builds
unchanged on Lean v4.36.0 / Mathlib v4.36.0, so it is now a root of the `MagicSquareOfSquares` library of
`lean_rank/lakefile.toml`, and this file states the six exclusions for the main library's own
`MagicSquareOfSquares.IsMSS a b c d e f g h i T` (the squares `a² … i²` in rows `a² b² c² / d² e² f² / g² h² i²`
form a magic square with sum `T`), with the nine entries distinct in the form the main library uses for the open
statement (`openStatement_iff_exists_isMSS`, `Statement.lean`):
`[a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup`.

"Edge ratio `m : n`" is in the sense of `CensusClosure.lean`: `(b + d)·n = (h + f)·m`.

The only new content is the bridge: the fields of `IsMagic` are the eight hypotheses `r₁ … ad`, and distinct
entries give the two hypotheses `b² ≠ d²`, `b² ≠ f²` of `ratio_ne_m_n` (`bd_bf_of_nodup`). The exclusions
themselves are `LeanRank.CensusClosure.ratio_ne_4_3`, `ratio_ne_7_1`, `ratio_ne_12_5`, `ratio_ne_24_7`,
`ratio_ne_17_7`, `ratio_ne_31_17`, unchanged. `control_needs_distinct`: without distinctness the statement is
false (a square of `±1`'s has `b + d = h + f = 0`).
-/

namespace LeanRank.CensusClosureMSS

open MagicSquareOfSquares

/-- Distinct entries give the two hypotheses of `CensusClosure.ratio_ne_m_n`. -/
theorem bd_bf_of_nodup {a b c d e f g h i : ℤ}
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    b ^ 2 ≠ d ^ 2 ∧ b ^ 2 ≠ f ^ 2 := by
  simp only [List.nodup_cons, List.mem_cons, List.not_mem_nil, not_or, or_false,
    List.nodup_nil, and_true] at hd
  obtain ⟨-, ⟨-, hbd, -, hbf, -⟩, -⟩ := hd
  exact ⟨hbd, hbf⟩

/-- The eight fields of `IsMSS`, as the hypotheses `r₁ … ad` of `CensusClosure.lean`. -/
theorem lines {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T) :
    a ^ 2 + b ^ 2 + c ^ 2 = T ∧ d ^ 2 + e ^ 2 + f ^ 2 = T ∧ g ^ 2 + h ^ 2 + i ^ 2 = T ∧
      a ^ 2 + d ^ 2 + g ^ 2 = T ∧ b ^ 2 + e ^ 2 + h ^ 2 = T ∧ c ^ 2 + f ^ 2 + i ^ 2 = T ∧
      a ^ 2 + e ^ 2 + i ^ 2 = T ∧ c ^ 2 + e ^ 2 + g ^ 2 = T := by
  unfold IsMSS at hs
  exact ⟨hs.row₁, hs.row₂, hs.row₃, hs.col₁, hs.col₂, hs.col₃, hs.diag, hs.anti⟩

/-- **No magic square of distinct squares (main library's `IsMSS`) has edge ratio `4 : 3`.** -/
theorem no_mss_4_3 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    (b + d) * 3 ≠ (h + f) * 4 := by
  obtain ⟨r₁, r₂, r₃, c₁, c₂, c₃, dg, ad⟩ := lines hs
  obtain ⟨hbd, hbf⟩ := bd_bf_of_nodup hd
  exact CensusClosure.ratio_ne_4_3 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf

/-- **… nor `7 : 1`.** -/
theorem no_mss_7_1 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    (b + d) * 1 ≠ (h + f) * 7 := by
  obtain ⟨r₁, r₂, r₃, c₁, c₂, c₃, dg, ad⟩ := lines hs
  obtain ⟨hbd, hbf⟩ := bd_bf_of_nodup hd
  exact CensusClosure.ratio_ne_7_1 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf

/-- **… nor `12 : 5`.** -/
theorem no_mss_12_5 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    (b + d) * 5 ≠ (h + f) * 12 := by
  obtain ⟨r₁, r₂, r₃, c₁, c₂, c₃, dg, ad⟩ := lines hs
  obtain ⟨hbd, hbf⟩ := bd_bf_of_nodup hd
  exact CensusClosure.ratio_ne_12_5 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf

/-- **… nor `24 : 7`.** -/
theorem no_mss_24_7 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    (b + d) * 7 ≠ (h + f) * 24 := by
  obtain ⟨r₁, r₂, r₃, c₁, c₂, c₃, dg, ad⟩ := lines hs
  obtain ⟨hbd, hbf⟩ := bd_bf_of_nodup hd
  exact CensusClosure.ratio_ne_24_7 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf

/-- **… nor `17 : 7`** (the partner of `12 : 5`). -/
theorem no_mss_17_7 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    (b + d) * 7 ≠ (h + f) * 17 := by
  obtain ⟨r₁, r₂, r₃, c₁, c₂, c₃, dg, ad⟩ := lines hs
  obtain ⟨hbd, hbf⟩ := bd_bf_of_nodup hd
  exact CensusClosure.ratio_ne_17_7 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf

/-- **… nor `31 : 17`** (the partner of `24 : 7`). -/
theorem no_mss_31_17 {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hd : [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup) :
    (b + d) * 17 ≠ (h + f) * 31 := by
  obtain ⟨r₁, r₂, r₃, c₁, c₂, c₃, dg, ad⟩ := lines hs
  obtain ⟨hbd, hbf⟩ := bd_bf_of_nodup hd
  exact CensusClosure.ratio_ne_31_17 r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf

/-- The six ratios of `CensusClosure.lean`. -/
def censusClosureRatios : List (ℤ × ℤ) := [(4, 3), (7, 1), (12, 5), (24, 7), (17, 7), (31, 17)]

/-- **The six together**: no `IsMSS a b c d e f g h i T` with nine distinct entries has edge ratio
`(b + d) : (h + f) = m : n` for `(m, n)` one of `(4,3), (7,1), (12,5), (24,7), (17,7), (31,17)`. -/
theorem no_mss_census : ∀ p ∈ censusClosureRatios, ∀ a b c d e f g h i T : ℤ,
    IsMSS a b c d e f g h i T → [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup →
      (b + d) * p.2 ≠ (h + f) * p.1 := by
  intro p hp a b c d e f g h i T hs hd
  simp only [censusClosureRatios, List.mem_cons, List.not_mem_nil, or_false] at hp
  rcases hp with rfl | rfl | rfl | rfl | rfl | rfl
  · exact no_mss_4_3 hs hd
  · exact no_mss_7_1 hs hd
  · exact no_mss_12_5 hs hd
  · exact no_mss_24_7 hs hd
  · exact no_mss_17_7 hs hd
  · exact no_mss_31_17 hs hd

/-- **Control: distinctness is needed.** The square of `±1`'s with `b = h = 1`, `d = f = −1` is an `IsMSS`
(sum 3) whose edge sums are `b + d = h + f = 0`, so `(b + d)·3 = (h + f)·4`. -/
theorem control_needs_distinct :
    IsMSS 1 1 1 (-1) 1 (-1) 1 1 1 3 ∧ ((1 : ℤ) + -1) * 3 = (1 + -1) * 4 := by
  refine ⟨?_, by norm_num⟩
  unfold IsMSS
  constructor <;> norm_num

/-- … and so the statement of `no_mss_4_3` without the distinctness hypothesis is false. -/
theorem control_distinct_not_vacuous :
    ¬ ∀ a b c d e f g h i T : ℤ, IsMSS a b c d e f g h i T → (b + d) * 3 ≠ (h + f) * 4 := by
  intro H
  exact H 1 1 1 (-1) 1 (-1) 1 1 1 3 control_needs_distinct.1 control_needs_distinct.2

end LeanRank.CensusClosureMSS

#print axioms LeanRank.CensusClosureMSS.bd_bf_of_nodup
#print axioms LeanRank.CensusClosureMSS.no_mss_4_3
#print axioms LeanRank.CensusClosureMSS.no_mss_7_1
#print axioms LeanRank.CensusClosureMSS.no_mss_12_5
#print axioms LeanRank.CensusClosureMSS.no_mss_24_7
#print axioms LeanRank.CensusClosureMSS.no_mss_17_7
#print axioms LeanRank.CensusClosureMSS.no_mss_31_17
#print axioms LeanRank.CensusClosureMSS.no_mss_census
#print axioms LeanRank.CensusClosureMSS.control_distinct_not_vacuous
