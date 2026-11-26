import LeanRank.CensusClosure

/-!
# Negative test for `scripts/CensusClosureAxioms.lean` (must FAIL)

Written by an AI model; not read by a person.  The audit command of `CensusClosureAxioms.lean`,
after one theorem with `sorry` in the namespace `LeanRank.CensusClosure`.
**Expected: exactly 1 error (the audit's "non-standard axioms"), and the `sorry` warning.**
-/

theorem LeanRank.CensusClosure.forged_by_sorry : (1 : ℕ) = 2 := sorry

open Lean Elab Command

run_cmd do
  let env ← getEnv
  let decls := env.constants.fold (init := (#[] : Array Name)) fun acc c _ =>
    if (`LeanRank.CensusClosure).isPrefixOf c then acc.push c else acc
  let std : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut bad : Array (Name × Name) := #[]
  let mut used : Array Name := #[]
  for c in decls do
    let axs ← Lean.collectAxioms c
    for a in axs do
      unless used.contains a do used := used.push a
      unless std.contains a do bad := bad.push (c, a)
  if decls.size < 100 then throwError m!"only {decls.size} declarations found"
  unless bad.isEmpty do throwError m!"non-standard axioms: {bad}"
  logInfo m!"{decls.size} declarations in LeanRank.CensusClosure; axioms used: {used}"
