import MagicSquareOfSquares

/-!
The theorems that `scripts/Audit.lean` counts but does not print: the auxiliary theorems that Lean
generates (internal details: equation lemmas, `proof_` terms of tactic blocks, and the like).
Run with `lake env lean scripts/AuditInternal.lean`; output in `scripts/audit_internal.txt`.
Together with the names that `Audit.lean` prints, this accounts for its count `N`.
-/

open Lean Elab Command

run_cmd do
  let env ← getEnv
  let mut shown : Array Name := #[]
  let mut internal : Array (String × Name) := #[]
  for (n, ci) in env.constants.toList do
    match env.getModuleIdxFor? n, ci with
    | some idx, .thmInfo _ =>
      let modName := env.header.moduleNames[idx.toNat]!
      if (`MagicSquareOfSquares).isPrefixOf modName then
        if n.isInternalDetail then internal := internal.push (modName.toString, n)
        else shown := shown.push n
    | _, _ => pure ()
  let sorted := internal.qsort (fun a b => a.1 < b.1 || (a.1 == b.1 && a.2.toString < b.2.toString))
  for (m, n) in sorted do
    logInfo m!"{m} {n}"
  logInfo m!"AUDIT INTERNAL: {shown.size} theorems printed by Audit.lean, {internal.size} auxiliary, {shown.size + internal.size} in all"
