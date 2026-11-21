import MagicSquareOfSquares

/-!
Axiom audit. Run with `lake env lean scripts/Audit.lean`.

Lists every theorem declared in this project, collects every axiom that any of them depends
on, and fails if that set contains anything other than `propext`, `Classical.choice` and
`Quot.sound`.
-/

open Lean Elab Command

namespace Audit

structure State where
  visited : NameSet := {}
  axioms : Array Name := #[]

/-- Walk the dependencies of `c`, recording axioms. The visited set is shared between calls. -/
partial def collect (env : Environment) (c : Name) : StateM State Unit := do
  let collectExpr (e : Expr) : StateM State Unit := e.getUsedConstants.forM (collect env)
  unless (← get).visited.contains c do
    modify fun s => { s with visited := s.visited.insert c }
    match env.find? c with
    | some (.axiomInfo _) => modify fun s => { s with axioms := s.axioms.push c }
    | some (.defnInfo v) => collectExpr v.type *> collectExpr v.value
    | some (.thmInfo v) => collectExpr v.type *> collectExpr v.value
    | some (.opaqueInfo v) => collectExpr v.type *> collectExpr v.value
    | some (.quotInfo _) => pure ()
    | some (.ctorInfo v) => collectExpr v.type
    | some (.recInfo v) => collectExpr v.type
    | some (.inductInfo v) => collectExpr v.type *> v.ctors.forM (collect env)
    | none => pure ()

end Audit

run_cmd do
  let env ← getEnv
  let names := env.constants.fold (init := (#[] : Array Name)) fun acc n ci =>
    match env.getModuleIdxFor? n with
    | some idx =>
      let modName := env.header.moduleNames[idx.toNat]!
      if (`MagicSquareOfSquares).isPrefixOf modName then
        match ci with
        | .thmInfo _ => acc.push n
        | _ => acc
      else acc
    | none => acc
  let sorted := names.qsort (fun a b => a.toString < b.toString)
  let (_, st) := (sorted.forM (Audit.collect env)).run {}
  for n in sorted do
    unless n.isInternalDetail do
      logInfo m!"{n}"
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let bad := st.axioms.filter (fun a => !allowed.contains a)
  logInfo m!"axioms used: {st.axioms.qsort (fun a b => a.toString < b.toString)}"
  if bad.isEmpty then
    logInfo m!"AUDIT OK: {names.size} theorems (including auxiliary ones), standard axioms only"
  else
    throwError "AUDIT FAILED: non-standard axioms {bad}"
