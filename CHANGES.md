# Changes

- 0.3.0 (2026-11-21). Lean: `G1Centre.lean` (the centre and `m² + n²`), `G1Relations.lean` (algebra of the Demjanenko–Manin relations), matrix forms in `G1Statement.lean`. Computed: a third curve (`dm3.py`) closes 6 more pairs, 118 in all; every closure rerun with a wider margin and recorded in `scripts/g1_full/dm_evidence.txt`. Result 1b rewritten: scope (modulo `p` only), a definition of density, the counts by script.
- 0.2.0 (2026-11-21, commit b9d82f9). Lean: `G1Norm.lean`, `G1Statement.lean`; `x² = 2y²` without the reals. Computed: 112 pairs closed by the Demjanenko–Manin argument. CI checks the generated files against their generators.
- Commit fff15d8. Lean: `G1Size.lean`, `G1SizeCheck.lean` (an edge root exceeds 500). Its one-line change to `G1Full.lean` removed the keyword `private` from the lemma `sqFive_of_solution`, so that `G1Size.lean` can use it. No definition changed.
