# referee_z13: checks of review z13 of `PAIRING_FORMULA_2.md` (the second note on the pairing on `E_08`, `E_LU`)

Written by an AI model (reviewer z13), not read by a person. Experiment, not counted; tier D; nothing here changes the
census. The note reviewed is a separate AI research session's `PAIRING_FORMULA_2.md` (not in this repository at the
time of the review); the review itself is the reviewer's report. PARI 2.17 (cypari2), Python 3; the repository's
`rank_second.py` (`full_selmer2`, `DFull`) through `../conductor-t10/t10_pairct.py` (`DFullX`, `safe_pair`).

| Script | What it does | Output | CI |
|---|---|---|---|
| `z13_formula.py` | the reviewer's own implementation of formula (2), `S′` and `Θ` (3) of the note, from its text (the note's scripts are not imported) | — | — |
| `gen_jobs.py` | the job list (seed 20261127): 378 curves of 197 ratios outside the pilot; random `m ≤ 20000` and `m ≤ 2000`, both odd, `A` odd with `v₂(C ∓ A) = 5`, square and twice-square type, all small ratios `m ≤ 18`, and 16 controls outside H2 | `jobs.txt` | `python3 -I gen_jobs.py \| diff - jobs.txt` (1 s) |
| `z13_indep.py` | per curve: `Sel²` by `full_selmer2` (seed 1), the pairing by `DFullX` for every `a = (1, δ, δ)` in `Sel²` (at most 128; above that a basis and 10 more) against a basis of `Sel²`, extended bilinearly to all of `Sel²`; compared with formula (2) on every `(a, b)`, `S′` with `{b_r = 1}`, `Θ` with the pairing on `S′ × S′` | `indep_a.jsonl`, `indep_b.jsonl`, logs `indep_[ab].log` | no (about 25 minutes on 2 shared cores); `rerun_z13.py` reruns a subset |
| `rerun_z13.py` | reruns `z13_indep.run` on every 25th job and the first of each category, compares every field but the time | `rerun_z13.txt` | yes (89 s) |
| `analyse_z13.py` | the counts of the review from the two `.jsonl` files, and side checks (e1)–(e5) | `analyse_z13.txt` | yes (30 s) |
| `h2_rerun_z13.py` | revision 3 follow-up: the five curves of the note's `h2_check.py` (which imports a module not supplied), recomputed with `z13_indep.run` | `h2_rerun_z13.txt` | yes (54 s) |
| `claims_z13.py` | the note's pilot numbers per line and per distinct case (sections 4–6) | `claims_z13.txt` | yes (1 s) |

Rerun a check: `python3 -I X.py 2>/dev/null | grep -v 'PARI stack' | diff - X.txt` from any folder.
