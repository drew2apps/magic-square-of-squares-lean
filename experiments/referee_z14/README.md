# referee_z14: the Cassels–Tate pairing on the mixed pair curves (review of untrusted explorations 3, 4 §5, 5 §1)

Written by an AI model (reviewer z14), not read by a person. Experiment; nothing here changes a count. The report is
`reviews/review_z14.md` (outside the repository).

| file | what | run (CPU on a loaded 2-core machine) |
|---|---|---|
| `e3_cassels.py`, `e3_sel.py`, `e3_rank0.py` | explore3's `cassels.py`, `sel.py`, `rank0.py`, untrusted, copied with only the `sys.path` lines changed | — |
| `ratios.py` → `ratios.txt` | the 134 test ratios (seed 14), none in the pilot, none a pilot partner | `python3 -I ratios.py \| diff - ratios.txt`, 1 s |
| `mixtest.py` → `mix_[012].jsonl` | tasks 1(a)–(c) on the mixed curve of each ratio against `rank_second`/`t10_pairct.DFullX` | `mixtest.py OUT i0 i1`, about 8 s CPU per curve; records carry a `time` field |
| `theoremB_z14.py` → `thB_[01].jsonl` | EXPLORATION_5 Theorem B against DFull on E_mix (two labellings), E_08, E_LU | `theoremB_z14.py OUT i0 i1`, about 20 s CPU per ratio |
| `rank0_z14.py` → `rank0_z14.jsonl` | EXPLORATION_3 §4 recomputed with `full_selmer2` + DFull | about 2 s CPU per ratio |
| `e3_rank0_perratio.py` → `e3_rank0_perratio.jsonl` | explore3's own `bound()` per ratio on the same three sets | 11 s CPU |
| `e3_rank0_rerun.txt` | `e3_rank0.py 3000 150 nd` and `3000 80 sq` as run here (the nd part reproduces `rank0.txt`, the sq part does not) | 15 s |
| `mix_analyse.py` → `mix_analyse.txt` | the summary of all JSONL files (deterministic) | `python3 -I mix_analyse.py \| diff - mix_analyse.txt`, < 1 s |
| `symbolic.py` → `symbolic.txt` | sympy: §1 identities, the tangent lemma, §2 (discriminants, fibres, j, 4-torsion, twist, corner forms) | `python3 -I symbolic.py \| diff - symbolic.txt`, 3 s CPU |
| `section_z14.py` → `section_z14.txt` | psi(iota t0) has infinite order on E_mix (5 specializations) | 1 s |
| `kappa_z14.py` → `kappa_z14.txt` | EXPLORATION_4 §5's constants kappa on 20 non-pilot ratios via `t10_pairct.PairData` | about 1 min |
| `conic_test.py` → `conic_test.txt` | explore3's `conic_pairing` (function text verbatim from task23.py) against DFull on the 134 curves | about 2 min |
| `redei_trivial.py` → `redei_trivial.txt` | how many of mixtest's soluble Redei pairs have d = beta (replays mixtest's draws) | about 1 min |
| `claims_z14.py` → `claims_z14.txt` | distinct-case counts behind the notes' numbers (pilot and explore3 logs) | `python3 -I claims_z14.py \| diff - claims_z14.txt`, 1 s |
| `probe.py` | first probe (Pair models on new ratios, DFull timing); prints timings | — |

`claims_z14.py` reads explore3's logs from `/home/claude/untrusted/explore3/explore3/` (not in the repository).
