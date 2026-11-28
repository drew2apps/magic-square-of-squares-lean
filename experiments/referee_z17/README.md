# Review z17: roadmap I (parts 1 and 2)

Written by an AI model (reviewer z17), not read by a person. Tier D, not counted; nothing here changes the census.
The review itself is the conductor's `reviews/review_z17.md`. These are its checks, run under PARI 2.17.4 (cypari2 2.2.9), sympy 1.14.0.
Each script prints its output and `# matches` against the committed .txt when run without `--write`.

| Script | What it checks | Time |
|---|---|---|
| `z17_i2check.py` | I2's residual file recounted (32,474 tuples, even sizes, killing primes good and <= 229); 50 ratios (seed 17; 20 of them from G2's lines) recomputed with conductor-t10's `residual()` | 30 s |
| `z17_i8count.py` | I8's checkpoint numbers recounted from `i8_first360.jsonl` | < 1 s |
| `z17_i8rerun.py`, `run_i8rerun.sh` | conductor-t10's analysis from scratch (its own residual) on the 3 "neither" ratios and 15 others (seed 17), one process per ratio -> `z17_i8rerun.jsonl` | 21 min |
| `z17_i8compare.py` | that rerun against `i8_first360.jsonl` | < 1 s |
| `z17_neither.py a`..`e` | the "neither" closures: the DM curves equal I8's E+/E-; PARI's ranks; dm_pari rerun; exact rank bounds (H7's step, rank_indep); PARI's points of infinite order checked exactly | a, b, c, e: 10 s; d: 1.5-4 min |

`neither/` holds I8's `epm.py`, `ri.py` and `dm_pari.py` outputs rerun here (identical to the committed ones). `ci/` holds the proposed CI lines' outputs and times (`times.txt`: under load, about 11 on 2 cores; `hm_clean_tree.txt`: through CI's `hm` helper on this tree).
