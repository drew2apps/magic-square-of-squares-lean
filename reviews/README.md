# Review reports

Written by AI models (Claude, made by Anthropic); **none of these is a reading by a person.** An outside review of
commit `24d6ef6` asked whether the reviews that `CHANGES.md` and `LEDGER.md` mention ("reviewed by a separate AI
agent") actually took place, since no report was committed. These are the reports, committed on 26 November 2026
(UTC) as they were written, unedited.

How they were made: each piece of work was done by a build agent in its own git branch; a separate agent, which had not
seen the work being produced, was then given the diff (not the build agent's summary) and asked to act as a hostile
referee: restate the claims, recompute what it could with its own code, try to forge certificates, and list defects.
Before each publication, one more agent reviewed the whole integrated diff (`referee_*.md`). The conductor applied the
fixes; `CHANGES.md` says which. Commit hashes inside the reports are of local branches, not of `main`; paths are those of
the machines the agents ran on.

| File | What it reviews |
|---|---|
| `review_s1.md` | `ClassHensel.lean` |
| `review_s3.md` | `indep_rank.py`, `indep_class.py`, `SECOND_DESCENT_CHAIN.md`, `SecondDescentIdentities.lean` |
| `review_s4.md` | `exact_last_pair.py`, `rank_indep*.py`, `rank_indep_heights.py` |
| `review_s5.md` | `SIGMA_TRANSFER.md`, `scripts/g1_full/sigma_*` |
| `review_s6.md` | `UNIFORM_RANK.md` |
| `referee_int1.md`, `referee_int2.md`, `referee_fix3.md`, `referee_run2.md` | integrated diffs before publication, 25 November 2026 |
| `review_t1.md` | `lball_exact.py` (exact `L(E,1)` for the curve of (440,279), (719,161)) |
| `review_t2.md` | `rank_second.py` |
| `review_t3.md` | `TwoDescentMap.lean`, `TwoDescentCertificates.lean` |
| `review_t4.md` | `SELMER_SQUARE.md`, `selmer_square*.py` |
| `review_t5.md` | `scripts/open_routes/routes.py` (the label claim of `SECOND_DESCENT.md` §6) |
| `review_t7.md` | `lean_rank/` (rank ≥ 2 on the 254 twice-square curves) |
| `review_t8.md` | `SQUARE_SECTION.md`, `square_section.py` |
| `review_v2.md` | rank 0 in Lean through Stoll's `x − T` map (`RankZero*.lean`) |
| `review_v3.md` | `LemmaB`, `Theorem2Family`, `Theorem1Family`, `LemmaDEquiv`, `Theorem1Order`, `FamilyForge` |
| `review_v4.md` | the `rank_second.py` fix and `rank_fisher.py` |
| `referee_run3.md` | the integrated diff later published, after its fixes, as `391e693` |
| `review_w1.md` | the rank 0 of `E1` over `ℚ(√2)` (`Zsqrt2`, `E1RankGlue`, `E1Descent`, `E1RankZero`) |
| `review_w2.md` | `TorsionE1.lean` |
| `review_w3.md` | `CensusClosure.lean` |
| `review_w4.md` | Fisher's `Θ₄` for `E⁺` of (72,65) |
| `review_w6.md` | `lball_exact53.py` |
| `review_w7.md` | the 19 census `L`-value closures |
| `review_bundle.md` | the `experiments/` bundle (census for `1000 < m ≤ 9000`, (627,364) over `ℚ(i)`) |
| `referee_run4.md` | the integrated diff published as `24d6ef6` |
| `review_x1.md` | run 5: `EminusModel.lean`, `UniformRankAll.lean`, `CensusClosureMSS.lean`, `TwoDescentCertForge.lean`, `UniformRankAllForge.lean` |
| `review_x2.md` | run 5: `torsion_exact.py` and the torsion step of the 19 `L`-value closures, the `lball` negative tests, `rank_fisher_theta4_order.py`, `kctp_isog_exact.py` |
| `review_x3.md` | run 5: the Fable bundle (`experiments/census-wide/g4k/`, `rank1/`, `open-ratio/ctp/`); its scripts are in `experiments/referee_x3/` |
| `referee_run5.md` | the integrated diff of run 5, before publication as release 0.6.0 (its fixes are in `CHANGES.md`) |

Work that `CHANGES.md` marks "not reviewed" (for example the (627,364) scripts `kctp_*` of run 4) has no report here. A few earlier reviews, among them those of `ClassPadic.lean` and of `0b34613`, were not saved.
The outside reviews that Drew obtained from other tools are not included.
