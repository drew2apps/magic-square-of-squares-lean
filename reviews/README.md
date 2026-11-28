# Review reports

Written by AI models (Claude, made by Anthropic); **none of these is a reading by a person.** An outside review of
commit `24d6ef6` asked whether the reviews that `CHANGES.md` and `LEDGER.md` mention ("reviewed by a separate AI
agent") actually took place, since no report was committed. These are the reports, committed on 26 November 2026
(UTC; runs 6 to 10's on 27 November, run 11's on 28 November) as they were written, unedited apart from a provenance note at the top of some.

How they were made: each piece of work was done by a build agent in its own git branch; a separate agent, which had not
seen the work being produced, was then given the diff (not the build agent's summary) and asked to act as a hostile
referee: restate the claims, recompute what it could with its own code, try to forge certificates, and list defects.
Before each publication, one more agent reviewed the whole integrated diff (`referee_*.md`; run 8's was not committed, and release 0.7.0 had a panel of three, `panel_p1.md` to `panel_p3.md`). The conductor applied the
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
| `review_x4.md` | run 6: the Fable bundle, part 2 (`experiments/census-wide/table/`, `lvalue/`, `CENSUS11000.md` and the `9000 < m ≤ 11000` records); its scripts are in `experiments/referee_x4/` |
| `review_y5.md` | run 6: `Θ₄` by the general module `rank_fisher_theta4_gen.py` on `E⁺` of (920,129), (760,39) and the `D1` curve of (631,71) |
| `review_x6.md` | run 6: the Roadmap G bundle (`experiments/roadmap-g/`: `exact_local.py` and the stress test to `m = 10⁹`, the kill sets and test 6 of `COMPLETE.md` §7, the degenerate types for `11000 < m ≤ 13000`, the point search on the square-type fibres); its scripts are in `experiments/referee_x6/` |
| `referee_run6.md` | the integrated diff of run 6 (after release 0.6.0), before publication (its fixes are in `CHANGES.md`) |
| `review_z7.md` | run 7: Fable's roadmap H, H7 (rank 1 for `g4·K` by a full-2-torsion isogenous curve; the `E04` bounds of tier B), and `scripts/second_descent/rank0_ct.py`; its scripts are in `experiments/referee_z7/` |
| `review_z8.md` | run 7: roadmap H, H1 (phantoms) and H3 (local images at good primes); its scripts are in `experiments/referee_z8/` |
| `review_z9.md` | run 7: roadmap H, H2 (survivor hunt), H4, H6 ((627,364), (991,263)) and G3's block `13000 < m ≤ 14000`; its scripts are in `experiments/referee_z9/` |
| `referee_run7.md` | the integrated diff of run 7 (after release 0.6.0), before publication (its fixes are in `CHANGES.md`) |
| `review_z10.md` | run 8: `Θ₄` on `E⁺` of (945,248) (`rank_fisher_theta4_gen.py --side=E'`, the conic finder `rank_fisher_theta4_norm.py`) |
| `review_z11.md` | run 9: `experiments/conductor-t10/` (the Cassels–Tate pairing on the pair curves against the good-prime removals); its scripts are in `experiments/referee_z11/` |
| `review_z12.md` | run 10: `experiments/research-pairing/` (a separate AI session's closed formula for the Cassels–Tate pairing on `E_08`, `E_LU`, tested against `rank_second` on full `Sel²`); its scripts are in `experiments/referee_z12/` (paths in the review are the reviewer's scratch folder) |
| `review_z13.md` | run 11: the second pairing note (`experiments/research-pairing/note2/`), revision 2 and revision 3; its scripts are in `experiments/referee_z13/` |
| `review_z14.md` | run 11: explorations 3 and 5 §1 (the mixed curves, the Rédei form, Theorem B) and exploration 4 §5; `experiments/referee_z14/` |
| `review_z15.md` | run 11: exploration 4 §1–§4 (the pairing on the 937 ratios, the ghost shape, cosets, root numbers); `experiments/referee_z15/` |
| `review_z16.md` | run 11: exploration 5 §2–§3 (the centre in the bad-prime test, the ghost laws); `experiments/referee_z16/` |
| `review_z17.md` | run 11: roadmap I, parts 1–2 (`experiments/roadmap-i/`); `experiments/referee_z17/` |
| `panel_p1.md`, `panel_p2.md`, `panel_p3.md` | the release panel of 0.7.0: the maths referee (the new written proofs), the claims editor, the reproducibility and Lean-statement auditor |
| `referee_run12.md` | run 12 (the changes after the release candidate `c8b677a`: every parameter in Lean, the `covers.py` tolerance, the text after outside reviews), before publication; it could not build the Lean (its fixes are in `CHANGES.md`, *Twelfth run*; its scratch folder is not committed) |
| `referee_run9.md` | the integrated diff of run 9 (after release 0.6.0), before publication (its fixes are in `CHANGES.md`) |
| `referee_run10.md` | the integrated diff of run 10 (after release 0.6.0), before publication (its fixes are in `CHANGES.md`; its scratch folder is not committed) |

Work that `CHANGES.md` marks "not reviewed" (for example the (627,364) scripts `kctp_*` of run 4) has no report here. A few earlier reviews, among them those of `ClassPadic.lean` and of `0b34613`, were not saved.
The outside reviews that Drew obtained from other tools are not included.
