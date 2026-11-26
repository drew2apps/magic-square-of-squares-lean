# Reproducing and checking a copy

For a reader who wants to rerun things, possibly on a machine without network access. `LEDGER.md` says what each result rests on; this page says how to rerun it.

## 1. Is this copy the committed one?

```
python3 scripts/verify_manifest.py
```

`MANIFEST.sha256` lists the SHA-256 of every file. The command needs only Python 3 and ends with `MANIFEST OK` when every file has the listed hash and no file is unlisted (build directories `.lake` and Python caches are ignored). CI runs it on every push, so no result file can change without the manifest changing in the same commit. The archived releases at Zenodo (doi:10.5281/zenodo.24435351) are copies of tagged commits; the last tag is `v0.5.0`, and everything after it is unreleased work (`CHANGES.md`, first section).

## 2. Pinned versions

| Component | Version | Where it is pinned |
|---|---|---|
| Lean | `leanprover/lean4:v4.33.1` | `lean-toolchain` |
| Mathlib and its dependencies | commit `0df444a360eaa60ab8c11dca51a86af692955474` and the revisions of `lake-manifest.json` | `lake-manifest.json` |
| Python | 3.13 | `.github/workflows/ci.yml` |
| PARI (through `cypari2`), passagemath (Sage libraries, eclib) | `cypari2==2.2.9` (PARI 2.17.4), `passagemath-*==10.8.14` | `scripts/g1_full/requirements.txt`; `passagemath-eclib==10.8.14` is installed in the job `four-points` for `eclib_check.py` |
| sympy, python-flint, numpy | `sympy==1.14.0`, `python-flint==0.9.0`, `numpy==2.5.5` | `.github/workflows/ci.yml` |
| GitHub Actions | `actions/checkout`, `actions/setup-python`, `leanprover/lean-action`, each by commit hash | `.github/workflows/ci.yml` |
| C compiler | `gcc` of the runner image `ubuntu-latest`, not pinned | — |

Not pinned and not reproducible here: **Magma**. Its outputs are transcripts of the public calculator (`*_raw.txt`), with the version printed in each transcript; CI checks that the committed summaries follow from the transcripts, not the transcripts themselves. Python package hashes are not pinned, and only the top-level packages have pinned versions; `scripts/g1_full/requirements-freeze.txt` records what pip resolved for the others on 24 November 2026.

## 3. Without network access

On a machine with network access, once:

```
git clone https://github.com/drew2apps/magic-square-of-squares-lean && cd magic-square-of-squares-lean
lake exe cache get          # Mathlib build products for the pinned commit (or build Mathlib from source: several hours)
pip download -d wheels -r scripts/g1_full/requirements.txt sympy==1.14.0 python-flint==0.9.0 numpy==2.5.5
```

Copy the directory (with `.lake` and `wheels`) and the Lean toolchain (`~/.elan`) to the other machine. There:

```
python3 scripts/verify_manifest.py
pip install --no-index --find-links wheels -r scripts/g1_full/requirements.txt sympy==1.14.0 python-flint==0.9.0 numpy==2.5.5
lake build                                        # every Lean proof, about 20 minutes with the Mathlib cache
lake env lean scripts/Audit.lean                  # axioms of every theorem; ends with AUDIT OK
lake env lean scripts/CoverDescentChecks.lean     # the count 306 + 119 = 425, and forged certificates refused
```

## 4. The computations

The executable specification is `.github/workflows/ci.yml`: every job is plain shell, uses no network after `pip install`, and compares its output with committed files by `diff`. To rerun a job by hand, run the commands of its `run:` block from the root of the repository.

| Job | What it reruns | Time on GitHub's runner |
|---|---|---|
| `build` | Lean; the axiom audit; generated Lean files; the C searches on small ranges | about 20 minutes |
| `rank` (every push) | The sub-project `lean_rank/`: Lean v4.36.0, Mathlib v4.36.0 (compiled files from Mathlib's cache) and M. Stoll's `EllipticCurves` 2164fb5 (compiled from source); every `#print axioms` line of `LeanRank/*.lean` (file and line) answered in the output of `lake build`, the three standard axioms only, no `sorry` (147 lines at run 4); six files that must compile (`RankChecks`, `E1RankChecks`, `CensusClosureChecks`, `CensusClosureAxioms`, and the negative controls `RankZeroForge`, `RankZeroCensusForge`); eight files of forged claims, each of which must fail at exactly the listed lines (`RankForge`, `FamilyForge`, `LemmaDForge`, `LemmaDBridgeForge`, `E1DescentForge`, `TorsionE1Forge`, `CensusClosureForge`, `CensusClosureAxiomsForge`) | 6 minutes at run 3 on GitHub; the run-4 modules add about 15 minutes on a loaded 8-core machine (`E1Descent` alone about 8); the job's limit is 150 minutes |
| `quick` | The manifest; the transcript checkers; counts; a small rerun of every descent program; Proposition 1 of the rational distance note; the exact checks behind `lean_rank` (the `ℚ(√2)` points of Lemma D, the 2-adic classes of the certificate of `RankZero.lean`) and a regeneration of `RankZeroCensus.lean` by its generator, which must refuse `(15, 8)`; since run 4 the exact checks of `E1(ℚ(√2))` (`e1_rank0_check.py`, `torsion_e1_check.py`), a regeneration of `CensusClosure.lean` with `census_closure_check.py`, `selmer_square.py --max 300`, `review_t8_rank1.py`, and `lball_exact53.py` on three curves | about 15 to 20 minutes (13 at `78659d1`; limit 45) |
| `four-points`, `descents`, `last-pair` | The four-point tier and the two descent experiments in full, stated samples of the rest | up to about an hour each (`descents` about 90 minutes since run 4, with `lball_exact53.py` on 51 curves and the `Θ₄` and Fisher scripts; limit 150) |
| `scan-5000` | The four-class descent on all 20,423 ratios to 5000; both rational-distance scans and the exhaustion on six lines; the third implementation on the 157 ratios that ranks leave | about an hour |
| `scan-15000` | The criterion and the four-class descent for `5000 < m ≤ 15000` | about 50 minutes (48 at `c122d93`) |
| `twist-search` | The C search of twists | minutes |
| `square-section` (since 26 November 2026) | `square_section.py --mmax 300` (`SQUARE_SECTION.md`), against its two outputs | about 5 minutes (7 on a loaded 2-core container) |

The long jobs run on a pull request, on a tag, by hand, and on a push whose commit message contains `[full ci]`; on other pushes they are skipped, and a green mark then says less. `OUTPUTS.md` lists the committed outputs that no job reruns, with the command and time of each (the wide runs, the scan beyond 1100, the Magma transcripts and others); the jobs `build`, `quick`, `four-points`, `descents`, `last-pair`, `square-section`, `scan-15000` and `scan-5000` end with `git diff --exit-code`.

## 5. What a rerun does not show

A rerun shows that the programs give the committed outputs with the pinned versions. It does not show that the programs compute what the prose says, that the prose proofs are right, or that PARI, Sage and Magma are right. `LEDGER.md` lists, tier by tier, what has to be trusted, and `AUDIT_REQUEST.md` lists the questions for a reader.
