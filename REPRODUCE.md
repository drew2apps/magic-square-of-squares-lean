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
| `rank` (every push) | The sub-project `lean_rank/`: Lean v4.36.0, Mathlib v4.36.0 (compiled files from Mathlib's cache) and M. Stoll's `EllipticCurves` 2164fb5 (compiled from source); the axioms of every theorem named by a `#print axioms` line in its files (76 lines, the three standard axioms only, no `sorryAx`); `scripts/RankChecks.lean` and two files of negative controls that must compile (`RankZeroForge`, `RankZeroCensusForge`); four files of forged claims that must fail with 4, 8, 5 and 2 errors (`RankForge`, `FamilyForge`, `LemmaDForge`, `LemmaDBridgeForge`) | 3 minutes at run 2; not yet measured with the modules of run 3 (the job's limit is 150 minutes) |
| `quick` | The manifest; the transcript checkers; counts; a small rerun of every descent program; Proposition 1 of the rational distance note; the exact checks behind `lean_rank` (the `ℚ(√2)` points of Lemma D, the 2-adic classes of the certificate of `RankZero.lean`) and a regeneration of `RankZeroCensus.lean` by its generator, which must refuse `(15, 8)` | a few minutes |
| `four-points`, `descents`, `last-pair` | The four-point tier and the two descent experiments in full, stated samples of the rest | up to about an hour each |
| `scan-5000` | The four-class descent on all 20,423 ratios to 5000; both rational-distance scans and the exhaustion on six lines; the third implementation on the 157 ratios that ranks leave | about an hour |
| `twist-search` | The C search of twists | minutes |

The long jobs run on a pull request, on a tag, by hand, and on a push whose commit message contains `[full ci]`; on other pushes they are skipped, and a green mark then says less. `OUTPUTS.md` lists the committed outputs that no job reruns, with the command and time of each (the wide runs, the scan beyond 1100, the Magma transcripts and others); the jobs `four-points`, `descents`, `last-pair` and `scan-15000` (and `quick`) end with `git diff --exit-code`.

## 5. What a rerun does not show

A rerun shows that the programs give the committed outputs with the pinned versions. It does not show that the programs compute what the prose says, that the prose proofs are right, or that PARI, Sage and Magma are right. `LEDGER.md` lists, tier by tier, what has to be trusted, and `AUDIT_REQUEST.md` lists the questions for a reader.
