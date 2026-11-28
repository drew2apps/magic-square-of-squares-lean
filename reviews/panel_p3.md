# Panel P3: reproducibility / supply chain / formal statements, release 0.7.0

Written by an AI model (panel P3), not read by a person.

> Committed as written apart from this note. It reviewed the release candidate of 0.7.0 (commit c857655 of the integration branch); its fixes are in `CHANGES.md` ("After the release panel"). Paths under `/home/claude/` are the panel's scratch folders.

Reviewed: /home/claude/conductor/int11 (branch run11_int11, HEAD c857655), against 0.6.0 = dabf7b4.
Clone used: /home/claude/conductor/panel_p3/clone. Findings appended as they are made.

## 1. CI (findings as made)

- **[blocker] ci.yml:743 (job `research`) fails on most runs: `grep -q` on a live pipe under `pipefail`.**
  `python3 -I t11_theta.py --negative 2>/dev/null | grep -q '^negative tests: 9 of 9 as expected'`: the matched line is
  the second-last line; `t11_negative.py:75` then prints `# matches`, grep has already exited, Python gets EPIPE on the
  final flush and exits 120, and `pipefail` makes the subshell fail, so `set -e` stops the job. Measured from a clean
  clone: 5 of 6 runs through the research helpers' harness exited 120, and 4 of 4 under `bash -c 'set -o pipefail; ...'`
  (`PIPESTATUS` = `120 0`, stderr "BrokenPipeError ... flushing sys.stdout"). It is a false red, not a false green, but
  the job cannot be expected to pass. Fix: `hm experiments/conductor-t11 -I t11_theta.py --negative` (the script prints
  `# matches` last and compares the whole output with `t11_theta_negative.txt`, which the present line does not check),
  or write to a file and grep the file (the pattern already used at ci.yml:220-222, 276).
- [ok] Lines 1-14 of ci.yml: actions pinned by full commit hash (checkout fbc6f39 v5.1.0, setup-python ece7cb0 v6.3.0,
  lean-action ba3fdfe v1). YAML parses (`yaml.safe_load`): 11 jobs; timeouts build none (GitHub default 360 min), rank 150,
  quick 45, four-points 90, descents 150, research 90, last-pair 120, twist-search 30, square-section 30, scan-15000 120,
  scan-5000 150.
- [ok] **The `research` helpers fail closed** (ci.yml:714 `d`, :716 `hm`; tested in a second clone, the function
  bodies taken verbatim from ci.yml by line number, `/tmp/` replaced by a private folder because other panels share
  this machine, under `set -e -o pipefail`; after each case a sentinel `echo` tells whether the job went on):
  `d` with the expected file changed in one character, with a line added, with the last line dropped, missing, inside a
  `for` loop, and with a script that prints the exact expected text but exits 1: exit 1 each time, job stopped. `hm` with
  `i1_laws.txt` changed in one character, last line dropped, missing (the script prints `# DIFFERS`); with
  `i8_compare50.txt` missing (the script prints no verdict); tampered inside the z17 `for` loop; a script printing
  `# DIFFERS` then `# matches`; a script printing `# matches` then exiting 3: exit 1 each time. Controls (untampered
  `kcurve_exact`, `i1_laws`, `i8 --compare50`, a fake script printing `# matches`): exit 0 and the job goes on. The older
  `hm` of `descents` (ci.yml:646, no `|| exit 1`, relies on `set -e`): tampered `i1_laws.txt` stops the job too.
  One property to know: `d` with an empty expected file and a script that prints nothing fails (`grep -v` exits 1);
  no expected file of the research job is empty, so this is only a trap for later lines.
- [ok] The run-10 z11 loop (ci.yml:672-675): tested with `recount.txt` tampered (exit 1, `== recount: stderr` printed)
  and with `recount.py` made to print its full output and then exit 2 (exit 1, the forced stderr line shown). Good.
- [nit] ci.yml:713, the comment of `d` says `d DIR SCRIPT [ARGS] < ... the file given by OUT=...`; the signature is
  `d DIR OUTFILE PYTHON-ARGS...` and there is no `OUT=` variable.
- [ok] `hm` scripts of `research` (i1_laws, i1_negative, i2_residual, i3_scan, i4_exchange, i5_rank, i8_pairct, the z17
  scripts, t11_theta): each compares its whole stdout text with its committed file (`open(f).read() == s`) and writes
  only under `--write`; those that print a verdict only `if os.path.exists(f)` print none when the file is missing,
  which `hm` treats as failure (tested above).
- [fix] ci.yml:911-966, job `scan-5000` does not end with `git diff --exit-code` (it has one for `class_good_table.txt`
  only, line 945; the last line is `pt.py scan`). REPRODUCE.md:63 says it does. Add the line (or correct the text).
  `twist-search` and `rank` have none (they write only /tmp and untracked logs; fine). `git diff --exit-code` does not
  see new untracked files; no new CI line writes one in the repository except `rank`'s logs (checked: after the spot runs
  `git status --porcelain` showed only ignored `__pycache__/`).
- [fix] **`descents` time.** Not measurable here (no full run; about 2 h). Estimate from the committed time notes:
  about 95 min at 0.6.0 plus the rows added since (runs 6-10, summed from OUTPUTS.md's time column) about 27 min,
  so about 2 h of the 150-minute limit (about 20 % margin), measured on a loaded 2-core container. The repository has
  no record of a *green* full `descents` run since 0.6.0: CHANGES.md:19, :29 record the failures at 1eb18a9 (run 6,
  spot.py glob order) and f0d7edc (run 9, h1_mine.py stopped after 4 lines, cause not found). So the "about 2 hours"
  of REPRODUCE.md:56 is an estimate. The release rule (ci.yml:11, "a release is made only from a commit with a full green
  run") must be met by the final manifest commit with `[full ci]`; if `descents` takes more than about 130 min there,
  move the run-6-10 tier-D blocks (ci.yml:582-688) to `research` or a job of their own. No step-level
  `timeout-minutes` anywhere (a hung step costs the whole job limit).
- [nit] Pins: `python-version: '3.13'` resolves to the newest 3.13.x on each run (REPRODUCE says "3.13", correct);
  `ubuntu-latest`, gcc, and transitive pip dependencies unpinned (documented in REPRODUCE.md:23-25). The `build` job has
  no `timeout-minutes` (default 360). Lean: lean-toolchain v4.33.1, lean_rank/lean-toolchain v4.36.0, Mathlib by
  manifest commit: as documented.
- [process note] P3's dry run of ci.yml (to map file references) briefly truncated the shared files `/tmp/hm.out`,
  `/tmp/d.err`, `/tmp/seven.out` between about 17:39 and 17:50 UTC (machine clock), and started `/tmp/seven_brute 1300`
  and four `/tmp/twist_find 400` binaries for about a minute before they were killed. If another panel saw an
  unexplained `hm`/`d` failure in that window, rerun it. (Other panels' CI reruns on this machine also share these
  /tmp names with each other.)

## 2. Coverage of the 1,204 files added since dabf7b4

Method: every `run:` block of ci.yml executed in a dry-run shell (python3, diff, cmp, xargs, compiled binaries replaced
by loggers; cat/head/tail/sed/grep passed through and logged), so that loops, `$(basename ...)`, `${s#*:}` and the
cwd of every command are expanded by bash itself; every argument resolved against the logged cwd; `hm` scripts mapped to
the files they compare by the string literals next to their `# matches` test. Then every added data file (781 of the
1,204 are not .py/.md/.sh/.lean/.c) matched against OUTPUTS.md rows (row prefix + name, globs, "outputs of `x`").

- Regenerated in full and compared by `diff`/`cmp`: 142 files. Compared in full by an `hm` script: 135. Named only as
  an input of a CI command: 15 (`t10_*` lists, the ten `conductor-t13/search_*.json`, `altmodels.json`).
  Never named by CI: 489.
- **[fix] OUTPUTS.md:297 claims class A (job `research`) for `conductor-t13/search_*_H2000_torsion.txt`, but CI
  regenerates only three of the four: `search_E_alt_H2000_torsion.txt` is never rerun** (ci.yml:747-749 have `E`, `E3`,
  `E3 --models=altmodels.json`; no `--curve=E --models=altmodels.json`). Add
  `d $T search_E_alt_H2000_torsion.txt search.py --models=altmodels.json --curve=E --H=2000 --classes=torsion --skipbase`
  (check the flags against NOTES.md) or narrow the row.
- [fix] Never regenerated and not in OUTPUTS.md under any row (by name, glob, folder row or "outputs of"):
  - `experiments/conductor-t13/covers.json` (written by `covers.py --write`, read by `search.py` and `certify.py`; CI
    compares only `covers.txt`), `hints_t13.json` (written by `certify.py --write`), the ten `search_*.json` read by
    `certify.py` (the `*_coset*.txt` row covers the .txt only), `altmodels.json`, `logs/*.log` (8). Suggest class D/B
    rows; better, run `covers.py --write` (and the three torsion `search.py ... --write`) in CI and let the final
    `git diff --exit-code` compare the JSON too.
  - `experiments/genus3/genus3_points_hints.json` (read by `genus3_points.py`; class D row missing).
  - `experiments/referee_z15/`: `ghost_order.txt`, `notA_order.txt`, `notA_order_mid.txt`, `notA_order_rev.txt`,
    `s80_a.txt`, `s80_b.txt`, `sample80_seed20261127.txt`, `t10_times.txt`, `session_data/` (4 files, copies of the
    session's records).
  - `experiments/referee_z16/`: `ratios937.txt`, `ratios_open.txt`, `ratios_samples.txt` (inputs), `time*.txt` (3),
    `data/` (9 files incl. `SHA256SUMS`, `SHA256SUMS.sess`, copies of the session's records), `sess/` (2).
  - `experiments/referee_z17/`: `ci/` (18 files: the reviewer's transcripts of the proposed CI lines, `.out`/`.err`),
    `neither/` (6), `z17_i8rerun.err`, `z17_i8rerun_times.txt`.
  - `experiments/roadmap-i/i8/neither/` (8), `i8/i8_first360_summary.txt`; `i2/` lists (`i1_classkept_list.txt`,
    `i2_i1_classkept_list_0of1.txt`, `i2_i2_needgood_list_{0,1}of2.txt`); `i4/` non-chunk records (`dm_*`, `r0_*`,
    `sd_*`, `stage_*`, `fam_11000_30000.json`: about 30 files; the row's "`i4` chunks" is loose).
  - `experiments/referee_z14/e3_rank0_rerun.txt`; `experiments/referee_z9/g3_samples_in/` (3 inputs);
    `experiments/roadmap-g/g3/done_11000`, `done_12000` (row OUTPUTS.md:241 names only `done_13000`).
  - timing notes `times.txt` in conductor-t11, genus3, roadmap-g/g2, g3, g4, roadmap-h/h1-h7 (6), roadmap-i/i1, i2, i3,
    i5, i8 (about 20; only r6557's and g1's are listed). One class-C row "timing notes (`times.txt`, `*_times.txt`)"
    would do.
  Everything else of the 489 is covered by a row (B 232, C 20, D 8, A rows that name a B companion).
- [ok] No `__pycache__`, `.pyc`, `.part`, editor or OS files added. 36.1 MB added in all; largest: 3.29 MB
  `census-wide/table/CENSUS_TABLE_1000_11000.tsv`, 2.52/2.51/2.50 MB `roadmap-g/g4/g4_{1,2,0}_3.json`, 1.99 MB
  `referee_z12/full0.jsonl`, 1.37 MB `full1.jsonl`, 1.06 MB `roadmap-i/i2/i2_residual.jsonl`, then 0.82 MB and below
  (30 largest listed in P3's notes; none over 3.3 MB). Duplicated copies: `out_15000_{0,1}.jsonl` three times
  (explore4, referee_z15/session_data, referee_z16/data), `g2_15000_all.txt` twice (about 1.7 MB) [nit].
- [ok] No credentials or tokens (regex scan for key/secret/password/token, PEM headers, ghp_/sk-/AKIA patterns: one hit,
  the word "token" in prose of reviews/referee_run10.md). Absolute home paths: 17 review reports name the reviewers'
  scratch folders (explained by reviews/README); outside reviews/: `referee_z14/claims_z14.py:81` and its README read
  `/home/claude/untrusted/explore3/explore3/` (documented in OUTPUTS.md:289, class B), `genus3/session/NOTE_2026-11-27.md:1`
  (provenance line), `referee_z17/ci/i1_negative.err:2` (a traceback path). None matters for a CI rerun.
- **[fix] Two CI-run scripts put a world-writable path ahead of the repository on `sys.path`:**
  `experiments/research-pairing/note2/closed_formula.py:6` and `note2/shape.py:2`, `sys.path.insert(0,'/tmp/t10')`
  (`shape.py` then does `from parse import sqf`). On GitHub's runner /tmp/t10 does not exist, so CI is unaffected; on a
  reader's machine where /tmp/t10 holds the session's old modules, the research job's note2 lines would import those
  instead of `experiments/research-pairing/parse.py`. Remove the two inserts with a comment, as review z14 did for its
  copies (`referee_z14/e3_rank0.py:4`, `e3_sel.py:3`, `e3_cassels.py:10`). (Run 6's census-wide and open-ratio scripts
  carry `/tmp/claude-0/...` and `/mnt/project-files/...` paths too; they are class B/C and only `py_compile`d.)
- [ok] **MANIFEST.sha256 verifier logic** (scripts/verify_manifest.py): walks the tree (skipping `.git`, `.lake`,
  `__pycache__`, `*.pyc`, itself), and fails on a hash mismatch, a listed file missing, a present file not listed, or a
  listed file not tracked by git; sound for its stated purpose (it cannot detect a manifest rewritten together with the
  files: the tag/Zenodo DOI is the anchor). At HEAD it fails as expected (1677 listed, 2222 tracked). No tracked name
  needs git's quoting (0 non-ASCII names), so the `git ls-files` comparison has no false mismatches.

## 3. Negative tests for the new tier-C/D results (LEDGER rows added since dabf7b4)

With a forged-input, mutated-rule or refusal test run in CI: Θ₄ general (`rank_fisher_theta4_gen.py --negative` on four
curves, `_order.py --negative`), Θ₄ on (945,248) (`--negative`), the census table 1000-11000 (`census_table.py
--negative`, `lvalue_stream.py --negative[-small]`, `referee_x4/recount.py --negative`), roadmap G1-G4 (`g1_stress`,
`g2_check`, `g3_table`, `g4_search` `--negative`; `referee_x6/local_mutation.py`), roadmap H (h1, h2_control, h3_tables,
h3_weil, h4, h6_ranks, h6_richelot, h7 `--negative`), rank0_ct (`--negative`), t10 (`--negative`, z11's `forge_t10`,
`forge_mine`), genus-3 points and the K-curve (`genus3_points.py --negative`, `kcurve_exact.py --negative`), I5's
Fisher pairings (`t11_theta.py --negative`, but see the blocker: the CI line is broken), roadmap I (i1, i2, i3, i4, i5,
i8 `--negative`). The closed pairing formula (note2, tier C) has a natural sensitivity test: review z13's controls
outside hypothesis H2 (`referee_z13/gen_jobs.py:56`), where formula (2) does disagree, so the comparison can fail.

Without one (where it matters):
- [fix] **Theorem B** (tier C, `explore5/` §1; LEDGER "Theorem B" row): CI reruns only z14's symbolic checks
  (`referee_z14/symbolic.py` itself says some of its controls "cannot fail"); the 454-pair comparison with
  `rank_second` (`theoremB_z14`, class B) is not in CI and has no mutated-formula run (e.g. one local symbol flipped
  must give disagreements). Same for **the Rédei-symbol law on the mixed curves** (tier D, `explore3/`, z14
  `mixtest.py`): agreement counts only.
- [fix] **(627,364) over ℚ(i)** (tier D, `conductor-t13/`): the claim is a negative search result; its only sensitivity
  evidence, `controls.py` (rank 2 certified on 7 control curves from the pipeline's own points), is class B and not in CI,
  and `certify.py` has no forged-hit test (a hit off the curve or in the wrong class must be refused). A cheap CI
  control (one control curve at small H) would make the torsion-search lines mean more than "same output".
- [nit] **The nine open ratios** (tier D, `exchange/open_pairs.py`): no positive control that the square test sees
  squares (e.g. at `t = 1` on a square-type open ratio three forms are squares). The 160 `F4` hits show it is not
  constantly false, so this is minor.
- [ok, by nature] Genus-3 splitting (tier C): ring/isogeny identities with a trace control; a mutation test adds little.
  The centre in the bad-prime test (z16): planted degenerate solutions and p-adic points must not be refused (0
  refused), which is the direction that matters for an exclusion rule.

## 4. Lean: Exchange.lean, ExchangeCensus.lean (read, not built: no toolchain here)

Definitions used: `nrm (x, y) = x² + y²`; `P0 = (mv, nu)`, `P8 = (nv, mu)`, `PL = ((m+n)(u−v), (m−n)(u+v))`,
`PU = ((m−n)(u−v), (m+n)(u+v))` (CoverDescent.lean:234-239); `FullSquares m n t` for `m n : ℕ`, `t : ℚ`: the five
rationals `F0 = m² + n²t²`, `F4 = (m² + n²)(t² + 1)/2`, `F8 = n² + m²t²`, `L = F4 − 2mnt`, `U = F4 + 2mnt` are squares in
ℚ (G1.lean:36, G1Full.lean:29).

Exchange.lean, in plain maths (all for arbitrary integers `m, n, u, v` unless said):
1. `P0_exchange`: `P0(u,v,m,n) = swap P0(m,n,u,v)`; `PL_exchange`: same for `PL`; `P8_exchange`, `PU_exchange`:
   `P8(u,v,m,n) = P8(m,n,u,v)`, `PU(u,v,m,n) = PU(m,n,u,v)`. (Checked by hand: e.g. `PL(u,v,m,n) = ((u+v)(m−n), (u−v)(m+n))`.)
2. `nrm_P0_exchange`, `nrm_P8_exchange`, `nrm_PL_exchange`, `nrm_PU_exchange`: the four norms are invariant.
3. `nrm_PL_add_nrm_PU`: `nrm PL + nrm PU = 4(m²+n²)(u²+v²)`; `nrm_P0_add_nrm_P8`: `nrm P0 + nrm P8 = (m²+n²)(u²+v²)`.
4. `fiveSquaresAt_exchange`: with `FiveSquaresAt(m,n,u,v)` := the four norms and `2(m²+n²)(u²+v²)` are squares in ℤ,
   `FiveSquaresAt(u,v,m,n) ↔ FiveSquaresAt(m,n,u,v)`. No hypotheses.
5. `G1.isSquare_mul_sq_iff`: for `a, c ∈ ℚ`, `c ≠ 0`: `a c²` is a square iff `a` is.
   `G1.isSquare_congr_of_mul_sq`: `c, d ≠ 0`, `a c² = b d²` ⇒ (`a` square ⇔ `b` square).
6. `G1.fullSquares_exchange`: for `m, n, u, v ∈ ℕ` with `n ≠ 0`, `v ≠ 0`:
   `FullSquares m n (u/v) ⇔ FullSquares u v (m/n)`. Each form at `(m,n,u/v)` times `v²` equals the same form at
   `(u,v,m/n)` times `n²` (checked by hand for all five; F0↔F0, F4↔F4, F8↔F8, L↔L, U↔U).
   The hypotheses are needed (Lean's `x/0 = 0`; computed: `(m,n,u,v) = (1,2,0,0)` gives `FullSquares 1 2 0` false
   (`F4 = 5/2`) and `FullSquares 0 0 (1/2)` true; 24 such failures with entries ≤ 3 when `n = 0` or `v = 0`) and
   satisfiable.
   Instance with both sides true: `(m,n,u,v) = (4,3,1,1)` and `(1,1,4,3)` (computed exactly).

ExchangeCensus.lean:
7. `CoverDescent.mem_exceptions_of_fullSquares`: for `m, n, u, v ∈ ℕ` with `n ≠ 0`, `0 < v < u ≤ 1000`,
   `gcd(u, v) = 1` and `FullSquares m n (u/v)`: `(u, v) ∈ CoverDescent.exceptions` (= `exceptions1000` filtered by
   `∉ ratios`; 425 − 119 = 306 pairs, `exceptions_length`). Proof: exchange, then `no_solution_le_1000` and
   `no_fullSquares_of_mem_ratios` at the ratio `u : v` with `t = m/n`.
   **Not vacuous**: exhaustive search over all coprime `0 < v < u ≤ 1000` with `m:n ∈ {1:1, 0:1, 2:1, 3:1, 3:2, 4:3,
   5:4}` finds 306 instances, `(1,1,u,v)` for 179 pairs and `(0,1,u,v)` for 127 pairs, which are exactly the 306
   pairs of `exceptions` (each attained; all with `m/n ∈ {0, 1}`, i.e. exchanged trivial points, as NOTES.md says).
   So the hypotheses are satisfiable and the conclusion is sharp on the census list; the theorem says nothing new about
   non-trivial solutions beyond the census, as the docs say.
- [ok] Proofs use `ext/simp/ring/field_simp/linarith/decide_eq_true` and the census theorems; no `native_decide`,
  `sorry`, axiom or `implemented_by` in either file.
- [ok] Docs: LEDGER.md:42 (tier A row), README.md:172, G1_FULL_SYSTEM.md:52, experiments/exchange/NOTES.md state
  `FullSquares m n (u/v) ↔ FullSquares u v (m/n)` for naturals with `n, v ≠ 0` and the census corollary with
  `t = u/v > 1` only, `t ↦ −t`, `t ↦ 1/t` as prose: this matches the Lean.
- [nit] LEDGER.md:42 and G1_FULL_SYSTEM.md:52 state the corollary as "`u : v` coprime and `v < u ≤ 1000`" without
  `0 < v`; read literally, `(u, v) = (1, 0)` satisfies "coprime, v < u" and `FullSquares 1 1 (1/0)` = `FullSquares 1 1 0`
  holds, while `(1, 0) ∉ exceptions`. The Lean has `0 < v`; README's "`u/v > 1`" implies it. Add "`v ≥ 1`".
- [ok] Audit files consistent: 961 → 976 `theorem` lines at the start of a line in `MagicSquareOfSquares/*.lean`
  (14 in Exchange.lean, 1 in ExchangeCensus.lean) = scripts/declared_theorems.txt (976); scripts/audit_output.txt adds
  exactly the 15 new names (964 names + 2 lines = 966 lines); scripts/audit_internal.txt adds 7 `_simp_1_*`
  auxiliaries of `fullSquares_exchange` (710 → 717; summary line "964 printed, 717 auxiliary, 1681 in all");
  audit_line.txt 1659 → 1681 = 1659 + 15 + 7. `MagicSquareOfSquares.lean:43-44` imports both files, so `Audit.lean`
  (unchanged since 0.6.0) sees them.
- [fix, process] The Exchange files have been built only "in a container before CI, no log committed" (CHANGES.md,
  tenth run); LEDGER.md:42 says "checked by CI's build". No green `build` on a commit containing them is recorded in
  the repository (the run-10 commits' CI outcome is not recorded). The tier-A row should not ship until the `[full ci]`
  run of the release commit is green; consider committing that run's audit line or run URL.
- [fix] OUTPUTS.md:295 lists `conductor-t11/t11_theta_negative.txt` as class A in `research`, but CI never compares it:
  ci.yml:743 greps one count line (and fails, see the blocker). The `hm ... --negative` fix makes the row true.

## 5. REPRODUCE.md against ci.yml

- [ok] The job table has the new `research` row (about 30 minutes; limit 90, as ci.yml:700) and says the long jobs,
  `research` included, are skipped on ordinary pushes (ci.yml:698). Pins table: `research` installs
  `requirements.txt` + `sympy==1.14.0`, both listed; numpy, scipy, mpmath, which `conductor-t13/covers.py`,
  `search.py` and `rank_second.py` import, arrive transitively with passagemath (scipy 1.18.1 is required by
  passagemath-schemes; numpy by passagemath-symbolics/flint), so the job does not depend on this container's extras;
  their versions are those of `requirements-freeze.txt` only by luck of resolution (documented at REPRODUCE.md:25).
- [fix] REPRODUCE.md:49 "run the commands of its `run:` block from the root of the repository": `quick` (lines
  296-319), `descents` (546-688) and all of `research` (`W="$GITHUB_WORKSPACE"`, ci.yml:712) need
  `GITHUB_WORKSPACE` set to the root; unset, `research` stops at its first `d` with `cd: /experiments/...`. Add
  `export GITHUB_WORKSPACE=$PWD`. Also: the helpers write `/tmp/d.err`, `/tmp/hm.out`, so two reruns on one machine
  interfere (this machine has three panels doing it).
- [fix] REPRODUCE.md:63 lists `scan-5000` among the jobs that end with `git diff --exit-code`; it does not (above).
  OUTPUTS.md:6-7 gives the older list (no `research`, `square-section`) [nit].
- [nit] REPRODUCE.md:11 "the last tag is `v0.6.0`" will be stale once v0.7.0 is tagged. REPRODUCE.md:34, 41 (offline
  install) omit `passagemath-eclib==10.8.14`, which `four-points` installs (pre-existing).
- [ok] PYTHONPATH: only the note2 lines need it; ci.yml:719 exports it and `note2/NOTES.md:89` gives the relative form.
  With the `/tmp/t10` inserts removed (section 2) nothing outside the repository can shadow these imports.
- [ok] Files not committed: every `research` line P3 ran from a clean clone found its inputs (section 6). Scripts
  that need uncommitted files are class B/C and say so (`referee_z14/claims_z14.py`, the explore3-5 scripts with
  `/tmp/rn`, `/tmp/mix`, `note2/h2_check.py`). `referee_x6/samples/samples.sh:9` and `referee_z9/g3_samples.sh:11`
  hold a `/tmp/claude-0/...` path only as a `sed` pattern to rewrite the session's path to the repository's; fine.

## 6. Spot runs from the clean clone (2-core container shared with two other panels, load 1.3-1.8)

Each line run through the research helpers taken verbatim from ci.yml (private /tmp), from
`/home/claude/conductor/panel_p3/clone` at c857655, Python 3.13.18, cypari2 2.2.9 / PARI 2.17.4, sympy 1.14.0,
passagemath 10.8.14:

| CI line | Result | Wall time |
|---|---|---|
| quick, ci.yml:151 `open_pairs.py \| diff` | pass | 0.03 s |
| research :753 `hm $R/i1 i1_laws.py` | pass | 0.7 s |
| research :753 `hm $R/i1 i1_laws.py --negative` | pass | 0.4 s |
| research :753 `hm $R/i1 i1_negative.py` | pass | 1.8 s |
| research :754 `hm $R/i2 i2_residual.py --negative` | pass | 1.6 s |
| research :756 `hm $R/i4 i4_exchange.py --negative` | pass | 0.5 s |
| research :757 `hm $R/i5 i5_rank.py --negative` | pass | 1.0 s |
| research :758 `hm $R/i8 i8_pairct.py --compare50` | pass | 3.4 s |
| research :758 `hm $R/i8 i8_pairct.py --negative` | pass | 9.6 s |
| research :755 `hm $R/i3 i3_scan.py --negative` | pass | 15.0 s |
| research :759 z17 loop (i8count, i8compare, i2check) | pass | 3.6 s |
| research :760 z17_neither a b c e | pass | 4.6 s |
| research :735 `kcurve_exact` and `--negative` | pass | 1.2 s |
| research :731 `genus3_split` | pass | 3.6 s |
| research :746 `covers.py` | pass | 23.6 s |
| research :750 `certify.py` (ten JSON files) | pass | 2.8 s |
| research :728 z16 loop (four scripts) | pass | 24.1 s |
| research :721 note2 `two`, `two2` (with the PYTHONPATH export) | pass | 2.0 s |
| research :743 `t11_theta.py --negative \| grep -q` | **fail** (exit 120, EPIPE; 5 of 6 runs) | 5.5 s |
| proposed: `d $T search_E_alt_H2000_torsion.txt search.py --models=altmodels.json --curve=E --H=2000 --classes=torsion --skipbase` | pass | 3 s |

After the runs `git status --porcelain` in the clone: only ignored `__pycache__/` directories.

## Bottom line

**Ready after fixes.** One blocker, a one-line CI change: ci.yml:743 makes `research` fail on most runs (EPIPE from
`grep -q` on a live pipe), so the full green run the release rule needs cannot happen as committed. The helpers `d`
and `hm` fail closed in every tampering case tried. The Lean statements are what the docs say and are not vacuous: all
306 census pairs occur as instances, through the trivial points. The audit counts agree (976 / 964 + 717 = 1681).
The rest is documentation and coverage: one OUTPUTS row claims a file CI does not rerun (the proposed line passes in
3 s); one claims a file CI only greps; there are OUTPUTS gaps; `scan-5000` has no final `git diff`;
`GITHUB_WORKSPACE` is needed for hand reruns; two `/tmp/t10` sys.path inserts sit in CI-run scripts; and some tier-C/D
results lack mutated-rule tests.
