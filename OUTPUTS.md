# Committed outputs that CI did not regenerate

Written 25 November 2026 (UTC) by an AI model (Claude); not read by a person.

**What CI does.** The jobs of `.github/workflows/ci.yml` rerun programs and compare what they print with committed files
by `diff` (`REPRODUCE.md`, section 4). The jobs `build`, `quick`, `four-points`, `descents`, `last-pair` and `scan-15000`
end with `git diff --exit-code` (the first had it before), so a program that rewrites a committed file instead of
printing it fails the job. **What CI does not do.** It does not run Magma (its answers are transcripts of the public
calculator) or the long runs of class B below, and it cannot check the notes of class C that were written by hand. A
rerun shows that a program gives the committed output with the pinned versions; it does not show that the program
computes what the prose says.

**This list.** The 78 committed files under `scripts/` whose name did not occur anywhere in `ci.yml` at commit `6a54ef4`
(`rank_indep_all.txt`, also missing then, was added to the job `quick` separately), and two more that this search by
name missed: `g1_full/magma_check.txt`, because a file of the same name in `second_descent` is in `ci.yml` (no other
name is shared in this way), and `second_descent/magma_last_pair_raw.txt`, whose name occurs in `ci.yml` only in a
comment (no other file is named only there). Files that `ci.yml` names only as inputs of a command are not listed.
For each file: the class, the command (run from the directory of the file unless another is given), the CI job that
now regenerates it, and the time.

- **A**: regenerable in a few minutes; the CI job named reruns it and compares it byte for byte (after removing PARI's
  line "PARI stack size set to ..." or header lines that the program does not print), unless the row says otherwise.
- **B**: regenerable with the command shown, but not put under CI: too long (time measured or estimated), or, where the
  row says so, not reproducible byte for byte or assembled by hand from the program's output.
- **C**: a transcript of an external tool (the public Magma calculator), or a note of a run written by hand, or a summary
  whose program is not committed. CI cannot regenerate it; the checker that reads it, if any, is named.
- **D**: an input or hint file read by a program; what reads it is named.

Times were measured on 25 November 2026 on a shared 2-core container (load 3 to 11 from other jobs), with Python 3.13,
cypari2 2.2.9 (PARI 2.17.4), sympy 1.14.0, passagemath 10.8.14; GitHub's runner should be faster.

| File | Class | Command | CI job | Time |
|---|---|---|---|---|
| `cover_descent/structure/class_selmer_notes.txt` | C | hand-run with functions of `lean_cert.py` and `class_good.py` (the two short programs are not kept, as its header says); no checker | none |  |
| `cover_descent/structure/xcheck_15000.txt` | A | `(python3 class_selmer_linalg.py 5100 5100 --from 5000; python3 class_selmer_linalg.py 15000 15000 --from 14900)` | scan-15000 | 75 s |
| `cover_descent/tuple_bruteforce_15000_sample.txt` | B | `tail -n +2 tuple_bruteforce_15000_sample_in.txt \| xargs -P 2 -I{} sh -c 'echo "{} $(python3 tuple_bruteforce.py magic {} \| tail -1)"' \| LC_ALL=C sort` (the file is in that order) | none | about 25 s per ratio, 219 ratios: about 1.5 CPU-hours |
| `cover_descent/tuple_bruteforce_15000_sample_in.txt` | D | input of the line above; its lines after the header are `awk 'NR % 40 == 1' structure/class_selmer_15000_list.txt \| cut -d' ' -f1`, which CI now checks | quick (that check) | instant |
| `g1_full/dm4_cross.jsonl` | A | already regenerated: `dm4.py` on the 52 pairs of `dm4_pairs.txt`, compared by `check_evidence.py --dm4` in the fields that do not depend on the basis (not byte for byte: generators and timings can differ); `check_evidence.py` (job `build`) also checks the committed records | four-points (already) | part of that job |
| `g1_full/dm4_log.jsonl` | A | as `dm4_cross.jsonl` (the first 9 pairs of `dm4_pairs.txt`) | four-points (already) | part of that job |
| `g1_full/dm4_degenerate.txt` | A | `python3 dm4_degenerate.py 4,3 15,8 21,20 72,65 425,168 920,129`, against the file without its first line (the command) | four-points | 27 s |
| `g1_full/dm4_validate.txt` | A | `(python3 dm4_validate.py 3,2 2 7/3 11/5 101/17 1009/331 20011/7919; python3 dm4_validate.py 72,65 2 7/3 101/17 1009/331)` (the file names neither pair; found by rerunning) | four-points | 47 s |
| `g1_full/dm_closures.txt` | B | a table assembled from the runs of `dm.py`, `dm2.py`, `dm3.py` (Result 4) and `dm4.py`; no program prints it. `check_evidence.py` (job `build`) checks it against `dm_evidence.txt`, `dm_log.jsonl` and `open_pairs_1000.txt` | none | the runs: about 19 CPU-minutes recorded (sum of the `sec` fields of `dm_log.jsonl`) |
| `g1_full/dm_evidence.txt` | B | as `dm_closures.txt`: assembled from `dm.py - m,n`, `dm2.py m,n`, `dm3.py m,n` on the run pairs (needs passagemath-eclib); checked by `check_evidence.py` | none | as above |
| `g1_full/dm_log.jsonl` | B | the JSON lines of `dm.py`, `dm2.py`, `dm3.py`, one per line of `dm_evidence.txt` (a partner line repeats its run); not byte for byte (field `sec` is a timing); checked by `check_evidence.py` | none | about 19 CPU-minutes recorded |
| `g1_full/dm_sigma_24.jsonl` | A | `python3 dm_sigma.py --effort=4 --magma=magma_sigma2_raw.txt,magma_sigma3_raw.txt` on its 28 pairs in file order (the flags of the job `four-points`, which reruns 9 of them), compared without the field `sec` (a timing). **The committed file is stale in 3 of 56 records**: in two the field `upper` ends with " (magma_sigma2_raw.txt)", which the program no longer prints, and (609,200) E- has rank bounds [2, 4] where the rerun, which now reads Magma's bound 2 in `magma_sigma3_raw.txt`, gives [2, 2] (status "rank not 1" either way; `check_sigma.py` on the rerun reports this one record). So no `diff` line was added; regenerate the file first | none (stale; see the command) | 379 s here under load |
| `g1_full/dm_sigma_306.jsonl` | B | `python3 dm_sigma.py` on its 358 pairs in file order; the flags used are not recorded. A rerun without flags (two processes, 1059 s here) agrees with 641 of the 716 records apart from the timing `sec`; 26 differ only in the generator `G` (PARI can return `-G` or `G + T`; the job `four-points` compares those up to this with `sigma_exact.py --same-group`), 32 have a Magma bound in `magma_sigma2_raw.txt` or `magma_sigma3_raw.txt` (so the file was made with `--magma=...`), and 17 differ in PARI's rank bounds without one. Not byte for byte, so no `diff` line; the job `four-points` reruns 2 of the pairs and compares fields with `check_sigma.py` | none | 1059 s here in two processes (about 5 CPU-minutes recorded) |
| `g1_full/eclib_check.txt` | A | `python3 eclib_check.py`, against the file without its two comment lines; needs passagemath-eclib, which `requirements.txt` does not list, so the job installs `passagemath-eclib==10.8.14` | four-points | 69 s |
| `g1_full/genus2/frobenius.txt` | A | `python3 frobenius.py`, against the file without its first line | quick | 1 s |
| `g1_full/genus2/g2.txt` | C | answers of the public Magma calculator to `magma_g2.py` (network); no checker | none |  |
| `g1_full/genus2/g2_grh.txt` | C | answers to `magma_g2_grh.py` (Magma, network); no checker | none |  |
| `g1_full/genus2/g2_grh_raw.txt` | C | verbatim Magma requests and responses of `magma_g2_grh.py`; no checker | none |  |
| `g1_full/genus2/g2_raw.txt` | C | verbatim Magma requests and responses of `magma_g2.py`; no checker | none |  |
| `g1_full/height_certify.txt` | A | `python3 height_certify.py` | four-points | 29 s |
| `g1_full/height_precision.txt` | A | `python3 height_precision.py` | quick | 1 s |
| `g1_full/local_all_census.txt` | B | `python3 local_all_census.py S 4` for S = 0, 1, 2, 3 (each writes `local_all_census_S.json` in the current directory and prints a line `DONE ...`; the table of the file was assembled from those four lines by hand) | none | not measured (304,191 pairs, PARI) |
| `g1_full/magma_e0l.txt` | C | bounds read by hand from `magma_e0l_raw.txt` (Magma); `dm4.py --upper` cites it; no checker | none |  |
| `g1_full/magma_e0l_raw.txt` | C | verbatim Magma transcript; no checker | none |  |
| `g1_full/magma_open.txt` | C | bounds read from `magma_open_raw.txt` (Magma); no checker | none |  |
| `g1_full/magma_open_raw.txt` | C | verbatim Magma transcript; no checker | none |  |
| `g1_full/magma_points_raw.txt` | C | verbatim Magma requests for points on `E0L`; the point used for (740,539) is checked on the curve by `dm4.py`; no checker of the file | none |  |
| `g1_full/magma_raw.txt` | C | verbatim transcript of `magma_check.py` (Magma, network), for 42 of the 88 curves; no program reads it (`check_evidence.py` reads the parsed answers in `g1_full/magma_check.txt`) | none |  |
| `g1_full/magma_check.txt` | C | answers of the public Magma calculator to `magma_check.py` (network), parsed; verbatim for 42 of the 88 curves in `magma_raw.txt`; read by `check_evidence.py` (job `build`). Not among the 78: a file of the same name in `second_descent` is in `ci.yml` | build (as input) |  |
| `g1_full/magma_sigma.json` | D | the curves sent by `magma_sigma.py` (Magma) | none |  |
| `g1_full/magma_sigma2.json` | D | the curves sent by `magma_sigma2.py`; read by `rank_indep.py` (job `quick`) | quick (as input) |  |
| `g1_full/magma_sigma3.json` | D | read by `rank_indep.py` (job `quick`) | quick (as input) |  |
| `g1_full/magma_sigma_raw.txt` | C | verbatim Magma transcript written by `magma_sigma.py`; no checker | none |  |
| `g1_full/open_pairs_1000.txt` | D | the 157 pairs left by `local.py 1000` and `global.py` (no program prints the file in this form); read by `check_evidence.py` (job `build`) | build (as input) |  |
| `g1_full/rank_indep_control.txt` | A | `python3 rank_indep_control.py 30` (the default count is 40; the committed file has 30 curves, and the first 30 of a run with 40 are the same) | four-points | 41 s with 40 |
| `g1_full/rank_indep_control_census.txt` | A | `python3 rank_indep_control.py --census` | four-points | 41 s |
| `g1_full/rank_indep_hints.json` | D | hints read by `rank_indep.py` and `rank_indep_heights.py` (job `quick`); written only with `--write` | quick (as input) |  |
| `g1_full/requirements-freeze.txt` | D | informational `pip freeze`; no program reads it (`REPRODUCE.md` cites it) | none |  |
| `g1_full/sigma_generic_pari.txt` | A | `python3 sigma_generic.py --pari \| grep -v "PARI stack"` (a guide only; the first 66 lines are those of `sigma_generic.txt`) | four-points | 107 s |
| `gauss_descent/magma_twocover_bounded_raw.txt` | C | verbatim Magma transcript; read by `twocover_check.py --bounded` | quick (checker) |  |
| `gauss_descent/magma_twocover_raw.txt` | C | verbatim Magma transcript; read by `twocover_check.py` | quick, descents (checker) |  |
| `gauss_descent/real_points_681.txt` | A | `python3 real_points.py 200 $(cat pairs_notrivial_1000.txt)`, in the background while the job runs the rest (the job already ran the same program on the 157 ratios) | descents | 323 s here under load |
| `gauss_descent/scan_1000_5000.txt` | B | `python3 scan.py M0 M1 ../g1_full/modp_classes.txt` for (M0, M1) = (1000, 2000), (2000, 2600), (2600, 4000), (4000, 5000), with `sed 's/  \[[0-9]* s\]$//'` on both sides (the lines end with timings); the job `descents` reruns `scan.py 1000 1100` | none | 139 + 164 + 743 + 919 s recorded (about 33 CPU-minutes) |
| `gauss_descent/small_all_60.txt` | A | `python3 small_all.py 60` | descents | 159 s |
| `gauss_descent/subsets_34.txt` | A | `python3 subsets.py $(cat pairs34.txt)`, in the background while the job runs the rest | descents | 669 s (one core, under load) |
| `hill/ap_exact_2000.txt` | B | `gcc -O2 -Wall -Wextra -Werror -fopenmp -o ap_exact ap_exact.c; ./ap_exact 2000` | none | about 40 minutes on two cores (provenance file) |
| `hill/ap_exact_2000_provenance.txt` | C | a note written by hand on the run of `ap_exact_2000.txt` (compiler, flags, hash, time); no checker | none |  |
| `hill/ap_relaxed_1000.txt` | B | `python3 ap_relaxed.py 1000 0 0.29; python3 ap_relaxed.py 1000 0.29 1`, against the file without its first line (written by hand); the job `quick` runs `ap_relaxed.py 250` | none | the first part alone ran for more than 6 CPU-minutes here and was stopped; not finished |
| `perfect_triangle/pt_scan_300_summary.txt` | B | `python3 pt.py scan 1 300 0.867`, of which the file keeps the last line and the SHA-256 (the 23,755 lines are not committed) | none | about 16 minutes (the file) |
| `perfect_triangle/pt_tors.txt` | A | `python3 pt_tors.py 2>/dev/null \| grep -v "PARI stack"` | quick | 0.3 s |
| `search/seven_brute_1300.txt` | A | `gcc -O2 -Wall -Wextra -Werror -o /tmp/seven_brute seven_brute.c && /tmp/seven_brute 1300 \| LC_ALL=C sort -u`, against the file without its two comment lines; the job `build` runs it to 430 | descents | about 3 minutes (the file); 697 s here under load |
| `second_descent/exact_last_pair_hints.json` | D | hints read by `exact_last_pair.py` (job `quick`); written only with `--write` | quick (as input) |  |
| `second_descent/exact_last_pair_route.txt` | A | `python3 exact_last_pair_route.py 187,84 55826 model=...`, the same for `440,279 231518` and `420,341 120238`, with the minimal models (PARI's `ellminimalmodel`; for (440,279) the model of `saturation_440_279.txt`), because `weierstrass()` is not integral there; the program prints each line twice (as it goes, then all again), so CI removes the repeat with `awk '!seen[$0]++'`; it exits with 1 for the last two on purpose | last-pair | 66 s for the three |
| `second_descent/exact_last_pair_route_hints.json` | D | hints read by `exact_last_pair_route.py`; written only with `--write` | last-pair (as input) |  |
| `second_descent/indep_rank_440_279.txt` | A | `python3 indep_rank.py --ainvs 0,-491040,0,-13400146081,0 440,279 \| grep -v "PARI stack"` | quick | 0.6 s |
| `second_descent/kclose_open4.txt` | C | a summary written from a run of `python3 kclose.py 440,279 719,161 627,364 991,263 4,3`; the committed program prints JSON lines instead. A rerun (48 s) gives the same facts: 210 quartics and no closure for the first four, closures by the same six quartics for (4,3). No checker | none | 48 s (the run) |
| `second_descent/kdist_627_364.txt` | A | `python3 kdist.py 627,364` | descents | 22 s |
| `second_descent/lcheck.txt` | B | `python3 lcheck.py 1e14 second_123.jsonl` (the earlier L-value check with PARI's `lfun`, superseded by `lproof.py` and `lball.py`, which the job `descents` reruns; the limit 1e14 is inferred from the largest conductor listed) | none | not measured |
| `second_descent/lproof_440_279.txt` | A | `python3 lproof.py curve 0,-491040,0,-13400146081,0 '(440,279), curve y^2 = x(x - 719^2)(x + 161^2), the twist by -1 of the Frey curve:' \| grep -v "PARI stack"` | quick | 3 s |
| `second_descent/magma_last_raw.txt` | C | verbatim Magma transcript, quoted by hand in the text; no checker | none |  |
| `second_descent/magma_last_pair_raw.txt` | C | verbatim Magma transcript for LAST_PAIR.md, quoted by hand; no checker. Not among the 78: its name occurs in `ci.yml`, but only in a comment | none |  |
| `second_descent/magma_second_ranks_raw.txt` | C | verbatim Magma transcript of `magma_second_ranks.py`; no program reads it (`indep_rank.py` reads the requests, `magma_second_ranks.json`) | none |  |
| `second_descent/magma_second_raw.txt` | C | verbatim Magma transcript; read by `magma_check.py` | quick, descents (checker) |  |
| `second_descent/rank_effort_wide.txt` | B | `python3 rank_effort.py two_map_scope_wide.txt` (PARI `ellrank` with effort 3 and 6, up to 150 s per call; the text notes that one line of the file was not printed by the committed program and was deleted) | none | up to about an hour (27 calls of up to 150 s); not measured |
| `second_descent/rank_pairs_wide.txt` | B | `python3 rank_pairs.py two_map_scope_wide.txt` (effort 6 on 80 curves; "slow" by its docstring) | none | not measured |
| `second_descent/two_map_scope.txt` | A | `python3 two_map_scope.py pairs123.txt \| grep -v "PARI stack"` | last-pair | 37 s |
| `second_descent/two_maps_independent.txt` | A | `python3 two_maps_independent.py`, compared after sorting the dictionary of each line: the program iterates over a set of strings, so the order of the entries depends on Python's hash seed (the counts do not) | descents | 17 s |
| `second_descent/wide_2000_2500.jsonl` | B | `python3 wide.py 2000 2500 ../g1_full/modp_classes.txt` (the job `descents` reruns `wide.py 1000 1020`) | none | about 44 s per ratio here (the 5 ratios of the job `descents` took 218 s); with the partners some 2 hours here |
| `second_descent/wide_2000_2500_lvalue.jsonl` | B | `python3 lvalue_close.py LIMIT` on the open ratios of `wide_2000_2500.jsonl` (LIMIT is not recorded; the job `descents` uses `1e13` for `lvalue_open.jsonl`) | none | not measured |
| `second_descent/wide_2000_2500_lvalue_partners.jsonl` | B | the same on their partners | none | not measured |
| `second_descent/wide_2000_2500_summary.txt` | A | `python3 wide.py --summary wide_2000_2500.jsonl` | quick | instant |
| `second_descent/wide_2000_2500_tally.txt` | C | written by a program that is not committed; read by `wide_open.py` (job `quick`, the lines "not closed"). A reconstruction (not committed) gives all its lines from `wide_2000_2500.jsonl` and the two L-value files | quick (as input) |  |
| `second_descent/wide_2500_3000.jsonl` | B | `python3 wide.py 2500 3000 ../g1_full/modp_classes.txt` (the job `descents` reruns `wide.py 1000 1020`) | none | about 44 s per ratio here (the 5 ratios of the job `descents` took 218 s); with the partners some 2 hours here |
| `second_descent/wide_2500_3000_lvalue.jsonl` | B | `python3 lvalue_close.py LIMIT` on the open ratios of `wide_2500_3000.jsonl` (LIMIT is not recorded; the job `descents` uses `1e13` for `lvalue_open.jsonl`) | none | not measured |
| `second_descent/wide_2500_3000_lvalue_partners.jsonl` | B | the same on their partners | none | not measured |
| `second_descent/wide_2500_3000_summary.txt` | A | `python3 wide.py --summary wide_2500_3000.jsonl` | quick | instant |
| `second_descent/wide_2500_3000_tally.txt` | C | written by a program that is not committed; read by `wide_open.py` (job `quick`, the lines "not closed"). A reconstruction (not committed) gives all its lines from `wide_2500_3000.jsonl` and the two L-value files | quick (as input) |  |
| `second_descent/wide_lvalue.jsonl` | B | as `wide_2000_2500_lvalue.jsonl`, for `wide_1000_2000.jsonl` | none | not measured |
| `second_descent/wide_lvalue_partners.jsonl` | B | the same on the partners | none | not measured |
| `second_descent/wide_tally.txt` | C | as `wide_2000_2500_tally.txt`, for `wide_1000_2000.jsonl` | quick (as input) |  |

**Added on 26 November 2026 (UTC), after the reviews of `391e693`** (written by an AI model; not read by a person). The
outputs of `SELMER_SQUARE.md` and `SQUARE_SECTION.md` committed on 25 November 2026, which this list did not have (all
under `scripts/g1_full/`; the times are from a 2-core container under load). The `--negative` outputs of both, and
`selmer_square_rule_5000.txt`, were already in `quick`.

| File | Class | Command | CI job | Time |
|---|---|---|---|---|
| `selmer_square_5000.tsv` | A for its first 108 rows, B for the rest | `python3 selmer_square.py --max 5000 --write` writes it (with `--jobs 2` on two cores); `python3 selmer_square.py --max 300 --dump-tsv` prints the 108 rows of `m ≤ 300`, compared with the first 109 lines of the file (the rows are in the order of `m`) | quick (`m ≤ 300`) | 10 s (`m ≤ 300`); 26 minutes in full |
| `selmer_square_5000.json` | A for its first 108 records, B for the rest | as above; `--dump-json` prints the records of `m ≤ 300`, compared as data with the first 108 records | quick (`m ≤ 300`) | 10 s (`m ≤ 300`) |
| `selmer_square_5000.txt` | B | the summary of the full run above (5 lines); the counts of `m ≤ 5000` are not rerun | none | 26 minutes |
| `selmer_square_rule_5000.tsv` | A | `python3 selmer_square_rule.py --max 5000 --tsv-to PATH` (the same run as the `.txt`, which the job already compared) | quick | 85 s with the `.txt` |
| `selmer_square_rank2.txt` | A | `python3 selmer_square_rank2.py` (reads the hints and checks every point exactly; no PARI); `--negative` checked by its last line | quick | 1 s |
| `selmer_square_rank2_hints.json` | D | points found by PARI's `ellrank` (`selmer_square_rank2.py --find --write`, a finder only); read by `selmer_square_rank2.py` | quick (as input) |  |
| `selmer_square_pari.txt` | B | `python3 selmer_square_pari.py --write` (PARI's `ellrank` on 120 random curves, a guide only; its counts depend on PARI's internal choices, so it is not compared in CI) | none | not measured |
| `square_section.txt` | A up to the line `Part 2: ...`; B for the two control lines | `python3 square_section.py --mmax 300 --bound 60 --controls --write` writes it; the job runs it without `--controls` (control C2 counts the points PARI's `ellrank` returns) and compares the file up to the line `Part 2: ...` | square-section | about 5 minutes |
| `square_section.tsv` | A | the same run, `--tsv-to PATH` | square-section | (same run) |
| `square_section_hints.json` | D | points of class 2 found by PARI's `ellrank` (`square_section.py --find-hints --write`, a finder only); read and checked exactly by `square_section.py` | square-section (as input) |  |
| `review_t8_class2.txt` | A | `python3 review_t8_rank1.py --class2 300` | quick | 13 s |
| `review_t8_negative.txt` | A | `python3 review_t8_rank1.py --negative` | quick | 30 s |
| `review_t8_rank1.txt` | B | three runs of `review_t8_rank1.py --curves ... --brute --crosscheck` (the commands are its comment lines) and a note of a PARI run; the lines end with timings; one curve alone took 1365 s | none | about 35 minutes |
| `review_t8_alphaK.txt` | B | `review_t8_rank1.py --curves` on nine curves (its first line); the lines end with timings | none | about 5 minutes (the timings in the file add up to 299 s) |

**Counts.** Of the 80 files, 26 are of class A: 23 regenerated by the lines added here, 2 already regenerated before (the records of `dm4.py`), and 1 not put under CI because the committed file is stale (`dm_sigma_24.jsonl`); 21 are of class B, 24 of class C, 9 of class D.

**How the lines were checked.** Every line added to `ci.yml` was run here from the directory the job is in at that
point (each job's added lines together with all of its top-level `cd` lines, taken from the YAML), and ended with
status 0; then `git diff --exit-code` found no changed file. A forged byte in `pt_tors.txt` made its line fail and then
made `git diff --exit-code` fail; a changed count in `two_maps_independent.txt` and a changed bound in
`exact_last_pair_route.txt` were refused by their comparisons. The jobs were not run as a whole: for the lines that
were already there, the claim that none of them writes a committed file rests on reading the programs (they print to
standard output; the files they write need `--write`, or are named on the command line and are under `/tmp`).

**Not done.** The class B runs were not repeated, except `dm_sigma_306.jsonl` (described in its row), one ratio of
`tuple_bruteforce_15000_sample.txt`, and part of `ap_relaxed_1000.txt` (stopped). The class C summaries were not
checked against their transcripts, except the three wide tallies (by a reconstruction that is not committed).

**Added on 26 November 2026 (UTC), run 4** (written by an AI model; not read by a person). Outputs committed with
run 4 (and run 3's Cassels–Tate and Fisher files, committed with it). Times from a 2-core container.

| File | Class | Command | CI job | Time |
|---|---|---|---|---|
| `g1_full/rank_fisher_theta4.txt`, `_k4.txt`, `_k6.txt`, `_recheck.txt`, `_control.txt`, `_negative.txt` | A | `python3 rank_fisher_theta4[_k4,_k6,_recheck,_control].py`, `python3 rank_fisher_theta4.py --negative` | descents | about 70 s together |
| `g1_full/rank_fisher.txt`, `rank_fisher5.txt`, `rank_fisher5_local.txt`, `rank_second.txt` and their `_negative.txt` | A | `python3 <script>.py` and `--negative` | descents | about 2.5 minutes together |
| `g1_full/e1_rank0_check.txt`, `torsion_e1_check.txt` | A | `python3 <script>.py` | quick | 1 s each |
| `second_descent/kctp.txt`, `kctp_controls.txt`, `kctp_negative.txt`, `kctp_search.txt` (run 3) | B | `python3 kctp.py` (and `--controls`, `--negative`); Sage and PARI | none | not measured here |
| `second_descent/kctp_recheck.txt`, `_controls.txt`, `_negative.txt`; `kctp_isog.txt`; `kctp_isodescent.txt`, `_negative.txt`; `kctp_points.txt` | B | `python3 kctp_recheck.py` (`--controls`, `--negative`), `kctp_isog.py`, `kctp_isodescent.py` (`--negative`), `kctp_points.py`; Sage and PARI | none | one of them about 12 CPU-minutes |
| `second_descent/lball_exact53.txt` | A for lines 1–51, B for lines 52–53 | `python3 lball_exact53.py --maxN 3e10 --nopari` (lines 1–51); lines 52–53, (41,1) and (751,289), by `--apchunk` in 2–3 processes and `--bigmerge` (the commands are in the script's docstring) | quick (lines 1–3), descents (1–51) | about 25 minutes for 1–51 on one core; about 20 minutes each for 52, 53 on 2–3 cores |
| `second_descent/lball_exact53_negative.txt`, `_controls.txt` | A | `python3 lball_exact53.py --negative`, `--controls` | quick | about 25 s together |
| `second_descent/lball_exact53_aptest.txt` | A | `python3 lball_exact53.py --aptest` | descents | about 4 minutes |
| `second_descent/lball_exact53.jsonl`, `_timing.txt` | B | `python3 lball_exact53.py --merge <the run outputs> --write`; timings | none | |
| `experiments/` (every `.json`, `.jsonl`, `.txt`, `.log` under `census2000/`, `census-wide/`, `open-ratio/`) | B | the drivers `census-wide/run7k.sh`, `run9k.sh`, `pipe9k.py`, the scripts named in each write-up; paths as in `experiments/README.md`; PARI 2.19.1 | none | hours per range |
| `g1_full/census_lvalue_exact.txt`, `_negative.txt` | A for lines 1–11 and the refusals; B for lines 12–19 | `python3 census_lvalue_exact.py --pairs '…'` (lines 1–11), `--negative`, `--merge census_lvalue_exact.jsonl` (the isogeny lines); lines 12–19 by `--pairs` one at a time, (700,529) by `census_lvalue_big.py` in three chunks (commands in its docstring) | quick (lines 1–2, refusals, merge), descents (lines 1–11) | about 13 minutes for 1–11; 3 to 6 minutes each for 12–18; about 23 minutes for (700,529) on three cores |
| `g1_full/census_lvalue_identify.txt`, `.json`; `census_lvalue_exact.jsonl`, `_timing.txt`; `census_lvalue_big_700_529.jsonl` | B | `python3 census_lvalue_identify.py --write` (PARI; the `.json` holds timings); the records of the runs above | none | about 9 minutes for the identification |
