# Referee report: fix round 3 (repo2, 8c0a0a3..HEAD = 0cea08a + 8784fbc)

Status: complete

## Checks done (all from /home/claude/review/repo2 at 8784fbc, read-only)

- `python3 scripts/degenerate_squares.py 1000 | diff - scripts/degenerate_squares.txt`: passes (0.2 s).
- `python3 scripts/verify_manifest.py`: MANIFEST OK, 433 listed, 433 tracked.
- Mutation test (copies in /tmp/ref3): 6 mutations (swapped corner type 1, sign of h type 1, roots m<->n type 2,
  corners swapped type 2, v changed type 2, all-equal square type 2) each give 179 or 127 FAIL lines; the CI diff
  would fail. Independent recount: 179 + 127 = 306 coprime n < m <= 1000.
- Explicit squares of BARRIER.md checked by hand against PAPER.md Theorem 1 (b+d = mu, h+f = nu, b-d = nv,
  f-h = mv, 8e^2 = S(u^2+v^2)) and the grid a b c / d e f / g h i: type 1 (u = v = 2) all four equations and the
  centre equation hold, b = f, d = -h, edge ratio 2m : 2n; type 2 (u = 0, v = 2) the equations hold, b = -d,
  f = -h; corners agree with 2a^2 = f^2 + h^2, 2c^2 = d^2 + h^2, 2g^2 = b^2 + f^2, 2i^2 = b^2 + d^2; every line
  sums to 3k^2.
- CI line: job `quick`, repo root (before `cd scripts/gauss_descent`), under `set -e -o pipefail`. Correct.
- PAPER.md: now byte-identical to v0.5.0 apart from the added header; 313 = count of `^theorem |^private theorem `
  at v0.5.0 (and scripts/declared_theorems.txt there); HEAD has 568. "24 open" matches v0.5.0 section E.
- G1Sigma.lean: diff 8c05df4..HEAD is inside the module docstring only; no statement changed. The math of the new
  sentence is right (lowest terms suffices; involution gives the reverse inequality).
- Published tree: `git ls-remote` shows main = dcdb895; cloned to /tmp/ref3/pub (431 files, no .pyc, .gitignore
  `.lake/` only). HEAD tree = dcdb895 tree + exactly 14 changed/added files; the three magma_sigma files are
  byte-identical to the published ones; overlaying HEAD on dcdb895 gives MANIFEST OK (433/433). NOTE: a patch made
  with `git diff 8c05df4 HEAD` does NOT apply to dcdb895 (MANIFEST.sha256 hunk fails: 8c05df4's manifest is not
  the published one). The patch must be made against dcdb895 (or by copying the 14 files). I did not see the
  actual patch file.
- greps: "558" only in CHANGES.md:62 (historical entry, correct then); "M1" only in SECOND_DESCENT.md:90 and
  scripts (unrelated factor names); "25 November": none; "not uniformly": none.

## Findings

F1 (date, fix 5 is wrong by the repo's own rule). LEDGER.md:3 says "Dates in this repository are in UTC". BARRIER
was written at 2026-11-25 02:47 UTC (/tmp/barrier.patch header) and 0cea08a is 2026-11-25 04:23 UTC, so by that rule
the original "25 November" was right and BARRIER.md:3 "24 November 2026" is wrong. Same pre-existing drift in the
published text: 78659d1 (00:52 UTC 25 Nov), e58a514 (02:36 UTC), dcdb895 (03:27 UTC) are all labelled
"24 November 2026" (CHANGES.md:37-40, "9 open since 24 November 2026" in PAPER.md:5, README, LEDGER). Fix: either
change LEDGER.md:3 to state the convention actually used (owner's local date, UTC-8, from 24 Nov on) or re-date;
do not leave both.

F2 (provenance). BARRIER.md:3 "after an outside review": the repo always says whether a review was by a person or
an AI model (e.g. CHANGES.md:38 "neither is a reading by a person"); here it does not, and removing "(M1)" left the
review unidentifiable. Say "after a review by an AI model (not a person)" or whatever it was.

F3 (dangling reference). BARRIER.md:29 "the 306 pairs of Theorem 7" and BARRIER.md:46 "Theorem 7 removes them":
the file cites "Theorem 1 of PAPER.md" two lines before, but PAPER.md has Theorems 1-6 only; Theorem 7 is
CLASS.md:43 (`IsMSS.ratio_norm_sq_or_two_sq_of_le_1000`). Write "Theorem 7 of CLASS.md".

F4 (type 2 vs Theorem 1's hypothesis; minor). For type 2, b + d = h + f = 0, which Theorem 1 of PAPER.md excludes
(h + f != 0) and Theorem 7 excludes too. The equations hold, and the point does lie on C(m : n) of
G1_FULL_SYSTEM.md:57-59 (n(b+d) = m(h+f) and m(b-d) = n(f-h)) and on the five-form curve at t = 0, so the
Consequence is right; but BARRIER.md:24-25 says "they satisfy ... of Theorem 1" and BARRIER.md:31 "the variety of
solutions" is never defined. Say which variety (C(m : n), equivalently the five forms in t) and that for type 2 the
square is the point t = 0, tied to m : n by b - d = 2n, f - h = 2m (as the script docstring already says).
BARRIER.md:37 "the three edge ratios" of the type-2 square: the middle-line ratio is 0 : 0 there.

F5 (nit). G1Sigma.lean:20 "H(X') <= (|A| + 4)·H(X)": `A` is not defined in this header (the map is written with
c1); write |c1| (an integer, as `moebius_height` takes integers).

Notes (not issues): the full-name internal-audit diff (ci.yml:65, from 8c0a0a3) has not been dry-run here (no Mathlib
build available); audit_internal.txt was regenerated in c0024ed and no Lean statement changed since, so it should
pass. second_exact.py | diff and sigma_twice.py | diff both pass locally (sympy 1.14.0). Publishing: the change
touches Lean, CI and descent code, so the commit needs [full ci]; build the patch against dcdb895, not 8c05df4.

## Verdict: FIX-FIRST (docs only: F1, F3 needed; F2, F4, F5 cheap). Fixes 1-4 and 6-9 are correct; fix 5 is half
right (M1 removed) and half wrong by LEDGER.md:3.
