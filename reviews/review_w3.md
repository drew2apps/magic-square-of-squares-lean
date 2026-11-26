# Hostile referee review: run4_w3 (CensusClosure.lean)

Started 01:12 UTC. Reviewer has not seen the production of this work. Complete at 01:25 UTC.

## Verdict
**HOLDS.** (01:25 UTC.) `ratio_ne_4_3`, `ratio_ne_7_1`, `ratio_ne_12_5`, `ratio_ne_24_7`, `ratio_ne_17_7`,
`ratio_ne_31_17` say exactly: for integers a..i, T whose squares a² b² c² / d² e² f² / g² h² i² form a 3×3 magic
square with sum T (the eight fields of the main library's `IsMagic`, i.e. `IsMSS`), with b² ≠ d² and b² ≠ f², one has
(b+d)·n ≠ (h+f)·m. That is the census's ratio (README.md:78, G1Bridge.lean:123), and the two inequations follow from
distinct entries (conversely, with Fermat's four-squares theorem, they imply distinct entries). Every link of the chain
was checked by reading and re-derivation; the torsion lists, F_11 counts and τ cases agree with an independent PARI
computation; the generator reproduces the file byte for byte; one compile of a reviewer file importing the module
(EXIT 0) confirms the statements, standard axioms, and the sign and orientation corollaries. The only external
dependency is run 3's rank 0 (`finrank_m_n`, kernel-checked, not re-elaborated by me).

**FIX-FIRST: none.**

**FIX** (none affects the mathematics):
1. Integration (conductor; w3 was not allowed to do it): MANIFEST.sha256 lacks the 10 new files, so
   `scripts/verify_manifest.py` fails and the `quick` CI job will be red; the generator diff, Python check, Checks,
   Forge, Axioms, AxiomsForge are not in ci.yml yet (proposed in w3.md:109-118, :151).
2. Forge in CI: use the existing loop (`CensusClosureForge:13`, anchored `^scripts/…: error`) and pin the 13 error
   lines, instead of the bare `grep -c ': error'` of w3.md:114.
3. Wording: every public sentence must keep the hypothesis "b² ≠ d², b² ≠ f² (true for nine distinct entries)";
   degenerate magic squares of squares with genuine ratio 4:3, 12:5, 24:7 exist (τ = ±1), so "admit no magic square of
   squares" unqualified is false. w3.md:185-186 "so nothing is vacuous" is overstated: the hypotheses (minus the
   ratio) are satisfiable iff the open problem has a solution; what is shown is that each inequation is necessary.

**NIT:** (a) the |b+d| : |h+f| sign convention and the reversed ratio n:m are one-line corollaries (reviewer R1, R2
compile) but are not stated; (b) the necessity of b² ≠ d² (reviewer R3) is not among the checks; (c) docstring
CensusClosure.lean:31-33 names the partner pair the other way round from `ratio_ne_partner`'s variables;
(d) CensusClosureChecks.lean:13-15: τ = 0 is attained for all four pairs, τ = ±1 not for (7,1);
(e) four identical `instance : Fact (Nat.Prime 11)` (:356/:466/:560/:670); (f) pre-existing run-3 doc error
`Y = m³n³y` (RankZeroCensus.lean:10, gen_rank_zero_certs.py:135) should be `Y = m²n²y` (w3 found and proposed it);
(g) proposed Result-3 sentence (w3.md:93-94) "degenerate … is exactly t ∈ {0, ±1}" omits t = ∞ (b = d), which the
Lean hypothesis b² ≠ d² also excludes; (h) optional: docstring :35 could add the converse (the two inequations imply
nine distinct entries, by Fermat), so that the theorem reads as the distinct-entry case exactly.

**Not checked by me:** run 3's rank-0 certificates (trusted as kernel-checked); re-elaboration of CensusClosure.lean
itself (my compile imports its olean); the Forge/Checks/Axioms files were read, not run; the PARI scan
(census_closure_scan.*) and the census counts of G1_FULL_SYSTEM were not re-derived beyond membership of the six pairs
in `exceptions200` (G1FullCertificates.lean:127) and the PARI ranks of the six curves.

## Findings (written as found)

### 01:17 First pass over lean_rank/LeanRank/CensusClosure.lean (834 lines, read in full)

**P1. Hypotheses of `ratio_ne_m_n` in plain maths (CensusClosure.lean:448-453 and the five siblings).**
Integers a..i, T with a²+b²+c² = d²+e²+f² = g²+h²+i² = a²+d²+g² = b²+e²+h² = c²+f²+i² = a²+e²+i² = c²+e²+g² = T,
b² ≠ d², b² ≠ f²; conclusion (b+d)·n ≠ (h+f)·m. The square is
```
a² b² c²
d² e² f²
g² h² i²
```
with roots a..i (any signs). The eight equations are literally the eight fields of `IsMagic` applied to
(a², ..., i²) (Defs.lean:31-45: row₁..anti), i.e. `IsMSS a b c d e f g h i T`, in the same order. Checked field by field.
The ratio is the one of README.md:78 and G1Bridge.lean:22/123 ((b+d):(h+f) with edge roots, same letters), and of
G1_FULL_SYSTEM "The system" (b = mA+nB, d = mA−nB, h = nA−mB, f = nA+mB; (b+d):(h+f) = 2mA:2nA = m:n). Consistent.

Degenerate cases: h+f = 0 with (b+d)n = (h+f)m forces b = −d, so b² = d² (handled, :272-278); b = d gives b² = d²
(:279). So under the theorem's hypotheses τ = (h+f)/(b−d) is always defined and nonzero-denominator. Sound.

**Strength of b² ≠ d², b² ≠ f² (reviewer's own argument, not in the repo).** In the standard parametrisation of a
3×3 magic square (centre E, a = E+x, i = E−x, c = E+y, g = E−y, b = E−x−y, h = E+x+y, d = E−x+y, f = E+x−y), b² ≠ d²
⟺ y ≠ 0 and b² ≠ f² ⟺ x ≠ 0. Any coincidence of two of the nine entries is a linear relation y = λx with
λ ∈ {±1, ±2, ±1/2} (or x = 0, y = 0), and then the nine entries lie in an arithmetic progression of 5 or 7 terms with
nonzero difference; four squares in arithmetic progression are impossible (Fermat/Euler). So, for a magic square of
squares, **b² ≠ d² ∧ b² ≠ f² is equivalent to all nine entries distinct**: the theorem covers exactly the magic squares
of nine distinct squares. The file's docstring (:35) states only the easy direction ("hold when the nine entries are
distinct"), which is what is needed. Consequence for vacuity: the hypotheses minus the ratio are satisfiable iff a
magic square of nine distinct squares exists (the open problem). So "nothing is vacuous" in w3.md:186 is an
overstatement; what the checks file shows is that each of the two inequations is needed (see P5).

Sign conventions: README defines the ratio as |b+d| : |h+f|. The Lean statement excludes only (b+d)n = (h+f)m; the
case (b+d)n = −(h+f)m follows by applying the theorem to (−b, −d) (squares unchanged), and the reversed ratio n:m by
applying it to the 180°-rotated square (i h g / f e d / c b a sends (b,d,h,f) to (h,f,b,d); the two inequations become
h² ≠ f², h² ≠ d², again implied by distinctness). Neither corollary is stated in Lean. NIT (same convention as the
main library's `IsMSS.ratio_ne_of_le_200`, which also needs n < m).

### 01:20 P2-P4 (read + independent computation)

**P2 `mss_forms`/`mss_point` (:173-255).** Re-derived by hand: from c₂, r₂: b²+h² = d²+f²; from the lines, 2a² = f²+h²,
2i² = b²+d² (corner = mean of the two far edges). With (b+d)n = (h+f)m and h+f ≠ 0: m(b−d) = n(f−h). Put
A = (h+f)/2n, B = (b−d)/2n: b = mA+nB, d = mA−nB, h = nA−mB, f = nA+mB, a² = n²A²+m²B², i² = m²A²+n²B². With τ = A/B,
r0 = a/B, r8 = i/B: r0² = m²+n²τ² (= f0(w), w = τ²), r8² = n²+m²τ² (= f8(w)). X = m²n²τ² gives X+m⁴ = m²·f0,
X+n⁴ = n²·f8, so X(X+m⁴)(X+n⁴) = m⁴n⁴·w·f0·f8 = (m²n²·τ·r0·r8)². **So Y = m²n²·y with y = τ·r0·r8, y² = w·f0·f8**;
the file is consistent with this everywhere (:20, :28, :226, :235-236, and the constants 144, 49, 3600, 28224 = m²n² in
`tau_m_n`). `curve e₂ e₃` is y² = x(x−e₂)(x−e₃) (RankZeroGenCurve.lean:24, :49), so `curve (−m⁴) (−n⁴)` is the
intended curve. Degenerate m, n: n > 0 is a hypothesis; m = 0 or m = n make the curve singular but the statement is
`Equation`, not a point, so it stays true; no hidden ℕ subtraction. Sound for every m ≥ 0 and n > 0
(n > 0 is needed for the division by 2n). `mss_point` is not used by the headline chain: `ratio_ne_of_tau` calls `mss_forms` and the
`tau_m_n` lemmas rebuild the point themselves. Harmless.

The wrong normalisation is in run 3, not here: RankZeroCensus.lean:10 says `Y = m³n³y`; it must be `m²n²y`
(e.g. (4,3), w = 1: f0 = f8 = 25, y = 25, Y = 3600 = 144·25, not 1728·25). Comment only (no theorem uses it); w3
found it and proposes the fix (w3.md:21, :106). NIT (pre-existing).

**P3 torsion (:61-121, per-pair blocks).** Stoll's lemma (EllipticCurves/ReductionAtPrime.lean:365, rev 2164fb5,
read in /home/claude/ext/EllipticCurves) needs: an integral model W₀ with `hE : W₀.map = E`, good reduction as the
instance `[(redCurve v W₀).IsElliptic]`, `[E.IsElliptic]`, `[CharZero K]`, p prime, `(p) ∈ v`, `(p) ∉ v^(p−1)`
(ramification e ≤ p−2), P torsion, red P = 0 ⇒ P = 0. `red_injective` (:80-92) supplies exactly these: W₀ a b =
[0, a+b, 0, ab, 0] with `map_W₀`; good reduction from `redCurve_isElliptic` (11 ∤ 16(ab(a−b))², by `norm_num`;
I recomputed: 256,81,175 / 2401,1,2400 / 20736,625 / 331776,2401 are all prime to 11 mod 11);
`Ideal.mem_span_singleton_self`; hram = 11^10 ∤ 11. **11-torsion is not excluded separately and need not be**:
e = 1 ≤ p−2 = 9 makes the formal group E₁(ℚ₁₁) torsion-free, so reduction is injective on the whole torsion
subgroup, 11-part included (Silverman VII.3.1). Torsion of every point from rank 0 (`finrank_m_n`, run 3, same
`curve (−m⁴) (−n⁴)` terms, RankZeroCensus.lean:74/112/152/191) + Stoll's Mordell–Weil (`fg_point_of_numberField`) +
`Module.finrank_eq_zero_iff_isTorsion` (needs Module.Finite, supplied). Then #E(ℚ) ≤ #Ẽ(𝔽₁₁) = N (decide +kernel),
and N distinct listed points ⇒ E(ℚ) = the list (`bijective_of_nat_card_le`). Correct shape.

Independent computation (PARI 2.17 via cypari2, /tmp/rw3/check.py, /tmp/rw3/pts.py): torsion Z/8×Z/2, Z/4×Z/2,
Z/8×Z/2, Z/8×Z/2 and ellrank [0,0] for (4,3), (7,1), (12,5), (24,7); #Ẽ(𝔽₁₁) = 16, 8, 16, 16 by brute force; the
torsion group generated in PARI equals the Lean list `pts_m_n` as a set, for all four. Abscissae /m²n² (= w):
negative, 0, 1, or 1/6, 6; 2/15, 15/2; 3/28, 28/3 (non-squares): matches every `tau_nonsq` bracket (Q, N, k) used
(:432-443, :636-647, :746-757). For (7,1) w = 1 is on the curve but τ = ±1 does not make f0 = f8 = 50 squares; the
lemma's conclusion {0, ±1} is a safe over-approximation. Generator rerun (python3 -I, /tmp/rw3/gen): output
byte-identical to the committed file (sha256 1cf37183…); `census_closure_check.py`: ALL OK, 20/20 mutations refused.

**P4 partner (:316-349, :777-805).** Symmetry: reflection in the middle column, then root h ↦ −h:
c b a / f e d / i (−h) g. It preserves the eight line equations (rows reversed, columns 1↔3, diagonals swapped; signs
irrelevant for squares). New edges (b', d', h', f') = (b, f, −h, d), so (b'+d') : (h'+f') = (b+f) : (d−h) =
(m+n)(A+B) : (m−n)(A+B) = (m+n) : (m−n); A+B = 0 is harmless (both sides 0; the contradiction then comes from the
partner theorem's own inequations). b'² ≠ d'² is b² ≠ f², b'² ≠ f'² is b² ≠ d², passed crosswise (:349 `hbf hbd`).
17:7 → 24:10 = 12:5 (hMN: 24·5 = 10·12), 31:17 → 48:14 = 24:7 (48·7 = 14·24); both targets are proved directly.
PARI: (17,7) has rank 0 on its own w·f0·f8, (31,17) rank 1 — so the partner route is genuinely needed for 31:17.
Sound. Docstring :31 names the partner M:N = (m+n):(m−n) of a listed pair, while the theorem's M:N is the listed
pair and m:n the new one; the relation is an involution up to scaling, so this is only confusing. NIT.

### 01:23 One compile (allowed): scripts/ReviewW3.lean in the w3 slice, `lean4.sh w3 check`, EXIT 0, 21 s, 290 MB
Preconditions checked: 0 lean processes in `rhreplay`; slice and container copies of CensusClosure.lean have sha256
1cf37183… (= committed = regenerated); container olean 01:08:46 is newer than the source copy 01:07:58.
The file (imports LeanRank.CensusClosure; reviewer-written, not in the repo):
- R1 `(b+d)·3 ≠ −(h+f)·4` from `ratio_ne_4_3` applied to roots (−b, −d): compiles. So the README's |b+d| : |h+f|
  convention is covered (by a one-line corollary not in the repo).
- R2 `(b+d)·4 ≠ (h+f)·3` (ratio 3:4) with h² ≠ f², h² ≠ d², from `ratio_ne_4_3` on the half-turn i h g / f e d / c b a:
  compiles. So reversed orientation is covered.
- R3 the square 49 1 25 / 1 25 49 / 25 49 1 with roots 7 1 5 / −1 5 7 / 5 −7 1 satisfies the eight equations, the
  4:3 relation (0 = 0), b² ≠ f², and b² = d²: compiles. So `b² ≠ d²` is necessary too (the repo's checks only show
  `b² ≠ f²` is necessary).
- `#check` output of `ratio_ne_4_3`, `ratio_ne_17_7`, `ratio_ne_31_17` matches the plain statement in P1 exactly
  (∀ {a..i T : ℤ}, eight equations → b² ≠ d² → b² ≠ f² → (b+d)·n ≠ (h+f)·m).
- `#print axioms` of the six `ratio_ne_*`, the four `points_*`, `mss_point`: `[propext, Classical.choice, Quot.sound]`.
  (Caveat: an import trusts the olean; it does not re-run the kernel on CensusClosure itself. The olean was built by
  `lake build` from this exact source per the timestamps.)

### 01:24 P5 vacuity, forge, integration

**Vacuity.** The hypotheses minus the ratio are satisfiable iff a magic square of nine distinct squares exists (P1),
so non-vacuity cannot be shown and is not the right test. The right test is that each inequation is necessary:
`b² ≠ f²` by the repo's degenerate square at 4:3 (Checks file, Forge F11), `b² ≠ d²` by my R3 (not in the repo). All
six ratios are trivial-point pairs (m²+n² = 25, 2·25, 169, 2·169, 625, 2·625), so degenerate squares exist at or
through every one of them; the theorems are false without the inequations, and correct with them.

**Forge (scripts/CensusClosureForge.lean, 13 claims), read, not run.** F1-F3 counts/model, F4 bad prime (the
`Fact (Nat.Prime 5)` instance exists globally from RankZeroPoints.lean:62, so F4 fails for the intended reason,
norm_num on 5 | 16(…)²), F5 the ramification condition at 2, F6 point off curve, F7 duplicate, F8 missing affine point,
F9 bound 8, F10 square passed as non-square, F11 degenerate square (needs b² ≠ f²), F12 τ = −1 dropped, F13 wrong
partner. These test the certificate layer well. Not tested: the necessity of b² ≠ d² (R3); a wrong abscissa
normalisation in `mss_point` (e.g. X = mnτ²); a curve without rank 0 (cannot be forged without `sorry`). The proposed
CI line (w3.md:114) only counts `: error` lines; the existing CI loop (ci.yml:105-109) also anchors the file name, and
RankForge pins line numbers (ci.yml:110). FIX: add `CensusClosureForge:13` to that loop and pin the 13 lines.

**Integration (not a defect of w3: BRIEF_COMMON4.md:9-11 forbids it to edit MANIFEST/CI).**
`python3 -I scripts/verify_manifest.py` on the worktree: MANIFEST FAILED, 10 files present but not listed (all the
new files). The `quick` job (ci.yml, first step) will be red on run4_w3 until the conductor runs `--write`.
The generator diff, Python check, Checks, Forge, Axioms and AxiomsForge are not in CI yet (w3.md:109-118, :151 propose
the lines). The rank job's `lake build` picks up CensusClosure via the glob; its axiom-line threshold (≥ 76) passes.

**Docs.** The checks docstring (CensusClosureChecks.lean:13-15) says τ = 0 is attained "for (7,1)"; it is attained for
all four pairs (f0 = m², f8 = n²), and τ = ±1 is attained for (4,3), (12,5), (24,7) but not for (7,1) (f0 = f8 = 50).
NIT. Four identical `instance : Fact (Nat.Prime 11)` (:356, :466, :560, :670): harmless generator artefact. NIT.
w3.md:185-186 "Hypotheses are satisfiable …, so nothing is vacuous": overstated (see above). FIX (report wording).

### 01:25 Close
Reviewer artefacts outside the repo: /tmp/rw3/check.py, /tmp/rw3/pts.py (PARI), /tmp/rw3/gen/cc.lean (regenerated
file), and the reviewer compile file /Users/drew/rh-replay/mss/slices/w3/lean_rank/scripts/ReviewW3.lean (Mac slice,
not in the git worktree; can be deleted). The worktree /home/claude/conductor/w3 was not modified (git status clean).
Extra worked example for NIT (b): at 12:5 the degenerate square with roots 13 17 7 / 7 13 17 / 17 −7 13
(entries 169 289 49 / 49 169 289 / 289 49 169, sum 507) has (b+d):(h+f) = 24:10 = 12:5 with h+f ≠ 0, b ≠ d, b = f.
