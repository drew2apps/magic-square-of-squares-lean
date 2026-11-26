# review_t1 — hostile referee of `sess_t1` (`lball_exact.py`, exact inputs for L(E,1) of the (440,279)/(719,161) curve)

Written by an AI model (reviewer, did not write the code); not read by a person. Worktree
`/home/claude/conductor/t1`, head 68f4436, base 8214be3. Scratch under `/tmp/rev_t1/`. Worktree not edited.

## Log

- Read `lball_exact.py` (792 lines) in full, the six committed outputs, `lball.py`, SECOND_DESCENT.md section 6, t1.md.
- `--all` (no `--write`), run in the worktree: output = `cat lball_exact.txt lball_exact_negative.txt` byte for byte
  (diff empty); 76 s wall, 20 s CPU. `git status --porcelain --ignored` identical before and after (only the
  pre-existing ignored `__pycache__/`); mtimes of the committed .txt unchanged. Every `open(..., 'w')` in the file is
  behind `if write:` (lines 614, 657, 717, 758, 782).
- `--census second_123.jsonl` (8 s), `--rowtest 3000` (1 s), `--tatetest 4000` (12 s): each byte-identical to the
  committed file; git status unchanged afterwards. Census: 89 rows, 17 `w = None` (all "not implemented at [2]"),
  72 determined, all agree with PARI; every Delta factored.
- `--nopari` under an import blocker (`/tmp/rev_t1/nopari.py`: cypari2, sage, flint, mpmath, sympy, gmpy2 refused at
  import): completes, exit 0, no such module loaded, and its output equals the first 24 lines of `lball_exact.txt`
  (everything but the PARI section). So the certified path is Python integers + numpy int64 only (sympy is not even
  reached: all prime factors of Delta are < 10^3).

### 1. Conductor (recomputed by hand and by three other programs)
- By hand (`/tmp/rev_t1/code/hand.py`, no code of the author): A = 719², B = 161², A + B = 2·521²,
  Delta = 16A²B²(A+B)² = 2⁶·7⁴·23⁴·521⁴·719⁴; v₂(c4, c6, Delta) = (4, 10, 6), c4' ≡ 3 (16), v₂(j) = 6 (pot. good).
  Minimal at 2: v₂(Delta) = 6 < 12, and Tate's run uses only x → x+1 (u = 1, Delta unchanged). Tate at 2:
  singular point of y² = x³ + x mod 2 is (1,0); after x → x+1, v₂(a4, a6) = (1, 6), v₂(b2) = 2, v₂(b8) = 2 < 3:
  additive, not II (4 | a6), type III, f₂ = v(Delta) − 1 = 5, c₂ = 2 (Ogg: 6 + 1 − 2). Odd p: v_p(c4) = 0, I4;
  split iff −A (p = 7, 23: −1 non-square, non-split), A (521), B (719) is a square: split. N = 2⁵·7·23·521·719 =
  1,929,934,048. Agrees with the script's trace line for line.
- Sage (passagemath 10.8.14): `conductor(algorithm='all')` (PARI, eclib/mwrank, Sage's generic Tate agree) =
  1929934048; generic local data III/5/2, I4 ×4 with c = 2,2,4,4. `is_minimal()` is False only because the model is
  not the *reduced* minimal model; `minimal_model()` has the same Delta, `is_p_minimal(2)` True. Not an issue.

### 2. L-series, tail, interval arithmetic
- Formula re-derived: Lambda(s) = ∫₀^∞ f(it/√N) t^{s−1} dt = (√N/2π)^s Γ(s) L(s); with f(i/(t√N)) = w t² f(it/√N)
  (below), splitting at t = y gives L(E,1) = Σ aₙ/n (e^{−2πny/√N} + w e^{−2πn/(y√N)}); at y = 1, w = +1:
  2 Σ aₙ/n e^{−2πn/√N}. Code `evaluate` (517–542) matches, including the y = 6/5 / 5/6 pairing.
- Tail: |a_{p^k}| ≤ (k+1)p^{k/2} (Hasse) at good p, ≤ 1 at bad p; d(n) ≤ 2√n (pair divisors), so |aₙ| ≤ 2n.
  Σ_{n>M}|aₙ|/n qⁿ ≤ 2q^{M+1}/(1−q); Σ_{n>M}|aₙ|qⁿ ≤ 2q^{M+1}((M+1) − Mq)/(1−q)². Code 487–491 uses an upper
  bound for q^{M+1}, a lower bound for 1 − q and (M+1) − M·q_lo: directions correct. Valid upper bound.
- Floats only in `terms_needed` (choice of M) and never in an end point; `lower_str`/`upper_str` round floor/ceil
  in integers. √N: isqrt(N·2^{2P}) = s with s ≤ √N·2^P < s+1, used outward in `rate`. π: Machin, per-term floor
  error and truncation bounded, 20 guard bits. exp(−t): alternating partial sums ending on a negative (lower, at
  t_hi) / positive (upper, at t_lo) term, valid since t < 1 makes the terms decrease.
- Primitives tested against mpmath at 400 bits (`prim_check.py`): π enclosed (width 1 ulp); `exp_neg` on 3000 random
  t ∈ [0,1): 0 failures (excess width ≤ 60 ulp at 2^−192); `rate` and q enclosed for y = 1, 6/5, 5/6; `series` on
  random coefficients of both signs encloses the true finite sum.
- Independent enclosure (`iv_check.py`, mpmath `iv` outward rounding, 130 bits, my own tail, a_n from PARI):
  y = 1: L(E,1) ∈ [1.1049249737488, 1.1049249743032] (tail 1.39e−10 per sum); y = 6/5:
  [1.1049249739169, 1.1049249741351]; F ∈ [3130.0553606, 3130.0553633]. Same as the script's
  ([1.1049249737, 1.1049249744] / [1.1049249739, 1.1049249742] / [3130.055, 3130.056]).
- a_p: the script's own table for all 26,175 primes ≤ 302,194 (= π(302194)) equals PARI `ellap` on every prime, and
  aₙ = PARI `ellan` for all n ≤ 302,194 (`ap_check.py`); the script itself only cross-checks n ≤ 20,000.

### 3. Root number via the Fricke involution
- Precisely: f the newform of level N (modularity + Carayol), trivial character, rational coefficients, so
  f|W_N = εf, ε = ±1 (Atkin–Lehner), i.e. f(−1/(Nτ)) = εNτ²f(τ). At τ = it/√N: f(i/(t√N)) = −ε t² f(it/√N);
  with w := −ε this gives Lambda(s) = wLambda(2 − s). At t = 1 (τ = i/√N, the fixed point of W_N):
  F = f(i/√N) = wF, so w = −1 ⇒ F = 0. F enclosed in [3130.0553606, 3130.0553633] (my enclosure; the tail is
  1.3e−6) excludes 0 by a factor 10⁹: w = +1 for the sign of the functional equation of L(f,s) = L(E,s), which is
  exactly what the L-formula needs. Consistent with the local route: w = −w₂·(+1)(+1)(−1)(−1) ⇒ w₂ = −1.
- Caveat (not a bug): F ≠ 0 has no power against a wrong N (F is nonzero for N/2 too); the route needs N exact,
  which comes from Tate (confirmed four ways above). The docstring says so (line 37).
- Rizzo's table: not load-bearing (v₂(c6) = 10, both sign readings give w₂ = −1; and the Fricke route alone
  suffices). arXiv:1612.03095 (footnote 13) lists misprints in Rizzo's Table III — second line should read
  (0, 0, ≥0), and (2, 3, 1) has type I₂* — neither is the row used here. Not read from the printed table by me.

### 4. Negative tests
- The ten refusals are genuine: each forged pair of enclosures is as narrow as the true one and disjoint by
  ≥ 10⁻⁶ (a_719) to 0.16 (w = −1); the true inputs pass at both tail levels (control). I recomputed the a_719 = −1
  case independently with mpmath iv: [1.1024691552, 1.1024691558] vs [1.1024680428, 1.1024680430] — same as the
  committed line. The guard refuses because the functional equation fails, not because of a tail or assert.
- Reach of the guard (`neg_probe.py`): ±2 or sign changes of a_p are refused at p = 101, 1009, 10007, 30011; at
  p = 100003, a_p + 2 is NOT refused (sign flip is). So the guard does not protect a_p beyond ~3·10⁴; there the
  a_p rest on the point count (now checked against PARI for every p ≤ 302,194). As stated by the author.

### 5. PARI in the certified path
- None: PARI is imported only in `pari_crosscheck` (562, after the asserts at 605–606), `local_only`, `census`,
  `tatetest`, `rowtest`; the blocked-import run proves it. Indirect PARI influence: the sign reading of the Rizzo row
  was chosen by agreement with PARI (irrelevant for this curve); the hard-coded `assert N == 1929934048` (606) is a
  check against the old PARI value, not an input.

### Findings
FIX-FIRST (merge, not maths)
- MANIFEST.sha256 does not list the six new files: `python3 scripts/verify_manifest.py` at 68f4436 prints
  "MANIFEST FAILED: 473 files listed; git: 479 tracked files", so the CI `quick` job fails at its first check. The
  author flagged it; the hashes in t1.md match `sha256sum` of the files.
- `lball_exact.py` is not in `.github/workflows/ci.yml`: no CI job regenerates its outputs. A doc must not say CI
  reruns it until the proposed lines are added (numpy==2.5.5 is already pinned elsewhere in ci.yml).
NIT
- lball_exact.py:51–52: "Kolyvagin ... and modularity (or by Kato's theorem)": Kato also needs modularity; write
  "modularity, then Kolyvagin (with Gross–Zagier, BFH / Murty–Murty) or Kato".
- lball_exact.py:34–37, 55: the trusted list should name Atkin–Lehner (f is a W_N-eigenform, ε = ±1, w = −ε).
- lball_exact.py:605: the certified run also asserts `w_local == 1`, so it depends (as an assert) on the Rizzo row,
  which the proof does not need; say in the output that the Fricke sign is the proof and the table a corroboration.
- lball_exact.py:562/577: cross-check aₙ only to 20,000; ellan to Mmax costs ~1 s (all agree, checked here).
- lball_exact.py:616–630 / t1.md "ten forged inputs refused": say the forged a_p are at p ≤ 719 and the guard stops
  seeing a ±2 change of a_p near p ≈ 10⁵.
- lball_exact.py:101–103: sympy is called on the cofactor even when trial division has already shown it prime
  (d² > m); harmless (finder only).
- lball_exact.py:465: the ±4 ulp margin in `exp_neg` is unexplained (harmless; the partial sums already bracket).

### 6. What a doc may claim
- Holds: N, aₙ (p ≤ 302,194), w = +1 and the enclosure, all computed in exact/integer interval arithmetic with no
  PARI in the certified path; independently reproduced here (hand Tate, Sage/eclib, PARI ellap/ellan, mpmath iv).
- Still trusted: modularity (Wiles; BCDT), Carayol (level = conductor), Atkin–Lehner, Kolyvagin + Gross–Zagier +
  BFH/Murty–Murty (or Kato); Python integers and numpy int64; one implementation of Tate's algorithm (cross-checked).
  For the *closure* of (440,279) and (719,161), all of the second-descent chain other than step 7: that this E is
  the right curve (`second_exact.py`), the torsion and equation step, the class step, the partner map — not reviewed
  here. Rank 0 also has the separate algebraic route (Magma's 2-power isogeny descent, transcript only).
- One sentence a doc may claim: "For E: y² = x(x − 719²)(x + 161²), `lball_exact.py` computes without PARI the
  conductor 1,929,934,048 = 2⁵·7·23·521·719 (Tate's algorithm), the aₙ (point counts) and the sign +1 of the
  functional equation (f(i/√N) ≠ 0), and encloses L(E,1) in [1.104924973, 1.104924975] in integer interval
  arithmetic; with modularity, Carayol and Kolyvagin (or Kato), E(ℚ) has rank 0 — one program, written by an AI
  model, reproduced by an independent review, not read by a person."
