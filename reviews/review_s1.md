# Review s1: ClassHensel.lean (commit ad737f3), started Wed Nov 25 05:45:37 UTC 2026

## Scope and method
Source review only (no Lean toolchain). Read ClassHensel.lean in full, the definitions it uses
(CoverDescent: nrm/P0/P8/PL/PU/Q/gm/cj/cjIf; ClassLocal: GSqM, ClassAt, exists_sq_modEq_pow;
ClassTable: Passes, Row, passes_of_table), COMPLETE.md s.5. Ran scripts/class_hensel_witness.py
(python3 -I, exit 0, "all checks passed"; its output matches the 2x20 integers in the Lean file
exactly) and an independent checker of my own (/tmp/s1rev/indep_check.py, does not import the
builder's script).

## 1. Statements (plain mathematics)
- classAt_norms_all (L177): p odd prime, K >= 5, ClassAt(p^K, p^3, d,m,n,u,v) => for every k there
  are w_j with N_j == w_j^2 (mod p^k) and p^3 not| w_j (j = 0,8,L,U). No hypothesis on d. STANDS.
- gsqM_lift (L461) (+ _split L280, _inert L378): p odd prime, K >= 18 (split) / K >= 10 (inert),
  p^2 not| d, p^17 not| N(q), d*q a Gaussian square mod p^K => d*q a Gaussian square mod p^k for
  every k. STANDS.
- classAt_all_of_classAt (L475): p odd prime, K >= 18, p^2 not| d, ClassAt(p^K, p^3, ...) =>
  ClassAt(p^k, p^3, ...) for every k, same (u,v). Same q = p^3 on both sides, same M = p^power
  convention as ClassAt. Primitivity not assumed (not needed). STANDS.
- passes_all_of_passes (L488): Passes p 17 d m n (= ClassAt at p^(17+1), primitive point) =>
  Passes p k d m n for all k, at the same point (k+1 used correctly). STANDS.
- exists_point_of_table (L509): hypotheses of passes_of_table => ONE primitive (u,v) with
  ClassAt(p^k, p^3) for all k. p^2 not| d is derived from Row (not_sq_dvd_of_row, correct for both
  disjuncts). STANDS.
- classAt_seven(_all), classAt_seventeen(_all): concrete instances; non-vacuity witnesses. STANDS
  (numerically verified, see 4).
None is vacuous: hypotheses are satisfied by the controls (d = 21, 34; p | d).

What it gives toward Proposition C: everything except (a) "an integer x that is a square mod p^k for
every k (with bounded root valuation / x != 0) is a square in Z_p (resp. Z_p[i])" -- immediate
(compactness of Z_p, or the same Hensel argument), norms are nonzero since v(N_j) <= 4; and (b)
Lemma 1 over Z_p (section 2, by hand) to call d "the class". So Prop C's arithmetic content is in
Lean; the Z_p wording and Lemma 1 over Z_p remain by hand. The commit message title "Proposition C
... in Lean" should say "modulo every power of p" (the file docstring L26-27 does say this).

## 2. Kernel-check hygiene
grep over ClassHensel.lean and its whole import closure (CoverDescent, ClassLocal, ClassConic,
ClassOdd, ClassFive, ClassTable): no sorry/admit, native_decide, implemented_by, extern, unsafe,
opaque, axiom, set_option, ofReduceBool, macros/elabs. `decide +kernel` in ClassHensel only on
closed concrete congruences/non-divisibilities of the controls (L531-573) and `¬ p^2 ∣ 21/34`
(L579, L585): all true and non-trivial-in-form (verified numerically), not statements that are true
by definition. Classical reasoning = by_cases/by_contra only (standard). lakefile has no weakening
options. I could NOT verify the build claim: no .lake/olean on this machine; tree clean at ad737f3.

## 3. Mathematics (all checked by hand)
- sq_lift (L133): y = p^e c, e < j, 2e <= K-1; x = p^(2e)(c^2 + p^(K-2e)s); Hensel on the unit part;
  root p^e x' with p^j not| it. Correct.  j = 3 needs K >= 5 (2j <= K+1).
- norm_not_dvd (L164): p^3 not| w => p^5 not| w^2 => (K>=5) p^5 not| N, i.e. v(N) <= 4, and
  v(N) = 2v(w) even (inside sq_lift). Correct.
- nrm_Q_not_dvd (L193): 5+5<=9+1, 9+5<=13+1, 13+5<=17+1: v(N(Q_f)) <= 16. Correct.
- split (L247-309): case is "-1 a square mod p" (equivalent to p = 1 mod 4; fine). v(X+iota Y) <= 16
  since it divides N(Q) mod p^K; with p^2 not| d, v(d(X+iota Y)) <= 17 < 18 <= K, so a root mod p^K
  has v <= 8 (< 9, j = 9 needs K >= 17); lift each component; recombine with h = 1/2:
  A = h(x1+x2), B = -iota h(x1-x2) gives A +- iota B = x1, x2 (checked). iota is re-chosen per k
  (mod p^(K+k+1)) - legitimate. Correct.
- inert (L313-454): p not| (x+iy) => p not| x^2+y^2 (else -1 a square). d Q_f not == 0 mod p^K
  because p^19 not| d^2 N(Q_f) (3+17 <= 19+1) and 2K >= 19; root content r with 2r < K; Newton
  step with alpha = (2(x^2+y^2))^-1 mod p, correction p^(k+1) alpha (t1+it2)(x-iy). Correct;
  K >= 10 suffices there (Gaussian valuation of d Q_f <= 9).
- Why 18: v(d * component) <= 1 + 16 = 17 < K is what makes d Q_f nonzero mod p^K so the root's
  valuation is controlled. The file does not claim 18 is sharp, and I do not claim it is either.
- Builder's d = p^19 remark (L24-26): right. At any ClassAt point, d = p^19 gives d Q_f == 0 ==
  0^2 mod p^18, but d Q_f = p^19 * (unit or component of valuation 0) has odd valuation and is not a
  Gaussian square mod p^20 (if v(z) >= 10, z^2 == 0; if v(z) <= 9, v(z^2) <= 18 < 19).

## 4. Positive controls
Script: exit 0, "all checks passed"; its 40 integers equal the Lean file's. Independent check
(/tmp/s1rev/indep_check.py): all 4 norm congruences and 8x2 Gaussian congruences hold mod p^18;
p^3 not| w_j; each Gaussian witness fits ONLY its own flag, in the order (fff),(fft),(ftf),... which
is the rcases order of `rintro ⟨_|_, _|_, _|_⟩`. Non-triviality: p | d (21 = 7*3, 34 = 17*2), both
in Row (7 = 7 mod 8, 17 = 1 mod 8); norm valuations [2,0,0,0], root valuations [1,0,0,0]; Gaussian
content of d*Q_f = 2 for all 8 f, root content 1; at 17 both split components have valuation 2.
So genuinely non-trivial (not unit-only), but mild: no control reaches v(w_j) = 2, unbalanced split
components, or v(d*Q_f) near 17. Fine for non-vacuity; the general theorem is what matters.
Doc error: ClassHensel.lean:523-524 says both controls have "N_0 = p^2 * 50"; true at 7
(2450 = 49*50), false at 17 (N_0 = 83810 = 17^2 * 290).
Script negative tests: the mutated root and the p^3*w root are refused by the congruences; the
"d * p^2" test is refused trivially by check()'s own `d % p^2 != 0` line (script L103), so it tests
nothing about the congruences. Script is not run by CI.

## 5. Draft negative control (reports/s1_ClassHenselNeg_DRAFT_unchecked.lean)
Claim checked by my own computation (/tmp/s1rev/neg_check.py): for 3:1 at (1,1) the pairs are
(3,1),(1,3),(0,4),(0,8), norms 10,10,16,64; 40613041^2 == 10 mod 3^18, 27 not| 40613041, 4, 8;
d = 3^19 kills all eight d*Q_f mod 3^18, so ClassAt(3^18) holds with a = b = 0. Q_(fff) = (0,-320),
3 not| 320, so d*Q = (0, -320*3^19); not a Gaussian square mod 3^20 (exhaustive analogues
(0, 3^3 c) mod 3^4 and (0, 3^5 c) mod 3^6: no root). no_gsq_odd's argument (2ab = p^19 * unit
forces v(a)+v(b) = 19, so v(a) != v(b) and v(a^2-b^2) = 2 min <= 18 < 20) is correct for every odd
p. The claim is RIGHT. Lean text read line by line; plausible (rcases `-` clears are deferred, the
rw targets match the GSqM unfolding) but it must be compiled before use. It shows only that some
bound on d is needed; it does not test sharpness of p^2 not| d or of K >= 18.

## 6. Integration (would break CI or leave the file unchecked)
- MagicSquareOfSquares.lean:35 imports ClassTable but not ClassHensel. Default `lake build`
  (defaultTargets = MagicSquareOfSquares) and scripts/Audit.lean / AuditInternal.lean
  (`import MagicSquareOfSquares`) never compile or axiom-audit ClassHensel; only the explicit
  target the builder used does. Add the import; regenerate scripts/audit_output.txt,
  scripts/audit_line.txt, scripts/audit_internal.txt.
- scripts/declared_theorems.txt = 568, but `cat MagicSquareOfSquares/*.lean | grep -c '^theorem
  \|^private theorem '` is now 597 (+29 from ClassHensel): CI build job (ci.yml, the `test` line in
  "Generated Lean files match their generators") fails as committed.
- MANIFEST.sha256 does not list ClassHensel.lean or class_hensel_witness.py (manifest was already
  stale at base 8c0a0a3: 32 unlisted, ~10 hash diffs; regenerate at publish).
- Docs not updated in this slice: COMPLETE.md:129, :131, :133, :135 still say Proposition C is by
  hand / "Not in Lean"; README.md Lean table (~L165) and CHANGES.md have no ClassHensel entry.
  Correct wording: "modulo every power of p, at one integer point, in Lean
  (classAt_all_of_classAt, exists_point_of_table); by hand: square mod every p^k at a nonzero
  integer => square in Z_p / Z_p[i], and Lemma 1 over Z_p".
- Commit title "Proposition C ... in Lean" overclaims slightly (the X_p / class wording is not in
  Lean); the file docstring (L26-27) is accurate.

## Verdicts
classAt_norms_all STANDS; gsqM_lift (+split, +inert) STANDS; classAt_all_of_classAt STANDS;
passes_all_of_passes STANDS; exists_point_of_table STANDS; classAt_seven(_all),
classAt_seventeen(_all) STAND (conditional, like all of the above, on the unverified build claim).
Draft negative control: claim RIGHT, Lean unchecked.
