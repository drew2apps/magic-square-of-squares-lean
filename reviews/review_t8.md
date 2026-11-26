# Hostile review of sess_t8 (SQUARE_SECTION.md, scripts/g1_full/square_section.py, reports/t8.md)

Reviewer: slice v5, run 3. Written by an AI model (Claude); not read by a person. My code:
`v5/scripts/g1_full/review_t8_rank1.py` (sympy + numpy + the standard library; imports only `sigma_exact.model`;
written after reading square_section.py's header docstring and function names only, sharing no code with it).
Outputs: `review_t8_rank1.txt` (the six curves, with `--brute --crosscheck`), `review_t8_class2.txt` (class 2 on all
108 curves m ≤ 300), `review_t8_alphaK.txt` (nine more curves), `review_t8_negative.txt` (6 forged/control inputs).

## Verdict
**Confirmed.** Rank exactly 1 for E⁺(15,8), E⁻(15,8), E⁺(56,33), E⁻(60,11), E⁻(99,20), E⁺(180,19), and C₂(ℚ) = ∅ (the
class 2 of Theorem S1 in Ш(E′/ℚ)[φ̂]) on exactly the 37 curves t8 lists, are reproduced by an independent
implementation. No mathematical error found in the covering, in the restriction on λ, or in the rank formula. Tier: C
for the listed curves (a written argument whose every computation is an exact script), as t8 says. Minor wording
points below.

## 1. The covering, re-derived
E: y² = x(x² + ax + b), a = −4a₀, b = 4a₀² + κ² (checked from `sigma_exact.model` on every curve). For a positive
squarefree d | b with d = α² + β² (every odd prime of b is ≡ 1 mod 4), C_d: N² = dM⁴ + aM²e² + (b/d)e⁴, and
multiplying by d: dN² = X₁² + X₂², X₁ = dM² − 2a₀e², X₂ = κe² (since b − a²/4 = κ²). A point has e ≠ 0 and N ≠ 0;
take gcd(M, e) = 1. z = (X₁ + iX₂)/(α + iβ) has N(z) = N², so z/N has norm 1 and z = λ·w² with λ ∈ ℚ*, w = s + it
(Hilbert 90), s, t coprime integers. Writing q₁ + iq₂ = (α + iβ)(s + it)²:
  X₂ = λq₂ ⇒ F₁ := λκ·q₂(s,t) = (κe)²;   dM² = λq₁ + 2a₀e² = λ(q₁ + 2a₀q₂/κ) ⇒ F₂ := dλκ(κq₁ + 2a₀q₂) = (dκM)².
This is t8's D_λ exactly (for d = 2 t8 also uses the Pythagorean form λ(s² − t², 2st); equivalent).

**Restriction on λ.** Let p be odd, p ∤ dκ. A common prime of X₁, X₂ divides κ or e; if e, it divides dM², hence d
(gcd(M, e) = 1). So p ∤ gcd(X₁, X₂). Inert p: v_p(X₁ + iX₂) = 0 = v_p(λ) + 0 + 2·min(v_p s, v_p t) ⇒ v_p(λ) = 0.
Split p = ππ̄: X₁ + iX₂ is prime to one of π, π̄, s + it is prime to one of them, and the valuation identity at that
prime gives v_p(λ) = 0 (t8 states "≡ (mod 2)"; it is in fact exactly 0). So λ ∈ ±{squarefree divisors of 2dκ} modulo
squares, as t8 uses. Remark: D_λ ≅ D_{−λ} by (s, t) ↦ (t, −s) (as (t − is)² = −(s + it)²), which is why t8's
survivors always come in pairs ±λ; harmless (it doubles the work).

**Rank formula.** α(E(ℚ)) ⊇ K = ⟨α(T), α(R)⟩ = {1, HH⁻, 2H, 2H⁻}; a group meeting a coset meets all of it, so
excluding one element per non-trivial coset of K in Sel(E) gives α(E(ℚ)) = K; with Sel(E′) = {±1} (and −1 = α′(T′)),
2^r = 4·2/4 (Silverman–Tate III.6). Correct.

**Good primes.** For odd p dividing neither the contents nor disc(F₁)disc(F₂)Res(F₁,F₂), e² = F₁, M² = F₂ is a smooth
complete intersection of genus 1 over 𝔽_p (Jacobian criterion: a singular point needs a double root or a common root),
so it has (√p − 1)² > 0 points and they lift. My script computes that product for every cover and divides out the
tested primes; nothing was ever left over (all primes are among those of 2·m·n·(m²−n²)·k·H·H⁻).

## 2. The local test, re-implemented
- **Method S** (all odd p): forms divided by their content at each level; non-root residues have constant class
  (parity of the content valuation, Legendre symbol); only roots mod p are refined; for p > 10⁵ a structural
  certificate (content of odd valuation; a form ≡ c·h² with (c/p) = −1; a product of two forms ≡ c·h² with
  (c/p) = −1) is required before "insoluble" is accepted.
- **Method T** (p = 2, and odd p < 60 as a cross-check of S): Taylor-coefficient disc refinement on both charts.
  S and T agreed on every cover at every odd p < 60 (assertion in the run), and on 300 random pairs of forms.
- **A third, structurally different certificate** for every cover declared insoluble (one element per coset): pruned
  brute force mod p^K (the residues (1 : t) mod p^K and (p s : 1) at which every form is a square mod p^K, lifted level
  by level; empty ⇒ no ℚ_p-point), at the place found, else at the others. All certified: 1088, 1088, 15488, 1152, 2304, 1280
  covers on E⁺(15,8), E⁻(15,8), E⁺(56,33), E⁻(60,11), E⁻(99,20), E⁺(180,19), 0 uncertified (numbers per curve in `review_t8_rank1.txt`). (For E⁺(15,8) I also certified every element of
  every coset: 5696 covers, 0 uncertified.)
- Real place: one exact rational point per arc between the real roots (exact isolating intervals).

## 3. Results
- **Selmer groups**, recomputed with the same local code on the quartics C_d (every d | b, d | b′ squarefree, both
  signs; places ℝ and p | 2bb′): equal to `selmer_square_5000.json` on all six curves (|Sel(E)| = 16, 16, 64, 16, 16, 16;
  Sel(E′) = {±1} on all six).
- **Every element of every non-trivial coset** of K (not one per coset) has all its covers D_λ insoluble somewhere,
  on all six curves, so the cosets are internally consistent (a coset with one excluded element must have all
  excluded). Obstruction places for E⁺(15,8), class 2: {2, 3, 7, 13, 17, 37}, the same set t8 prints.
- **Rank exactly 1 on the six curves: confirmed.** PARI `ellrank` (guide only) gives [1, 1] and root number −1 for
  E⁺(180,19), E⁻(60,11), E⁻(99,20) as well (t8 had not run the sixth).
- **Class 2 on all 108 square-type curves m ≤ 300** (`--class2 300`): for every curve my number of covers and of
  surviving covers equals t8's (108/108); in particular the 37 curves with 0 survivors are exactly t8's "sha" list.
- **The three points of class 2** (E⁻(4,3) X = −1968, E⁻(12,5) X = 11640720, E⁻(24,7) X = 237904800): on the curve,
  x = X + 2S³ of class 2, checked exactly. (Not said by t8: on these three curves this gives α(E(ℚ)) ⊇ ⟨HH⁻, 2H, 2⟩,
  so rank ≥ 2 there.)
- **α(E(ℚ)) = K on 9 of t8's other 15 curves** (E⁺(4,3), E⁺(21,20), E^±(77,36), E⁺(91,60), E⁻(143,24), E⁺(165,52),
  E⁻(180,19), E⁻(195,28); `review_t8_alphaK.txt`, without the brute-force pass): every coset excluded, Sel groups equal
  to the JSON, |Sel(E′)| = 4, 4, 4, 4, 8, 4, 4, 4, 4 as t8 says; so 1 ≤ rank ≤ log₂|Sel(E′)| there. Not rerun: E⁻(56,33),
  E⁻(91,60), E⁻(140,51), E⁻(252,115), E⁻(253,204), E⁻(275,252).
- A cheap next step t8 did not take: on these 15 curves rank = log₂|α′(E′(ℚ))|, so the same second descent on the
  coverings of E′ (classes of Sel(E′) outside {±1}) could prove rank 1 on more of them, or find the points of rank 2.

## 4. Points on the text (none changes a result)
1. §2 "p inert gives v_p(gcd(X₁, X₂)) ≡ v_p(λ) (mod 2)": true but weaker than what holds (v_p(λ) = 0); fine.
2. "a non-trivial element of the Tate–Shafarevich group (of the 2-isogenous curve, in the kernel of the isogeny to E)":
   correct; clearer as "a non-zero element of Ш(E′/ℚ)[φ̂], φ̂: E′ → E the dual isogeny" (C_d are φ̂-coverings).
3. "No curve with s = 2 has C₂(ℚ) = ∅ (that would contradict the finiteness of Ш, by the parity remark)": precisely,
   C₂(ℚ) = ∅ with s = 2 forces α(E(ℚ)) = K, rank 1 and s even, which the 2-parity theorem allows only if
   Ш(E)[2^∞] is infinite. State it that way.
4. "rank ≥ 2 fails at these six curves ... the analogue of Theorem 2 of UNIFORM_RANK.md is false for the square type":
   correct. Note what is not claimed (t8 says it): nothing about the generic rank of the family; the six fibres could
   be exceptional fibres of a section.
5. "PARI's ellrank ... also gives [1, 1] on the first five; the sixth was not run": now run (above), [1, 1].
6. Status line "tier C ... for the listed curves only": appropriate. The rank-1 statement also rests on the
   Selmer groups of SELMER_SQUARE.md, which this review recomputed independently for these six curves.

## 5. Suggested LEDGER/doc line (for the conductor)
"Square type, m ≤ 300 (SQUARE_SECTION.md, tier C for the listed curves; experiment, not counted): the class 2 of
Theorem S1 is not the image of a rational point on 37 of the 108 curves (it is a non-zero element of Ш(E′)[φ̂]);
rank E(ℚ) = 1 exactly for E⁺(15,8), E⁻(15,8), E⁺(56,33), E⁻(60,11), E⁻(99,20), E⁺(180,19), by a second descent on the
coverings of every non-trivial coset of K. So no analogue of Theorem 2 of UNIFORM_RANK.md holds at the square type.
One independent AI review re-derived the covers and re-implemented the local tests (three methods); same results on
all 108 curves for the class 2 and on the six rank-1 curves for every coset."
