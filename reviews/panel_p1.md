# Panel P1 (maths referee) -- release 0.7.0

Written by an AI model (panel P1), not read by a person.

> Committed as written apart from this note. It reviewed the release candidate of 0.7.0 (commit c857655 of the integration branch); its fixes are in `CHANGES.md` ("After the release panel"). Paths under `/home/claude/` are the panel's scratch folders.

Scope: tier-C written proofs new since 0.6.0: genus3/NOTES.md (t12), conductor-t11/NOTES.md (t11),
reviews/review_z16.md law (b), spot-check of tier-C labels for pairing formula (note2) and Theorem B (explore5).
Repository read at /home/claude/conductor/int11 (branch run11_int11), not edited. Scratch: /home/claude/conductor/panel_p1/.

Findings are appended as the review proceeds.

## Item 1: experiments/genus3/NOTES.md (build agent t12)

Runs (python3 -I, this machine): genus3_split.py, genus3_points.py, genus3_points.py --negative, genus3_isog.py,
genus3_mult.py, kcurve_exact.py all print their committed outputs exactly (stdout differs only by PARI's
"stack size" line where PARI is loaded; run_genus3.sh filters it). genus3_rank_upper.py not rerun (5 min, tier D).

**Theorem A (splitting): holds.** Read line by line and re-derived with my own sympy script
(panel_p1/scripts/p1_genus3_check.py, written from the prose, not from genus3_split.py):
- F0F8 = x(Az+B), LU = 4x(Cz+D), AD − BC = −(m⁸+14m⁴n⁴+n⁸), Y² = f/(4t⁴) = (Az+B)(Cz+D); (t²−t⁻²)² = z²−4,
  (t±1/t)² = z±2: all true over Q(m,n).
- the sign table of h_i under σ, τ (t²−t⁻²: +,−; t+1/t: −,+; t−1/t: −,−) gives exactly K1 = <σ, ιτ>,
  K2 = <τ, ισ>, K3 = <στ, ισ>; the group-ring sum is 3 + σ + τ + στ + 2(ιτ + ισ + ιστ) and with ι* = −1 it is
  4 − Σ_{G'} g; the character trivial on K_i with χ(ι) = −1 is unique for each i (so each E_i is the χ_i-part, and
  the genus-2 quotients X/<ιτ>, X/<ισ>, X/<ιστ> are ~ E1×E3, E2×E3, E1×E2 as the note says).
- disc_t f = 2^48 (mn(m²−n²)S q1 q2)^8, the E_i discriminants as printed: genus 3 / genus 1 exactly when
  mn(m²−n²) ≠ 0 over Q (S, q1q2 > 0 for real (m,n) ≠ 0).
What the argument actually needs (all present, all standard): (a) π_H^*π_{H*} = Σ_{h∈H} h^* on Pic⁰ (true on
divisors, ramification included); (b) ι^* = −1 on J of a hyperelliptic curve; (c) J(genus 0) = 0, i.e. the
conic X/G' (index count 8/2 = 4 = |G'| is correct); (d) a homomorphism with finite kernel between abelian
varieties of equal dimension is an isogeny. The full Kani–Rosen theorem is not needed; the note's remark that this
is Kani–Rosen for the Klein group G' (J(X) ~ J(X/σ) × J(X/τ) × J(X/στ), J(X/G') = 0) is also correct
(ε_1 + 2ε_{G'} = ε_σ + ε_τ + ε_στ). Defined over Q: yes, φ_i are over Q and each E_i has a rational point
(z = ±2 on E1, ∞ on E2, E3). The condition mn(m²−n²) ≠ 0 is exactly right.
- [note] Notation: in NOTES.md §A the same symbol `φ_i*` is used for push-forward (in Φ) and pull-back (in
  Σ φ_i*), and "π_H*π_H*"; write φ_{i*} and φ_i^*. Cosmetic, the proof is unaffected.

**(B) points of infinite order on the 27 curves: holds** (method read: point on curve exactly, nP ≠ O for
n ≤ 12, Mazur; model change to y² = x(x−e2)(x−e3) checked symbolically — I re-derived the E1 constant
c = 4(2A+B)(2C+D) = 16S²(m²−n²)² by hand; output reproduced). The rank table is tier D (rank_second).

**(B′) "each of E1, E2, E3, F0·F8, L·U occurs exactly once in J(X5)": holds as a tier-D statement given the
quoted decomposition.** I re-counted the decomposition: Riemann–Hurwitz gives g(X5) = 49; C_I have genus
1 (10), 2 (10), 3 (5), 4 (1); the even/palindromic C_I are 6 of the genus-2 ones, all 5 of genus 3 and the genus-4
one, giving 29 genus-1 and 10 genus-2 pieces (29 + 20 = 49), as the note says. The upper bound uses only "isogenous
⇒ equal a_p / L_p factor", correct. For Demjanenko–Manin the right count is rank Hom_Q(J(X5), E) =
multiplicity (End_Q E = Z), so multiplicity 1 with rank ≥ 1 is indeed no route.
- [note] "or any curve between it and the solutions": correct if it means quotients of X5 (curves Z with
  X5 → Z); a curve that covers X5 can carry more copies of E. Say which.

**(D) the corner curve over Q(i) is not a Q-curve: holds; the criterion is a complete argument.**
- Recomputed independently (panel_p1/scripts/p1_kcurve.py: sympy over Z[i], Cremona's I, J, own valuation
  code): v(j) = −8,−8 at 5; −4,−4 at 29; {−4, 3} at the two primes above 37 and {−4, 3} above 181 (my π for 37 is
  1+6i = i·(6−i), so my "3 at π, −4 at π̄" is the note's "−4 at 6+i"): agrees.
- The criterion is complete for Q̄-isogeny, so all twists (quadratic or not; ramified or not) are covered: a twist
  of E^c is Q̄-isomorphic to E^c and has the same j; potentially good reduction is read off j alone; an isogeny
  over a finite L extends to L', good reduction is an isogeny invariant over a fixed field (Serre–Tate), and
  v_𝔓'(j) = e·v_𝔭(j) < 0. CM is excluded by non-integrality of j (and the criterion does not need non-CM).
  The note's criticism of the session's "unramified twist" argument is right.
- [note] "L(E/K,s) is not a product of two L-functions over Q": true in the sense that Ind_K^Q ρ_{E,ℓ} is
  irreducible (Mackey: ρ_E ≇ ρ_E^c, since E is not K̄-isogenous to any twist of E^c) — worth one clause. Also
  every elliptic curve over Q(i) is modular by Caraiani–Newton (arXiv:2301.10509: all elliptic curves over
  Q(√−d), d = 1, 2, 3, 5; abstract checked by web search, the paper not read), so the note's sentence about Ribet's
  route is correct but should not be read as "E is not known to be modular"; what is missing is a rank-0-from-L
  theorem over imaginary quadratic fields, not modularity.
- [fix] LEDGER: NOTES.md calls §D a tier-C written proof, but LEDGER has §D only inside the tier-D row (line 126)
  and the tier-C genus-3 row (line 87) has only the splitting. Either add the criterion to the tier-C row or
  call §D tier D in NOTES.md.

**LEDGER row (line 87) strength: right.** "isogenous over Q ... at every pair with mn(m²−n²) ≠ 0 ... Kani–Rosen" is
what is proved. [note] "Kani–Rosen" in the row is more than is used; "(π^*π_* = Σ h^*, ι^* = −1)" would be the
honest input list.

## Item 2: experiments/conductor-t11/NOTES.md (build agent t11)

Runs: `python3 -I t11_theta.py --negative` (9 of 9 as expected, `# matches`), `t11_theta.py "Sig+(22057,2351)"`
(`# matches`, 30 s) and `t11_theta.py "E1(8249,3720)"` (`# matches`, 112 s). Fisher's paper: the PDF could not be
downloaded (proxy 403); Sections 2 and 3 read through WebFetch summaries (as the repository did before), so this
check has the same gap. What the summaries give: S_1 = S^(φ)(E/K), S'_1 = S^(φ̂)(E'/K); S_{2n} = image of
S^(2^n)(E'/K), S_{2n+1} = image of S^(2^nφ)(E/K); E'(K)/φE(K) ⊆ every S_m; Thm 2.1: m odd, Θ_m alternating on
S_m × S_m (and Θ'_m on S'_m), kernel S_{m+1}; m even, Θ_m: S_m × S'_m, kernels S_{m+1}, S'_{m+1}; Lemma 2.2:
rank ≤ dim S_{2n−1} + dim S'_{2n} − dim E(K)[φ] − dim E'(K)[φ̂]; Section 3: for x in S^(ψ)(E'/K) and y in
S^(φ̂)(E'/K) (both torsors under E'), ⟨x, y⟩ = Σ_v (f(P_v), y)_v with f a pushout, div f = 2b, b ↔ the generator of
E'(K)[φ̂]; Θ_m(ξ, η) = ⟨ξ1, η⟩ for any lift.

**Bookkeeping: holds.**
- "Θ₂ = 0 automatically for type A": holds. The images of rational points lie in S_3 and S'_3 (the kernels), and for
  type A one side's S_2 consists of images of points only, so Θ_2 vanishes on it.
- "Θ₃ with one Lemma 7.1 pushout, z on the same side as y": holds, and I re-derived why. For m = 3, Θ_3(ξ, η) =
  ⟨ξ1, η⟩ with ξ1 in S^(2φ)(E), a torsor under E; Section 3 (with E, E' swapped) pairs it with η in S^(φ)(E), the
  first Selmer group of the 2-isogeny out of the same curve, i.e. the same side; the pushout is attached to the
  only rational 2-torsion point of Jac(C3) (Q(E[2]) = Q(i)), and both Galois-stable cone pairs give the same class
  b (in an E[2]-torsor of four theta classes, D1 − D2 = D3 − D4), which is why f/f† = c·h². The local term is
  constant in P_v because the values f(P_v) form a coset of the image of the other side's Kummer map, the exact
  annihilator of this side's under the Hilbert symbol; so the point-independence check is a real test (the
  negative test shows it).
- [note] The older `scripts/g1_full/rank_fisher_note.md`, §"Formulas used" 3, states the Section 3 formula with "y
  in the first Selmer group of the other side" for x in Sel(m) of any m. That is right for even m only (what that
  note used); for odd m it is the same side, as t11 uses and as rank_fisher.py's docstring (6) implicitly does.
  Worth one clause there, since t11 now relies on the odd case.
- The two bounds: **E1(8249,3720)**: S'_2 = <known, 13278577, 755496193> (dim 3), S_2 = <known (2), 113,
  2991753233> (dim 4); the only non-zero entry is Θ_2(13278577, 2991753233), found from both sides
  (λ = −7, 7 and λ = −704698813, −55369215681829); rank 1, so dim S_3 = 3, dim S'_3 = 2; Θ_3, Θ'_3 are alternating
  with one open class, hence 0; rank ≤ dim S_3 + dim S'_4 − 2 = 3 + 2 − 2 = 3. **Sig+(22057,2351)**: S'_2 = S'_1
  is all known, so Θ_2 = 0 and S_3 = S_1 (dim 5, basis −1, 14, 274, 339, 108146; −1 = image of T'); computed rows
  Θ_3(14, ·) = (0,0,0,1,0) and Θ_3(−274, ·) = 0 (the row of −274, equal to that of 274 since −1 is in the kernel);
  kernel ⊇ <−1, 274> gives rank ≤ 3, alternating gives rank even, Θ_3(14, 339) ≠ 0 gives rank ≥ 2: rank 2,
  dim S_4 = 3, and Lemma 2.2 with E, E' swapped, n = 2: dim S'_3 + dim S_4 − 2 = 2 + 3 − 2 = 3. **Both 5 → 3 are
  correct consequences of the stated pairing values.**
- Local terms: places S = {R} ∪ {p < 47} ∪ bad primes ∪ the pencil-discriminant primes (factored completely by
  PARI, every factor proved prime) ∪ primes where f0 ∈ span(QA, QB) mod p (gcd of 3×3 minors, exact) ∪ primes of
  the z. The written argument outside S is complete: good reduction of C3', f0 not identically 0 on the reduction,
  ≥ p + 1 − 2√p > 8 points for p ≥ 47, Hensel, units. I tested the repository's Hilbert symbol (used for every term)
  against PARI's `hilbert` on 2100 random pairs at p = 2, 3, 5, 7, 23, 113, 2351: no mismatch
  (panel_p1/scripts/p1_hilbert.py).
- The replaced `choose_Xi_lin`: harmless. The set {f(P_v)} at a place is a full coset of the local image, so
  closing sampled classes with the certified local image is the correct set, and in any case C4 is re-checked
  locally soluble at every place. It is used only for Θ_4 on E1(8249,3720), whose value 0 lowers nothing, so no
  bound rests on it.

**Gaps (holds with gaps):**
- [note] Fisher read through tool summaries (as before), including by me.
- [note] (closed by my control below) Sig+(22057,2351): the single non-zero entry Θ_3(14, 339) = −1 that carries the bound comes from one
  level-3 lift and one pushout (terms −1 at 2, 23, 113); its controls on the same C3 (Θ_3(14,14) = 0,
  Θ_3(14,−1) = 0, point-independence) are good but not independent of that lift. "Symmetry (Θ₃(14, 274) =
  Θ₃(274, 14))" compares two zeros. The natural independent check, Θ_3(339, 14) from a lift of 339 (or a second
  lift/the other cone pair of 14), was not done. Same for E1(8249,3720) is not needed: there the non-zero Θ_2 entry
  is found from both sides.
- [fix] NOTES.md bottom line: "E1(4440,2231): Θ'₃(2, ·) = 0, so ... no bound below 2 is possible" does not match
  its own table ("5 (3 at best)"), and both are statements about levels ≤ 4 only: with 2 in the kernel and the
  E-side S'_3 of dim 5, rank Θ'_3 ≤ 2, so the level-4 bound is ≥ 3; Θ'_5 (not in Section 7) could still lower it.
  Say "no bound below 3 at level 4".
- [fix] LEDGER row (line 88): "on E1(10512,145)− and E1(4440,2231) the deciding pairing is 0" is right for
  (10512,145)− (Θ'_3 has one possible entry and it is 0) but overstated for (4440,2231): only the row of the class 2
  of Θ'_3 was computed; Θ'_3 itself may have rank 2 (which would lower 5 to 3); the other rows need fields with
  regulator ≈ 895 or more. Suggested: "on E1(4440,2231) the row of the class 2 is 0, so the bound cannot go below 3
  at this level".
- LEDGER row strength otherwise right ("not reviewed by a separate agent", Fisher sections, unchanged module with a
  replaced sampler — "corrected choice of Ξ" reads as if the original were wrong; it was a sampler that missed a
  class; [note]).

## Item 3: law (b) (reviews/review_z16.md, Task 4; note2 Lemmas 1–2)

**Holds, for 4 | mn and for m, n both odd, exactly as LEDGER (line 125) states it.** Re-derived line by line.
- Setting: E_08 = y² = x(x − A²)(x − C²), (A, B, C) = (m²−n², 2mn, m²+n²) (halved when m, n both odd; x ↦ 4x is
  an isomorphism and 4 is a square, so the classes do not change). κ = (1, β, β), β = sqf(2mn).
- Torsion images, computed: T_r = (0,0) ↦ (A²C², −A², −C²) = (1, −1, −1); T_m = (A²,0) ↦ (A², −A²B², −B²) =
  (1, −1, −1); T_* = (C²,0) ↦ (1, 1, 1) (indeed T_* = 2P, since C² and C² − A² = B² are squares). So the torsion
  classes are {(1,1,1), (1,−1,−1)} at every place, as z16 says.
- Odd p | A: I re-derived Lemma 2 (x ≡ 0: b_* = −1; x ≡ C²: x, x − A² unit squares, b_* = 1; x a unit off C², or
  x non-integral: b_* = 1), so the local image has b_* ∈ {±1}, and (1, ±1, ±1) are torsion images: κ is in the
  image iff ±β is a square mod p. Since m ≡ ±n (mod p), β ≡ ±2n²: the condition is "2 or −2 is a square mod p",
  i.e. p ≢ 5 (mod 8). gcd(m − n, m + n) = 1 for opposite parity, so this is exactly "no prime ≡ 5 (mod 8) divides
  m² − n²".
- Odd p | B: Lemma 1 gives image ⊆ {(1, γ, γ)}, a set of 4 classes, and |E(Q_p)/2E(Q_p)| = 4, so equality; κ is in.
  Odd p | C: 2mn ≡ 2ιn² = ((1+ι)n)², a square. Good primes: β is a unit there. ∞: β > 0.
- p = 2, opposite parity, 4 | mn: C − A = 2n², C + A = 2m², the even one of m, n is ≡ 0 mod 4, so 32 divides one
  of them: H2. I re-derived Lemma 1 at 2 (H2 gives v₂(B²) ≥ 6; if v₂(x − A²) ≤ v₂(B²) − 3 the ratio
  (x − C²)/(x − A²) ≡ 1 mod 8; otherwise x ≡ A² ≡ 1 mod 8; non-integral x falls in the first case), so the image is
  all 8 classes (1, γ, γ) and κ is in.
- m, n both odd: A' = (m² − n²)/2 ≡ 0 (mod 4), C' odd, v₂(β) = 1. I re-derived Lemma 2 at 2 case by case
  (v(u) ≤ −2; v(u) = −1: classes 1 and 5; u odd; u ≡ 2 mod 4: 3 against 4·{1,5}; u ≡ 0 mod 4: −δ square), which
  needs A ≡ 0 mod 4 (true here; note2 states "A even" but the proof uses 4 | A). ±β is not a square in Q₂: κ ∉ Sel².
- [note] z16 writes "given note2's Lemmas 1, 2 and Cassels". Law (b) is a statement about local images only;
  Cassels' theorem is not used. z16 says it re-derived the Lemmas "at odd primes"; the 2-adic halves (Lemma 1 at 2
  for 4 | mn, Lemma 2 at 2 for both odd) are the ones that carry the two parity cases; I have now re-derived both,
  so the proof is checked at every place.
- [note] note2's Lemma 2 at p = 2 is stated for "A even" but proved (and true as used) for A ≡ 0 (mod 4).
- LEDGER strength: right ("proved ... when 4 divides mn ... both odd proved too; tested only when mn ≡ 2 (mod 4)").

## Item 4: tier-C labels (pairing formula; Theorem B)

**Pairing formula (note2, z13): tier C is justified; the input list in LEDGER (line 85) is nearly complete.** The
proof of (2) uses: Cassels' formula (1) with the component assignment (z12's derivation; z13 matched it to FSS 2010
§3/§5/Thm 5.2 in revision 3, with the "common factor l_u l_v l_w" correction); Cassels' local constancy (cases (i),
(ii) use one point each); FSS Thm 5.2 (Cassels = Cassels–Tate); the sizes of local images (4, 8, 2) and the
unramified image at good odd primes; Lemmas 1, 2 (proved); H2 (a hypothesis). All named in the row except:
- [fix] note2/NOTES.md "Tier" paragraph still says the assignment is "not yet matched to a numbered statement in
  Cassels 1998 or Fisher–Schaefer–Stoll 2010", while its own "Revision 3" section and LEDGER say z13 matched it to
  FSS §3/§5/Thm 5.2. Make them agree (and say Cassels 1998 itself was opened by nobody; FSS was fetched by z13).
- [note] Tier E (line 144) says Cassels 1998 was "Read by AI models"; note2 says the session could not open it. If
  no agent read Cassels 1998, say so there.
- [note] "standard sizes of local images" should also name "the local image at a good odd prime is the unramified
  subgroup" (z13's T4), which the good-prime step uses.

**Theorem B (explore5 §1, z14): the proof holds (I re-checked Steps 1–3: L = −N(z − q)/2 on u² − dv² = c;
N(ε − 1) = −4dy²/N(η); norm compatibility and Galois equivariance of Hilbert symbols; Hilbert reciprocity over Q
and over K), so tier C "given Cassels' theorem" is justified; the LEDGER row (line 86) lists the inputs
incompletely.**
- [fix] The row says only "given Cassels' formula". It also uses: Cassels' local constancy ("γ_w ... at any point
  of the covering"); FSS Thm 5.2 if "pairing" means the Cassels–Tate pairing (as in row 85); the component
  assignment in the weak form "the conic omitting u1 goes with b_1" (the assignment between b_2, b_3 is irrelevant
  since both are β, but if b_1 = 1 were attached to another conic the formula would change); rational points on the
  two conics (Hasse–Minkowski, since a ∈ Sel²); and standard local class field theory (Hilbert 90 in K ⊗ Q_v, norm
  compatibility, Galois equivariance, Hilbert reciprocity over K). Copy the list of row 85.
- LEDGER wording "the Rédei form ... follows ... up to local factors at 2, ∞ and the primes of d, β, AB": z14 adds
  the primes of z0 (the denominator of π), where π is not a unit; row 86 omits z0. [fix] add "and of z0".

## Added control for item 2 (closes the Sig+(22057,2351) gap above, as far as the same code can)

In a scratch copy of the needed folders (panel_p1/repo/, the repository untouched), `p1_row339.py` builds a new
level-3 lift of the class −339 (λ = 13, ξ = 3173523211, C3 locally soluble at R and 19 primes), finds a conic point
for its U pair over Q(√339) with the t11 finder (5 s), and computes the row with `t11_lib.theta3`:
Θ_3(−339, ·) on (−1, 14, 274, 339, 108146) = (+1, −1, +1, +1, −1). So Θ_3(−339, 14) = −1 = Θ_3(14, −339) (the
bound-carrying entry, now from a different lift, a different cone pair and a different quadratic field);
Θ_3(−339, 274) = +1 (agrees with 274 in the kernel); Θ_3(−339, ±339) = +1 (alternating). New entry:
Θ_3(339, 108146) ≠ 0. The 4×4 matrix on (14, 274, 339, 108146) is then [[0,0,1,0],[0,0,0,0],[1,0,0,1],[0,0,1,0]],
rank 2, kernel S_4 = <−1, 274, 14·108146> (14·108146 = 4·378511; −378511 is among the listed open cosets). Fully
consistent with "Θ₃ has rank exactly 2, dim S_4 = 3". It is the same implementation, so this checks the value, not
the code.

## Bottom line (panel P1)

No blocker. All four items hold mathematically; the gaps are of labelling and of recorded inputs.
1. genus3 (t12): Theorem A (J(X) ~ E1 × E2 × E3 over Q iff-condition mn(m²−n²) ≠ 0) holds, re-derived
   independently; §D's j-valuation criterion is a complete argument (Q̄-isogeny, so all twists; CM excluded), the
   valuations re-computed independently; (B), (B′) hold as stated (tier D). Every script reproduces.
2. conductor-t11 (t11): the reading of Fisher §§2–3 (Θ_odd same side, Θ_even opposite sides, Lemma 2.2) is right;
   Θ₂ = 0 for type A is right; both 5 → 3 follow from the stated values; the values reproduce, and an independent
   lift of −339 (my control) confirms the one entry that carries Sig+(22057,2351)'s bound and gives the full rank-2
   matrix. One LEDGER sentence overstates E1(4440,2231).
3. law (b) (z16): holds for 4 | mn and for m, n both odd; checked at every place including both 2-adic lemmas.
4. Tier C for the pairing formula and Theorem B is justified; Theorem B's LEDGER row lists its inputs too briefly.

Time used: about 1 h 40 min. Not done: genus3_rank_upper.py (5 min, tier D) not rerun; Fisher's PDF not read
(download refused by the proxy; summaries only); no re-implementation of the t11 pushout code.
