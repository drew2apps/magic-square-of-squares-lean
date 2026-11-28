# The genus-3 quotient `y² = F0·F8·L·U` and the corner K-curve of (627,364): a separate session's note, checked

Written by an AI model (build agent t12), not read by a person. Not reviewed by a separate agent.
**Experiment, not counted; tier D of `LEDGER.md`, with the written proofs of sections A and D in tier C.**
Nothing here changes a count or closes a ratio.

The input is `session/NOTE_2026-11-27.md` (a separate AI session, untrusted; sections 2 and 3 are checked here,
sections 1, 4, 5 elsewhere) with its scripts, committed unchanged apart from a banner line in `session/`. Their four
outputs reproduce byte for byte under PARI 2.17.4 (CI's version). These notes take precedence over the session's note.

## Bottom line

| Claim of the note | Verdict | Where |
|---|---|---|
| (A) `J(y² = F0F8LU) ~ E1 × E2 × E3` over ℚ (checked there by traces at 12 primes for (627,364)) | **Holds exactly**, at every rational pair with `mn(m² − n²) ≠ 0` (written proof below, every identity by sympy) | `genus3_split.py` |
| (B) no `E_i` has rank 0 at the nine open ratios, "so no new closure" | **Holds exactly** (a point of infinite order on each of the 27 curves). The ranks in the note's table are PARI's: exact here are rank = 1 for 11 curves, rank ≤ 2 for 5 more; rank 2 is proved for none | `genus3_points.py`, `genus3_rank_upper.py` |
| (B′) "no rank-1 class occurs twice, so Demjanenko–Manin gets nothing" | **Holds exactly for `E1, E2, E3` at all nine ratios and for `F0·F8`, `L·U`**: each occurs exactly once in the Jacobian of the curve of the five forms (on the quoted Kani–Rosen splitting). The note's list was not the whole Jacobian | `genus3_mult.py` |
| (C) the quotients of (627,364) form 15 isogeny classes (equal traces at 40 primes) | **Holds exactly**: explicit isomorphisms and one explicit 2-isogeny inside the groups, a differing `a_p` across them. The "twist pairs" are twists only up to isogeny for `F4·L`, `F4·U` | `genus3_isog.py` |
| (D) the corner curve `L1·L3·M1·M2` of (627,364) over ℚ(i) is not a ℚ-curve | **Holds exactly**, but the note's argument has a gap (it excludes only twists of the conjugate "by a character unramified there"). Closed here by the `j`-invariant alone | `kcurve_exact.py` |

## A. The splitting (tier C: written proof; the identities computed)

Notation of `FOUR_POINTS.md` §1: `F0 = m² + n²t²`, `F8 = n² + m²t²`, `L = GL`, `U = GU`, `S = m² + n²`; `x = t²`,
`z = x + 1/x`; `A = m²n²`, `B = m⁴ + n⁴`, `C = S²`, `D = 2S² − 16m²n²`;
`E1: w² = (z² − 4)(Az + B)(Cz + D)`, `E2: w² = (z + 2)(Az + B)(Cz + D)`, `E3: w² = (z − 2)(Az + B)(Cz + D)`.

**Theorem A.** Let `m, n ∈ ℚ` with `mn(m² − n²) ≠ 0`, and `X` the smooth projective curve `y² = f(t)`,
`f = F0·F8·L·U`. Then `X` has genus 3, each `E_i` has genus 1 and a rational point, and
`Φ = (φ1*, φ2*, φ3*): J(X) → E1 × E2 × E3` is an isogeny over ℚ, where, with `Y = y/(2t²)`,
`φ1 = (z, Y(t² − t⁻²))`, `φ2 = (z, Y(t + t⁻¹))`, `φ3 = (z, Y(t − t⁻¹))`. Moreover `(Σ φ_i*)∘Φ = [4]`.

*Proof.* (Every identity is checked by `genus3_split.py` over `ℚ(m, n)`, items (1)–(7) of its output.)
`f` has degree 8 and is squarefree: its four quadratic factors have discriminants `−4m²n²`, `−4m²n²`,
`−16(m² − n²)²` (twice) and pairwise resultants that are products of `m, n, m ± n, S` and the two quartic factors of
`m⁸ + 14m⁴n⁴ + n⁸`, which is positive for real `(m, n) ≠ 0` (7). So `X` is hyperelliptic of genus 3; let `ι` be the
hyperelliptic involution. `f(−t) = f(t)` and `t⁸f(1/t) = f(t)` (1), so `σ(t, y) = (−t, y)` and
`τ(t, y) = (1/t, y/t⁴)` are automorphisms over ℚ; with `ι` they generate a group `(ℤ/2)³` of eight distinct
automorphisms (2). For a finite group `H` of automorphisms and `π_H: X → X/H`, `π_H*π_H* = Σ_{h∈H} h*` on
`J(X) = Pic⁰(X)`, and `ι* = −1` on `J(X)` (`P + ιP` is the hyperelliptic class).

(i) `F0F8 = x(Az + B)`, `LU = 4x(Cz + D)` (1), so `Y² = (Az + B)(Cz + D)` on `X`, and `z, Y` are fixed by
`G′ = ⟨σ, τ⟩` (3). `z` has degree 8 on `X` and 2 on the conic `Y² = (Az + B)(Cz + D)` (smooth: `AD − BC =
−(m⁸ + 14m⁴n⁴ + n⁸) ≠ 0`), so `ℚ(z, Y)` has index 4 = `|G′|` in `ℚ(X)` and equals `ℚ(X)^{G′}`: `X/G′` has genus 0
and `Σ_{g∈G′} g* = 0` on `J(X)`.

(ii) `φ_i` lands on `E_i` and is invariant exactly under `K1 = ⟨σ, ιτ⟩`, `K2 = ⟨τ, ισ⟩`, `K3 = ⟨στ, ισ⟩` (4); it has
degree `8/2 = 4 = |K_i|` (5), so `ℚ(E_i) = ℚ(X)^{K_i}` and `φ_i*φ_i* = Σ_{g∈K_i} g*`. The polynomials of `E_i` in
`z` are squarefree of degrees 4, 3, 3 under the same condition (7), so `g(E_i) = 1`; `E1` has the rational points
`z = ±2`, `E2`, `E3` the point at infinity.

(iii) In the group ring, `Σ_i Σ_{g∈K_i} g = 3 + σ + τ + στ + 2ιτ + 2ισ + 2ιστ` (6); with `ι* = −1` this acts as
`4 − Σ_{g∈G′} g* = [4]` on `J(X)`. So `(Σ φ_i*)∘Φ = [4]`, `ker Φ ⊆ J(X)[4]` is finite, and since both sides have
dimension 3, `Φ` is an isogeny; it is defined over ℚ because the `φ_i` are. ∎

This is the Kani–Rosen relation for the Klein group `G′` (`J(X) × J(X/G′)² ~ J(X/σ) × J(X/τ) × J(X/στ)`), in the
form needed here. The intermediate curves are 2-isogenous to the `E_i`: `X/σ` is `[F0F8LU|x=t²]` of the note, mapped
onto `E1` by `(x, y) ↦ (x + 1/x, y(x − 1/x)/(2x))` (`genus3_isog.py`); `X/τ → E2` and `X/στ → E3` by
`(w, Y) ↦ (w² ∓ 2, wY)` with `w = t ± 1/t`.

*Control (not a proof):* `a_p(X) = Σ a_p(E_i)` at all 892 (ratio, prime) pairs with `5 ≤ p < 600` good for the nine
ratios, and `#X(F_{p²})` agrees at 87 pairs with `p < 60`, all counts in pure Python.

*Consequences, exact once (B) holds:* `rank J(X)(ℚ) = Σ rank E_i(ℚ) ≥ 3 = g(X)`, so Chabauty–Coleman does not apply
to `X`; the genus-2 quotients `X/⟨ιτ⟩`, `X/⟨ισ⟩`, `X/⟨ιστ⟩` have Jacobians `~ E_i × E_j` (by the same character
count) of rank `≥ 2 =` genus, so not to them either.

## B. Ranks of the 27 curves (tier D, computed)

`genus3_points.py` uses its own integral models `y² = x(x − e2)(x − e3)` (full rational 2-torsion; the change of
variables from the note's `E_i` is checked with sympy for general `m, n`). For each curve a point found by PARI's
`ellrank` (a finder only; the points are stored in `genus3_points_hints.json`, so the output is deterministic) is
checked exactly: on the curve; `nP ≠ O` for `n ≤ 12`, so of infinite order by Mazur (quoted); a multiple `2P` or `3P`
with non-integral `x` (Nagell–Lutz, quoted, a second check); and its image under the 2-descent map together with the
torsion spans dimension 3 of `E(ℚ)/2E(ℚ)`. The torsion is proved complete (a closed list whose order is
`gcd #E(F_q)`, `17 ≤ q < 400`).

| ratio | E1 | E2 | E3 |
|---|---|---|---|
| (195,28) | **1** | 1–3 (PARI 1–3) | **1** |
| (223,167) | 1–3 (PARI 1–3) | **1** | **1** |
| (440,279) | **1** | **1** | 1–2 (PARI 2) |
| (609,200) | 1–4 (PARI 0–4 alone; the note: 2–4 over the class) | 1–4 (note: 2) | 1–2 (PARI 2) |
| (627,364) | 1–3 (PARI 1–3) | 1–3 (PARI 1–3) | **1** |
| (719,161) | **1** | **1** | 1–2 (PARI 2) |
| (809,409) | 1–4 (note: 2) | 1–4 (note: 2–4) | 1–2 (PARI 2) |
| (861,620) | 1–2 (PARI 2) | **1** | 1–3 (PARI 1–3) |
| (991,263) | 1–3 (PARI 1–3) | 1–3 (PARI 1–3) | **1** |

Bold: rank exactly 1. The lower bound 1 is exact everywhere. The upper bounds are exact from
`genus3_rank_upper.py`: the complete 2-descent with the Cassels–Tate pairing of `scripts/g1_full/rank_second.py`
(imported unchanged; its theory quoted as there), on the curve and on the curves with full 2-torsion 2-isogenous to it.
PARI's own bounds (in parentheses) use the whole isogeny class and are sometimes better ((609,200) E2,
(809,409) E1); they are PARI's.

**Rank 2 is not proved for any curve.** On the nine curves of root number +1 (PARI's `ellrootno`), every point found,
also on the other curves of the isogeny class and mapped back by PARI's dual isogenies (images checked exactly), is
dependent on the first one modulo torsion (exact `F_3` test through `E(F_q)[3]` at 60 primes). PARI's lower bound 2
there comes from parity, not from points. For the five curves with exact upper bound 2 the rank is 1 or 2; it is 2
if `Ш[2^∞]` is finite (root number +1, the 2-parity theorem of Dokchitser–Dokchitser, quoted).

## B′. Demjanenko–Manin: multiplicities in the Jacobian of the curve of the five forms (tier D, computed)

`X5`, the fibre product of the five double covers `y_i² = F_i(t)`, has genus 49. Quoted: `J(X5) ~ ∏_I J(C_I)` over the
31 subsets `I` with `|I| ≥ 2`, `C_I: y² = ∏_{i∈I} F_i` (Kani–Rosen for `(ℤ/2)⁵`), and the further splitting of each
`C_I` that is even (`y² = h(x)`, `y² = xh(x)`) or palindromic (by `w = t + 1/t`). `genus3_mult.py` lists the resulting
29 pieces of genus 1 and 10 of genus 2 (the docstring has the list) and checks the list against direct point counts:
`L_p(C_I)` equals the product of the `L_p` of its pieces for the eight `C_I` of genus 2 to 4 at `p = 17, 23` (counts
over `F_{p^k}`, `k ≤ 4`, Zech logarithms built in the script). Then, exactly: a genus-1 piece with one differing `a_p`
is not isogenous to `E`; a genus-2 piece contains no factor `E` if at one good prime `1 − a_pT + pT²` does not divide
its `L_p` (true for all ten pieces and all classes, at some `p ≤ 41`).

Result at each of the nine ratios: `E1`, `E2`, `E3`, `F0·F8` and `L·U` each occur **exactly once** in `J(X5)` (at least
once by the explicit maps). So the number of independent maps from `X5` (or any curve between it and the solutions)
to any `E_i` is 1, and Demjanenko–Manin needs rank 0: no route through the `E_i`. The other classes occur at most twice
(the class of `F0·L`, the four-point curve, four times). Computed here: the multiplicity bounds. Reasoned: the
decomposition theorem, and "D–M needs more independent maps than the rank". Chabauty on `X` and its genus-2
quotients: section A. Not looked at: elliptic Chabauty over number fields, and the classes that occur twice.

The note's list of quotients of (627,364) (`isogeny_classes.py`) is not the full decomposition: it misses the four
copies of the genus-2 curve `C_{F0,F4,L}`, the genus-2 parts of `C_{F0,F4,L,U}`, `C_{F0,F8,F4,L}` and the two genus-2
parts of `C_I` with all five forms, and the genus-1 pieces `y² = F0F8F4(x)` and `y² = g4(w)` (each has the `a_p` of a listed curve or of
its partner at every good `p < 600`; whether they are isogenous is not needed for the upper bound and not proved). Its sentence "no rank-1 class occurs twice" was therefore a
statement about the listed quotients; with `genus3_mult.py` it holds for the whole Jacobian, for the classes above.

## C. The isogeny classes of (627,364) (tier D, computed)

`genus3_isog.py`, for the 25 curves of the note's `isogeny_classes.py`:
- inside the 15 groups, explicit isomorphisms `x ↦ −x, 1/x, −1/x` (identities of binary quartics, sympy, general
  `m, n`) and the degree-2 map `[F0F8LU|x=t²] → E1` above (a nonconstant map of genus-1 curves induces an isogeny of
  the same degree, quoted);
- across the groups, all 105 pairs differ in `a_p` at some prime `≤ 41` (isogenous curves have equal `a_p`, quoted);
- "twist pairs": `L·U` is isomorphic to the twist of `F0·F8` by `−1` (invariants `I, J` of the quartics);
  `F4·L` and `F4·U` are **2-isogenous to the twist of `F0·F4` (or `F8·F4`) by `−2`**, not twists of them (Vélu
  2-isogenies from the twist, exact, until a curve isomorphic to the target appears).

So the grouping into 15 classes holds exactly. The ranks in the note's table of classes are PARI's.

## D. The corner curve of (627,364) is not a ℚ-curve (tier C: written proof; the valuations computed)

`E/K`, `K = ℚ(i)`, the Jacobian of `y² = L1·L3·M1·M2 = (m − ns)(n + ms)(ps − p̄)(p̄s + p)`, `p = m + in` (labels of
`scripts/second_descent/kclose.py`). `kcurve_exact.py` computes `j` exactly from the invariants of the quartic
(`j = 6912I³/(4I³ − J²)`, the same as PARI's) and its valuations at the two primes above each split prime:

| `p` | `v(j)` at `π`, `π̄` | conductor exponents (PARI, a cross-check) |
|---|---|---|
| 37 | −4, 3 | 1, 0 |
| 181 | 3, −4 | 0, 1 |
| 23490241 | −4, 3 | 1, 0 |
| 395249137273 | 3, −4 | 0, 1 |

(the four primes proved prime by sympy's deterministic test below `2⁶⁴`; 5 and 29 are symmetric, 263 and 991 inert).

**Criterion.** If `v_𝔭(j(E)) < 0 ≤ v_𝔭̄(j(E))` for a prime `𝔭` of `K`, then `E` is not isogenous over `K̄` to its
conjugate `E^c`, so `E` is not a ℚ-curve (and not CM).
*Proof.* `j(E^c) = c(j(E))`, so `v_𝔭(j(E^c)) = v_𝔭̄(j(E)) ≥ 0`: `E^c` has potentially good reduction at `𝔭`. If
`μ: E → E^c` were an isogeny over a finite `L/K`, take `L′/L` finite with `E^c` of good reduction at a prime `𝔓′`
above `𝔭`; then `E`, isogenous to it over `L′`, has good reduction at `𝔓′` (Néron–Ogg–Shafarevich; Serre–Tate;
Silverman AEC VII.7.2, quoted), so `v_{𝔓′}(j(E)) ≥ 0`, against `v_𝔭(j(E)) < 0`. CM `j`-invariants are algebraic
integers. ∎

**On the note's argument.** "Not isogenous to its conjugate (nor to a twist of it by a character unramified there)"
leaves out the twists ramified at those primes. A ℚ-curve only needs a `K̄`-isogeny, which for non-CM `E` is a
`K`-isogeny to some quadratic twist of `E^c`, ramified twists included. With the data it is closed (a twist of a
curve with potentially good reduction keeps `v(j) ≥ 0`, so it never has conductor exponent 1), and the criterion
above needs neither the conductor nor twists. Reproduced: `session/kcurve_not_qcurve.txt` (same `j`, conjugate `j`
outside the six `j`-invariants of the 2-power class, conductor norm `≈ 5.7·10³⁹`).

What this does and does not give: `E` is not isogenous to a base change from ℚ, so `L(E/K, s)` is not a product of two
`L`-functions over ℚ, and Ribet's GL₂-type route does not apply. The note's further sentences (no theorem "`L ≠ 0`
gives rank 0" for such curves; about `10²⁰` terms) are its reading of the literature and an estimate, not checked.

## What rests on what

- **sympy only** (no PARI): A (`genus3_split.py`), C (`genus3_isog.py`), B′ (`genus3_mult.py`), the criterion of D.
- **PARI as a finder**, every output checked exactly: the points and torsion points of B (`ellrank`, `ellisomat`,
  `ellisogenyapply`, `elltors`); factorisations and conic points inside `rank_second.py`.
- **PARI's results, reported as PARI's**: rank intervals, root numbers, the conductor of D (cross-check only).
- **Quoted theory**: Mazur, Nagell–Lutz, the 2-descent map, the Cassels–Tate pairing (as in `rank_second.py`),
  Néron–Ogg–Shafarevich, `π*π* = Σ h*`, Kani–Rosen, Dokchitser–Dokchitser 2-parity (only for the conditional remark).

## Reproduce (times from `times.txt`, one core, other jobs running)

```
python3 -I genus3_split.py                 #  8 s   genus3_split.txt
python3 -I genus3_points.py                # 13 s   genus3_points.txt (reads genus3_points_hints.json)
python3 -I genus3_points.py --negative     #  2 s   9 forged inputs refused
python3 -I genus3_isog.py                  #  8 s   genus3_isog_627_364.txt
python3 -I genus3_mult.py                  # 53 s   genus3_mult.txt
python3 -I genus3_rank_upper.py            #  5 min genus3_rank_upper.txt (rank_second.py, MSS_REPO or ../..)
python3 -I kcurve_exact.py                 #  1 s   kcurve_exact.txt
python3 -I kcurve_exact.py --negative      #  2 s   5 forged inputs refused
cd session; python3 -I genus3_ranks.py     # 34 s   (and the three other session scripts, 5 s together)
bash run_genus3.sh                         # all of the above against the committed outputs; times.txt
```
Each script prints exactly its committed output (add `--write` to rewrite it); PARI's stack message goes to stderr.
