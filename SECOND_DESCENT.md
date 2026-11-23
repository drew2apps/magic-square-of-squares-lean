# A second descent for the ratios with a degenerate solution (experiment)

**Status.** An experiment made after release 0.4.0 and included, as an experiment, in release 0.5.0. Nothing here is counted in the census or in the README numbers. It rests on the first descent of `GAUSS_DESCENT.md`, which is itself unreviewed, on rank bounds computed by PARI, and on torsion and polynomial roots over `ℚ(i)` computed by Sage. Nobody other than the author of the scripts (an AI model) has read the argument. Treat every number below as a claim to be checked.

**What it would give, if both descents are right.** Of the 26 ratios (`m ≤ 1000`) that the methods with elliptic curves leave, 2 have no degenerate solution; they are eliminated by the first descent and, independently, by Magma's two-cover descent, on which the census now counts them as closed (`GAUSS_DESCENT.md`, section 4). Of the other 24, the open ones, this second step closes 22: 20 with rank bounds from descents, and (440,279) with its partner (719,161) with a rank 0 that is read from an `L`-value (section 6). Two are left, one pair of partners: **(627,364) and (991,263)**. The smallest open `m` would be 627.

## 1. Where the first descent stops

Notation as in `GAUSS_DESCENT.md`: ratio `m:n`, `t = u/v` in lowest terms, `S = m² + n²`, `D = m² − n²`, and over `ℤ[i]`

    α0 = mv + i·nu,   α8 = nv + i·mu,   γ4 = u + iv,   γL = (m+in)u − i(m−in)v,   γU = (m+in)u + i(m−in)v.

A solution gives a tuple of five classes in `K*/K*²`, `K = ℚ(i)`, one for each form. A ratio has a degenerate solution exactly when `S` is a square (`t = ±1`) or twice a square (`t = 0, ∞`). These are real points of the system, so no congruence can remove their tuples. For all 123 such ratios in the list of 157, the first descent leaves exactly two tuples and keeps the tuples of both degenerate solutions (`scripts/gauss_descent/gauss_descent_157.txt`). The two degenerate solutions have different tuples in all 123 cases (`degenerate_tuples.txt`), so the two tuples left are exactly theirs.

So the first descent says: **every solution has, form by form, the same class as one of the two degenerate solutions.** The second step uses this.

## 2. The curves

Put `s = −i·t`. The values of the forms and of their conjugates, divided by `v`, are

    L1 = m − ns = α0/v        L2 = m + ns = ᾱ0/v        L3 = n + ms = ᾱ8/v        L4 = n − ms = α8/v
    L5 = 1 − s,  L6 = 1 + s           (L5·L6 = N(γ4)/v²)
    q₊ = S s² − 2D s + S = −γL·γ̄U/v²           q₋ = S s² + 2D s + S = −γ̄L·γU/v²

All of them have rational coefficients in `s`. Take a quartic `Q(s)` that is a product of four of the `L`, or of one `q` and two `L`. At a solution `t0`, with `s0 = −i·t0`, the class of `Q(s0)` in `K*/K*²` is a product of the classes of the tuple and of their conjugates (`v⁴` is a square). By section 1 it equals the class of `Q` at one of the two degenerate solutions. If those two classes are equal, and equal to the class of a rational number `d`, then

    d·Q(s0) is a square in K,

so `(s0, Y)` is a `K`-point of the curve `C: Y² = d·Q(s)`. This curve is defined over `ℚ`, and it has a rational point, a rational root `r` of `Q`. Let `E` be its Jacobian in Weierstrass form with `x = A/(s − r)` (`A` a rational constant), and let `P ∈ E(K)` be the point of the solution.

**The main curve is a Frey curve.** For `Q_b = L1·L3·L5·L6 = (m − ns)(n + ms)(1 − s²)` the Jacobian is, for any `m, n`,

    E_b:  y² = x·(x + (m+n)²)·(x − (m−n)²),        and its twist by −1 is   y² = x·(x − (m+n)²)·(x + (m−n)²)

(the map is `x = (n² − m²)(s + 1)/(s − 1)`, `y = 2(n² − m²)·Y/(s − 1)²`; `quartic_b_frey` in Lean, and checked in Sage for several pairs). When `m² + n² = k²` the class is `d = 1`, because `Q_b(−i) = ((1 − i)(m + in))²` (`quartic_b_degenerate` in Lean), and `(m+n)² + (m−n)² = 2k²`: `E_b` is the Frey curve of the three squares `(m−n)², k², (m+n)²` in arithmetic progression. So for a ratio with `m² + n²` a square the main criterion reads: *if the first descent leaves only the two degenerate tuples, and one of these two curves has rank 0 over `ℚ`, the ratio has no solution.* This is the case for 50 of the 74 such ratios.

**What these curves are.** Let `W` be the curve of the triples `(t, Y, Y′)` with `Y² = d·Q(−it)`, `Y′² = d·Q(it)`. It is defined over `ℚ` (conjugation exchanges `Y` and `Y′`), it has genus 3, and a solution gives a rational point of it. Its Jacobian is isogenous over `ℚ` to `E × E⁽⁻¹⁾ × Jac(D)`, where `D: Z² = Q(−it)·Q(it)` is one of the classical curves in `t` (for `Q_b` it is `Z² = F0·F8`). The census with ranks used `D`. The two factors `E` and `E⁽⁻¹⁾` are new quotients of the curve of solutions; they become visible only after the first descent has fixed the twist `d`.

Up to `s → −s` there are 24 such quartics. For 12 of them the class is not rational or differs at the two degenerate solutions; they are not used.

## 3. Rank 0 on one side

Let `σ` be complex conjugation. `E` is defined over `ℚ`, so `σ` acts on `E(K)`, and `P + σP ∈ E(ℚ)`, while `P − σP` lies in the part of `E(K)` on which `σ` acts as `−1`, which is `E⁽⁻¹⁾(ℚ)`, the twist by `−1`.

- If `E(ℚ)` has rank 0, then `P + σP` is torsion: `σP = −P + T`.
- If `E⁽⁻¹⁾(ℚ)` has rank 0, then `P − σP` is torsion: `σP = P + T`.

In both cases `T` is a torsion point of `E(K)`. Now `t0` is rational, so `σ(s0) = −s0`, and `s` is a function of `x` alone. So

    x(P + T) = g(x(P))   or   x(P − T) = g(x(P)),      g(x) = −A·x / (2r·x + A),

for some torsion point `T` of `E(K)`. For each `T` this is a polynomial equation for `x(P)` (from the addition formula, squared to remove `y`). If the polynomial is not identically zero it has finitely many roots in `K`. Each root gives `s0 = r + A/x`, `t0 = i·s0`; the rational `t0` are tested in the five-form system directly. A superset of candidates is harmless.

For some curves the polynomial is identically zero for one torsion point `T`. Then a solution with that `T` is not constrained at all, so the curve cannot close the ratio, whatever the other torsion points give. In the computation this happens exactly for the three quartics that are even in `s` (`L1L2L3L4`, `L1L2L5L6`, `L3L4L5L6`): there `s → −s` is an automorphism of the curve, of the form `P → P + T1` or `P → −P + T2`, and the equation holds for every `P`. Those curves are the ones already defined over `ℚ` in the variable `t`. The program goes through all torsion points of every curve and records the identities; it does not stop at the first.

If no rational `t0` passes the five-form test except the degenerate ones, the ratio has no solution. A ratio is also closed when its partner `(m+n):(m−n)` (halved if both are even) is closed; a solution for one gives a solution for the other (`partner_forms` in Lean).

## 4. Result

Run on the 123 ratios of the list of 157 that have a degenerate solution: 74 with `S` a square, 49 with `S` twice a square (`second_123.jsonl`, `second_summary.txt`). All 24 curves were run on every ratio, also after a curve had closed it.

- Closed by a rank bound 0 from descent: 89 directly (60 and 29), **101 of 123** with partners. With the four closures from an `L`-value of section 6 added: 93 directly, **108 of 123** with partners (`tally.txt`, which is the count to quote).
- The first curve that closed, in the order of the family: `L1·L3·L5·L6` in 59 cases (`Y² = (m − ns)(n + ms)(1 − s²)`, `d = 1` for the square type), `q₊·L1·L4` in 17, `q₊·L2·L3` in 13. In all, 127 curves close a ratio: 57 ratios are closed by one curve, 26 by two, 6 by three, and they agree.
- The equation is an identity for a torsion point on 110 curves, all of them even quartics; no even quartic reaches the equation step without one. So no closure was missed by stopping early.
- The summand of rank 0 was `E(ℚ)` in 38 cases and `E⁽⁻¹⁾(ℚ)` in 51.
- The only rational `t` returned by the equations, on all 127 curves, were `0`, `1`, `−1`. None is a non-degenerate solution.
- Not closed by descent bounds: 22 ratios (15 after the `L`-values). On the curves of each, and of the partner when it is in the list, no usable summand has rank bound 0.

Against the census. 99 of the 123 were already closed there (rank-1 method or the four-point tier); the second descent closes 81 of these again by a different route, and fails on 18. Of the **24 open ratios with a degenerate solution it closes 20**:

    91,60  151,31  195,28  209,120  223,167  329,89  420,341  572,555  608,105  609,200
    680,111  713,503  741,580  760,39  761,79  791,569  799,721  809,409  861,620  945,248

and leaves 4 at this stage: `440,279`, `719,161`, `627,364`, `991,263`. For these all 24 curves were computed for the ratio and for its partner. On the usable ones the 2-descent bounds of the summands are between 1 and 4. Two even quartics of (719,161) and of (991,263) have a summand of rank 0, but for an even quartic the equation is an identity (section 3). Section 6 closes the first two of the four in another way.

## 5. Checks made

- **The degenerate solutions come back.** For each of the 60 closed ratios of square type, the equations return `t = 1` and `t = −1`. These are not put in by hand; they are points with `σP = ±P + T`, so the equation step finds points that are known to be there.
- **Real points that are not degenerate** (`planted.py`, `planted.txt`). For any coprime `m > n` and any rational `t0` with `F0(t0)·F8(t0)` a square (a rational point of the classical curve `Y² = (m² + n²t²)(n² + m²t²)`), the number `Q_b(−i·t0)` has a square norm, so it is `d` times a square in `K` for a rational `d`, and `(−i·t0, Y)` is a `K`-point of `Y² = d·Q_b(s)`: the situation of sections 2 and 3 with a point that exists, and without the first descent. The argument predicts that such a curve has no summand of rank 0, unless `t0` is in the finite list. For `m ≤ 30` and `t0 = u/v`, `u, v ≤ 60`, there are 164 such points on 82 curves: **none has a summand of rank 0.** Negative control: with `d` replaced by `3d`, `5d`, `7d` (the point is then not on the curve) a summand of rank 0 appears for 96, 94 and 76 of the 164 points, and the solver does not return `t0`, as it should not. So the curve, the class `d` and the rank bounds, taken together, see real points.
- **The polynomial.** For a random `x1 ∈ K`, the product `(x(P+T) − g(x1))·(x(P−T) − g(x1))`, computed with Sage's group law over `K(√f(x1))`, equals the polynomial of the script divided by the stated denominator (`second_control.py`, four ratios, all torsion points not of order 2).
- **Rank 0, a second program.** For all 89 closing curves Magma gives the upper bound 0 as well (public calculator, version 2.29-10; 37 requests, verbatim in `magma_second_raw.txt`): `RankBounds` returns `0 0` for 80 of them, and `TwoPowerIsogenyDescentRankBound` returns 0 for 31, among them the 9 for which `RankBounds` did not finish within the calculator's 60 seconds or returned the weaker `0 2` (two curves). `magma_check.py` checks that every closing curve is in the transcripts with exactly the coefficients of `second_123.jsonl` and with the bound 0; CI runs it. Magma itself cannot be run in CI.
- **Rank 0 from `L(E,1) > 0`, with an error bound** (`lproof.py`, `lproof.txt`). For 53 of the 89 closing curves (conductor up to `10¹⁴`) the series `L(E,1) = 2·Σ aₙ/n·qⁿ`, `q = exp(−2π/√N)`, is summed with exact `aₙ` until the tail, bounded with `|aₙ| ≤ 2n`, is below `10⁻¹⁸`. All 53 values are positive (the smallest is 0.027), and `L(E,1)/Ω` is an integer to nine decimals. With modularity and Kolyvagin's theorem this gives rank 0 without any descent. It is a numerical computation with an explicit error bound, not interval arithmetic with directed rounding: the tail is bounded rigorously, and the rounding error is bounded under a stated assumption on the relative error of PARI's real operations (derived in the docstring of `lproof.py`; below `10⁻⁴⁰` here). What it trusts: that assumption, and PARI's conductor, root number and `aₙ`. CI recomputes all 53. This covers 16 of the 20 ratios above, directly or through the partner. (`lcheck.txt` is the earlier version with PARI's `lfun` and no error bound; it gives the same values.) The four not covered are (572,555), (609,200) with its partner (809,409), and (861,620): their curves have conductors of `10²⁴` to `10³⁰`. The same 53 values have since been enclosed in ball arithmetic (`lball.py`, `lball.txt`; Arb through python-flint, 128 bits, every operation enclosed, tail below `10⁻⁹`): all 53 lower ends are positive, and the enclosures contain the values of `lproof.py`. That removes the assumption on PARI's rounding; what is still trusted is PARI's conductor, root number and `aₙ`, Arb, modularity and Kolyvagin's theorem.
- **Root numbers.** All 127 rank-0 curves have root number `+1`, as they must.
- **No conflict with the census.** No candidate is a solution, on ratios where the census already says there is none.

## 6. Beyond the 2-descent: a descent over `ℚ(i)`, and rank 0 from `L`-values

The 24 quartics use the factors of the forms in pairs that are defined over `ℚ`. Over `K` there are ten linear factors (`L1, …, L6` and the four factors `M1, …, M4` of `q₊`, `q₋`), hence 210 quartics, and for each the classes at the two degenerate solutions give two curves `δ·Y² = Q_J(s)` over `K`. A solution is a `K`-point of one of them. These curves have all their points of order 2 over `K`, and if both have rank 0 over `K` the solution is a torsion point, and there are finitely many to test. No equation is needed.

Their ranks cannot be computed with a program for curves over `ℚ`. `kdescent.py` is a 2-descent over `ℚ(i)` written for this: the 2-Selmer group as the `S`-units modulo squares whose classes lie in the local images at all places of `S`, with the local images generated by actual local points until they have their known order (4 at an odd place, 16 at `1 + i`), and the local classes taken from the labels of the first descent. Test (`kdescent_test.py`, `kdescent_test.txt`): on 295 curves `y² = x(x + a)(x + b)` over `ℚ`, where the rank over `K` is the sum of the ranks of the curve and of its twist by `−1` and PARI knows both, the Selmer bound over `K` is never below the rank; it equals the rank for 267 curves and exceeds it by 2 for 28 (an even excess, as the Cassels pairing requires).

Result (`kclose.py`, `kclose_open4.txt`): for (440,279), (719,161), (627,364), (991,263), none of the 210 quartics has both curves of Selmer rank 0; the least bounds are 2 and 3. For (4,3), as a check, six quartics close the ratio again. So this closes none of the four. It also shows a limit of the tool: for curves that come from `ℚ` the 2-Selmer group over `K` can be larger than over `ℚ` (for (15,8) the bound over `K` is 2 although both ranks over `ℚ` are 0), so a 2-descent over `K` alone is weaker than what PARI does over `ℚ` with isogenous curves.

**Rank 0 that the 2-descent cannot see** (`lvalue_close.py`, `lvalue_open.jsonl`, `tally.txt`). A 2-descent bounds the rank by the Selmer rank, which also counts the 2-part of the Tate–Shafarevich group. For (440,279) the twist by `−1` of the Frey curve,

    y² = x·(x − 719²)·(x + 161²),        conductor 1,929,934,048, root number +1,

has 2-descent bound 2 over its whole isogeny class, so it was not used in section 4. But `L(E,1) = 1.104924974026 ± 3·10⁻²¹`, summed over 402,483 terms with the error bound of `lproof.py` (`lproof_440_279.txt`); in ball arithmetic, with no assumption on rounding, `1.104924973 ≤ L(E,1) ≤ 1.104924976` (`lball.py`, `lball_440_279.txt`), and `L(E,1)/Ω = 256.000000000`. So the rank is 0, and what the 2-descent sees is a Tate–Shafarevich group. With this summand the equation step returns only `t = 0, ±1`: **(440,279) is closed, and with it its partner (719,161).** The same happens for three ratios that the census had already closed, (612,35), (660,259) and (920,129). In all, 108 of the 123 ratios are closed, and of the 24 open ones only (627,364) and (991,263) are left. These four closures rest on an `L`-value. For the curve of (440,279) there is also an algebraic confirmation: Magma's `TwoPowerIsogenyDescentRankBound` gives 0 for two of its three points of order 2, after five steps (`magma_last_raw.txt`).

For the last pair this does not help. Every curve of the family with root number `+1` either has rank exactly 2 (PARI finds the points; for the two with a small conductor `L(E,1) = 0` to 12 decimals), or, for `q₊·L1·L4` and `q₊·L2·L3`, Magma's `TwoPowerIsogenyDescentRankBound` still gives 2 after five steps (`magma_last_raw.txt`).

Where the last pair could go next, not done:

- **Rank 1 with two maps.** `s → 1/s` (that is `t → −1/t`) exchanges the quartics `q₊·L1·L3` and `q₊·L2·L4` up to the factor `−1`, so `E(q₊L2L4) ≅ E(q₊L1L3)⁽⁻¹⁾`. For (627,364) the ranks are 1 and 2 for the first and 2 and 1 for the second. So the same curve over `ℚ`, of rank at most 1, receives two maps from the curve of solutions: `P + σP` from the first quartic and `P − σP` from the second. Two independent maps into a curve of rank 1 is the situation of the rank-1 method of the census (Demjanenko–Manin); it needs a generator and explicit height bounds for these two maps, which have not been worked out.
- **Chabauty on `W`**, when the three ranks add up to less than 3.

**Where the last pair stands (23 November 2026, `kdist.py`, `kdist_627_364.txt`).** To see whether a curve over `ℚ(i)` of rank at most 1 exists, on which an elliptic-curve Chabauty argument would apply (rank below the degree of `ℚ(i)`), the 2-Selmer bounds over `ℚ(i)` of all 210 quartics were listed for (627,364), taking for each quartic the larger of the bounds of its curves: 2 for 4 quartics, 3 for 8, 4 for 18, and 5 to 11 for the other 180. None is 0 or 1. The four with bound 2 are `L1·L3·M1·M2`, `L1·L3·M3·M4`, `L2·L4·M1·M2`, `L2·L4·M3·M4`, that is `q₊·L1·L3`, `q₋·L1·L3`, `q₊·L2·L4`, `q₋·L2·L4`: curves defined over `ℚ`, whose bound 2 is the sum of the ranks 1 and 1 (or 2 and 0 with a Tate–Shafarevich group) of the curve and its twist. So, as far as 2-descent sees, every route through a single curve of genus 1 is used up, over `ℚ` and over `ℚ(i)`. What is left, in the order in which we would try it:

1. *Rank 1 with two maps* (the first item above). It uses machinery that the census already has (`dm.py`, the height bounds), on a curve over `ℚ` of rank 1. The two maps come from quartics exchanged by `s → 1/s`, so they should have the same degree; then the height argument forces `P + σP = ±(P′ − σP′) + T` for all solutions above a height bound, which is an equation in `t` of the kind solved in section 3. Needed: a generator of the rank-1 curve (its coefficients have about 20 digits, so finding it is the first risk), and explicit bounds for the difference between the canonical height and the degree times the height of `t` for the two maps.
2. *Rank 0 hidden behind a Tate–Shafarevich group* on one of the curves with bound 2, as for (440,279): needs `L(E,1)` for a conductor far beyond what `lproof.py` can sum, or a 4-descent, which the public calculator does not finish.
3. *Chabauty or a Mordell–Weil sieve on a quotient of genus 2 or 3* whose Jacobian has rank below its genus: the quotients whose Jacobian splits into the elliptic curves above all have rank equal to or above the genus here, so this needs a quotient that has not been looked at.

**A first result on route 1 (23 November 2026; `mw_sieve.py`, `mw_sieve_627_364.txt`, `fibre.py`, `fibre_627_364.txt`).** For (627,364), with `d = 521266`:

- `E₁`, the Jacobian of `Y² = d·q₊·L1·L3`, has rank exactly 1 (PARI and Magma: bounds 1 and 1), torsion `{O, T}`, and a point `G` of canonical height 20.007 that generates it modulo torsion (`saturation.py`). The feared difficulty, finding a generator on a curve of conductor about `2·10⁴²`, does not arise, for a reason: **`G` is the point of the degenerate solution**, `G = P + σP` at `t = ±1`.
- The Jacobian of `Y² = d·q₊·L2·L4` is the twist of `E₁` by `−1`. So a solution `t` gives two points of `E₁(ℚ)`: `A = P + σP` from the first curve and `B = R − σR` from the second, `A = ±(aG + εT)`, `B = ±(bG + ε′T)` with integers `a, b ≥ 0`.
- The two maps are independent: modulo 10009 and 20089, among 1,208 and 2,420 values of `t` only 8 and 4 satisfy a relation `A = ±B + torsion` for some choice of signs, which is what chance gives.
- *Sieve.* For a prime `p ≡ 1 (mod 4)` the pairs `(a, b)` modulo the order of `G` in `E₁(𝔽_p)` that can come from some `t ∈ ℙ¹(𝔽_p)` form a small set (1,632 of 695,556 pairs for `p = 10009`). With the four primes 10009, 20089, 30013, 40009, **the only pair with `0 ≤ a, b ≤ 596,309` that survives is `(1, 1)`**, the pair of the degenerate solution, which has to survive.
- *The fibre over `(1, ·)` and `(·, 1)`.* `A = ±G` or `±G + T` means `σP = c − P` for one of four constants `c`, a polynomial equation for `x(P)` over `ℚ(i)` of the kind solved in section 3; the same for `B`. Its only rational solutions are `t = ±1`.

*A larger range* (`mw_sieve_big.py`, `mw_sieve_big_627_364.txt`): with the six primes 10009, 20089, 30013, 40009, 50021, 60013 the only pair with `0 ≤ a, b < 59,631,000` that survives is again `(1, 1)`.

So, if the two descents and these programs are right: **a non-degenerate solution for (627,364) has `a` or `b` at least 59,631,000**, that is, one of the points `A`, `B` has canonical height above `7·10¹⁶` (PARI's normalisation).

**From the sieve to a proposed proof (`LAST_PAIR.md`).** The pairs `(A, B)` of solutions are rational points of a curve `C ⊂ E₁ × E₁`, the image of the curve of solutions. For an irreducible curve of genus at least 2 in `E²`, with `E` without complex multiplication, Checcoli, Veneziano and Viada (*The explicit Mordell conjecture for families of curves*, Forum Math. Sigma 7, 2019, Theorem 1.2) bound the height of every point of rank at most 1:

    ĥ(P) ≤ C₁·h₂(C)·deg C + C₂(E)·(deg C)² + C₃(E),        C₁ = 72.251.

`LAST_PAIR.md` writes the whole argument as a theorem with seven lemmas and gives every number; `closure_estimate.py` computes them. In short:

- *The curve.* The curve of solutions has two components, each a cover of the line of `t` with group `(ℤ/2)³`, branched at 8 points, of **genus 9** (the first version of this section said 16 branch points and genus 25; that was wrong, and a reviewer found it). The four differentials `dt/yᵢ` belong to four different characters of the group, so no combination `uA + vB` is constant, and `C` is not in a translate of a proper algebraic subgroup.
- *Degree and height.* `deg C ≤ 96` and `h₂(C) ≤ 7.78·10⁵` (from the points above roots of unity and Zhang's inequality; logarithms rounded up).
- *The bound.* Every point of `C` with coordinates in `E₁(ℚ)` has `ĥ ≤ 6.03·10⁹`, that is `a, b ≤ 14,174`. The sieve leaves only `(1, 1)` below 59,631,000 (already below 596,309 with four primes), and the fibre over `(1, 1)` is `t = ±1`.
- *Tests on real values.* The degenerate solution; a planted solution that is not degenerate, on twisted systems, which survives the whole sieve with exactly its computed pair `(2, 1)` (`planted_control.py`); the whole route on all 24 ratios of square type where `E₁` has rank 1 (`route_all.py`), which gives the known answer on the 23 closed by other means; `G` against the height comparison of the paper.
- *Second programs.* For `E₁`, Magma gives rank bounds 1 and 1, torsion of order 2, the same height of `G`, and `G` saturated at the primes up to 200 (`magma_last_pair_raw.txt`); the statements of the paper were read twice with an automated reader, with the same result.

*For every ratio, and twelve more.* The construction does not depend on the ratio (`LAST_PAIR.md`, section 9; `lemma2_general.py` checks it as identities in `m, n, d, t`): the four polynomials are `d·A·B, d·A′·B′, d·A·B′, d·A′·B` with `A·A′ = GL·GU/4` and `B·B′ = F0·F8`, the curve of solutions always has genus 9, and it is an unramified double cover of the classical curve `u² = F0·F8, v² = GL·GU` of genus 5, on whose rational points the class `d` is always a rational number. Run on the ratios left with `1000 < m ≤ 3000`, the route closes 12 more of square type and, with partners, 18 of the 67 (`LAST_PAIR.md`, section 10; `route_wide.txt`); 49 are left there.

**If every step is right, (627,364) and its partner (991,263) have no solution. Nobody has read the argument, and the census is not changed.** The order in which to check it is in section 7 of `LAST_PAIR.md`; the first item is the three statements taken from the paper, the second the reduction argument of the sieve.

*How far the route reaches* (`two_map_scope.py`, `two_map_scope.txt`). For all 74 ratios of square type in the list of 123, the point `G = P + σP` of the degenerate solution has infinite order on `E₁` (so this curve never has rank 0, which is why section 4 never used it on this side), and `E₁` has rank exactly 1 for 24 of them, (627,364) and (440,279) among them. For those 24 the sieve applies as it stands, with a generator that costs nothing to find, and it has been run on all of them (`route_all.txt`; `LAST_PAIR.md`, section 6).

What it rests on, beyond the two descents: PARI's rank bound for `E₁`; that `G` generates `E₁(ℚ)` modulo torsion (Sage's lower bound for the height, and divisibility by the primes up to 19); the reduction argument of `mw_sieve.py` (points of `C₁`, `C₂` over `𝔽_p × 𝔽_p` with every choice of sign of the square roots, the points with `s = r` or `s = ∞` included); and the identification of the anti-invariant points of `E₂(ℚ(i))` with `E₁(ℚ)`. Nobody has read any of it.

**Beyond the census** (`wide.py`; `wide_1000_2000.jsonl`, `wide_2000_2500.jsonl`, `wide_2500_3000.jsonl` and their summaries and tallies). Both descents were run on the ratios with a degenerate solution and `1000 < m ≤ 3000` that avoid the 484 classes: 604 ratios. The first descent leaves exactly the two degenerate tuples for all 604, and for the 169 partners that were needed. The second descent, with rank bounds from PARI only, closes 516. Rank 0 from an `L`-value closes 21 more, for the ratio or for its partner. 67 are left (32 up to 2000, 21 up to 2500, 14 up to 3000). The only rational `t` returned are again `0, ±1`. With the scan of `GAUSS_DESCENT.md` this means: if both descents are right, a magic square of squares with a ratio of height at most 3000 has that ratio in a list of 69 (the 2 left below 1000 and these 67). The rank bounds of these runs were not checked with Magma.

## 7. What could be wrong, and what is not done

- **Everything depends on the first descent.** If a tuple was wrongly removed there, the class `d` here is wrong and the conclusion fails. `GAUSS_DESCENT.md` lists what is unproved in that step.
- **Rank bounds** are PARI's (`ellrank`, effort 1, least bound over the isogeny class). Magma agrees for all 89 closing curves (transcripts of the public calculator; not reproducible in CI); `L(E,1) > 0` is a third route for 53 of them. The other 38 curves that also close a ratio in the run of all curves were checked with PARI only.
- **Torsion of `E(K)` and roots over `K`** are Sage's. No certificate is stored.
- **The argument in sections 2 and 3** was written by the same AI model that wrote the scripts. The steps to check: that the class of `Q(s0)` is determined by the tuple; that `C` and `E` are as stated; the equation for `x(P)`; the special points (`P = O`, `x = 0`, `2r·x + A = 0`). The same holds for `kdescent.py`, which is new code for a standard method.
- **Not done:** the two routes at the end of section 6. The six unclosed ratios whose partner has `m > 1000` were not run on the partner. In Lean there are only the identities (`GaussDescent.lean`: the factorisations of `Q_b` and `q₊` and the square at `s = −i`); no rank, no torsion and no closure.
- **The family of 24 quartics is a choice, not a claim.** Each curve is used on its own: a curve that closes a ratio closes it whatever the other curves do. Nothing depends on the family being complete; other curves could close more.
- **Known technique.** Passing to covers and to points over a quadratic field is standard (covering techniques of Bruin, and of Flynn and Wetherell; named from memory, not checked). No claim that the method is new, and the literature on this system was not searched for it.

## 8. Reproduce

    pip install -r scripts/g1_full/requirements.txt
    cd scripts/second_descent
    python3 second_descent.py --all $(cat pairs123.txt) > out.jsonl    # about 10 minutes on one core
    diff out.jsonl second_123.jsonl
    python3 second_descent.py --summary second_123.jsonl           # second_summary.txt
    python3 second_control.py                                      # second_control.txt
    python3 degenerate_tuples.py pairs123.txt                      # degenerate_tuples.txt
    python3 planted.py 30 60; python3 planted.py 30 60 3           # planted.txt (with the negative controls 3, 5, 7)
    python3 lproof.py 1e14 second_123.jsonl                        # lproof.txt, about 4 minutes, 3 GB of memory
    python3 lball.py 1e14 second_123.jsonl                         # lball.txt: the same in ball arithmetic (needs python-flint)
    python3 lvalue_close.py 1e13 $(the ratios left open)           # lvalue_open.jsonl; then tally.py
    python3 magma_check.py                                         # the Magma transcripts against second_123.jsonl
    python3 kdescent_test.py 1 300                                 # the 2-descent over Q(i) against PARI, kdescent_test.txt
    python3 wide.py 1000 2000 ../g1_full/modp_classes.txt          # beyond the census, about 70 minutes; wide_1000_2000.jsonl
    python3 kclose.py 440,279 719,161 627,364 991,263              # section 6

The full CI run (job `descents`: on a pull request, a tag, by hand, or with `[full ci]` in the commit message; see the README) runs all of these and compares each output with the committed file; the job `quick`, on every push, reruns three ratios and five `L`-values and the checkers. The job `last-pair` of the full run reruns everything in section 8 of `LAST_PAIR.md` (`saturation.py`, `closure_estimate.py`, `mw_sieve.py`, `mw_sieve_big.py`, `fibre.py`, `planted_control.py`). Committed outputs that no CI job regenerates: `kclose_open4.txt`, `kdist_627_364.txt`, `two_map_scope.txt`, `two_maps_independent.txt`; the wide runs beyond the first five ratios (`wide_*.jsonl`, their `L`-value files and tallies); `magma_last_raw.txt` and `magma_last_pair_raw.txt`, which have no checker and are quoted by hand; in `scripts/gauss_descent`, `subsets_34.txt`, `real_points_681.txt` and the scan beyond `m = 1100`. They are data, reproducible with the commands given, and an error in them would not be caught automatically. Magma itself is not run in CI: the requests were sent to the public calculator by hand, and `magma_check.py` checks that the transcripts contain every closing curve with the answer 0.
