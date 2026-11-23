# A second descent for the ratios with a degenerate solution (experiment)

**Status.** An experiment made after release 0.4.0. Nothing here is counted in the census, in the README numbers or in any release. It rests on the first descent of `GAUSS_DESCENT.md`, which is itself unreviewed, on rank bounds computed by PARI, and on torsion and polynomial roots over `ℚ(i)` computed by Sage. Nobody other than the author of the scripts (an AI model) has read the argument. Treat every number below as a claim to be checked.

**What it would give, if both descents are right.** Of the 26 ratios of the census (`m ≤ 1000`) that are open, 2 have no degenerate solution and are eliminated by the first descent (and, independently, by Magma's two-cover descent; `GAUSS_DESCENT.md`, section 4). Of the other 24, this second step closes 22: 20 with rank bounds from descents, and (440,279) with its partner (719,161) with a rank 0 that is read from an `L`-value (section 6). Two are left, one pair of partners: **(627,364) and (991,263)**. The smallest open `m` would be 627.

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

- Closed directly: 89 (60 and 29). With partners: **101 of 123**.
- The first curve that closed, in the order of the family: `L1·L3·L5·L6` in 59 cases (`Y² = (m − ns)(n + ms)(1 − s²)`, `d = 1` for the square type), `q₊·L1·L4` in 17, `q₊·L2·L3` in 13. In all, 127 curves close a ratio: 57 ratios are closed by one curve, 26 by two, 6 by three, and they agree.
- The equation is an identity for a torsion point on 110 curves, all of them even quartics; no even quartic reaches the equation step without one. So no closure was missed by stopping early.
- The summand of rank 0 was `E(ℚ)` in 38 cases and `E⁽⁻¹⁾(ℚ)` in 51.
- The only rational `t` returned by the equations, on all 127 curves, were `0`, `1`, `−1`. None is a non-degenerate solution.
- Not closed: 22 ratios. On the curves of each, and of the partner when it is in the list, no usable summand has rank bound 0.

Against the census. 99 of the 123 were already closed there (rank-1 method or the four-point tier); the second descent closes 81 of these again by a different route, and fails on 18. Of the **24 open ratios with a degenerate solution it closes 20**:

    91,60  151,31  195,28  209,120  223,167  329,89  420,341  572,555  608,105  609,200
    680,111  713,503  741,580  760,39  761,79  791,569  799,721  809,409  861,620  945,248

and leaves 4 at this stage: `440,279`, `719,161`, `627,364`, `991,263`. For these all 24 curves were computed for the ratio and for its partner. On the usable ones the 2-descent bounds of the summands are between 1 and 4. Two even quartics of (719,161) and of (991,263) have a summand of rank 0, but for an even quartic the equation is an identity (section 3). Section 6 closes the first two of the four in another way.

## 5. Checks made

- **The degenerate solutions come back.** For each of the 60 closed ratios of square type, the equations return `t = 1` and `t = −1`. These are not put in by hand; they are points with `σP = ±P + T`, so the equation step finds points that are known to be there.
- **Real points that are not degenerate** (`planted.py`, `planted.txt`). For any coprime `m > n` and any rational `t0` with `F0(t0)·F8(t0)` a square (a rational point of the classical curve `Y² = (m² + n²t²)(n² + m²t²)`), the number `Q_b(−i·t0)` has a square norm, so it is `d` times a square in `K` for a rational `d`, and `(−i·t0, Y)` is a `K`-point of `Y² = d·Q_b(s)`: the situation of sections 2 and 3 with a point that exists, and without the first descent. The argument predicts that such a curve has no summand of rank 0, unless `t0` is in the finite list. For `m ≤ 30` and `t0 = u/v`, `u, v ≤ 60`, there are 164 such points on 82 curves: **none has a summand of rank 0.** Negative control: with `d` replaced by `3d`, `5d`, `7d` (the point is then not on the curve) a summand of rank 0 appears for 96, 94 and 76 of the 164 points, and the solver does not return `t0`, as it should not. So the curve, the class `d` and the rank bounds, taken together, see real points.
- **The polynomial.** For a random `x1 ∈ K`, the product `(x(P+T) − g(x1))·(x(P−T) − g(x1))`, computed with Sage's group law over `K(√f(x1))`, equals the polynomial of the script divided by the stated denominator (`second_control.py`, four ratios, all torsion points not of order 2).
- **Rank 0, a second program.** For all 89 closing curves Magma gives the upper bound 0 as well (public calculator, version 2.29-10; 37 requests, verbatim in `magma_second_raw.txt`): `RankBounds` returns `0 0` for 80 of them, and `TwoPowerIsogenyDescentRankBound` returns 0 for 31, among them the 9 for which `RankBounds` did not finish within the calculator's 60 seconds or returned the weaker `0 2` (two curves). `magma_check.py` checks that every closing curve is in the transcripts with exactly the coefficients of `second_123.jsonl` and with the bound 0; CI runs it. Magma itself cannot be run in CI.
- **Rank 0 from `L(E,1) > 0`, with an error bound** (`lproof.py`, `lproof.txt`). For 53 of the 89 closing curves (conductor up to `10¹⁴`) the series `L(E,1) = 2·Σ aₙ/n·qⁿ`, `q = exp(−2π/√N)`, is summed with exact `aₙ` until the tail, bounded with `|aₙ| ≤ 2n`, is below `10⁻¹⁸`. All 53 values are positive (the smallest is 0.027), and `L(E,1)/Ω` is an integer to nine decimals. With modularity and Kolyvagin's theorem this gives rank 0 without any descent. What it trusts: PARI's conductor, root number and `aₙ`. CI recomputes all 53. This covers 16 of the 20 ratios above, directly or through the partner. (`lcheck.txt` is the earlier version with PARI's `lfun` and no error bound; it gives the same values.) The four not covered are (572,555), (609,200) with its partner (809,409), and (861,620): their curves have conductors of `10²⁴` to `10³⁰`.
- **Root numbers.** All 127 rank-0 curves have root number `+1`, as they must.
- **No conflict with the census.** No candidate is a solution, on ratios where the census already says there is none.

## 6. Beyond the 2-descent: a descent over `ℚ(i)`, and rank 0 from `L`-values

The 24 quartics use the factors of the forms in pairs that are defined over `ℚ`. Over `K` there are ten linear factors (`L1, …, L6` and the four factors `M1, …, M4` of `q₊`, `q₋`), hence 210 quartics, and for each the classes at the two degenerate solutions give two curves `δ·Y² = Q_J(s)` over `K`. A solution is a `K`-point of one of them. These curves have all their points of order 2 over `K`, and if both have rank 0 over `K` the solution is a torsion point, and there are finitely many to test. No equation is needed.

Their ranks cannot be computed with a program for curves over `ℚ`. `kdescent.py` is a 2-descent over `ℚ(i)` written for this: the 2-Selmer group as the `S`-units modulo squares whose classes lie in the local images at all places of `S`, with the local images generated by actual local points until they have their known order (4 at an odd place, 16 at `1 + i`), and the local classes taken from the labels of the first descent. Test (`kdescent_test.py`, `kdescent_test.txt`): on 295 curves `y² = x(x + a)(x + b)` over `ℚ`, where the rank over `K` is the sum of the ranks of the curve and of its twist by `−1` and PARI knows both, the Selmer bound over `K` is never below the rank; it equals the rank for 267 curves and exceeds it by 2 for 28 (an even excess, as the Cassels pairing requires).

Result (`kclose.py`, `kclose_open4.txt`): for (440,279), (719,161), (627,364), (991,263), none of the 210 quartics has both curves of Selmer rank 0; the least bounds are 2 and 3. For (4,3), as a check, six quartics close the ratio again. So this closes none of the four. It also shows a limit of the tool: for curves that come from `ℚ` the 2-Selmer group over `K` can be larger than over `ℚ` (for (15,8) the bound over `K` is 2 although both ranks over `ℚ` are 0), so a 2-descent over `K` alone is weaker than what PARI does over `ℚ` with isogenous curves.

**Rank 0 that the 2-descent cannot see** (`lvalue_close.py`, `lvalue_open.jsonl`, `tally.txt`). A 2-descent bounds the rank by the Selmer rank, which also counts the 2-part of the Tate–Shafarevich group. For (440,279) the twist by `−1` of the Frey curve,

    y² = x·(x − 719²)·(x + 161²),        conductor 1,929,934,048, root number +1,

has 2-descent bound 2 over its whole isogeny class, so it was not used in section 4. But `L(E,1) = 1.104924974026 ± 3·10⁻²¹`, summed over 402,483 terms with the error bound of `lproof.py` (`lproof_440_279.txt`), and `L(E,1)/Ω = 256.000000000`. So the rank is 0, and what the 2-descent sees is a Tate–Shafarevich group. With this summand the equation step returns only `t = 0, ±1`: **(440,279) is closed, and with it its partner (719,161).** The same happens for three ratios that the census had already closed, (612,35), (660,259) and (920,129). In all, 108 of the 123 ratios are closed, and of the 24 open ones only (627,364) and (991,263) are left. These four closures rest on an `L`-value. For the curve of (440,279) there is also an algebraic confirmation: Magma's `TwoPowerIsogenyDescentRankBound` gives 0 for two of its three points of order 2, after five steps (`magma_last_raw.txt`).

For the last pair this does not help. Every curve of the family with root number `+1` either has rank exactly 2 (PARI finds the points; for the two with a small conductor `L(E,1) = 0` to 12 decimals), or, for `q₊·L1·L4` and `q₊·L2·L3`, Magma's `TwoPowerIsogenyDescentRankBound` still gives 2 after five steps (`magma_last_raw.txt`).

Where the last pair could go next, not done:

- **Rank 1 with two maps.** `s → 1/s` (that is `t → −1/t`) exchanges the quartics `q₊·L1·L3` and `q₊·L2·L4` up to the factor `−1`, so `E(q₊L2L4) ≅ E(q₊L1L3)⁽⁻¹⁾`. For (627,364) the ranks are 1 and 2 for the first and 2 and 1 for the second. So the same curve over `ℚ`, of rank at most 1, receives two maps from the curve of solutions: `P + σP` from the first quartic and `P − σP` from the second. Two independent maps into a curve of rank 1 is the situation of the rank-1 method of the census (Demjanenko–Manin); it needs a generator and explicit height bounds for these two maps, which have not been worked out.
- **Chabauty on `W`**, when the three ranks add up to less than 3.

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
    python3 lvalue_close.py 1e13 $(the ratios left open)           # lvalue_open.jsonl; then tally.py
    python3 magma_check.py                                         # the Magma transcripts against second_123.jsonl
    python3 kdescent_test.py 1 300                                 # the 2-descent over Q(i) against PARI, kdescent_test.txt
    python3 wide.py 1000 2000 ../g1_full/modp_classes.txt          # beyond the census, about 70 minutes; wide_1000_2000.jsonl
    python3 kclose.py 440,279 719,161 627,364 991,263              # section 6

CI runs all of these on every commit, in a job of its own, and compares each output with the committed file. Magma itself is not run in CI: the requests were sent to the public calculator by hand, and `magma_check.py` checks that the transcripts contain every closing curve with the answer 0.
