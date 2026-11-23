# A second descent for the ratios with a degenerate solution (experiment)

**Status.** An experiment made after release 0.4.0. Nothing here is counted in the census, in the README numbers or in any release. It rests on the first descent of `GAUSS_DESCENT.md`, which is itself unreviewed, on rank bounds computed by PARI, and on torsion and polynomial roots over `ℚ(i)` computed by Sage. Nobody other than the author of the scripts (an AI model) has read the argument. Treat every number below as a claim to be checked.

**What it would give, if both descents are right.** Of the 26 ratios of the census (`m ≤ 1000`) that are open, 2 have no degenerate solution and are eliminated by the first descent. Of the other 24, this second step closes 20. Four are left, which are two pairs of partners: (440,279) with (719,161), and (627,364) with (991,263). The smallest open `m` would be 440.

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

Up to `s → −s` there are 24 such quartics. For 12 of them the class is not rational or differs at the two degenerate solutions; they are not used.

## 3. Rank 0 on one side

Let `σ` be complex conjugation. `E` is defined over `ℚ`, so `σ` acts on `E(K)`, and `P + σP ∈ E(ℚ)`, while `P − σP` lies in the part of `E(K)` on which `σ` acts as `−1`, which is `E⁽⁻¹⁾(ℚ)`, the twist by `−1`.

- If `E(ℚ)` has rank 0, then `P + σP` is torsion: `σP = −P + T`.
- If `E⁽⁻¹⁾(ℚ)` has rank 0, then `P − σP` is torsion: `σP = P + T`.

In both cases `T` is a torsion point of `E(K)`. Now `t0` is rational, so `σ(s0) = −s0`, and `s` is a function of `x` alone. So

    x(P + T) = g(x(P))   or   x(P − T) = g(x(P)),      g(x) = −A·x / (2r·x + A),

for some torsion point `T` of `E(K)`. For each `T` this is a polynomial equation for `x(P)` (from the addition formula, squared to remove `y`). If the polynomial is not identically zero it has finitely many roots in `K`. Each root gives `s0 = r + A/x`, `t0 = i·s0`; the rational `t0` are tested in the five-form system directly. A superset of candidates is harmless.

The polynomial is identically zero when `Q` is even in `s`: then `s → −s` is an automorphism of the curve, and the equation says nothing. Those curves are the ones already defined over `ℚ` in the variable `t`.

If no rational `t0` passes the five-form test except the degenerate ones, the ratio has no solution. A ratio is also closed when its partner `(m+n):(m−n)` (halved if both are even) is closed; a solution for one gives a solution for the other (`partner_forms` in Lean).

## 4. Result

Run on the 123 ratios of the list of 157 that have a degenerate solution: 74 with `S` a square, 49 with `S` twice a square (`second_123.jsonl`, `second_summary.txt`).

- Closed directly: 89 (60 and 29). With partners: **101 of 123**.
- The curve that closed: `L1·L3·L5·L6` in 59 cases (`Y² = (m − ns)(n + ms)(1 − s²)`, `d = 1` for the square type), `q₊·L1·L4` in 17, `q₊·L2·L3` in 13.
- The summand of rank 0 was `E(ℚ)` in 38 cases and `E⁽⁻¹⁾(ℚ)` in 51.
- The only rational `t` returned by the equations were `0`, `1`, `−1`. None is a non-degenerate solution.
- Not closed: 22 ratios. On the curves of each, and of the partner when it is in the list, no usable summand has rank bound 0.

Against the census. 99 of the 123 were already closed there (rank-1 method or the four-point tier); the second descent closes 81 of these again by a different route, and fails on 18. Of the **24 open ratios with a degenerate solution it closes 20**:

    91,60  151,31  195,28  209,120  223,167  329,89  420,341  572,555  608,105  609,200
    680,111  713,503  741,580  760,39  761,79  791,569  799,721  809,409  861,620  945,248

and leaves 4: `440,279`, `719,161`, `627,364`, `991,263`. For these all 24 curves were computed for the ratio and for its partner. On the usable ones the rank bounds of the summands are between 1 and 4. Two even quartics of (719,161) and of (991,263) have a summand of rank 0, but for an even quartic the equation is an identity (section 3).

## 5. Checks made

- **The degenerate solutions come back.** For each of the 60 closed ratios of square type, the equations return `t = 1` and `t = −1`. These are not put in by hand; they are points with `σP = ±P + T`, so the equation step finds points that are known to be there.
- **The polynomial.** For a random `x1 ∈ K`, the product `(x(P+T) − g(x1))·(x(P−T) − g(x1))`, computed with Sage's group law over `K(√f(x1))`, equals the polynomial of the script divided by the stated denominator (`second_control.py`, four ratios, all torsion points not of order 2).
- **Rank 0, second route.** For 53 of the 89 rank-0 curves (conductor up to `10¹⁴`), PARI's `L(E,1)` is not zero, and `L(E,1)/Ω` is an integer between 4 and 32768 to six decimals (`lcheck.txt`). With Kolyvagin's theorem that gives rank 0 without the 2-descent. This covers 16 of the 20 ratios above, directly or through the partner. The four not covered are (572,555), (609,200) with its partner (809,409), and (861,620): their curves have conductors of `10²⁴` to `10³⁰`, and they rest on PARI's `ellrank` alone.
- **Root numbers.** All 89 rank-0 curves have root number `+1`, as they must.
- **No conflict with the census.** No candidate is a solution, on ratios where the census already says there is none.

## 6. What could be wrong, and what is not done

- **Everything depends on the first descent.** If a tuple was wrongly removed there, the class `d` here is wrong and the conclusion fails. `GAUSS_DESCENT.md` lists what is unproved in that step.
- **Rank bounds** are PARI's (`ellrank`, effort 1, least bound over the isogeny class). No second program was run for the upper bounds; the `L`-values are a second route for 53 curves only.
- **Torsion of `E(K)` and roots over `K`** are Sage's. No certificate is stored.
- **The argument in sections 2 and 3** was written by the same AI model that wrote the scripts. The steps to check: that the class of `Q(s0)` is determined by the tuple; that `C` and `E` are as stated; the equation for `x(P)`; the special points (`P = O`, `x = 0`, `2r·x + A = 0`).
- **Not done:** a 2-descent over `ℚ(i)` for curves that are not defined over `ℚ` (products that use a single factor of `q₊`); that is where the four remaining ratios would have to go next. The six unclosed ratios whose partner has `m > 1000` were not run on the partner. Nothing is in Lean.
- **Known technique.** Passing to covers and to points over a quadratic field is standard (covering techniques of Bruin, and of Flynn and Wetherell; named from memory, not checked). No claim that the method is new, and the literature on this system was not searched for it.

## 7. Reproduce

    pip install -r scripts/g1_full/requirements.txt
    cd scripts/second_descent
    python3 second_descent.py $(cat pairs123.txt) > out.jsonl      # about 10 minutes on one core
    diff out.jsonl second_123.jsonl
    python3 second_descent.py --summary second_123.jsonl           # second_summary.txt
    python3 second_descent.py --all 440,279 719,161                # all 24 curves of a ratio
    python3 second_control.py                                      # second_control.txt
    python3 degenerate_tuples.py pairs123.txt                      # degenerate_tuples.txt
    python3 lcheck.py 1e14 second_123.jsonl                        # lcheck.txt, about 40 minutes

CI runs the first five on every commit (with `lcheck.py 3e9`).
