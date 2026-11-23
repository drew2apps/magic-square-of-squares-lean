# Descent over the Gaussian integers: an experiment

23 November 2026, after release 0.4.0. Written with Claude (an AI model made by Anthropic). **An experiment: nobody has read it, it is not in Lean, and it is not counted in the census.** If it is right, it does three things:

- it closes the two open ratios without a trivial point, (820,799) and (889,320), by congruences alone;
- it gives a second proof, with no elliptic curve, rank or height in it, for the other 32 ratios without a trivial point that the census closes with ranks (among them (724,665), which rests on Magma's rank bound, and (740,539), a four-point closure);
- for the 123 ratios with a trivial point it shows that a solution would have to look, in a precise sense, like the degenerate one.

The census would then read: 26 open becomes 24 open, all with a trivial point. It is not changed until this has been checked.

Programs: `scripts/gauss_descent/`. Everything is elementary: unique factorisation in `ℤ[i]` and congruences.

## 1. The system over `ℤ[i]`

A solution for the ratio `m : n` is a coprime integer pair `(u, v)`, `t0 = u/v`, at which the five forms are squares (`G1_FULL_SYSTEM.md`, "The system"). Each form is the norm of a linear form over `ℤ[i]`. With `S = m² + n²`:

| form | linear form over `ℤ[i]` | condition |
|---|---|---|
| `v²F0` | `α0 = mv + inu` | `N(α0)` is a square |
| `v²F8` | `α8 = nv + imu` | `N(α8)` is a square |
| `v²G4/(2S)` | `γ4 = u + iv` | `2S·N(γ4)` is a square |
| `v²GL/2` | `γL = (m + in)u − i(m − in)v` | `2·N(γL)` is a square |
| `v²GU/2` | `γU = (m + in)u + i(m − in)v` | `2·N(γU)` is a square |

`ℤ[i]` has unique factorisation, so each value is `δ_j·β_j²` with `δ_j` squarefree, determined up to squares of units, that is, up to sign.

**Lemma (where `δ_j` can live).** Let `ℓ` be one of the five forms, `c` the constant of its condition (`N(ℓ) ∈ c·ℚ²`), and `R` the resultant of `ℓ` and its conjugate `ℓ̄`. For coprime integers `u, v`:

- an inert prime `q` (`q ≡ 3 mod 4`) divides `δ` only if `q | R`;
- for a split prime `p = ππ̄`: if `v_p(c)` is odd, exactly one of `π`, `π̄` divides `δ`; if `v_p(c)` is even, either none or both, and both only if `p | R`;
- `1 + i` divides `δ` exactly when `v_2(c)` is odd.

*Proof.* `v_π(N(ℓ)) = v_π(ℓ) + v_π(ℓ̄)` and `v_π(ℓ̄) = v_π̄(ℓ)`. If `v_π(ℓ)` is odd and `v_p(c)` is even, then `v_π(ℓ̄)` is odd, so `π` divides `ℓ` and `ℓ̄`, hence `R` (the values of two linear forms at a coprime pair have their common divisors among the divisors of the resultant), and `π̄` divides `ℓ` to an odd power as well. The other cases are the same count. ∎

The resultants, as determinants of the coefficients of `(u, v)` in `ℓ` and `ℓ̄`: for `α0 = (in)u + mv` and `ᾱ0 = (−in)u + mv` it is `2imn`; for `α8` the same with `m`, `n` exchanged; for `γ4 = u + iv` and `u − iv` it is `−2i`; for `γL` and `γ̄L = (m − in)u + i(m + in)v` it is `i(m + in)² + i(m − in)² = 2i(m² − n²)`, and for `γU` its negative. So

- `δ_0, δ_8 ∈ { iᵃ·d : d | odd part of mn, squarefree }`,
- `δ_4 ∈ { iᵃ·(1 + i)ᵇ·∏ (π_p or π̄_p) }`, the product over the primes `p | S` with odd exponent, `b ≡ 1 + v_2(S) mod 2`,
- `δ_L, δ_U ∈ { iᵃ·(1 + i)·d : d | odd part of m² − n², squarefree }`.

Call these finite sets `D_0, …, D_U`. For (820,799) they have 32, 32, 4, 16, 16 elements: about a million tuples.

## 2. The local step

Let `p` be a rational prime and `v` a place of `ℚ(i)` above it. A solution gives a primitive pair `(u, v) ∈ ℤ_p²`, and for each `j` the class of `ℓ_j(u, v)` in `K_v*/K_v*²` equals the class of `δ_j`. So:

> a tuple `(δ_0, δ_8, δ_4, δ_L, δ_U)` can only come from a solution if, for every prime `p`, some primitive `(u, v) ∈ ℤ_p²` gives its five classes at all places above `p`.

For each prime the program computes a set `W_p` that contains all class tuples that occur (when a value vanishes to an order that the working precision does not resolve, every class is allowed for it, so `W_p` can only be too large), and keeps the tuples of `D_0 × … × D_U` that match some element of `W_p`. If no tuple is left after finitely many primes, the ratio has no solution. Classes: at an odd place, the parity of the valuation and the quadratic character of the unit part (for an inert prime the residue field is `F_{p²}`, and a unit is a square there exactly when its norm is a square modulo `p`). At `π = 1 + i`, the parity of the valuation and the class of the unit part modulo squares, read off modulo `π⁵`: *a unit `x` of `ℤ_2[i]` with `x ≡ 1 mod π⁵` is a square.* Proof: Hensel's lemma for `f(y) = y² − x` at `y = 1` needs `v_π(f(1)) > 2·v_π(f′(1)) = 2·v_π(2) = 4`. So two units are in the same class exactly when their quotient is congruent to a square modulo `π⁵`; the units modulo `π⁵` form a group of order 16, the squares among them a subgroup of order 2, and there are 8 classes, as there must be (`[ℚ_2(i) : ℚ_2] = 2` gives `2³` unit classes). The program works modulo `8 = π⁶` and takes as the label of `x` the smallest element of `x·{z², z² + 4 + 4i}`, where `4 + 4i` generates `π⁵` modulo 8.

*Why finite precision is safe.* The program knows `(u, v)` modulo `pᵏ`, hence each value `ℓ_j(u, v)` modulo `pᵏ`. A class is taken as known only if the valuation is below `k` at an odd place (then the leading term is the same for every lift), or `v_π + 5 ≤ 2k` at `π`. Otherwise the place is marked as undetermined, and an undetermined place matches every candidate. A point is refined (all `p` lifts are examined) until every place is determined or a depth limit is reached; at the limit the undetermined marks stay. So every class tuple of an actual point of `ℤ_p²` is matched by a tuple of `W_p`.

This is finer than testing the five forms modulo `p`, which only sees the norms. The ratios of the census all pass that test at every prime, by definition.

## 3. Result

`gauss_descent_157.txt` has one line per ratio of the list of 157.

- **All 34 ratios without a trivial point: no tuple is left**, after the primes up to 5 (7 ratios), 7 (1), 13 (10), 17 (9), 19 (2), 29 (2), 37 (3).
- **All 123 ratios with a trivial point: exactly 2 tuples are left** (primes up to 200 and the bad primes), and they are the tuples of the two degenerate solutions (`t0 = ±1` when `m² + n²` is a square, `t0 = 0, ∞` when it is twice a square).

For the two open ratios, the number of tuples left after each prime:

| ratio | start | 2 | 3 | 5 | 7 | 11 | 13 | 17 |
|---|---|---|---|---|---|---|---|---|
| (820,799) | 1,048,576 | 4096 | 1024 | 128 | 128 | 128 | 0 | |
| (889,320) | 1,048,576 | 4096 | 1024 | 512 | 384 | 384 | 128 | 0 |

**More ratios** (`gauss_descent_notrivial_1000.txt`, `small_all_60.txt`). Of the coprime pairs with `m ≤ 1000`, 987 avoid the 484 residue classes of `modp_classes.txt` (a superset of the 425 of the census, which also uses the prime 2 and the classes with `p | n`). 681 of them have `m² + n²` neither a square nor twice a square, and **for all 681 no tuple is left**, each time after primes up to 37 at most (353 at 2, 163 at 5, 3 at 7, 79 at 13, 43 at 17, 3 at 19, 31 at 29, 6 at 37). For every coprime pair with `m ≤ 60` (1101 pairs) the set is empty unless `m² + n²` is a square or twice a square, and then it has exactly 2 elements.

**Beyond the census** (`scan.py`, `scan_1000_5000.txt`; first implementation only). For `1000 < m ≤ 5000`, 20,117 coprime pairs avoid the 484 residue classes; 1,217 have a trivial point and were not run; **the other 18,900 all have an empty descent set**, each after primes up to 97 at most. With the 681 above: for `m ≤ 5000`, no ratio that avoids the 484 classes and has `m² + n²` neither a square nor twice a square is left. If the experiment is right, this is a statement about every magic square of squares with `m ≤ 5000`, obtained without any elliptic curve; the census with ranks stops at 1000.

So in every case computed, the descent set is the set of degenerate solutions. If that held for all `m, n`, the ratio `(b + d) : (h + f)` of a magic square of squares would always have `m² + n²` a square or twice a square, the exceptional case of Theorem 3′ of `PAPER.md`. We have no proof, and we do not expect a uniform one from this method: the sets `D_j` depend on the primes of `m`, `n`, `m² − n²`, and what the computation shows is the usual experience that such a descent set is empty unless there is a rational point. It is a statement ratio by ratio.

## 4. Checks made

- **A system with a solution is never excluded** (`gd_control.py`). For random `m, n, u, v`, let `c_j` be the squarefree part of `N(ℓ_j(u, v))`. The system "`N(ℓ_j) ∈ c_j·ℚ²`" has the solution `(u, v)`, and the Lemma gives its sets `D_j`. In 100 such systems (two runs) the tuple of the solution survived every prime, and 1 or 2 tuples were left in all. This tests the Lemma (a wrong support would lose the true tuple), the choice between `π` and `π̄`, and the agreement between the global and the local classes. The rule of the Lemma, applied to the real constants, gives exactly the sets of Section 1.
- **The degenerate solutions are never excluded**: in all 123 ratios with a trivial point their tuples are among those left.
- **A second implementation of the local step** (`gd_check.py`), without class labels and without the refinement by lifting: it decides for each `(u, v)` modulo `pᵏ` and each `δ` whether `ℓ_j(u, v)·δ` is a square, and counts an undecided case as possible. **It gives 0 tuples for all 34 ratios without a trivial point** (`gd_check_34.txt`) and 2 for (4,3) and (15,8), as the first. At the fixed precision `p²` it left 8 tuples for (320,289) and (609,31), where 17² divides `n` or the forms vanish to higher order modulo 17; at precision `17³` it leaves none. It shares with the first program the sets `D_j` and the five forms, so it does not check the Lemma.
- **Lower precision only enlarges the set**: with lifting depth 1, 2 and 5 the sets left are nested, as they must be (and still empty for (820,799), (889,320), (740,539)).

## 5. What could be wrong, and what is not done

- **A bug that a reader found within an hour of the first commit of this file** (`25c1253`): the function that compares a candidate with a local pattern treated an undetermined place as a mismatch instead of a wildcard, which is the unsafe direction, so the first program's empty sets were not proofs. Fixed (the comparison is now per place, with a self-test; `gauss_descent.py`), and everything was run again: no output changed, for the 157 ratios, the 681, or the pairs up to 60. Why it made no difference was not investigated; it was luck, not design. The second implementation never had this defect.

- Both implementations were written by the same AI model in one evening. An error in the shared ideas (the Lemma, the handling of the place above 2, the claim that a solution gives a primitive pair over `ℤ_p` with these classes) would not be caught by the second implementation. The control with known solutions is the main protection against that.
- Not in Lean. The statement for one ratio is a finite computation of the same kind as the table of 484 classes, so it could be formalised.
- The 2 tuples left for a ratio with a trivial point prove nothing about that ratio: other rational `t0` could give the same two patterns. They are not explained away: the next step there is a second descent, on the elliptic curves `δ·y² = ℓ_aℓ_bℓ_cℓ_d` over `ℚ(i)` that the twisted cover maps to. A curve of rank 0 over `ℚ(i)` among them would close the ratio. Not started: it needs a 2-descent over `ℚ(i)`, which the software at hand does not provide.

## 6. Reproduce

```
cd scripts/gauss_descent
python3 gauss_descent.py -v 820,799 889,320        # under a second each
python3 gauss_descent.py --pmax=200 $(cat pairs157.txt)   # all 157, about 10 minutes
python3 gd_control.py 1 60                         # the control with known solutions
python3 gd_check.py $(cat pairs34.txt)            # the second implementation on the 34 ratios, under 2 minutes
python3 small_all.py 60                            # every coprime pair with m <= 60
python3 scan.py 1000 2000                          # the ratios with 1000 < m <= 2000 that avoid the 484 classes, about 2 minutes
```

Needs Python with numpy only.
