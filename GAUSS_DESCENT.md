# Descent over the Gaussian integers: an experiment

23 November 2026, after release 0.4.0. Written with Claude (an AI model made by Anthropic). **An experiment: nobody has read it, only its starting identities are in Lean, and it is not counted in the census.** If it is right, it does three things:

- it closes the two open ratios without a trivial point, (820,799) and (889,320), by congruences alone (since confirmed by Magma's two-cover descent on a curve of genus 2 for each; section 4);
- it gives a second proof, with no elliptic curve, rank or height in it, for the other 32 ratios without a trivial point that the census closes with ranks (among them (724,665), which rests on Magma's rank bound, and (740,539), a four-point closure);
- for the 123 ratios with a trivial point it shows that a solution would have to look, in a precise sense, like the degenerate one.

The census would then read: 26 open becomes 24 open, all with a trivial point. It is not changed until this has been checked.

Programs: `scripts/gauss_descent/`. Everything is elementary: unique factorisation in `ℤ[i]` and congruences.

**What it is, in known terms.** The descent is a two-cover descent in the sense of Bruin and Stoll on the hyperelliptic curves `Z² = (product of some of the five forms)`, carried out for all five forms together and restricted to the classes at which each form is a square. The étale algebra is a product of copies of `ℚ(i)`, which is why everything can be done by hand in `ℤ[i]`. The method is not new; what is done here is its application to this system, and the organisation over `ℤ[i]`. Section 4 compares it with Magma's implementation of the two-cover descent.

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

*Proof, case by case.* Write `N(ℓ) = c·r²` with `r` rational, and let `ℓ = ℓ(u, v)`, `ℓ̄ = ℓ̄(u, v)`. Two facts are used throughout: `N(ℓ) = ℓ·ℓ̄`, and conjugation exchanges the valuations at `π` and `π̄`, so `v_π(ℓ̄) = v_π̄(ℓ)`.

*The resultant step.* If a Gaussian prime divides both `ℓ` and `ℓ̄` at a coprime pair `(u, v)`, it divides `R`. Indeed, with `au + bv = 1`, `R` is a combination of `ℓ` and `ℓ̄` with coefficients in `ℤ[i]`; the combinations are written out, and checked, in Lean (`support_alpha0`, `support_alpha8`, `support_gamma4`, `support_gammaL`, `support_gammaU` in `GaussDescent.lean`, with `R = 2mn`, `2mn`, `2`, `2(m² − n²)`, `2(m² − n²)` up to a unit).

*The prime above 2.* Let `π = 1 + i`. Then `π² = 2i`, so `v_π(2) = 2`, and `π̄ = 1 − i = −i·π`: conjugation fixes this place, so `v_π(ℓ̄) = v_π(ℓ)` and `v_π(N(ℓ)) = 2·v_π(ℓ)`. Also `N(ℓ)` is a rational integer, so `v_π(N(ℓ)) = 2·v_2(N(ℓ))`. Hence

    v_π(ℓ) = v_2(N(ℓ)) = v_2(c) + 2·v_2(r),

and `π` divides `δ` (that is, `v_π(ℓ)` is odd) exactly when `v_2(c)` is odd. No resultant enters. For `α0`, `α8` (`c = 1`): never. For `γL`, `γU` (`c = 2`): always. For `γ4` (`c = 2S`): exactly when `v_2(S)` is even; `v_2(S)` is 0 or 1 because `m`, `n` are coprime.

*An inert prime `q`.* `q` is its own conjugate, so `v_q(ℓ̄) = v_q(ℓ)`. If `v_q(ℓ)` is odd then `q` divides `ℓ` and `ℓ̄`, hence `R` by the resultant step.

*A split prime `p = π·π̄`.* `v_π(N(ℓ)) = v_p(N(ℓ)) = v_p(c) + 2·v_p(r)`, and `v_π(N(ℓ)) = v_π(ℓ) + v_π(ℓ̄) = v_π(ℓ) + v_π̄(ℓ)`. If `v_p(c)` is odd, exactly one of `v_π(ℓ)`, `v_π̄(ℓ)` is odd: exactly one of `π`, `π̄` divides `δ`. If `v_p(c)` is even, both are even or both are odd; if both are odd, `π` divides `ℓ` and (because `v_π̄(ℓ) = v_π(ℓ̄)` is odd) also `ℓ̄`, hence `R`, so `p | R`. ∎

*Units.* The units of `ℤ[i]` are `±1, ±i` and the squares among them are `±1`. So `δ` is determined up to sign, and there are two unit classes, `1` and `i`: the factor `iᵃ`, `a ∈ {0, 1}`, in the sets below. The Lemma says nothing about `a`; the local step decides it. (`2 = −i·π²`, so the class of `2` is the class of `i`.)

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

**Beyond the census** (`scan.py`, `scan_1000_5000.txt`, and one line per ratio in `scan_certificates.txt`: the ratio, the last prime used, the number of tuples left; first implementation, with the second on the 58 ratios that need a prime of 43 or more). For `1000 < m ≤ 5000`, 20,117 coprime pairs avoid the 484 residue classes; 1,217 have a trivial point and were not run; **the other 18,900 all have an empty descent set**, each after primes up to 97 at most. With the 681 above: for `m ≤ 5000`, no ratio that avoids the 484 classes and has `m² + n²` neither a square nor twice a square is left. If the experiment is right, this is a statement about every magic square of squares with `m ≤ 5000`, obtained without any elliptic curve; the census with ranks stops at 1000.

So in every case computed, the descent set is the set of degenerate solutions. If that held for all `m, n`, the ratio `(b + d) : (h + f)` of a magic square of squares would always have `m² + n²` a square or twice a square, the exceptional case of Theorem 3′ of `PAPER.md`. We have no proof, and we do not expect a uniform one from this method: the sets `D_j` depend on the primes of `m`, `n`, `m² − n²`, and what the computation shows is the usual experience that such a descent set is empty unless there is a rational point. It is a statement ratio by ratio.

The ratios with a trivial point are taken one step further in `SECOND_DESCENT.md`, which uses the two tuples left here.

## 4. Checks made

- **The Lemma on real points** (`real_points.py`, `real_points_157.txt`). Each single condition is a conic, so it has many rational points. For the 157 ratios and `|u|, |v| ≤ 300` there are 51,434 coprime pairs at which one of the five forms satisfies its condition. At every one of them exactly one element `d` of the program's set `D_j` makes `ℓ_j(u, v)·d` a square in `ℤ[i]`, tested with integer square roots only. The same holds at 29,656 points of the 681 ratios without a trivial point (`|u|, |v| ≤ 200`). At those points the class labels of `ℓ_j(u, v)` and of `d` were compared at every place above the primes up to 37, 874,378 comparisons, 51,434 of them at `1 + i`: no difference.
- **The local sets on real points** (`local_points.py`, `local_points_157.txt`). `W_p` must contain the class pattern of every primitive point, and every coprime integer pair is one. For the 157 ratios, 1,979,148 integer pairs (all residues modulo a power of `p` with random higher digits, random large pairs, and pairs at which a form is divisible by a high power of a place above `p`) were tested with exact classes against `W_p`, for the primes up to 37 and the bad primes up to 100: no pattern is missing. Together with the previous item this tests, link by link and on actual numbers, that a real solution could not be removed.
- **The class labels against PARI** (`hilbert_check.py`, `hilbert_check.txt`). Two numbers have the same label at the places above `p` exactly when their product is a local square, which holds exactly when its Hilbert symbols with a spanning set are all 1. For 27,000 random pairs of Gaussian integers (3,000 for each of `p = 2, 3, 5, 7, 11, 13, 17, 29, 37`) the labels of the program and PARI's `nfhilbert` agree in both directions, and the symbols see 16 classes at `1 + i` and 4 at every other place, as the theory says. This is the check of the rule "a unit congruent to 1 modulo `π⁵` is a square" that does not use the program's own reasoning.
- **Magma's two-cover descent gives the same numbers** (`subsets.py`, `subsets.txt`, `magma_twocover_raw.txt`, `twocover_check.py`). The descent can be restricted to a subset `T` of the five forms: keep only the tuples `(δ_j)`, `j ∈ T`. A solution gives a rational point of the hyperelliptic curve `Z² = ∏_{j∈T} f_j(t)` at which every `f_j` is a square, and the tuples left correspond to the elements of the fake 2-Selmer set of that curve (Bruin and Stoll's two-cover descent) at which each of the norms is a square. For three forms the curve has genus 2 and Magma's `TwoCoverDescent` computes that set on the public calculator. The two-cover descent works modulo `ℚ*`, so tuples that differ by the same rational factor in every form count once (`subsets.py` prints both numbers). On the curves for which the calculator finished within its limits (9 at the time of writing, for the ratios (4,3), (15,8), (59,56), (820,799), (889,320); `twocover_check.txt`), the two numbers agree every time: 2 where there is a trivial point, out of fake Selmer sets of 2, 8, 2, 2 and 11 elements, and 0 otherwise. In particular Magma finds the fake 2-Selmer set **empty** for `Z² = G4·GL·GU` with the ratio (820,799) and for `Z² = F0·F8·G4` with (889,320): these two curves have no rational point, so the two open ratios without a trivial point have no solution, **by Magma's two-cover descent alone**, with none of the programs of this repository involved beyond the definition of the five forms. All the quadratic factors define `ℚ(i)`, so no class group or unit computation of any size enters. For all 34 ratios without a trivial point some three of the five forms already have an empty set (`subsets_34.txt`: `F0·F8·G4` for 31 of them, `G4·GL·GU` for 29; for 8 ratios two forms suffice). This is an implementation of a published algorithm, written by other people, that shares nothing with the programs here except the five forms.
- **A system with a solution is never excluded** (`gd_control.py`). For random `m, n, u, v`, let `c_j` be the squarefree part of `N(ℓ_j(u, v))`. The system "`N(ℓ_j) ∈ c_j·ℚ²`" has the solution `(u, v)`, and the Lemma gives its sets `D_j`. In 100 such systems (`gd_control.py 1 100`, output in `gd_control.txt`) the tuple of the solution survived every prime, and 1 or 2 tuples were left in all. This tests the Lemma (a wrong support would lose the true tuple), the choice between `π` and `π̄`, and the agreement between the global and the local classes. The rule of the Lemma, applied to the real constants, gives exactly the sets of Section 1.
- **The degenerate solutions are never excluded**: in all 123 ratios with a trivial point their tuples are among those left.
- **A second implementation of the local step** (`gd_check.py`), without class labels and without the refinement by lifting: it decides for each `(u, v)` modulo `pᵏ` and each `δ` whether `ℓ_j(u, v)·δ` is a square, and counts an undecided case as possible. **It gives 0 tuples for all 34 ratios without a trivial point** (`gd_check_34.txt`) and 2 for (4,3) and (15,8), as the first. At the fixed precision `p²` it left 8 tuples for (320,289) and (609,31), where 17² divides `n` or the forms vanish to higher order modulo 17; at precision `17³` it leaves none. It shares with the first program the sets `D_j` and the five forms, so it does not check the Lemma.
- **The second implementation on the hardest ratios of the scan** (`gd_check_scan_high.txt`). Of the 18,900 ratios with `1000 < m ≤ 5000`, 58 need a prime of 43 or more (37 need 53, 9 need 61, 2 need 73, 2 need 97). The second implementation gives 0 tuples for all 58, with the same primes. For one of them, (3968,2633), it first left 4 tuples after 53, where the first has none. The cause was in the second implementation: at a split prime it returned "undecided" as soon as one of the two places was undecided, even if the other place already said "not a square". That errs on the safe side. It now lets a place that says no decide (the output for the 34 ratios is unchanged). For that ratio and prime the set allowed by the first implementation, 4,194,304 tuples, is exactly the set of tuples realised by integer points `(1, v)`, `(u, 1)` to depth `53³`.
- **Lower precision only enlarges the set**: with lifting depth 1, 2 and 5 the sets left are nested, as they must be (and still empty for (820,799), (889,320), (740,539)).

## 5. What could be wrong, and what is not done

- **A bug that a reader found within an hour of the first commit of this file** (`25c1253`): the function that compares a candidate with a local pattern treated an undetermined place as a mismatch instead of a wildcard, which is the unsafe direction, so the first program's empty sets were not proofs. Fixed (the comparison is now per place, with a self-test; `gauss_descent.py`), and everything was run again: no output changed, for the 157 ratios, the 681, or the pairs up to 60. Why it made no difference was not investigated; it was luck, not design. The second implementation never had this defect.

- Both implementations were written by the same AI model in one evening, and so were the tests on real points. An error in the shared ideas would not be caught by the second implementation. What protects against it: the Lemma is tested on about 81,000 real points and in the systems with known solutions, the labels are tested against PARI's Hilbert symbols, and the local sets against two million integer points. None of this is a proof read by a person, and there is still no implementation written by someone else.
- Not in Lean, except the identities it starts from (`GaussDescent.lean`: the five norms, and that a common divisor of a form and its conjugate divides `2mn`, `2` or `2(m² − n²)`). The statement for one ratio is a finite computation of the same kind as the table of 484 classes, so it could be formalised; that has not been done.
- The 2 tuples left for a ratio with a trivial point prove nothing about that ratio by themselves: other rational `t0` could give the same two patterns. `SECOND_DESCENT.md` takes them further, on the curves `d·Y² = Q(s)` that are defined over `ℚ`; the curves that are only defined over `ℚ(i)` would need a 2-descent over `ℚ(i)`, which has not been written.

## 6. Reproduce

```
cd scripts/gauss_descent
python3 gauss_descent.py -v 820,799 889,320        # under a second each
python3 gauss_descent.py --pmax=200 $(cat pairs157.txt)   # all 157, about 10 minutes
python3 gd_control.py 1 100                        # the control with known solutions, gd_control.txt, about 7 minutes
python3 gd_check.py $(cat pairs34.txt)            # the second implementation on the 34 ratios, under 2 minutes
python3 gd_check.py --pmax=100 $(cat pairs_scan_high.txt)   # and on the 58 hardest ratios of the scan
python3 small_all.py 60                            # every coprime pair with m <= 60
python3 scan.py 1000 2000                          # the ratios with 1000 < m <= 2000 that avoid the 484 classes, about 2 minutes
python3 real_points.py 300 $(cat pairs157.txt)    # the Lemma and the labels on real points, real_points_157.txt
python3 local_points.py $(cat pairs157.txt)       # the local sets against integer points, local_points_157.txt, about 3 minutes
python3 hilbert_check.py 3000                      # the labels against PARI's Hilbert symbols (needs cypari2), hilbert_check.txt
python3 subsets.py 4,3 15,8 59,56 820,799 889,320  # the descent on subsets of the forms, subsets.txt
python3 twocover_check.py                          # Magma's two-cover descent against subsets.txt
```

Needs Python with numpy only, except `hilbert_check.py`.
