# Ledger: what each claim of this repository rests on

One page for a reader who wants to know, for any sentence in this repository, which kind of sentence it is. Written by Claude (an AI model made by Anthropic), directed by Drew (`drew2apps`), on 24 November 2026, after a review of commit `9d59c67`. Dates in this repository are in UTC (a commit made on the evening of 23 November in the Americas is dated 24 November here). **No human mathematician has verified anything here**; the first column says what a reader would have to trust.

The problem is open and is not solved here. No novelty is claimed for the methods; they are classical (congruences, descent on norms, ranks of elliptic curves, Demjanenko–Manin, Mordell–Weil sieve). Where a statement may be new it is in tier C or D and says so.

## Tiers

| Tier | What it means | What has to be trusted |
|---|---|---|
| **A. Lean** | A theorem of `MagicSquareOfSquares/*.lean`, checked by the Lean kernel, standard axioms only, no `sorry`, no `native_decide` | The kernel and Mathlib; that the statement says what the prose says (the statements have been read by nobody outside) |
| **B. Census** | Counted in the census of ratios `m ≤ 1000` (268 + 122 + 9 + 2 closed, 24 open) | PARI, Sage with eclib, Magma; derivations written here and checked by nobody; for the 2: transcripts of the public Magma calculator, which CI can check but not regenerate |
| **C. Experiment with a written proof** | A statement with a proof in prose in this repository, not counted in the census | The proof, read by nobody; the programs it uses |
| **D. Experiment, computed** | The output of a program of this repository, reproducible, not counted | The program; the rules it implements (proved in prose in tier C or tested only) |
| **E. Quoted** | A statement of someone else's paper or program used as an input | The source, as read by an AI model; listed with what was and was not read |
| **W. Withdrawn** | Was claimed in an earlier commit and is not claimed any more | — |

## Tier A: Lean (kernel-checked)

| Claim | Where | Released |
|---|---|---|
| The nine necessary conditions of the Wikipedia list; reformulations of Robertson; published examples | `Basic`, `Primes`, `Residues`, `Primitive`, `Progressions`, `Curve`, `Elliptic`, `Triangles`, `Examples`, `Statement` | 0.5.0 |
| Conditions on the two ratios of edge roots; `m² + n²` divides `8e²`; the factorisation of the centre | `G1Norm`, `G1Centre`, `G1Statement` | 0.5.0 |
| The ratio `m : n` avoids 484 residue classes at 34 primes; for `m ≤ 1000` it is one of 425 pairs (`IsMSS.ratio_ne_of_le_1000`) | `G1Full`, `G1FullCertificates`, `G1Census1000`, `G1Bridge` | 0.5.0 |
| One edge root exceeds 500 | `G1Size`, `G1SizeCheck` | 0.5.0 |
| Algebraic identities behind the tiers B and C (no curve, no height, no rank is in Lean) | `G1Relations`, `G1FourPoints`, `GaussDescent`, `SevenSquares` | 0.5.0 |
| **For 119 of the 425 pairs there is no magic square of squares with that ratio; 306 pairs are left, and for all of them `m² + n²` is a square or twice a square.** So: in a magic square of squares with ratio `m : n`, coprime `n < m ≤ 1000`, `m² + n²` is a square or twice a square (`IsMSS.ratio_norm_sq_or_two_sq_of_le_1000`; `IsMSS.ratio_ne_of_coverDescent`, `IsMSS.ratio_ne_of_le_1000_descent`, `exceptions_degenerate`) | `CoverDescent`, `CoverDescentTable`; `LAST_PAIR.md` §15 | **not released** (after 0.5.0) |
| **The class of a solution, for every ratio** (coprime `m, n`; `CLASS.md`): it exists and its primes divide `2mn(m² − n²)` (`exists_class_prod`, no unique factorisation); its odd primes are `≡ ±1 (mod 8)` (E); it is even (F); it is a square modulo every odd prime dividing `m² + n²` (G) and modulo every odd prime dividing `m⁸ + 14m⁴n⁴ + n⁸` (O) (`exists_class_struct`) | `CoverDescent`; `LAST_PAIR.md` §17 | **not released** (after 0.5.0) |
| **The class of a solution is `∏ (w_j + x_j)`, up to squares and sign** (`class_iff`: Lemma 1 of `COMPLETE.md` over the integers, both directions); the same identity modulo any `M` (`gsqM_class`); **at an odd prime `p` dividing `m`, the classes of sections 4.1 and 4.3 (case `p² | m`) of `COMPLETE.md` pass modulo every power of `p`** (`pass_first_order`, `pass_odd_valuation`; with no hypothesis on a point: `unit_passes_of_three_mod_four`, `prime_mul_unit_passes_mod_eight`). This is part of Theorem C, not all of it: see `COMPLETE.md`, section 5 | `ClassLocal` | **not released** |
| **For 185 of the 306 ratios to 1000 with a degenerate solution, the class of every solution on the four forms is the one listed** (`structTable`, `class_of_mem_structTable`); **for 627 : 364 and 991 : 263 it is 521266**, from those theorems and the primes 5, 29, 37; no local computation. This is the input of Theorem A (tier C), up to the identification of the eight products with `Q₁`, `Q₂`, which is on paper | `CoverDescent`, `CoverDescentTable`; `LAST_PAIR.md` §17 | **not released** |

The 119 are exclusions with a proof in Lean, not new closures: all 119 were already closed in tier B (by ranks, Demjanenko–Manin, four points or Magma); what changes is what they rest on. The exact statements and their hypotheses (`h + f ≠ 0`, `b ≠ d`; entries not assumed distinct) are spelled out in `LAST_PAIR.md` §15. `scripts/CoverDescentChecks.lean` (run by CI, not part of the library) checks the count `306 + 119 = 425`, that the 119 are exactly the pairs with `m² + n²` neither a square nor twice a square, and that forged certificates of each kind are refused.

CI on every push: `lake build`, the axiom audit of every theorem (`scripts/Audit.lean`, compared with `scripts/audit_output.txt`), the generated Lean files against their generators, the count of declared theorems.

## Tier B: the census of ratios `m ≤ 1000` (release 0.5.0)

| Count | Method | Rests on |
|---|---|---|
| 268 | An elliptic curve of rank 0 | PARI and Sage rank computations (`scripts/g1_full/`) |
| 122 | Demjanenko–Manin on three curves of rank 1 | PARI, Sage, Magma (it confirms every upper rank bound; three pairs rest on it alone); derivations in `G1_FULL_SYSTEM.md`, unchecked |
| 9, counted apart | The same with four points on a curve of rank 2 or 3 | `FOUR_POINTS.md`: the explicit form was derived here; unchecked |
| 2, counted apart: (820,799), (889,320) | Magma's two-cover descent on a curve of genus 2 | Transcripts (`scripts/gauss_descent/magma_twocover_raw.txt`); since `9d59c67` both are also in tier A |
| **24 open** | — | All have a degenerate solution; (627,364) and (991,263) are among them |

The census changes only by a decision of the repository owner. Nothing in tiers C and D is counted in it.

## Tier C: experiments with a written proof (not counted, not released unless noted)

| Claim | Where | Independent checks that exist |
|---|---|---|
| Theorem D: a solution on the four forms has one squarefree class `d`, supported on the primes of `2mn(m² − n²)`, that passes every local test (Lemmas S1, S2, S3) | `LAST_PAIR.md` §11 | **The existence and support of the class are in Lean for every ratio** (`CoverDescent.exists_class_prod`, stated for all `m, n`, in integer coordinates, with no unique factorisation; tier A). What is not in Lean beyond the 119 ratios is the local computation that says which classes remain, (627,364) included, and the identification of the Lean statement with the prose form of S1 and S2 |
| Theorems E, F, G, O on the class (hand proofs) | `LAST_PAIR.md` §17; `CLASS.md` | **In Lean (tier A above).** Tests: the exhaustion on 3,043 ratios, brute force over all residues (`class_fg.py`; for O, modulo 13, 37, 61, 73), and the 306 ratios with a degenerate solution, where a candidate must pass and does (`class_struct.py`). Not read by a person. The observation behind O (only residues at the roots of the octic) is a reviewer's |
| The law at 5: if 5 divides one of `m, n, m ± n` exactly once, the class is a non-residue modulo 5 | `LAST_PAIR.md` §17 | **In Lean** (`class_mod_five`: a complete check modulo 25 by the kernel on pairs in normal form; in `exists_class_struct` and in `structOk`). Also `class_law5.py` (all residues modulo 25, in CI) and a proof by hand by another instance, not in the repository. Not read by a person |
| For 67 of the 119 ratios to 1000 without a degenerate solution, no solution on the four forms, from the theorems on the class alone; the class for 185 of the 306 others | `CLASS.md`, Theorems 6 and 6b | **In Lean** (`structZero`, `no_point_of_mem_structZero`, `structTable`, `class_of_mem_structTable`; the kernel checks primality, divisibility and Euler's criterion, no points). Experiments: not counted in the census |
| Theorem C: the set of classes of `p`-adic points at every bad prime, for every ratio; hence the five laws are all that the bad primes say about the class | `COMPLETE.md` | **A proof by hand written by an AI model; not read by a person. In Lean: the upper bounds for integer points, F′ (`class_mod_sixteen`), Lemma 1 over the integers and two of the four families of points of the lower bounds, modulo every power of `p` (`ClassLocal.lean`). Not in Lean: the conic and Hasse's bound, the case `p ∥ m` of the odd classes, the prime 5, the symmetries, the passage to `p`-adic points and the upper bounds there.** Checks: `class_local.py` (exact `p`-adic classes at 211,894 random points, none outside the table, every class seen; the points of the proof; in CI) |
| The same statement as tested before the proof | `CLASS.md` §3 | **Tested** (`class_complete.py`, `class_complete_random.py`, in CI with small bounds): 425 ratios with 92,752 candidates and 1,952 random ratios with `m ≤ 20000` and 965,504 candidates, bad primes below 400, points modulo small prime powers (weaker than `p`-adic); no candidate allowed by the laws is excluded there |
| The local rules R1 to R5 (when a `p`-adic class is determined) | `LAST_PAIR.md` §14 | Mutation tests (`tuple_mutation.py`); two exhaustion programs with no shared code; a third implementation of the one-class test written from the prose by a separate reviewer instance (`scripts/cover_descent/third/`), with witnesses: no class with a witness rejected at 3,024 pairs, and the proving inclusion at 3,006; a fourth, of the four-class test, by another instance, on 2,475 pairs (§16; `scripts/cover_descent/fourth/`; one pair where the exhaustion is the weaker side, named there). Rule R4's hint: a wrong hint is noticed by the comparison with the exhaustion on the cuboid system (8 of 24 pairs; all 24 if the hint is dropped), and by no test with witnesses (`mut_hint.py`) |
| Theorem A: for (627,364), `d = 521266`, the only solutions are the degenerate ones (seven lemmas; Mordell–Weil sieve and an explicit height bound) | `LAST_PAIR.md` §1–8 | Tests on planted points; rests on tier E (CVV constants, CPS bounds, Sage's height bound 0.045162) |
| Theorem C: the degenerate section has infinite order for every `τ ≠ 0, ±1` | `LAST_PAIR.md` §9 | Ten factored polynomials, checked by program |
| The four-point form of Demjanenko–Manin | `FOUR_POINTS.md` | Lean has the divisibility statements and two lattice inequalities only |
| Proposition 1 on the rational distance problem (the four classes are never equal) | `RATIONAL_DISTANCE.md` §1 | Every computation of the proof in exact arithmetic (`rd_prop1.py`). The statement is `n ≠ 1` in the notation of Si's preprint (tier E), which uses `n > 1`; we did not find a proof of it in the two sections we read |

## Tier D: experiments, computed (not counted, not released unless noted)

| Claim | Where | CI |
|---|---|---|
| Gaussian descent and second descent on the open ratios | `GAUSS_DESCENT.md`, `SECOND_DESCENT.md` | Samples rerun; `SECOND_DESCENT.md` §8 lists outputs that no job reruns |
| The class descent on all 157 ratios without a rank-0 closure, and the four-class descent | `LAST_PAIR.md` §11–12, `scripts/cover_descent/` | Rerun in full |
| Every ratio with `m ≤ 5000` avoiding the 484 classes: 19,581 with nothing left, 1,523 with exactly the degenerate tuples | `tuple_scan_5000.txt`, `tuple_trivial_5000.txt` | Rerun in full by the job `scan-5000` (long jobs run on `[full ci]`, tags and pull requests) |
| The proportion of ratios that the 484 classes leave: 0.2766% as a product over 34 primes, 0.2777% counted to 5000 | `LAST_PAIR.md` §12, `G1_FULL_SYSTEM.md` | Counts checked; the density statement itself is argued, not proved |
| Rational distance: 13,697 of 13,698 lines `y = r/s`, `s ≤ 300`, closed | `RATIONAL_DISTANCE.md` | Rerun in full by `scan-5000` |
| Perfect triangle (Heron triangle with three rational medians): all 10,604 lines to denominator 200 closed, with one second step for all lines | `PERFECT_TRIANGLE.md` | Small range rerun on every push; the full scan in `scan-5000` |
| The first descent removes non-exceptional ratios by accidents at small split primes, not by a law; `d = 2` passes every structural condition | `LAST_PAIR.md` §17 | `class_structure.py`, `class_two.py` rerun on a small range |
| Perfect cuboid: a two-step descent closes 274 of 411 open bases | repository `perfect-cuboid-census` | Its own CI |
| Searches (seven squares to centre root 10⁹; pairs of ratios to height 6000) | `README.md`, `scripts/search/` | Small ranges rerun; the source hash of the long run is checked |

## Tier E: quoted from others

| Statement | Source | Read how |
|---|---|---|
| Explicit Mordell bound (Theorem 1.2, Proposition 3.2, Theorem 2.2) | Checcoli, Veneziano, Viada, Forum Math. Sigma 7 (2019) | By an AI model, twice; normalisations not confirmed by a person |
| Bounds between naive and canonical height (constants 1.946, 2.14, Table 1) | Cremona, Prickett, Siksek (2006) | By an AI model |
| Lower bound 0.045162 for the canonical height | Sage, `height_function().min` | Its documentation; not recomputed independently |
| Two-cover descent | Magma (Bruin and Stoll's algorithm), public calculator | Transcripts committed; not regenerable |
| Corner classes and the rank requirement for the rational distance problem | Y. Si, doi:10.5281/zenodo.24212887, v1.1.0 | Sections 1 and 12 of the LaTeX sources, by an AI model; his corollary is not formalised |
| The list of 1,072 proven cuboid fibers | R. Peschmann, arXiv:2604.28072 | Lists compared, proofs not read |
| The surface of magic squares of squares | Auel and Singer, arXiv:2609.09351 | Statements and relevant sections, by an AI model |
| Hill's claimed proof | arXiv:2510.08286 | The last step only (`HILL_NOTE.md`); not verified, not refuted as a whole |

## Tier W: withdrawn (kept in `CHANGES.md`)

- A second step for the perfect cuboid with three factors (the class depended on the denominator); replaced by the step with four factors.
- The identification of a Jacobian by Sage in Proposition 1 of `RATIONAL_DISTANCE.md`; replaced by an explicit map and exact arithmetic.
- "The two local methods give the same sets" without qualification; they give the same final sets, not the same set at every prime.
- The first version of the seven-square search (two of six cases missing); rerun in full.

## What would raise a tier

- C → checked: a person reads S1/S2 (one page) and §14; or S1/S2 for all ratios go into Lean (the route of §15 already avoids them for 119).
- (627,364) and (991,263): a reading of Theorem A, in particular the application of Checcoli–Veneziano–Viada; a certified replacement for Sage's 0.045162.
- B → A for more ratios: nothing is left to do with the first descent for `m ≤ 1000`. For 150 of the 306 pairs with a degenerate solution the four theorems on the class fix the class (computed, `class_struct.py`; in Lean for the last pair only). For the 306 pairs of the degenerate types a class always survives, and Lean would need the second step (Theorem A) or ranks, neither of which is in Lean or in Mathlib.
- Integrity of a copy: `python3 scripts/verify_manifest.py` (see `REPRODUCE.md`).
