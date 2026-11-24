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
| **For 115 of the 425 pairs there is no magic square of squares with that ratio; 310 pairs are left**, 306 with `m² + n²` a square or twice a square and 4 others (`IsMSS.ratio_ne_of_coverDescent`, `IsMSS.ratio_ne_of_le_1000_descent`) | `CoverDescent`, `CoverDescentTable`; `LAST_PAIR.md` §15 | **not released** (after 0.5.0) |
| **The class of a solution, for every ratio** (coprime `m, n`): its primes divide `2mn(m² − n²)`; its odd primes are `≡ ±1 (mod 8)` (Theorem E); it is even (F); it is a square modulo every odd prime dividing `m² + n²` (G) (`exists_class_struct`; `mod_eight_of_dv`, `two_dvd_class`, `isSquare_class_mod`) | `CoverDescent`; `LAST_PAIR.md` §17 | **not released** (after 0.5.0) |

The 115 are exclusions with a proof in Lean, not new closures: all 115 were already closed in tier B. The exact statements and their hypotheses (`h + f ≠ 0`, `b ≠ d`; entries not assumed distinct) are spelled out in `LAST_PAIR.md` §15. `scripts/CoverDescentChecks.lean` (run by CI, not part of the library) checks the count `310 + 115 = 425`, the list of the four pairs without a certificate that are not of the degenerate types, and that forged certificates are refused.

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
| Theorem D: a solution on the four forms has one squarefree class `d`, supported on the primes of `2mn(m² − n²)`, that passes every local test (Lemmas S1, S2, S3) | `LAST_PAIR.md` §11 | **The existence and support of the class are in Lean for every ratio** (`CoverDescent.exists_class_prod`, stated for all `m, n`, in integer coordinates, with no unique factorisation; tier A). What is not in Lean beyond the 115 ratios is the local computation that says which classes remain, (627,364) included, and the identification of the Lean statement with the prose form of S1 and S2 |
| Theorems E, F, G on the class (hand proofs) | `LAST_PAIR.md` §17 | **In Lean (tier A above)**, each with a different proof route from the hand proof for F (a check of 4,096 cases modulo 8 by the kernel). Tests: the exhaustion on 3,043 ratios and a brute force over all residues (`class_fg.py`), no exception. Not read by a person |
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

- C → checked: a person reads S1/S2 (one page) and §14; or S1/S2 for all ratios go into Lean (the route of §15 already avoids them for 115).
- (627,364) and (991,263): a reading of Theorem A, in particular the application of Checcoli–Veneziano–Viada; a certified replacement for Sage's 0.045162.
- B → A for more ratios: for the 4 non-degenerate pairs without a certificate the descent with one class cannot conclude (§15); a checker for four classes in Lean would be needed. For the 306 others a class always survives, and Lean would need the second step (Theorem A) or ranks.
- Integrity of a copy: `python3 scripts/verify_manifest.py` (see `REPRODUCE.md`).
