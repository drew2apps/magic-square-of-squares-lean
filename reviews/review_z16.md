# Reviewer z16 report: EXPLORATION_5 §2 (centre form F4) and §3 (ghost data)

Written by an AI model (reviewer z16), not read by a person.

> Committed as written apart from this note. In this repository the reviewer's folder is `experiments/referee_z16/; the note reviewed is `experiments/research-pairing/explore5/` (sections 2 and 3)`. Paths under `/home/claude/` are the reviewer's scratch folders.

Worktree /home/claude/conductor/z16, branch run11_z16 (from run10). Scripts under experiments/referee_z16/.

## Progress log (appended as I go)

- Started. Worktree created. Read EXPLORATION_5 §2-3 and ghost/GHOST_NOTE.md. The session's code imports
  `pipeline.py` (copy in explore4/), which loads the repo's `exact_local` (as X, TD = tuple_descent) and
  explore4's `cassels.py`/`sel.py` (as /tmp/mix).

## Task 1: which stages impose the centre F4 (read from the repository at run10, 315ad1e)

Finding: the session's claim is **true of the stage-2 tuple test and everything built on it, but it is not a hidden
gap: the repository says so in four places**, and two other stages of the repository do impose F4.

Four forms only (F0, F8, GL, GU; no F4):
- `scripts/cover_descent/tuple_descent.py` `magic()` (docstring: "the four forms F0, F8, GL, GU") and
  `cover_descent.py` (one class of the product of the four); LAST_PAIR.md §11 says in its Setting:
  "(a magic square of squares with ratio m : n gives such a pair, with t = u/v; `G4` is not used)".
- Lean `CoverDescent.lean`: header "the four forms of the G1 family that are not `G4`"; `formsOk`/`formsOkE` test
  the four norms only; `exists_class_prod`, `exists_classes4`, `kills`, `killsE`, `recs4`, `entryOk4`.
- COMPLETE.md (Theorem C, C', the class criterion; `X_p` = "the four norms are nonzero squares"), `class_good.py`,
  `class_selmer.py`; LEDGER tier D row "**Without the fifth form**: for every coprime ratio with m <= 15000, F0, F8,
  GL, GU are squares together only when ..." (so the repository already labels this system as four-form).
- `experiments/roadmap-g/g1/exact_local.py`, `g2/g2_killsets.py` (the 937 = 101 of the 997 + 836 of the 8,739 with
  "bad-alone: k left"; 6,904 tuples, matches), `roadmap-h/h1` (`h1_scan.py`), `h2` (`h2_exhaust.py`, `h2_hunt.py`),
  `h3` (`h3_tables.py`), `h4` (`h4_split.py`), `conductor-t10/t10_pairct.py`, `research-pairing/` (six pair curves
  y^2 = F_j F_k, j, k in {0, 8, L, U}), and EXPLORATION_4's 937/ghost pipeline: all call `TD.magic`.

Five forms (F4 imposed):
- Stage 1: the 484 residue classes / the G1Full five-form test modulo p (`classify.c`, Lean `obstructsFull`),
  and `local_all.py` (p-adic solubility of the five forms, norms only, no classes) on all 304,191 pairs to 1000.
- **`scripts/gauss_descent/gauss_descent.py` (GAUSS_DESCENT.md, tier D): a descent on all five forms with Gaussian
  classes, `g4 = u + iv` with condition `2S·N(g4)` a square**, i.e. exactly the session's fifth form up to the constant
  (1+i)(m+in); Lean `GaussDescent.lean` has `support_gamma4`. It was run on the 157 census ratios and the scan
  1000 < m <= 5000 (`scan_1000_5000.txt`), with good primes. It was never run on the 937 (m <= 15000, bad primes only).
- The rank stages that use quotients containing F4 (E04 = F0·F4 and F8·F4, D1, g4·K; G1_FULL_SYSTEM.md Result 4,
  LEDGER tier B/D "F4's quotient"), the second descent, and the pair-of-ratios search (F1 = centre condition).

So the correct statement is: **"the four-class (rational-class) descent of tuple_descent.magic and all stage-2 work
built on it (G1, G2, H1-H4, T10, the pairing notes, EXPLORATION_4) does not impose F4; the five-form Gaussian descent
of GAUSS_DESCENT.md does, but was run only to m <= 5000 with good primes."** "No existing stage imposes F4"
(EXPLORATION_5 §2, "whether the census already imposes F4 elsewhere" left unchecked) is [fix]: too strong as worded.
Since with good primes the four-form test already empties all 9,736 (LEDGER tier D, G2 row), adding F4 changes no
count; it matters for the bad-primes-only question (what the pairings must explain), which is where §2 uses it.

## Task 2 (in progress): a finding while building p-adic witnesses

- At (3320,2303), a "ghost" ratio of EXPLORATION_4, `m² + n² = 2621 · 6229`, both primes `≡ 5 (mod 8)` to exponent 1.
  The five-form allowed set of `exact_local.allowed_at` is **empty** at 2621 and at 6229 (four-form set: 4 labels each).
  That is no descent at all: the five forms have no `ℚ_p` point there. It is the known condition **"every prime
  `≡ 5 (mod 8)` divides `m² + n²` to an even power" (Lean, tier A: `IsMSS.edge_norm`, `G1Primes.lean`;
  G1_FULL_SYSTEM.md Result 1b)**, which the stage-2 lists (the 997, the 8,739, hence the 937) were never filtered by,
  because they are lists of the four-form system. Hand check: at `p | m² + n²` with odd exponent, `F4 = 2C(u² + v²)`
  square forces `u ≡ ±ιv (mod p)`, and then `F0 ≡ 2m²v²`, a square only if 2 is, i.e. `p ≡ 1 (mod 8)`.
  So part of "F4 removes most of the 937" is a tier-A condition re-found; quantified below.

## Task 3: the counts, recomputed (`z16_five.py`, my own 40-line driver on the repository's `exact_local` and
`tuple_descent`, not the session's pipeline; `z16_summary.py`)

- All 937 ratios (taken from the repository's `g2_*_all.txt`, "bad-alone: k left"): four-form tuples 6,904, equal to
  EXPLORATION_4's tables at every ratio. **Five forms: 211 of 937 ratios keep a 5-tuple; 1,370 of 6,904 4-tuples
  extend; 2,108 5-tuples; pairing survivors that extend: only the 2 at (14903,9696)** (`(1,2,10414,49198)`,
  `(1,2,24599,5207)`, each with two d4). Identical to the session's `f4_all.jsonl` ratio by ratio. Confirmed.
- No pairing flag in EXPLORATION_4's tables is "outside Sel²" (f4_all.py counts that value as "not killed"; it never occurs).
- **318 of the 937 have a prime ≡ 5 (mod 8) to an odd power in m² + n²**, so are already excluded for a magic square by
  the Lean theorem of `G1Primes.lean` (tier A); the five-form test empties all 318 (a local insolubility, not a class
  effect). Among them are two of the five pairing-survivor ratios, (3320,2303) and its partner (5623,1017). **The new
  effect of F4 is 619 → 211 ratios and 5,194 → 1,370 tuples**, not 937 → 211.
- Support check: no 5-tuple left has a prime of even exponent of m² + n² in d4 (as the hand derivation predicts).

## The main finding: the 937 and the ghosts are four-form ratios; for a magic square almost all are already dead

`COMPLETE.md` §7 test 3 builds the 997 (and test 7 the 8,739) from *every* coprime ratio ("no list is used"; "the 484
classes of that list are about another curve"). So the 937 were never filtered by the five-form stage 1. Checking
them against two tier-A statements (`z16_ghost.py`; my membership test reproduces 987 = 306 + 681 ratios to 1000
outside the 484 classes, 0 of the 306 degenerate ratios inside, none of the 9 open ratios inside):

- **844 of the 937 lie in one of the 484 residue classes** (Lean `G1Full`/`obstructsFull`, LEDGER tier A row "The
  ratio m : n avoids 484 residue classes at 34 primes"); 318 have a prime ≡ 5 (mod 8) to an odd power in m² + n²
  (`G1Primes.lean`); 866 are excluded by one or the other.
- **All nine ghost ratios are in the 484 classes**: (14903,9696) and its partner (24599,5207) at 89, 103 and 151;
  (3320,2303), (5623,1017), (9113,3200), (12313,5913), (90384,58537), (126127,16889), (115513,18217) as well (the last
  three and the first pair also fail the 5 mod 8 condition). So for a magic square of squares none of the 18 ghost
  tuples is a candidate at all, and EXPLORATION_5's "next on this line: pair curves with the centre for (14903,9696)
  and its partner, the only ghost ratios left" has no object for the five-form problem: both ratios are excluded by
  congruences mod 89 (tier A). For the four-form problem (eight entries off the centre) F4 is not part of the system,
  so it cannot be added there. Either way the "ghost" question is a question about the four-form system only.
- The samples beyond m = 15000 are the same kind: 268 of 269 (biased), 14 of 14 (random), 38 of 38 (few primes) are
  in the 484 classes; the sampler (`sampler.py`) applies only a mod-8 pattern, not the 484 classes.

## Task 2: soundness of the five-form test

**Hand derivation.** Let a magic square of squares have ratio m : n (coprime, m > n > 0). Lean
`exists_coprime_squares_full` gives coprime integers u, V with F0 = m²V² + n²u², F8, GL, GU and
F4 = 2(m² + n²)(u² + V²) all squares; these are the norms of tuple_descent's γ0 = mv + inu, γ8 = nv − imu,
γL = (1+i)((m − in)u + i(m + in)v), γU, and γ4 = (1+i)(m+in)(u+iv) (I checked the five norm identities by hand
and in `z16_five.py`). (1) γ4 ≠ 0 has square norm, so by Lemma S1 γ4 = d4·z², d4 squarefree positive, z ∈ ℚ(i).
(2) If an odd p divides d4, Lemma S2 gives π | p dividing γ4 and its conjugate. With γ4 = A·u + B·v, A = c,
B = ci, c = (1+i)(m+in): B̄γ4 − B·conj(γ4) = R·u and Āγ4 − A·conj(γ4) = −R·v with R = AB̄ − ĀB = −2i|c|² =
−4i(m² + n²); with xu + yV = 1, π | R, so p | m² + n². Finer (not needed): an even-exponent prime of m² + n² never
divides d4 (π^e with e even from m + in, and u + iv is divisible by at most one of π, π̄), so the support is 2 and
the odd-exponent primes; the session's support (all primes of 2(m² + n²)) is a harmless superset, and in the 2,108
5-tuples left no d4 has an even-exponent prime. (3) At every p, (u, V) is a primitive point of ℤ_p² at which
d_j·γ_j = (d_j z_j)² is a square in ℚ_p ⊗ ℚ(i) for all five j, so the label tuple of (d0, …, d4) is the label tuple
of a local point; an allowed set that contains the labels of every primitive local point with five rational classes
contains it. Nothing in this uses J = 4. So the five-form test is sound exactly when the four-form one is, provided
the local code handles a form whose two coefficients share a factor (new here: at p | m² + n², p ≡ 1 mod 4, one
component of γ4 is p^e times a unit form; no four-form coefficient pair has a common odd factor).

**The code, line by line.** `f4_form.residual5`: γ4 = (c, c·i) as (u-coefficient, v-coefficient) — correct;
support {2} ∪ primes(C) — correct; same bad primes as `residual` — fine (any set of primes is sound); the full-set
skip and `X.join` are generic in J. `exact_local.ExactLocal._disc1`: when `va = vb` (both coefficients ≡ 0 mod p^v)
it factors out p^vb and labels (vb mod 2, χ(b′(x − r))) — correct for the common factor; `_disc3` is never reached
with a common factor (no prime ≡ 3 mod 4 divides m² + n²); `TD.Local` at 2 (`key2`) is per value and generic. The
unit closure (×3 at 2, flip at p ≡ 1 mod 4) is the same argument for five forms. I found no line that assumes four
forms. f4_all.py counts a pairing value "outside" as not killed, but no such value occurs in the data.

**Tests** (`z16_witness.py`, my own class code: Lemma S1 in integers for global points, valuations and Legendre
symbols at odd p, squares of units mod 8 at 2; seeds fixed):
- A. Planted rational points: every coprime degenerate ratio m ≤ 200 (60) and 40 random degenerate ratios with m up
  to ~10⁵ (both types), 200 rational points: 4,688 (point, prime) checks at every bad prime and the primes < 60,
  **0 refused**; and the end-to-end join (`z16_five.residual`, with the d4 support) keeps all 200 tuples.
- B. 46,255 ℚ_p points with five rational classes, random and near every root of every component (both charts,
  depths 1 to 8; the common factor of γ4 divided out first — my first version missed the F4 roots there), at all
  bad primes < 10⁷ and 3, 5, 7, 11, 13, 17, 29, 37, 41 for nine ratios (ghosts, (627,364), (465,416), (12,5), (7,1)):
  **0 outside the allowed set.** Primes with 0 points are exactly primes where the five-form allowed set is empty
  (2621, 6229 for (3320,2303); 17, 41, 109, 349 for (90384,58537); 7, 53 for (115513,18217)): 484-class primes and
  primes ≡ 5 mod 8 of m² + n².
- C. My class code against `tuple_descent.Local.cls` at all 46,255 points: 0 disagreements.
- Session's own self-test rerun (imports replaced, body unchanged): see below.
**Verdict: sound. No [blocker] found.** (Caveat: B samples points; it is not a proof of the allowed sets.)

- Session's self-test rerun (`sess/f4_selftest_rerun.py`: the session's body unchanged, imports pointed at the
  repository through `z16_five.py`): **output identical** to `f4_selftest.txt` (37,974 points with `TD.Local`,
  27,873 with `exact_local`, 0 failures). [note] It only uses bad primes below 200 of 2mn(m² − n²)(m² + n²), so it
  never tests a large prime of m² + n² (where γ4 has the common factor) nor an octic prime; my test B adds 12644473
  for (14903,9696) and (24599,5207) (327 and 302 points, 0 refused) and all bad primes below 10⁷ for nine ratios.

## Task 3, continued

- Samples beyond m = 15000 (my run on all 321 sample ratios): four-form tuples equal to the session's at every ratio;
  biased 46 of 269 keep a 5-tuple, 118 of 1,096 4-tuples extend, 6 pairing survivors, 0 extend; random 4 of 14;
  few-primes 4 of 38. **Confirmed.**
- Nine open census ratios (`out_open.jsonl`): the five-form test at the bad primes leaves exactly 2 tuples at each,
  and they are exactly the tuples of the two degenerate rational points (computed by Lemma S1 in integers); four forms
  left 4 at (609,200), (627,364), (809,409), (991,263) and 2 at the other five. **Confirmed.** [note] With good primes
  the four-form test already leaves exactly the degenerate tuples at these ratios (LAST_PAIR §12, e.g. (627,364)), so
  F4 here only replaces good primes; nothing new is excluded.
- Times (2 cores, PARI 2.17 via cypari2): 937 ratios 839 s of CPU (4 and 5 forms each; max 6.1 s at (9216,4225));
  321 sample ratios 190 s; open ratios 4 s.

## Task 4: §3, ghost data and laws

- "18 surviving tuples at 9 ratios, all of ghost shape": the 16 in the data files are S1 or S2 (8 and 8); the
  partner's 2 (from `t1a.txt`, not in a data file; I did not recompute the pairing) are S2 by the definition.
  Residual tuples at (24599,5207): 8 (my run). 321 sample ratios with 1,366 tuples: confirmed. The biased file has
  no closing "done" line and its last `tried` is 5,892: the caveat is accurate. [fix] The 9,736 ratios that §3 and
  `law_a.py` call "stage-1 ratios" are the class-keeping ratios of the four-form system (COMPLETE.md §7, built from all
  coprime ratios); they are not the census's stage 1 (the 484 classes), and 866 of the 937 fail the census's stage 1
  or the 5 mod 8 theorem. All 18 ghost tuples sit at ratios in the 484 classes.
- **Law (a)** (`z16_laws.py`, from EXPLORATION_4's tables and the repository's lists): S1 shape at 333 ratios, S2 at
  166, never both; odd-exponent condition 333/333 and 166/166; holds without such a tuple at 1,554 and 1,392; the
  primes not ±1 mod 8 that occur do so to even power and are all ≡ 3 mod 8 (94 and 46 ratios); S1 ratios all of
  opposite parity, S2 all both odd. **Reproduced exactly.** Status "necessary, tested" is right.
- **Law (b)**: my own Sel² local images (random ℚ_p points, exact rationals, image size checked = 4 or 8) on a fresh
  sample (seed 2016, 3,000 ratios m < 4000, C not a square): the rule agrees on 3,000 of 3,000 (2,019 opposite
  parity: in Sel² exactly when no prime ≡ 5 mod 8 divides m² − n², 629 in; 981 both odd: never). [note] The
  session's totals 1,994 / 1,006 are not in the shipped `law_b.txt` (it prints other splits; `law_b.json` is not in
  the bundle), so they could not be checked from its output; my sample replaces them.
- **Law (b), proof** (given note2's Lemmas 1 and 2, which I re-derived at odd primes): E_08: y² = x(x − A²)(x − C²),
  κ = (1, β, β), β = sqf(2mn). Torsion images are only {(1,1,1), (1,−1,−1)} at every place ((0,0) and (A²,0) have the
  same image). ∞: β > 0, κ is the class of x > C². Odd p ∤ ABC: unramified. Odd p | B: Lemma 1, the image is all of
  {(1,γ,γ)} (4 = 4). Odd p | A = (m − n)(m + n): Lemma 2 forces β ∈ ±□ in ℚ_p, and then κ is a torsion image;
  β ≡ ±2n² mod p, so the condition is "±2 is a square mod p", i.e. p ≢ 5 (mod 8). Odd p | C: m ≡ ιn, 2mn ≡ 2ιn² =
  ((1 + ι)n)², always a square. p = 2, opposite parity, 4 | mn: H2 holds (C ∓ A = 2n², 2m², 32 | the even one), and
  Lemma 1 gives the image {(1,γ,γ)} of size 8. Both odd: on the curve of the halved triple A′ is even and Lemma 2 at 2
  needs ±β ∈ □₂, impossible since v₂(β) = 1. **So law (b) is proved for 4 | mn and for m, n both odd (given note2's
  Lemmas 1, 2 and Cassels), and tested only for mn ≡ 2 (mod 4)** (H2 fails there; 316 of 316 such ratios in my sample
  in Sel² under the 5 mod 8 condition). All ghost ratios have 4 | mn or are both odd, so for them it is proved.
- Law (c), ranks (§3 ranks bullet), not checked (time). [note] (3320,2303) and (5623,1017), where PARI's [1,1] is
  quoted, are excluded for a magic square by tier-A congruences, so these ranks bear only on the four-form system.

## Claims editor, §2–§3

| Claim (EXPLORATION_5) | Verdict |
|---|---|
| §2 title, "the centre form is not in the stage-2 local test" | true of the four-class descent and all of G1/G2/H1–H4/T10/pairing; [fix] add that this is by design and documented (LAST_PAIR §11 "G4 is not used"; `CoverDescent.lean` header; LEDGER "Without the fifth form"), and that `gauss_descent.py` (GAUSS_DESCENT.md) is a five-form descent run to m ≤ 5000 |
| 937 → 211 ratios; 6,904 → 1,370; 5 → 1 pairing-survivor ratio (14903,9696) | confirmed (computed, my driver on the repository's code; same at every ratio as `f4_all.jsonl`) |
| [fix] missing context | 866 of the 937 are already excluded for a magic square by tier A (844 by the 484 classes, 318 by a prime ≡ 5 mod 8 to odd power in m² + n²); the F4 effect on the other 71 is 71 → 29 ratios, 688 → 186 4-tuples; none of the 71 has a pairing survivor |
| 269 → 46 sampled ratios; 6 survivors → 0 | confirmed; [note] 268 of the 269 are in the 484 classes |
| six degenerate ratios keep their solutions | confirmed, and for 100 degenerate ratios (200 rational points) |
| self-test 0 failures on 37,974 + 27,873 | confirmed (rerun identical); [note] small bad primes only |
| nine open ratios: exactly 2 tuples each; four forms left 4 at four of them | confirmed, and the 2 are the degenerate tuples; [note] good primes already did this (LAST_PAIR §12), nothing new is excluded |
| "If it does not, this is a free strengthening of stage 2: no new prime, no new curve" | [fix] the census does impose F4 (stage 1, G1Primes, the Gaussian descent); as a strengthening of the bad-primes test it is real but small for magic squares (71 → 29); no count changes (agree) |
| "Next: pair curves y² = F_j·F4 for (14903,9696) and its partner, the only ghost ratios left" | [fix] both ratios are in the 484 classes (mod 89, 103, 151): no magic square has them (tier A); the step has no object for the five-form problem, and F4 is not part of the four-form problem |
| §3 partner (24599,5207): 8 residual tuples, 2 survivors, ghost shape | 8 confirmed; survivors not recomputed; the two listed are S2 |
| 321 ratios, 1,366 tuples, 6 survivors at 3 new ratios; 18 at 9, all ghost shape; stopped by hand, ≥ 5,892 tried | confirmed from the data; all 9 ratios in the 484 classes, 5 of them also fail the 5 mod 8 theorem |
| law (a) 333/333, 166/166 (and 1,554, 1,392, 94, 46) | confirmed exactly |
| law (b) 1,994/1,994 and 1,006/1,006 | totals not in the shipped output; fresh sample 3,000/3,000; **proved** for 4 \| mn and for m, n both odd (given note2's Lemmas 1, 2), tested for mn ≡ 2 mod 4 |
| "stage-1 ratios" (law (a), sampler) | [fix] they are the class-keeping four-form ratios, not the census's stage 1 |
| law (c); ranks at the six ghost ratios | not checked |

## Proposed LEDGER wording (tier D, not released)

| **The centre added to the four-class descent at the bad primes** (run 11; EXPLORATION_5 §2, a separate AI research
session; reviewed by a separate AI agent, z16: holds, with context): with `γ4 = (1+i)(m+in)(u+iv)` (norm `F4`) as a
fifth form, its class supported on 2 and the primes of odd exponent in `m² + n²`, the same local code and bad primes,
211 of the 937 ratios `m ≤ 15000` whose four-class set survives the bad primes keep a 5-tuple (1,370 of 6,904 tuples
extend; of the pairing survivors only the two of (14903,9696)); at the nine open ratios exactly the two degenerate
tuples remain. The 937 are ratios of the four-form system: 866 of them, and all nine "ghost" ratios, are already
excluded for a magic square by tier A (844 by the 484 residue classes, 318 by a prime `≡ 5 (mod 8)` to an odd power in
`m² + n²`, `G1Primes.lean`); on the other 71 the centre leaves 29. No count changes | `experiments/referee_z16/`
(EXPLORATION_5 not committed) | `z16_witness.py` (planted degenerate points, 46,255 local points, 0 refused) |

## Status

Done within the time box. Not done: law (c); the ranks bullet; the pairing of the partner's 8 tuples; a proof of the
2-adic case mn ≡ 2 (mod 4) of law (b); an independent re-implementation of the exact local sets themselves (my
checks are witnesses from below and the repository's code from above).
