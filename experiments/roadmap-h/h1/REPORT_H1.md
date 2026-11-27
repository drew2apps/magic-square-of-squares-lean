# H1: the phantoms of the four-class descent

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Machine: PARI 2.19.1 with cypari2, not CI's PARI 2.17.4. Real sympy 1.14.0. PARI is a finder only (points, conic points, torsion points, factorisations), except where a line says "rank by PARI".

A **phantom** is a ratio whose four-class descent keeps a tuple beyond the tuples of its degenerate points t = 0, ±1, ∞ (the **known** tuples). The extra tuples are then classes that no known point explains. C is the curve where F_0, F_8, F_L, F_U are all squares.

## Bottom line

- **Phantoms are rare, occur only at degenerate ratios, and always come in partner pairs.** The scan of 20,009 random degenerate ratios with 10⁴ ≤ m ≤ 10⁹ finds 13 phantoms. With G1's (1109503,597671) that makes 14. The partner (m,n) → ((m+n)/g, (m−n)/g) of every one of the 14 is a phantom too, with the same numbers of known and extra tuples, so they form 12 pairs (two pairs have both sides in the scan). *(computed, exact)*
- **Their extra tuples are local everywhere.** Each phantom keeps its extra tuples through every bad prime, ℝ and good primes up to 20,000. By H3, good primes above 1153 never cut anything, so these tuples are in the full local four-class set. The descent's local test cannot remove them. *(computed, exact; H3 rests on Hasse–Weil, quoted)*
- **In 11 of the 12 pairs, every extra tuple is the class of no rational point of C.** On one side of each pair, either a rank-0 elliptic quotient of C forces every point of C over t ∈ {−1, 0, 1, ∞} (7 pairs, part 6), or the pair curves exclude each extra class (4 more pairs, parts 4 and 7). Then every rational point of C there has a known tuple, and the phantom is explained: it is a Selmer-type class with no point. *(computed, exact)* The step from one side of a pair to the other is now computed too: part 8 writes the partner isomorphism down, t′ = (t + 1)/(t − 1), and checks exactly that it carries known tuples to known tuples.
- **The partner isomorphism is explicit (part 8).** For all 12 pairs, in both directions, (u, v) → (−u − v, −u + v) carries the forms z_0, z_8, z_L, z_U of one side to z′_L, z′_U, z′_0, z′_8 of the other, times −1 or −4 (classes 1 or 2), as an identity of linear forms. No image of a known tuple is an extra tuple of the partner. *(computed, exact)*
- **1 pair stays open:** (11875080,6497849) / (18372929,5377231). On both sides the pair curves have 3 ≤ rank ≤ 4 and the extra classes are classes of rational points found on them, so they cannot decide; every rank-0 quotient tried has rank at least 1. *(open)*
- **(1109503,597671), the case that started H1, is closed twice.** Its pair curves E_0L and E_0U have rank exactly 3, and each extra tuple has a class that is not the class of any rational point of one of them (h1_pairs.txt). Independently, y² = w·f0·f8 has rank 0 exactly (h1_rank0_exact.txt). *(computed, exact)*

## Part 1: (1109503,597671) on its pair curves (h1_pairs.py, h1_pairs.txt)

- For each extra tuple and pair (j, k), the class of d_j d_k is compared with the classes of the known tuples. Only E_jk: y² = F_j F_k(t) with a new class can exclude anything.
- Each useful E_jk gets an integral model with full 2-torsion, a complete 2-descent and the Cassels–Tate pairing (the repository's rank_second.py), points from PARI as a finder (each checked on the curve), and the class x − e of the descent map matched to D_jk / D_jk(t = 0) on 8 or more points.
- **Result:** E_0L and E_0U have rank exactly 3 (dim Sel² = 7, pairing of rank 2, 3 points and the 2-torsion spanning dimension 5). The extra tuples (1109503, 597671, 2, 959) and (1195342, 2219006, 2, 959) are killed by E_0U; (1109503, 597671, 1918, 1) and (1195342, 2219006, 1918, 1) by E_0L. The two known tuples are not excluded, as they must not be.
- **A second test, when the rank bounds differ:** a class that is not the class of any element of the kernel of the Cassels–Tate pairing is excluded too (the "CT kernel test").
- **Negative controls (h1_pairs_negative.txt):** the pairing left out (not decided, as it must be); a wrong root (refused); a forged tuple (1109503, 597671, 6, 1) (excluded).

## Part 2: the scan (h1_scan.py, h1_scan_k13.txt to h1_scan_k29.txt)

1,177 random coprime degenerate ratios per dyadic range k = 13 to 29 (seed 3000 + k), 20,009 in all, each through G1's `one()` (class criterion, then the four-class descent with every odd prime exact and good primes to 150, 1500, 5000, 20000 while more than the known tuples is left).

| k | phantoms | time | | k | phantoms | time |
|---|---|---|---|---|---|---|
| 13 | 0 | 232 s | | 22 | 0 | 581 s |
| 14 | 0 | 201 s | | 23 | 1 | 1063 s |
| 15 | 0 | 252 s | | 24 | 0 | 666 s |
| 16 | 0 | 237 s | | 25 | 1 | 1451 s |
| 17 | 1 | 339 s | | 26 | 0 | 1084 s |
| 18 | 0 | 286 s | | 27 | 4 | 2618 s |
| 19 | 3 | 440 s | | 28 | 2 | 2123 s |
| 20 | 1 | 629 s | | 29 | 0 | 2972 s |
| 21 | 0 | 615 s | | | | |

- No ratio failed (a degenerate ratio always keeps its known tuples).
- 3 of the 13 keep a class of dimension 1 to 3; the other 10 have class dimension 0.
- **Negative control (h1_scan_negative.txt):** three ratios with the known tuples withheld are listed, as they must be, and (1109503,597671) comes out as a phantom.

## Part 3: the partners (h1_partners.py, h1_partners.txt)

- The partner map (m,n) → ((m+n)/g, (m−n)/g) gives isomorphic systems (CENSUS9000.md). G.one is run on the partner of each of the 14 phantoms.
- All 14 partners are phantoms, with the same numbers of known and extra tuples. The 12 pairs:

| ratio | partner |
|---|---|
| (204724,97893) | (302617,106831) |
| (413004,148603) | (561607,264401) |
| (756492,255635) | (1012127,500857) |
| (853587,255916) | (1109503,597671) |
| (1992873,999664) | (2992537,993209) |
| (11875080,6497849) | (18372929,5377231) |
| (32396425,3694248) | (36090673,28702177) |
| (175775580,15972419) | (191747999,159803161) |
| (189089400,156063271) | (345152671,33026129) |
| (231641256,72371617) | (304012873,159269639) |
| (234319635,115640732) | (349960367,118678903) |
| (518804652,31743595) | (550548247,487061057) |

## Part 4: the pair curves on every phantom of the scan (h1_allphantoms.py, h1_allphantoms.txt)

- Part 1's method, run on the 13 scan phantoms. PARI's points are pooled over the pair curves that share an integral model.
- **All extra tuples killed** at (561607,264401), (1012127,500857), (36090673,28702177) and (191747999,159803161).
- Elsewhere the pair curves decide nothing. Either the extra classes are classes of rational points found on E_jk (for example (349960367,118678903), h1_phantom.txt: rank 3 to 5, 10 points pooled), or the curve stops: "no base point t0 found" (32 times, mostly square type) or "the root of the descent map is not identified".
- A stopped curve decides nothing; a known tuple is never excluded (asserted).

## Part 5: rank-0 quotients of C, PARI's ranks (h1_rank0.py, h1_rank0.txt)

- With w = t², the repository's global.py lists six genus-1 quotients of C that do not use the fifth form g4: y² = w·f0·f8, f0·K, f8·K, w·K, w·f0·K, w·f8·K. If one has rank 0, every point of C lies over its torsion points.
- **Result:** 8 of the 14 phantoms have a rank-0 quotient (rank by PARI): w·f0·f8 at five, w·K at three. In every case the torsion points lie over t ∈ {−1, 0, 1, ∞} only, and their tuples are exactly the known ones.
- The g4 quotients (w·f0·g4 and the rest) are not quotients of C. At the two holdouts tried, (36090673,28702177) and (518804652,31743595), and their partners, all of them have PARI rank at least 1 anyway (a scratch probe, not in the bundle).

## Part 6: the same, with exact rank bounds (h1_rank0_exact.py, h1_rank0_exact.txt)

- **Curve A,** y² = w·f0·f8, as y² = x(x + m⁴)(x + n⁴) with x = m²n²w: complete 2-descent and Cassels–Tate pairing (rank_second.py), taken over up to 3 curves of its isogeny class with full rational 2-torsion, reached by 2-isogenies (Vélu). The rank is the same on isogenous curves (quoted, standard).
- **Curve B,** y² = w·K, as y² = x(x² + (2S² − 16m²n²)x + S⁴) with x = S²w: rank_indep.py's descent via 2-isogeny with its exact second descent. Its conic points come from rank_second.find_conic_point_any, because PARI's first point has N = 0 here.
- **Torsion:** PARI's elltors as a finder; each point checked exactly, and the count equal to the gcd of #E(F_p) over good odd p < 400, so the list is all of E(ℚ)_tors.
- **Result:** rank 0 exactly on curve A for one side of 7 of the 12 pairs: 413004/561607, 853587/1109503, 1992873/2992537, 175775580/191747999, 189089400/345152671, 231641256/304012873, 234319635/349960367. This agrees with PARI's 8 of 14 in part 5, pair by pair. The isogeny class matters: at (349960367,118678903) the first curve gives rank ≤ 2 and an isogenous one gives 0.
- On 4 isogenous curves the pairing's local search gave up ("no local point found at 2"); those curves give no bound and are reported as such.
- **Negative controls (h1_rank0_exact_negative.txt):** a withheld torsion point is refused; a forged "extra" equal to the t = 0 tuple is not excluded; a known tuple not over the torsion is refused.

## Part 7: the partners of the open pairs (h1_holdouts.py, h1_holdouts.txt)

- h1_holdouts.py reads the pairs left open by part 6 (h1_rank0_exact.txt), drops those already settled by part 4, and runs part 4's pair-curve analysis on the other side of each, the partner, with the known and extra tuples of h1_partners.txt. 3 partners: (302617,106831), (18372929,5377231), (550548247,487061057).
- **Result:** every extra tuple of (302617,106831) and of (550548247,487061057) is killed by E_0L, E_0U, E_8L and E_8U (the Cassels–Tate kernel test: the class is not the class of any element of the kernel of the pairing). At (18372929,5377231) the four pair curves share a model with 3 ≤ rank ≤ 4, and both extra classes are classes of rational points found, so nothing is decided.
- With part 4's kills at (1012127,500857) and (36090673,28702177), this settles 4 pairs beyond part 6's 7.
- **The tally per pair:**

| pair | settled by |
|---|---|
| (204724,97893) / (302617,106831) | pair curves at 302617 (part 7) |
| (413004,148603) / (561607,264401) | rank 0 at 561607 (part 6); pair curves at 561607 (part 4) |
| (756492,255635) / (1012127,500857) | pair curves at 1012127 (part 4) |
| (853587,255916) / (1109503,597671) | rank 0 at 1109503 (part 6); pair curves (part 1) |
| (1992873,999664) / (2992537,993209) | rank 0 at 1992873 (part 6) |
| (11875080,6497849) / (18372929,5377231) | **open** |
| (32396425,3694248) / (36090673,28702177) | pair curves at 36090673 (part 4) |
| (175775580,15972419) / (191747999,159803161) | rank 0 at 175775580 (part 6); pair curves at 191747999 (part 4) |
| (189089400,156063271) / (345152671,33026129) | rank 0 at 345152671 (part 6) |
| (231641256,72371617) / (304012873,159269639) | rank 0 at 304012873 (part 6) |
| (234319635,115640732) / (349960367,118678903) | rank 0 at 349960367 (part 6) |
| (518804652,31743595) / (550548247,487061057) | pair curves at 550548247 (part 7) |

- Where a pair is settled on its scan side by part 4 (1012127, 36090673, 561607, 191747999), that side needs no partner step; the partner follows by part 8.

## Part 8: the partner isomorphism (h1_iso.py, h1_iso.txt)

- The forms, with t = u/v: z_0 = mv + inu, z_8 = nv + imu, z_L = (m + n)(u − v) + i(m − n)(u + v), z_U = (m − n)(u − v) + i(m + n)(u + v) (checked against h1_pairs.forms).
- The script searches the integer matrices with entries in {−1, 0, 1} and the bijections of the four forms for which z′_s(j)(M(u, v)) = λ_j z_j(u, v) identically, in exact Gaussian arithmetic, with N(λ_j) a square.
- **Result:** for all 12 pairs and both directions, M = (u, v) → (−u − v, −u + v), that is t′ = (t + 1)/(t − 1), and s = (0, 8, L, U) → (L, U, 0, 8). The λ_j are −1 and −4 when g = 1, and their halves when g = 2, so the tuple of the image point is σ(T) times one class e for all four forms.
- **Transfer check:** for each known tuple T and each class e that keeps the image in the partner's tuple space, the image is not an extra tuple of the partner. Since every rational point has a tuple in known + extra, every rational point of the partner's curve then has a known tuple. Holds for all 12 pairs, both ways. The images of the known and the extra tuples also reach all of the partner's known and extra tuples (consistency).
- **Negative controls (h1_iso_negative.txt):** a forged partner (n′ + 2) has no isomorphism of this shape; a forged bijection (0, 8, L, U) → (L, U, 8, 0) fails the identity; an image of a known tuple planted among the partner's extra tuples is refused.

## What it rests on

- **Exact:** the descents, the pairings, the torsion counts, the tuples and every check of a point. No floating point in any certificate.
- **Quoted:** isogeny invariance of the rank; Hasse–Weil (through H3); the soundness of the four-class descent (every rational point of C has a tuple in known + extra).
- **The step between the two sides of a pair (part 8):** a rank-0 quotient on one side proves that every point of C there lies over t ∈ {−1, 0, 1, ∞}; the pair curves prove that every point has a known tuple. Part 8 carries either statement to the partner: a rational point of the partner's curve is the image of a rational point of C, its tuple is the image of a known tuple, and no such image is an extra tuple of the partner, so it is a known tuple. Computed (exact) for all 12 pairs; it no longer rests on CENSUS9000.md's description of the partners.
- **Rank by PARI:** part 5 only, consistent with PARI. Part 6 replaces it with exact bounds for the same 7 pairs.

## What failed

- The pair curves stop at most square-type phantoms ("no base point t0 found" within |t| ≤ 49).
- A first run of h1_rank0_exact stopped on every curve B ("conic point with N = 0"). Fixed by the repository's own find_conic_point_any. Curve B then gave no new rank 0.
- A bound through 2-isogeny descents on A's own models was useless (they all have full 2-torsion) and was removed before the committed run.
- An isogeny-class probe on (349960367,118678903) with PARI (scratch) was stopped; part 6 settles that ratio exactly.

- **New quotients for the open pair (scratch, not in the bundle).** The involutions t → 1/t and t → −1/t of C give further genus-1 quotients over s = t + 1/t or s = t − 1/t: E± from F_0F_8F_L, and two even quartics Q1: y² = (As² + B)(S²s² − 16A), Q2: y² = (As² + S²)(S²s² + 4B), with A = m²n², B = (m² − n²)², S = m² + n². None has rank 0 on the open pair: E± have PARI bounds [1,3] and [2,2] at (11875080,6497849); the exact 2-descent and pairing over the isogeny class of Q1 and Q2 give rank ≤ 4, 3 there and ≤ 2, 3 at (18372929,5377231), where PARI finds Q1's rank to be 2. A Mordell–Weil sieve on two pair curves at the same t is the next tool.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** The phantoms of the four-class descent (experiments/roadmap-h/h1). A scan of 20,009 random degenerate ratios up to m = 10⁹ finds 13 ratios whose descent keeps tuples beyond those of the degenerate points; with (1109503,597671) and the partners, they form 12 partner pairs, and the extra tuples pass every local test. In 11 pairs every extra tuple is shown to be the class of no rational point, by an exactly computed rank-0 quotient of the four-form curve or by the pair curves' descent and Cassels–Tate pairing; 1 pair, (11875080,6497849) / (18372929,5377231), remains open. The partner isomorphism t′ = (t + 1)/(t − 1) is written down and checked exactly, so the step from one side of a pair to the other is computed too. Computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-h/h1        # F = grep -v '^PARI stack size'; G1_DIR, MSS_REPO as in the scripts
python3 h1_pairs.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'          # about 113 s
python3 h1_scan.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'           # about 18 s
python3 h1_rank0_exact.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'    # about 4 s
python3 h1_iso.py 2>&1 | $F | tail -1 | grep -q '^# matches'                       # about 3 s
python3 h1_iso.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'            # about 1 s
# not for CI (times.txt): h1_pairs.py, h1_scan.py k (200 to 3000 s per range), h1_partners.py (429 s),
#   h1_allphantoms.py (2219 s), h1_rank0.py (1036 s), h1_rank0_exact.py (1783 s), h1_holdouts.py (1941 s)
```
