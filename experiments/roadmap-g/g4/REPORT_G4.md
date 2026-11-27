# G4: a search for a solution in the square-type fibres

> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds: "no solution" means none of the points listed gives one; 23 ratios had no curve searched; a search, not a certificate.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Run at 062c54f (0.6.0 release candidate); the programs used (scripts/g1_full/global.py, sigma_exact.py) are identical at 24d6ef6. Real sympy 1.14.0 is installed. Machine: PARI 2.19.1 and cypari2, not CI's PARI 2.17.4.

## Bottom line

- **No solution found.** On every coprime square-type ratio with m ≤ 20,000 (3,573 ratios, square_20000.json), no t from the points listed makes all five forms squares, and no non-degenerate t makes even four of them squares. *(computed (numerical): a search, not a certificate; every t found is tested exactly)*
- **The search is far from complete, and how far differs by curve.** It lists every point up to a height bound on the subgroup that PARI's points generate. Only 1,741 of the 3,573 ratios have a full search on E+ or E−, where PARI's rank bounds meet and it returns that many independent points. 1,473 more have only a partial one (fewer points than the rank bound). 359 have neither, and 23 have no curve searched at all.
- **The planted test passes: 15,724 of 15,724.** On every search of a curve that carries the known degenerate point t = ±1, that point is listed.

## How the search works (g4_search.py, g4_pass2.py)

A solution t of the five-form system (F0, F8, F4, L, U all squares; global.py's classify_t) gives a rational point on 14 curves:
- E+ and E−: Y² = (z ± 2) g(z) with z = t + 1/t (UNIFORM_RANK.md §1, sigma_exact.model(m, n, ±2));
- the 12 quotients y² = P(w), w = t², of global.py's curves().

For each ratio and curve:
1. PARI's ellrank gives rank bounds lo ≤ hi and some independent points.
2. Where lo = hi ≤ 3 and PARI returns lo independent points, those points are saturated at primes ≤ 100 (ellsaturation). Every point of canonical height ≤ 64 · max h(gens), torsion included, is listed by Fincke–Pohst (qfminim) on the height matrix. That is K = 8 times the largest generator in height.
3. Each point is mapped back to t (z² − 4 a square, or w a square). For each t, the five forms are evaluated in exact rational arithmetic and classify_t is called. Every t with four or five squares is recorded.

The second pass (g4_pass2.py) handles E+ and E− where lo = hi but PARI returned fewer than lo points. It reruns ellrank at effort 2, adds the point of t = ±1 when it is independent, and searches the subgroup it has in the same way. That search is marked partial: the missing direction is not searched.

## Results (g4_summary.txt)

```
square-type ratios: 3573 (m <= 19995); second pass (E+, E-) run on 2488 of them
SOLUTIONS: none
non-degenerate t with at least four of the five forms squares: none
ratios by their best search on E+ or E-: full 1741, partial 1473, none 359
ratios with no curve searched at all: 23 [(5561, 1200), (8911, 6600), (10125, 4532), (10961, 1680), (11543, 2376), (12852, 11245), (13837, 4284), (14111, 4080), (14985, 13832), (15120, 15089), (15561, 10480), (16095, 4352), (16095, 13328), (16109, 12420), (16815, 11152), (16975, 10608), (17139, 15620), (17363, 11316), (17711, 15840), (18688, 11055), (18721, 13920), (19000, 9849), (19599, 280)]
planted test: searches on a curve carrying the point of t = +-1: 15724; in which it is listed: 15724
per curve (full: rank bounds equal and a full set of PARI points, every point of height <= 64 max h(gens) listed;
  partial: rank bounds equal but fewer points, the subgroup they generate searched the same way):
  E+          partial (pass 2) 1506, skipped: rank bounds differ 1058, full 917, full (pass 2) 92
  E-          partial (pass 2) 1509, skipped: rank bounds differ 1019, full 941, full (pass 2) 103, skipped: rank above 3 1
  f0*K        skipped: too few points, not in pass 2 2118, skipped: rank bounds differ 996, full 435, skipped: rank above 3 24
  f0*g4*f8    skipped: too few points, not in pass 2 1689, full 1126, skipped: rank bounds differ 754, skipped: rank above 3 4
  f8*K        skipped: too few points, not in pass 2 2255, skipped: rank bounds differ 968, full 326, skipped: rank above 3 24
  g4*K        skipped: too few points, not in pass 2 1494, skipped: rank bounds differ 1300, full 774, skipped: rank above 3 5
  w*K         skipped: rank bounds differ 1853, skipped: too few points, not in pass 2 1183, full 537
  w*f0*K      skipped: too few points, not in pass 2 2252, skipped: rank bounds differ 964, full 332, skipped: rank above 3 25
  w*f0*f8     skipped: too few points, not in pass 2 1677, full 1394, skipped: rank bounds differ 502
  w*f0*g4     full 2654, skipped: too few points, not in pass 2 805, skipped: rank bounds differ 90, skipped: rank above 3 24
  w*f0*g4*f8  skipped: too few points, not in pass 2 1688, full 1126, skipped: rank bounds differ 754, skipped: rank above 3 5
  w*f8*K      skipped: too few points, not in pass 2 2100, skipped: rank bounds differ 993, full 450, skipped: rank above 3 30
  w*f8*g4     full 2656, skipped: too few points, not in pass 2 803, skipped: rank bounds differ 90, skipped: rank above 3 24
  w*g4*K      skipped: too few points, not in pass 2 1491, skipped: rank bounds differ 1300, full 777, skipped: rank above 3 5
```

"skipped: too few points, not in pass 2" on the quotient curves: the second pass was run on E+ and E− only.

Coverage on the open ratios of square type: the five of the census ≤ 1000, (7735,4968) from F0, and the three from G3:

| ratio | E+ | E− | quotient curves with a full search |
|---|---|---|---|
| (195,28) | not searched, bounds [0,4] | partial | 4 |
| (440,279) | not searched, bounds [1,3] | partial | 2 |
| (609,200) | partial | not searched, bounds [0,4] | 8 |
| (627,364) | not searched, bounds [1,3] | not searched, bounds [1,3] | 3 |
| (861,620) | partial | not searched, bounds [0,4] | 2 |
| (7735,4968) | full (pass 2) | partial | 5 |
| (11165,7548) | partial | partial | 2 |
| (12201,7360) | not searched, bounds [0,2] | partial | 4 |
| (12204,9853) | not searched, bounds [1,3] | full (pass 2) | 0 |

No point found on any of them gives a non-degenerate t with four squares. That says nothing about whether they are open.

## What holds and how it was checked

- **Exact where it counts:** each t is a Fraction; the five forms, squareness and classify_t are exact. A reported solution would be a real one.
- **Planted test:** on every search of a curve whose generators carry the point of t = ±1 (a degenerate point with all five forms squares on a square-type ratio), the search lists both t = 1 and t = −1. 15,724 of 15,724 searches, both passes. g4_search.py and g4_pass2.py stop with an error if one fails.
- **Negative mode (g4_search.py --negative, g4_negative.txt):** with classify_t told that every t with five squares is a solution, the search must report t = ±1 as SOLUTIONS on (12,5), (24,7) and (40,9). It does: 3 of 3. So the reporting path is live.
- **Input (gen_square.py):** square_20000.json is every coprime ratio with m² + n² a square and m ≤ 20,000, enumerated from (a² − b², 2ab). No congruence filter, so the census ratios m ≤ 1000 are in it too.

## What does not hold, and why this is a search, not a certificate

- **Floating point:** the height matrix and Fincke–Pohst bound are floating point (PARI's ellheightmatrix and qfminim). A point near the bound can be missed or listed. So the list "every point of height ≤ H" is numerical, not proved.
- **The rank is PARI's, not proved.** Where lo = hi, the rank is taken as lo. PARI's ellrank lower bound can exceed the number of points it returns (its documentation: r > #L), and then lo rests on the 2-part of Sha being of even rank. On those curves the rank is **conditional**, and the second pass searches only part of the group.
- **Generators:** PARI's, saturated only at primes ≤ 100. A group of index divisible by a larger prime would make the height bound miss points.
- **Coverage:** curves with rank bounds that differ, or rank above 3, are not searched (counts per curve above).
- **PARI version:** point counts, heights and which curves PARI decides depend on the PARI version and ellrank's randomness. The committed .json outputs are from PARI 2.19.1 and are not expected to match bit for bit under 2.17.4. That is why only the bookkeeping and the negative mode are proposed for CI.

## What failed

- The first pass left many curves with lo = hi but too few points (about half of E+ and E−). The second pass at effort 2 found the full set on only 195 of those 3,210 curves (92 on E+, 103 on E−); the other 3,015 are partial.

## Next step

A complete search on a curve needs the rank proved and the generators saturated at every prime up to the index bound. A height bound tied to how large a solution's t would have to be would turn "no point up to height H" into a statement. That is not attempted here.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** A point search on the square-type fibres (experiments/roadmap-g/g4): for the 3,573 coprime square-type ratios with m ≤ 20,000, every point up to 64 times the largest generator's height, on E± and the 12 quotient curves, on the subgroup that PARI's points generate, mapped to t and tested exactly. No solution, and no non-degenerate t with four of the five forms squares. A search, not a certificate: PARI rank bounds (conditional where the lower bound exceeds the points found), saturation at primes ≤ 100, floating-point heights; 1,741 ratios have a full search on E+ or E−. Computed (numerical). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-g/g4        # F = grep -v '^PARI stack size'
python3 g4_summary.py 2>&1 | $F | tail -1 | grep -q '^# matches'            # under 1 s
python3 g4_search.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'  # about 7 s
python3 gen_square.py 20000 | tail -1 | grep -q '^# matches'                # about 12 s
# not for CI: bash run_g4.sh i 3 (first pass, about 4.6 h per worker on three workers) and
#             bash run_g4_pass2.sh i 2 (second pass on E+, E-; times.txt)
```

g4_search.py finds the repository through MSS_REPO, or else `../../..`.
