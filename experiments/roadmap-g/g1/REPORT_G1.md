# G1: stress test of (R) on random ratios up to m = 10⁹

> **Read `NOTES.md` first.** This is the session's report, committed unchanged; where they differ, `NOTES.md` (after review x6) holds: `m² + n² = 2·(5³·7129)²`; "no prime left out" is false for G1; the four extra tuples of (1109503,597671) are locally soluble everywhere, and whether a rational point has one of these classes is open ("classes that no point explains" is not shown).

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Started from 24d6ef6, before 0.6.0 was on main. The runs finished at 062c54f, the release-candidate commit whose CHANGES.md starts with "- 0.6.0". None of the files G1 uses changed between the two commits: tuple_descent.py, tuple_bruteforce.py, structure/, g1_full/global.py and sigma_exact.py are identical. Real sympy 1.14.0 is installed (`pip install sympy==1.14.0`). It is used only where the repository's own modules import it; every factorisation in G1 is made by PARI, and each prime factor is proved prime by `isprime`. Machine: PARI 2.19.1 and cypari2, not CI's PARI 2.17.4.

## Bottom line

- **No survivor.** 340,000 random coprime ratios were drawn, 20,000 in each of the 17 dyadic ranges from m = 10⁴ to 10⁹. No non-degenerate ratio keeps a four-class tuple. *(computed, exact local tests; tier D as in the repository)*
- **About 1% of the ratios have local points everywhere** (3,133 of 340,000). Of those, **2.0% keep a class under the criterion** (63). The descent empties all 63 with good primes up to 150. The share kept doesn't grow with m: it runs 0% to 4.2% per range, with no trend.
- **One control fails, and the failure is real.** For the twice-square ratio (1109503,597671), the four-class descent leaves 6 tuples. It should leave only the 2 of its degenerate points (t = 0 and t = ∞). The 4 extra tuples survive:
  - every bad prime, tested exactly (the large ones too);
  - every good prime up to 20,000;
  - a 2-adic point found explicitly for each of them.
  This doesn't break (R), which is about non-degenerate ratios. It does show the descent can leave classes that no point explains. So a non-degenerate survivor is possible in principle, even though none turned up.

## Per range (g1_summary.txt)

| k | range of m | no local point | points everywhere | class kept | survivors | controls exact |
|---|---|---|---|---|---|---|
| 13 | [10000, 16384) | 19804 | 195 | 4 (2.05%) | none | 200/200 |
| 14 | [16384, 32768) | 19807 | 191 | 3 (1.57%) | none | 200/200 |
| 15 | [32768, 65536) | 19832 | 168 | 0 | none | 200/200 |
| 16 | [65536, 131072) | 19836 | 164 | 1 (0.61%) | none | 200/200 |
| 17 | [2¹⁷, 2¹⁸) | 19830 | 170 | 2 (1.18%) | none | 200/200 |
| 18 | [2¹⁸, 2¹⁹) | 19831 | 169 | 3 (1.78%) | none | 200/200 |
| 19 | [2¹⁹, 2²⁰) | 19818 | 182 | 4 (2.20%) | none | 200/200 |
| 20 | [2²⁰, 2²¹) | 19799 | 201 | 4 (1.99%) | none | **199/200** |
| 21 | [2²¹, 2²²) | 19806 | 194 | 3 (1.55%) | none | 200/200 |
| 22 | [2²², 2²³) | 19828 | 172 | 6 (3.49%) | none | 200/200 |
| 23 | [2²³, 2²⁴) | 19833 | 167 | 7 (4.19%) | none | 200/200 |
| 24 | [2²⁴, 2²⁵) | 19813 | 187 | 4 (2.14%) | none | 200/200 |
| 25 | [2²⁵, 2²⁶) | 19786 | 214 | 5 (2.34%) | none | 200/200 |
| 26 | [2²⁶, 2²⁷) | 19804 | 196 | 5 (2.55%) | none | 200/200 |
| 27 | [2²⁷, 2²⁸) | 19813 | 187 | 5 (2.67%) | none | 200/200 |
| 28 | [2²⁸, 2²⁹) | 19802 | 198 | 1 (0.51%) | none | 200/200 |
| 29 | [2²⁹, 10⁹] | 19822 | 178 | 6 (3.37%) | none | 200/200 |

Three of the random draws were degenerate (k = 13 and 14) and were set aside. Each ratio that keeps a class is listed in its g1_k<k>.txt with its class dimension. "Nothing is left" holds for every one of them at good primes ≤ 150.

## What holds and how it was checked

- **Pipeline (g1_stress.py), for each ratio:**
  - The repository's class criterion, unchanged: `class_selmer.empty`, then `class_selmer_linalg.system` / `solve`. The only change is that its factorisations come from PARI, with every factor proved prime (asserted). That covers 2mn(m²−n²)(m²+n²) and the octic.
  - For a ratio that keeps a class: the four-class descent of tuple_descent.py, with its local test replaced by the exact one of exact_local.py at odd primes. At p = 2 the repository's test (depth 12) is kept.
  - Primes used: every prime of the supports plus the good primes up to 150. While anything beyond the known tuples is left, it goes on to 1500, 5000 and then 20,000.
- **exact_local.py** decides the local condition at an odd prime exactly. It doesn't cut the line into p discs, so primes up to 10⁹ and beyond (the octic's) are handled. Each form is linear on a disc:
  - where it has no root mod p its class is constant;
  - where it has one, the generic residues give sign patterns of characters, settled by Weil's bound for p > 5·10⁶ and by listing them below that;
  - the root residues are refined;
  - common roots are handled exactly in ℤ[i].
- **Checks of exact_local.py** (test_exact_local.py, committed output test_exact_local.txt):
  - on 118 (ratio, prime) cases it gives the same set as the repository's test at increasing depth, never a smaller one;
  - 240,000 random rational points (61,856 alive) all land inside its sets;
  - descents agree with the repository's on 12 ratios, and degenerate ratios keep their points.
- **Controls:** 200 degenerate ratios per range (square and twice-square type, from Pythagorean parameters, seed 2000+k). 3,399 of 3,400 keep exactly their degenerate tuples.
  - (3728456,649617) and (984232377,743811136) needed good primes to 5000. The first run stopped at 1500 and showed them as failures. The driver now goes on to 5000 and 20,000, and k = 20, 21 and 29 were rerun with it.
  - The other ranges were made before that change. k = 13 and 15 rerun with the final code match byte for byte. The change can only alter ratios with tuples left at 1500, and those ranges have none.
- **The failing control, (1109503,597671):**
  - m² + n² = 2·(5·7129)², so it is twice-square type.
  - Known tuples: (1109503, 597671, 2, 1) at t = 0 and (1195342, 2219006, 2, 1) at t = ∞.
  - Extra tuples: each of those with (2, 959) or (1918, 1) in the last two places, where 959 = 7·137. So 6 tuples in all.
  - Still 6 after all bad primes (exact) and good primes to 20,000.
  - Each of the 6 has an explicit 2-adic point, (u, v) = (t, 1) or (1, 2t) with t < 2¹⁶ (g2_badalone_witness.witnessed_at_2).
  - The repository's own test at the odd bad primes up to 1871 keeps all 6.
  - tuple_bruteforce.py (primes ≤ 61 only) keeps 288 ⊇ these, as it must. It can't test the support primes 593, 1871 and 597671, so it isn't an independent check here.
- **Negative mode** (`--negative`, g1_negative.txt): 12 of 12 refused as they must be.
  - 9 degenerate ratios from three ranges, run with the degeneracy check off, are all reported as SURVIVORS.
  - 3 degenerate ratios whose known points are withheld still leave those tuples.
  So the survivor test isn't vacuous.

## What failed

- The first full run stopped at good primes ≤ 1500 and so reported 2 false control failures. Fixed as above.
- An O(p) search for √−1 inside the repository's `Local.__init__` made exact_local.py unusable at primes above about 10⁸. exact_local.py now sets those fields itself and takes √−1 from PARI. This doesn't change any output, and G1's ranges were not affected: they only reach such primes through supports, which are small enough.
- The extra tuples of (1109503,597671) are not explained. I haven't checked whether they come from a Selmer element that no rational point explains, or from points on the cover. *(open)*

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** Stress test of (R) (experiments/roadmap-g/g1): 340,000 random coprime ratios, 20,000 per dyadic range from m = 10⁴ to 10⁹, factorised by PARI with every factor proved prime. 3,133 have local points everywhere; 63 of those keep a class under the class criterion; the four-class descent, with an exact local test at odd primes (exact_local.py, no prime left out), empties all 63 with good primes ≤ 150. No survivor. Of 3,400 degenerate controls, 3,399 keep exactly the tuples of their degenerate points. The twice-square ratio (1109503,597671) keeps 4 more, which survive every bad prime and every good prime ≤ 20,000 and each have a 2-adic point: the four-class descent does not isolate the degenerate classes in general. Computed (exact), tier D; nothing is proved for all ratios. Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-g/g1        # F = grep -v '^PARI stack size'
python3 test_exact_local.py 2>&1 | $F | tail -1 | grep -q '^# matches'      # about 4 s
python3 g1_stress.py 13 20000 200 2>&1 | $F | grep -q '^# matches'           # about 17 s
python3 g1_stress.py 20 20000 200 2>&1 | $F | grep -q '^# matches'           # about 130 s (holds the failing control)
python3 g1_stress.py --negative 2>&1 | $F | grep -q '^# matches'             # about 3 s
python3 g1_summary.py | diff - g1_summary.txt                                 # under 1 s
# not for CI: bash run_g1.sh 13 ... 29   (all ranges; about 1.5 h on two cores; times.txt)
```

Timings (one core, shared): k = 13 to 19 take 15 to 40 s each, k = 22 takes 97 s, k = 25 takes 340 s, k = 28 takes 1036 s and k = 29 takes 891 s.

exact_local.py finds the repository through MSS_REPO, or else `../../..`.
