# H3: local images of the four-class cover at good primes

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted. Nothing here changes the census, which stops at m ≤ 1000.

Machine: PARI 2.19.1 with cypari2, not CI's PARI 2.17.4. PARI is used only for √−1 mod p inside G1's exact_local.py. Everything else is exact integer arithmetic.

## Bottom line

- **Good primes above 1153 never remove anything.** At a good prime p ≥ 1154, the four-class set allows all 16 label vectors, for every ratio. *(proved, given the Hasse–Weil bound, which is quoted; the argument is in h3_weil.py)*
  - Each twist y_j² = u^{e_j}F_j(t), j = 0, 8, L, U, is a smooth (ℤ/2)⁴-cover of ℙ¹ branched at 8 points, of genus 17, with good reduction at a good prime.
  - Hasse–Weil gives #C(F_p) ≥ p + 1 − 34√p, which is > 0 for p ≥ 1154. The threshold is checked as (p + 1)² > 1156p, in exact integers.
  - An F_p-point lifts by Hensel, so every vector is allowed.
- **Good primes p ≡ 3 mod 4 never remove anything either.** G2 showed that every divisor has the same label there.
- **So the full local four-class set of any ratio is a finite computation:** its bad primes, ℝ, and the 96 good primes p ≡ 1 mod 4 with 17 ≤ p ≤ 1153. Before this, a ratio left non-empty by good primes up to some cutoff (G1's control, to 20,000) was not known to be non-empty everywhere. Now that is settled. *(proved, given Hasse–Weil)*
- **Every prime in that range matters for some residue.** Losses occur at every p ≡ 1 mod 4 from 17 to 1153, the last one included (204 of its 1,140 good residues). In all, 17,224 of 49,764 good residues lose a vector. *(computed, exact)*
- **The shape of the sets (h3_sizes.py, h3_sizes.txt).** *(computed, exact, from h3_weil.json)*
  - Every set is closed under adding (1,1,1,1), the class of a non-square unit, so every size is even.
  - A good prime empties the set (no local point) only up to p = 137. Sets of at most 2 vectors occur up to 313, of at most 8 up to 617, of at most 10 up to 1069. Above 1069 a loss is one or two complementary pairs.
  - 15,984 of the 49,764 sets are not affine subspaces of F₂⁴. So the good-prime test is not linear algebra over F₂, and the four-class set need not be an affine space. A proof built on these tables has to handle that.
- **The mean size of the allowed set rises** from 2.7 of 16 at p = 17, to about 10 at p ≈ 190, about 14 at p ≈ 500, and 15.6 at p ≈ 1150. That is why good primes ≤ 181 have been enough in every case so far: the small ones remove most.

## The tables

- **h3_tables.py (p ≤ 197):** two methods must agree on every good residue r = m/n mod p.
  - **A:** G1's exact local test on tuple_descent.magic(r, 1).
  - **B:** an enumeration of ℙ¹(F_p) with the forms split over F_p, the root points handled by Hensel.
- **Results to 197 (h3_tables.txt):**
  - 5 and 13 have no good residue (always bad);
  - 36 of 1,848 good residues have an empty set (no local point at p);
  - full sets occur only at 137 (8 residues), 181 (16) and 193 (16).
- **h3_weil.py (p ≤ 1153):** method B in integer numpy arrays (residues mod p, no floating point).
  - It equals h3_tables.json at every p ≤ 197.
  - Method A agrees on all 47,916 good residues at 197 < p < 1154.
  - h3_weil.txt lists, per prime, the residues that lose something and the mean set size. h3_weil.json holds every set.

## Checks

- **h3_tables.py --check 2000 (h3_check.txt):** 2,000 random coprime (m, n) up to 10¹², each with a random good prime of the list. G1's exact test on the real forms gives the table's set at r = m/n, 0 disagree. This checks the reduction to the residue r.
- **h3_tables.py --negative (h3_negative.txt):** 3 forged entries refused.
- **h3_weil.py --negative (h3_weil_negative.txt):**
  - p = 1153 fails the threshold test, as it must;
  - forged sets at p = 229 and 409 (a lost vector put back) are refused by the exact test;
  - a bad residue is refused by the root assertion.

## What it rests on

- **Quoted, not proved here:** the Hasse–Weil bound; smoothness of the cover's reduction at a good prime (p odd, branch points distinct mod p, asserted in the code); Hensel's lemma.
- **Computed:** everything else, exactly.

## What failed

Nothing. The threshold is not sharp: losses still occur at 1153.

## Proposed LEDGER / CHANGES text

> **Experiment, not counted.** Local images of the four-class cover at good primes (experiments/roadmap-h/h3). At a good prime p ≥ 1154 every label vector is allowed: the twists are smooth genus-17 curves and Hasse–Weil gives them points. At good p ≡ 3 mod 4 nothing is lost either. So the four-class set of any ratio is determined by its bad primes, ℝ and the 96 good primes p ≡ 1 mod 4 up to 1153. Tables for those primes, computed two ways and checked on 2,000 random ratios, show losses at every one of them. Proved given Hasse–Weil (quoted); tables computed (exact). Written by an AI model; not read by a person.

## Proposed CI lines

```
cd experiments/roadmap-h/h3        # F = grep -v '^PARI stack size'; needs ../../roadmap-g/g1 (or G1_DIR)
python3 h3_tables.py 2>&1 | $F | tail -1 | grep -q '^# matches'              # about 17 s
python3 h3_tables.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'   # about 3 s
python3 h3_tables.py --check 2000 2>&1 | $F | tail -1 | grep -q '^# matches' # about 17 s
python3 h3_weil.py --negative 2>&1 | $F | tail -1 | grep -q '^# matches'     # about 2 s
python3 h3_sizes.py 2>&1 | $F | tail -1 | grep -q '^# matches'              # about 1 s
# not for CI: python3 h3_weil.py   (320 s, method A on every residue; times.txt)
```
