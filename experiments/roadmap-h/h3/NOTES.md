# H3: local images of the four-class cover at good primes (Roadmap H, run 7)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (z8, `reviews/review_z8.md`, its
own checks in `../../referee_z8/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.** Nothing
here changes the census, which stops at `m ≤ 1000`. `REPORT_H3.md` is the session's report, unchanged; the notes below
take precedence over it.

## What holds: the tables

For each of the **91** primes `p ≡ 1 (mod 4)` with `17 ≤ p ≤ 1153` (the report says 96) and each good residue
`r = m/n mod p`, the set of label vectors allowed at `p`, computed two ways (G1's exact local test, method A, and an
enumeration of `ℙ¹(F_p)`, method B) and checked against the exact test on 2,000 random ratios at primes up to 197.
Every one of these primes removes a vector at some residue (17,224 of 49,764 good residues); a set is empty only up to
`p = 137` (36 empty sets); 15,984 of the sets are not affine subspaces of `F₂⁴`; the other counts of the report's
"shape" paragraph hold (review z8 recomputed them from `h3_weil.json`, and its own local test equals the tables at 15
primes from 17 to 1153). Computed (exact arithmetic; method A is G1's `exact_local.py`, tested in run 6, not proved; method B is an enumeration of ℙ¹(F_p), lifted by Hensel's lemma, quoted); tier D.

## What is withdrawn: the Hasse–Weil step

The report's first three bottom-line bullets, its "So the full local four-class set of any ratio is a finite
computation …", the "What failed: Nothing", the STATEMENT / PROOF in the docstring of `h3_weil.py`, and line 2 of `h3_weil.txt` ("every good prime p >= 1154 keeps all 16 vectors") are **wrong**:

1. **Wrong curve.** The label of a point at a good `p ≡ 1 (mod 4)` is the class of the linear factor `l_j⁺` of
   `F_j = l_j⁺ l_j⁻`, at points where every `F_j` is a non-zero square. The twist `y_j² = u^{e_j} F_j(t)` with
   `e_j = 1` has `F_j(t)` a non-square, which is not an allowed point at all. The curve that carries a label vector is
   the `(ℤ/2)⁷`-cover of `ℙ¹` given by `√(l_j⁺ l_j⁻)` (four) and `√(u^{e_j + e_0} l_j⁺ l_0⁺)` (three), branched at 8
   points, of **genus 129**, not 17.
2. **Good primes above 1153 do remove vectors.** Review z8 found losses at every good `p ≡ 1 (mod 4)` it tried from
   1181 to 1249 (for example the ratio (5,1) at `p = 1181` loses the vectors 5 and 10, also by a brute force over
   `ℙ¹(ℤ/p²)` with the integer forms), and in full scans up to 5000 the last loss is at 3853
   (`../../referee_z8/`). Three methods agree (z8's own enumeration, G1's exact test, the brute force).
3. **Good primes `p ≡ 3 (mod 4)` can empty a set.** There the set is everything or nothing (G2: every candidate has the
   same label), everything exactly when some `t ∈ ℙ¹(F_p)` makes all four `F_j(t)` non-zero squares: always for
   `p ≡ 3 (mod 8)` (`t = 0` or `±1`), and for `p ≥ 1154` (Hasse–Weil on the genus-17 curve `y_j² = F_j(t)`, which is
   the right curve there). For some residues it is nothing at `p = 7, 23, 31, 47, 79, 103, 127, 151`, and never at
   `167 ≤ p ≤ 1151` (z8's enumeration; G1's exact test confirms empty sets, for example (2,1) at 7 and (14,1) at 151).

**A correct finite local test** (review z8's argument; Hasse–Weil quoted; not proved in Lean or in a script that
checks the argument): the bad primes, ℝ, the good primes `p ≡ 7 (mod 8)` below 1154 and the good primes `p ≡ 1 (mod 4)`
below 67,583. At a good `p ≡ 1 (mod 4)` with `p ≥ 67,583`, Hasse–Weil on the genus-129 curve, with the `8·64` points
over the branch points removed, gives a point with every label vector (`(p − 511)² > 66,564 p` in exact integers). The
largest good prime that removes a vector is not known; z8's full scans go to 5000.

## Rerunning

```
python3 h3_tables.py              # 17 s; --negative 1 s; --check 2000 5 s (primes <= 197 only)
python3 h3_weil.py --negative     # 1 s (it still tests the withdrawn threshold 1154: read it as a test of the arithmetic only)
python3 h3_sizes.py               # 1 s
python3 h3_weil.py                # 320 s, not in CI (method A on every residue)
```
Under PARI 2.17.4 (CI's version) every output above matches; the session used PARI 2.19.1.
