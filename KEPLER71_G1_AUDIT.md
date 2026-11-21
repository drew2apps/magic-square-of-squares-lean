# Audit of the G1 closures in Kepler71/magic-square-of-squares

[Kepler71/magic-square-of-squares](https://github.com/Kepler71/magic-square-of-squares) studies a two-parameter family of 3×3 magic squares called G1. For coprime `1 ≤ n < m ≤ 20` (127 pairs) it reports that 120 pairs cannot give a magic square of nine distinct squares: 62 because an elliptic curve has rank 0, and 58 because of congruences. Its handoff note says these 120 closures were inherited without a consolidated audit.

This is an independent check of those 120 closures. It was made from the statements in `RESULTS_FOR_CLAUDE_2026-09-12_G1.md` and `SUMMARY_2026-09-12.md` (read on 21 November 2026). The repository's code and logs were not read or run.

## The claim that is checked

The three diagonal cells of the family are

```
F0 = m² + n²t²,   F4 = (m² + n²)(t² + 1)/2,   F8 = n² + m²t².
```

A magic square of nine squares in the family needs a rational `t` for which all three are squares. With `t = U/V` in lowest terms this needs the integers `m²V² + n²U²`, `2(m² + n²)(U² + V²)`, `n²V² + m²U²` to be squares.

## Results

**Reproduced.** Every number below matches the Kepler71 report.

| Item | Kepler71 | This audit |
|---|---|---|
| Torsion of `E_B: Y² = X(X − (m²−n²)²)(X − (m²+n²)²)` | `Z/2 × Z/4` for all 127 pairs | same |
| Rank bounds from PARI `ellrank` | 62 × [0,0], 60 × [1,1], 4 × [2,2], 1 × [0,2] | same |
| Pairs without rank 0 that are excluded by a congruence | 58 | 58 |
| First obstruction among those 58, in the order 32, 25, 7, 23 | 32, 7, 15, 4 | 32, 7, 15, 4 |
| Pairs left open | (11,4), (13,8), (15,1), (15,8), (16,5), (19,5), (19,16) | the same seven |
| The 2-adic criterion (lemma 1) | stated for all coprime `m, n` | agrees with the scan on all 127 pairs |

The ranks were recomputed with the same software (PARI 2.17.2), so that row is a reproduction, not an independent proof. The congruence scan is independent code.

**New.**

1. **56 of the 62 rank-0 pairs are also excluded by a congruence.** So 114 of the 120 closures need no rank computation and no computer algebra system. Only six closures depend on the rank being 0: (4,3), (7,1), (12,5), (17,7), (20,1), (20,7).
2. This covers (13,3) and (19,11). For these two the report says that eclib gave only the bound [0,2] and that the closure rested on PARI alone. Both are excluded modulo 7.
3. The pair (12,7) is excluded modulo 47, and by no smaller modulus that was tested. The report scanned odd primes up to 31.
4. **All 114 congruence exclusions are checked by the Lean kernel.** `MagicSquareOfSquares/G1.lean` proves that a finite check modulo `M` rules out every rational `t`. `MagicSquareOfSquares/G1Certificates.lean` lists the 114 triples `(m, n, M)` and the kernel evaluates each check. The theorems depend only on `propext`, `Classical.choice` and `Quot.sound`.

The moduli used are 7 (68 pairs), 8 (23), 23 (8), 25 (9), 32 (5) and 47 (1).

## Limits

- This audit covers only the family G1 with `m, n ≤ 20`. It says nothing about the open problem.
- The formulas for `F0`, `F4`, `F8` are taken from the Kepler71 notes. The family itself, and the fact that it is magic, are not formalized here.
- The 13 pairs without a congruence exclusion were tested for all primes below 100, the primes dividing `2mn(m² − n²)(m² + n²)`, and some prime powers. No obstruction was found. That is not a proof of local solubility.
- The seven open pairs were not examined.

## Reproduce

```
python3 scripts/g1_audit.py            # scan; recomputes ranks if cypari2 is installed
lake build MagicSquareOfSquares.G1Certificates
```
