# H1: the phantoms of the four-class descent (Roadmap H, run 7)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (z8, `reviews/review_z8.md`, its
own checks in `../../referee_z8/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.** Nothing
here changes the census, which stops at `m ≤ 1000`. `REPORT_H1.md` is the session's report, unchanged; the notes below
take precedence over it.

## What it claims (as corrected)

- **The scan skipped primes.** The scan of 20,009 random degenerate ratios (`10⁴ ≤ m ≤ 10⁹`) ran G1's `one()`, which
  tests the primes of the supports (`2mn`, `2(m² − n²)`) and the primes up to 20,000 only, not the primes of `m² + n²` or
  of the octic above 20,000 (run 6's correction 2 to G1, `../../roadmap-g/g1/NOTES.md`). The report's "every bad prime"
  and "their extra tuples pass every local test" are therefore wrong as descriptions. The counts in the files are right
  for that test: 13 ratios keep tuples beyond those of `t = 0, ±1, ∞`; with (1109503,597671) and the partners, 12 pairs.
- **With every bad prime, 7 pairs are phantoms, not 12.** Review z8 tested the primes the scan left out: at an octic
  prime, every extra tuple is removed for 5 of the 12 pairs, by G1's `exact_local.py` and by a short exact argument of
  its own (two linear components of different forms meet modulo `p` with valuation 1, which forces
  `e_{j_a} + e_{j_b} = χ(κ)` on every allowed point). These 5 are not phantoms: (11875080,6497849)/(18372929,5377231)
  (the report's "open pair", removed at 5462461 and at 798775110477654781), (32396425,3694248)/(36090673,28702177),
  (175775580,15972419)/(191747999,159803161), (234319635,115640732)/(349960367,118678903),
  (518804652,31743595)/(550548247,487061057). The bundle's own `exact_local.descend_exact` with every bad prime gives the
  same: 14 ratios keep exactly the report's extra tuples, the other 10 keep only their known tuples
  (`../../referee_z8/full_descent.txt`).
- **The 7 true pairs** are (204724,97893)/(302617,106831), (413004,148603)/(561607,264401),
  (756492,255635)/(1012127,500857), (853587,255916)/(1109503,597671), (1992873,999664)/(2992537,993209),
  (189089400,156063271)/(345152671,33026129), (231641256,72371617)/(304012873,159269639). Their extra tuples are locally
  soluble everywhere: review z8 gave explicit local points at 2, at every odd bad prime (support primes, primes of
  `m² + n²` and of the octic) and at every good prime `p ≡ 1 (mod 4)` below 1200 and between 20,000 and 67,583 (its own enumeration; between 1200 and 20,000 G1's `exact_local.py`, as in the scan, and an unpublished check by the pre-publication referee agrees), and the
  genus-129 Hasse–Weil argument of `../h3/NOTES.md` covers the rest (quoted); at good `p ≡ 3 (mod 4)` these ratios keep
  everything, since they have rational points.
- **All 7 true pairs are settled:** every extra tuple is the class of no rational point of the four-form curve C. Five
  by an exactly computed rank-0 quotient `y² = w·f0·f8` (complete 2-descent and Cassels–Tate pairing over curves of its
  isogeny class; torsion by the gcd of `#E(F_p)`; the torsion points lie over `t ∈ {−1, 0, 1, ∞}` only): 413004/561607,
  853587/1109503, 1992873/2992537, 189089400/345152671, 231641256/304012873. Two by the pair curves and the
  "Cassels–Tate kernel test" (the image of `E_jk(ℚ)` lies in the kernel of the pairing on `Sel²`; a class outside that
  kernel is the class of no point): 204724/302617 at 302617, 756492/1012127 at 1012127. **No pair is open.** Review z8:
  the map `C → A`, the torsion fibres and the kernel test hold (the test compares classes modulo −1, which is weaker,
  hence sound); PARI 2.17.4 agrees on rank 0 and the torsion; the pairing rests on `rank_second.py`
  (`reviews/review_t2.md`).
- **The partner isomorphism** `t′ = (t + 1)/(t − 1)`, with `(m′, n′) = ((m + n)/g, (m − n)/g)` and the forms permuted
  `(0, 8, L, U) → (L, U, 0, 8)` up to the square constants −1, −4 (halves when `g = 2`), holds: derived by hand and checked
  as a polynomial identity by review z8; `h1_iso.py` checks it for the 12 pairs. An exclusion on either side settles a pair.
- Computed (exact), resting on the soundness of the complete 2-descent and of the Cassels–Tate pairing as `rank_second.py` computes them (Cassels 1998; Skorobogatov 6.2), Vélu's formula and isogeny invariance of the rank, the injectivity of torsion into `E(F_p)` at good odd `p` (quoted, tier E); the local tests on G1's `exact_local.py` (tested, not proved) and the repository's one-sided test at 2. Tier D.
- Computed (exact), resting on the soundness of the complete 2-descent and of the Cassels–Tate pairing as `rank_second.py` computes them (Cassels 1998; Skorobogatov 6.2), Vélu's formula and isogeny invariance of the rank, the injectivity of torsion into `E(F_p)` at good odd `p` (quoted, tier E); the local tests on G1's `exact_local.py` (tested, not proved) and the repository's one-sided test at 2. Tier D.
- Computed (exact), resting on the soundness of the complete 2-descent and of the Cassels–Tate pairing as `rank_second.py` computes them (Cassels 1998; Skorobogatov 6.2), Vélu's formula and isogeny invariance of the rank, the injectivity of torsion into `E(F_p)` at good odd `p` (quoted, tier E); the local tests on G1's `exact_local.py` (tested, not proved) and the repository's one-sided test at 2. Tier D.
- (1109503,597671) is settled twice on that ratio itself (pair curves `E_0L`, `E_0U` of rank exactly 3; `y² = w·f0·f8`
  of rank 0), as the report says.

## Corrections to `REPORT_H1.md` and the hand-off
- Bottom line: read "13 ratios keep extra tuples under G1's test; with every bad prime 7 of them, 7 partner pairs with
  (1109503,597671), all settled (5 by a rank-0 quotient, 2 by the pair curves and the Cassels–Tate pairing); the other 5
  pairs are removed by a local test at an octic prime. No pair is open."
- "Their extra tuples are local everywhere … By H3, good primes above 1153 never cut anything": withdrawn (`../h3/NOTES.md`);
  for the 7 true pairs see above.
- The "Next step" (a Mordell–Weil sieve for the open pair) and the scratch ideas aimed at that pair are moot: the pair
  is not a phantom.
- The docstrings of `h1_scan.py` and `h1_partners.py` ('the exact local test at every odd prime'; 'Since H3, good primes to 1153 are all that can cut, so good primes to 20000 is the full local test') are wrong in the same way: G1's `one()` skips the primes of `m² + n²` and of the octic above 20,000, and H3's bound is withdrawn.
- The scan sampled degenerate ratios only; REPORT_H1's "phantoms occur only at degenerate ratios" says nothing about the others (for those see H2 and G1: no survivor).
- The docstrings of `h1_scan.py` and `h1_partners.py` ('the exact local test at every odd prime'; 'Since H3, good primes to 1153 are all that can cut, so good primes to 20000 is the full local test') are wrong in the same way: G1's `one()` skips the primes of `m² + n²` and of the octic above 20,000, and H3's bound is withdrawn.
- The scan sampled degenerate ratios only; REPORT_H1's "phantoms occur only at degenerate ratios" says nothing about the others (for those see H2 and G1: no survivor).
- The docstrings of `h1_scan.py` and `h1_partners.py` ('the exact local test at every odd prime'; 'Since H3, good primes to 1153 are all that can cut, so good primes to 20000 is the full local test') are wrong in the same way: G1's `one()` skips the primes of `m² + n²` and of the octic above 20,000, and H3's bound is withdrawn.
- The scan sampled degenerate ratios only; REPORT_H1's "phantoms occur only at degenerate ratios" says nothing about the others (for those see H2 and G1: no survivor).
- A cheap fix (z8: 5 to 16 s per ratio) is to make `h1_scan.py` test every bad prime, as G2 does, and rerun; not done.

## Rerunning
```
python3 h1_pairs.py --negative          # 45 s under PARI 2.17.4
python3 h1_scan.py --negative           # 9 s
python3 h1_rank0_exact.py --negative    # 2 s
python3 h1_iso.py                       # 84 s
python3 h1_iso.py --negative            # 1 s
```
Each prints `# matches` (CI's PARI 2.17.4; the session used 2.19.1). The long runs (the scan, 200 to 3000 s per range;
the pair-curve and rank-0 runs, 1000 to 2200 s) are not in CI; review z8 reran single ratios of each and they match.
