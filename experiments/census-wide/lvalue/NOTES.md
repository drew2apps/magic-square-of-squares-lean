# lvalue: L(E,1) in exact interval arithmetic for (6557,324) and a rerun for (700,529) (run 5, F5)

**Experiment, not counted; tier D of `LEDGER.md`.** Scripts and outputs written by an AI model (Claude, the "Fable"
session of run 5); reviewed by a separate AI agent (`reviews/review_x4.md`). Not read by a person. Nothing here changes
the census or `census-wide/table/` (which keeps the numerical label for (6557,324), (6881,6233)).

**(6557,324).** For E = w·f0·f8 of (6557,324) (model [0, 1848516467034577, 0, 20370457151276745574584576, 0],
N = 11077451294972790 by Tate's algorithm), **L(E,1) ∈ [1.444085955623, 1.444085956081], conditional on PARI's
`ellap` at the good primes ≤ 1.21·10⁹** (`r6557/record.json`). The arithmetic is the repository's `lball_exact.py`
(192-bit fixed point with directed rounding, the exact tail bound from |aₙ| ≤ 2n, the guard: the enclosures at y = 1
and y = 6/5 meet and the Fricke residual contains 0); `lvalue_stream.py` restates `lball_exact.series` so that the
1.21·10⁹ terms can be streamed in segments, with aₙ from a numpy sieve. w = +1 because the Fricke sum F lies in
[−7802755.744, −7802755.743]. There is no floating point in the certificate: floats only choose the number of terms.
With L(E,1) ≠ 0, rank E(ℚ) = 0 (Kolyvagin; modularity), and the census-wide closure of (6557,324) and, through the
partner symmetry, (6881,6233) would then rest on PARI only through `ellap`: the review made the torsion step exact
(`../../referee_x4/tors6557.py`: E(ℚ)_tors = ℤ/2 × ℤ/8, torsion t = ±1, 0, ∞ only, `global.py`'s `analyse` closes).
This is a better footing than the numerical value 1.44409 of `../hard7k/` (same a_p source, but exact arithmetic and
tail), not an unconditional one.

What the condition and the checks are:
- **Condition:** the good aₚ are PARI's `ellap` on PARI's minimal model (Hasse asserted for each; PARI's conductor
  checked against Tate's). The repository's own certificate `lball_exact53.ap_one` agrees at 20,000 random primes
  (`r6557/apcheck.txt`) and, in the review, at 30,600 more (seed 6, including the 300 largest primes ≤ Mmax;
  `../../referee_x4/apsample.txt`). Removing the condition means `ap_one` at all 60,926,519 primes (about 24
  core-hours); not done.
- **One run:** the six exact sums (`r6557/series_A/B/C.json`) come from one run of `r6557/run.sh` on the session's
  machine (about 2 h on 4 cores; its prime and a_p arrays, about 465 MB, are not committed). `--merge` recomputes the record
  from these sums and `--negative` forges w = −1 on them; neither recomputes the sums. Not rerun in review.
- **Positive control:** `stream_selftest.py` — on the repository's curve (`lball_exact.CURVE`, N = 1929934048) the
  streamed aₙ equal `lball_exact.an_from_ap` and the six streamed sums equal `lball_exact.series` exactly, for one
  segment and for segments of 4093 (about 74 boundaries inside the sums).
- **Forged conductor:** `MSS_STREAM_TEST=1 lvalue_stream.py --negative-small` — N/2 and 2N fail the guard on that
  curve, the true N passes (the repository's forged-N controls, through the streamed code).
- **64-bit arithmetic** (review): n ≤ 1.21·10⁹ < 2³¹; every partial product in the sieve is a_m for a divisor m of n,
  so |·| ≤ d(n)√n < 2²⁶ given Hasse; overflow cannot occur and the |aₙ| < 2⁵⁰ assertion is only a sanity check. The
  sums themselves are Python integers.
- `r6557/primes.json`: 60,926,519 primes up to Mmax = 1,209,869,493 (PARI's `primepi` agrees, review).

**(700,529).** `lvalue_big.py` is the repository's `scripts/g1_full/census_lvalue_big.py` with the pair and curve read
from `MSS_PAIR`, `MSS_CURVE` (default (700,529), w·f0·f8). Its rerun on the session's machine (`r700_529/`: three a_p
chunks, merge, PARI check) gives `r700_529/p.jsonl`, equal to the committed `census_lvalue_big_700_529.jsonl` in every
field except the timings (checked in review). Not rerun in review (about 1.5 core-hours).

Files: `lvalue_stream.py` (stages `--primes`, `--ap`, `--series` write into the folder given; `--merge`, `--negative`,
`--negative-small`, `--apcheck` write only with `--write`), `lvalue_big.py`, `stream_selftest.py` (+ `selftest.txt`),
`test_lball/negative_small.txt`; `r6557/`: `run.sh`, `record.json` (= `merge.out`), `series_A/B/C.json`, `negative.txt`,
`apcheck.txt`, `primes.json`, `ap_4_*.json`, `times.txt`; `r700_529/`: `p.jsonl`, `merge.out`, `bigcheck.out`,
`p*.bin.json` (chunk statistics), `times.txt`.
