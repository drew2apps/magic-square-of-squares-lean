# REPORT_I2 (27 Nov 2026, about 23:45 UTC): what the bad primes leave

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review z17) holds.


Written by an AI model; not read by a person. Experiment, not counted. Nothing here changes the census. All 3,876 ratios that need a good prime are done. The shape run finished 3,838 of them; 38 were over its 5 GB memory cap. The residual file covers all 3,876.

## Bottom line

- **The set S left by the bad primes always has even size.** Every |S| is even, from 2 to 834, in all 3,876 ratios (residual file). In the shape run's 3,838: 1,547 have |S| = 2, 893 have 4, 393 have 8 and 266 have 6. *(computed, exact)*
- **70% are cosets.** 2,693 of the 3,838 sets are a coset of a subgroup of F₂ⁿ (tuples as vectors over the support primes), and 1,145 are not.
- **The killing good prime is small.** P, the first good prime that empties S, is 29 for 1,592 of the 3,876 ratios, 37 for 713, 41 for 491, 53 for 260 and 61 for 250 (residual file). The largest is 229, as in H2. *(computed, exact)*
- **Killing sets are short and lean on 2, 5 and 13.** An irredundant killing set (bad primes plus good primes up to P) has 3 to 5 primes for 2,557 of the 3,838. The primes in the most sets are 2 (3,115), 5 (2,482), 13 (1,935), 29 (1,905), 3 (1,471) and 17 (1,203).
- **At the killing prime the smallest separating label quotient never has 1 coordinate.** It has 2 coordinates for 2,477 ratios, 3 for 474 and 4 for 887: the quotient that kills at P always involves at least two of the four classes.
- **G2's 937 are all in the list.** My earlier "928" counted 9 of G2's 836 under H2's hunts, where they also occur. The list is 2,807 (H2, exhaustive) + 132 (H2's hunts only) + 836 + 101 (G2) = 3,876.

## The residual-tuple file (for the conductor)

`i2_residual.jsonl` has one JSON object per line for all 3,876 ratios of `i2_needgood_list.txt`, in its order. Its fields and their order are fixed, since the research session will use it as a test set:

`{"m": m, "n": n, "residual": [[d0, d8, dL, dU], ...], "kill": [p, ...], "killing_prime": P}`

- `residual` is S: the tuples allowed at every bad prime (the repository's exact local test, every bad prime; ℝ imposes nothing), sorted.
- `kill[j]` is the first good prime ≤ 5000 at which `residual[j]` is not allowed.
- `killing_prime` is the largest of them.
- There are 32,474 residual tuples. The most for one ratio is 834, every tuple is killed by a good prime at most 229, and no record has a status line or a tuple that survives every good prime up to 5000.
- **Cross-checked:** for the 3,838 ratios done by both runs, |S| and P agree with `i2_shape.py`, and I8 recomputed 50 records (seed 1) with conductor-t10's own `residual()`: 50 of 50 agree.
- **Negative controls** (`i2_residual.py --negative`, 3 s): a forged residual tuple, a forged smaller killing prime and a record with a tuple left out are all refused, and the true record of (15001,2081) is accepted. *(computed, exact)*

## What the columns of the shape run mean

`i2_shape.py` prints one line per ratio: `m,n |S| log2(space) coset|not-coset P kill=[p:k …] [S=…]`. P is the first good prime ≤ 5000 that empties S. In kill=[p:k], p runs over an irredundant killing set and k is the smallest number of label coordinates whose projection already separates; S is printed when |S| ≤ 8.

## What failed or is still open

- 38 ratios of the shape run hit its 5 GB cap (`MEMORY`). The residual file has them, so S and P exist for each; only their killing sets are missing.
- The prompt also asks for the shapes of all 30,434 + 2,131 class-keeping ratios, not only those that need good primes. For the ones that bad primes empty alone, only the killing sets apply. That run is not started.
- "The smallest quotient of the genus-129 cover" is tabulated as the number of label coordinates. Each coordinate set is a (ℤ/2)ᵏ subcover; naming its genus needs z8's description of the cover.
- Two runs hung for about two hours with PARI's threads waiting forever under the memory cap. Both scripts now set PARI's `nbthreads` to 1. The outputs don't depend on it.

## Files and checks

- `i2_shape.py`, `run_i2.sh`, `i2_residual.py`, `run_i2_residual.sh`, `i2_needgood_list.txt` (the input), `i1_classkept_list.txt` (I1's 29 class-keeping ratios).
- `i2_i2_needgood_list_0of2.txt` and `_1of2.txt` are the shape run's two halves (ratio j of the list goes to worker j mod 2). `i2_residual_partial.jsonl` is the part-1 file, now superseded.
- PARI 2.19.1 here; CI has 2.17.4. `i2_residual.py --negative` prints exactly the committed output under PARI 2.17.2 too (cypari2 2.2.3 in a venv); 2.17.4 itself was not available here. PARI is used for factorisations (each factor proved prime) and prime lists, and the outputs don't depend on its order.

```
python3 experiments/roadmap-i/i2/i2_residual.py --negative     # 3 s
```

## Proposed LEDGER/CHANGES text

> **Experiment, not counted.** I2: for the ratios of H2 and G2 that need a good prime, the set left by the bad primes always has even size, 70% are cosets, and the killing good prime is at most 229; at it, the separating label quotient always involves at least two classes. A machine-readable file gives each ratio's residual tuples and the good prime that kills each. Computed (exact). Written by an AI model; not read by a person.
