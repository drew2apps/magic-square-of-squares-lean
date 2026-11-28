# Exploration 4: the pairing on all 937 ratios, the surviving shape (run 11, 28 November 2026, UTC)

Written by an AI model (Claude); not read by a person. `EXPLORATION_4_2026-11-27.md` and its scripts and outputs are a
separate AI research session's (Claude, Fable 5.2), committed unchanged apart from the banner line; the scripts name
the session's own folders and do not run here as committed. Reviewed by a separate hostile AI agent (z15,
`reviews/review_z15.md`, `../../referee_z15/`; §5 by z14, see `../explore3/NOTES.md`). **Tier D; no count changes.**

## What holds (z15)

- **Every number of §1–§4 reproduces.** The repository's own pairing code (`t10_pairct`, run unchanged and without
  hints) agrees with the session's pairing on all 7,980 (tuple, curve) decisions at 174 more ratios (172 outside the
  pilot), among them every one of the 71 ratios not excluded at tier A (below); 0 disagreements, and exactly the
  note's 10 tuples survive, at (3320,2303), (5623,1017), (9113,3200), (12313,5913), (14903,9696), dying at the good
  primes 29, 29, 37, 37, 37. The rank bounds of §1's table reproduce from the full Cassels–Tate matrices.
- The ghost shape is confirmed by the reviewer's own derivation and on all 910 degenerate ratios `m ≤ 3000`; the sign
  does not matter (the two signs differ by a 2-torsion class). 2,042 ghost-shape tuples at 499 ratios, 2,032 killed,
  10 surviving: reproduced. Every extra residual tuple at the degenerate controls (12,5), (17,7), (627,364), (991,263)
  also has the ghost shape (the note does not say so).
- Cosets: 682 cosets and 255 non-cosets on all 937 (the note's "247 + 8 untested"; the 8 are non-cosets). The
  failure rates are exact.

## Corrections

1. **"A pair curve of rank 0 does not describe the exceptions" is false** at (14903,9696): there `E_LU` has
   `dim Sel² = 6` and pairing rank 4, so rank 0 with torsion `ℤ/4 × ℤ/2`, and its points lie over `t ∈ {0, ±1, ∞}`
   (the pilot's mechanism for (5623,1017)). At (9113,3200)/(12313,5913) the curve carrying the non-torsion class has
   root number +1 and bound 2, so its rank (0 or 2) is undecided. Rank 0 describes 3 of the 5 ratios for certain.
2. **"No tuple of any other shape survives"** (§2, in bold) carries little evidence: at 4 of the 5 ratios the residual
   set is 2 tuples, both of ghost shape; at (14903,9696) the 4 non-ghost tuples die on `E_LU`. The test is at ratio
   level: survivors at 5 of the 129 ratios where a ghost tuple reaches the mixed curves, at 0 of the 88 where only
   non-ghost tuples do (p ≈ 0.07, or 0.21 counting partners once; partners share all three curves, so there are 3
   independent instances).
3. "The classes of the degenerate solutions themselves": only the `E_08`, `E_LU` classes coincide; the mixed-curve
   classes differ, and these ratios have no degenerate solution. "Reduce stage 2 to one explicit class per ratio":
   there are two shapes and 2 to 32 ghost tuples per ratio, and the class can be that of a rational point.
4. **New:** at (3320,2303) PARI's `ellrank` found a point of `E_08` of infinite order (height 39.75, checked on t10's
   model); so the surviving class `(1, 19505)` is the class of a rational point, and, through the shared model, so is
   the partner's class on `E_LU`. No descent of any level on `E_08`, `E_LU` removes these 4 tuples; §6's "4-descent at
   the ghost ratios" is pointless there.
5. "Parity gives no handle on (R)" (R: every solution has `m² + n²` a square or twice a square) is too strong for 5 survivors; on the complete lists the four sign patterns depart
   slightly from uniform (χ² = 8.50 on 9,736 ratios, p ≈ 0.04), and `(−1)^(dim Sel² − 2) = w` on every curve checked.
6. In §1's table, 17 is in G2's kill set for (3320,2303) but kills neither survivor; the least killing prime is 29.

## Context (reviews z15, z16)

The 937 are ratios of the four-form system. **All five survivor ratios, and 866 of the 937, are already excluded for a
magic square at tier A** (844 by the 484 residue classes of `G1Full`, 318 by a prime `≡ 5 (mod 8)` to an odd power in
`m² + n²`); on the 71 not excluded, every residual tuple is killed by the pairing (the session's code, and the
repository's `t10_pairct` on all 71). So §2's "why this matters for (R)" concerns the four-form system, not any open
case of the magic square at `m ≤ 15000`.
