# t10: are the good-prime kills shadowed by a Cassels–Tate obstruction on a pair curve? (run 9)

Written by an AI model (Claude, build agent t10); not read by a person; reviewed by a separate AI agent (z11; `reviews/review_z11.md`), not by a person. **Experiment, not counted;
tier D of `LEDGER.md`.** Nothing here changes the census.

## Question

For a ratio `m : n` with `m² + n²` neither a square nor twice a square, the four-class descent (`LAST_PAIR.md`
§11–12) leaves, after every bad prime (G2: `../roadmap-g/g2/`), a non-empty set of tuples `τ = (d0, d8, dL, dU)` for
the 101 of the 997 ratios of `tuple_selmer_997.txt` (`g2_997_all.txt`; 291 ≤ m ≤ 5000) and the 836 of `g2_15000_need_good.txt`, and a small good prime `p ≡ 1 (mod 4)`
removes them. Is each such tuple also removed by a global, reciprocity-type mechanism: the Cassels–Tate pairing on
`Sel²` of a pair curve `C_jk: y² = F_j(t) F_k(t)`, `j, k ∈ {0, 8, L, U}`?

## Method (`t10_pairct.py`)

- Residual tuples exactly as G2 (`exact_local.py` at every odd bad prime, the repository's test at 2; ℝ imposes
  nothing); for each, the least good prime that removes it.
- Each of the six pair curves with a rational base point `t0`: `t0 = 0` for `08` and `LU`; `t0 = m/n, −m/n, n/m, −n/m`
  are points of `C_0L, C_0U, C_8L, C_8U` for every ratio (`F_0 F_L(m/n) = 4m²(m² − n²)²/n²`); a Möbius map moves
  `t0` to `s = 0`, Washington's map (as `h1_pairs.py`) gives `ψ: C_jk → E_jk`, an integral model with full rational
  2-torsion. The roots `e_r` (factorisation `{F_j, F_k}`) and `e*` (`{a_j a_k, conj}`) are identified exactly by the
  affine map between the roots `a(θ₁θ₂ + θ₃θ₄)` of the resolvent cubic and the `e`'s (exactly one permutation fits;
  asserted), and checked on the 2-torsion points and on every point found.
- A rational point with tuple `τ` has `F_j(t)` a square and `a_j a_k(t)` of class `D = d_j d_k`, so its Kummer image
  is `c(τ) = (x − e_r: κ_r, x − e*: ±κ* D)`, the constants `κ` read off every point found (all agree; asserted). The
  two signs differ by the class of `ψ(ι t0)`, a rational point (checked). So `c(τ) ∈ δ(E(ℚ))`, which lies in the kernel
  of the Cassels–Tate pairing on `Sel²`.
- **Pairing-killed** at `(j, k)`: `c(τ) ∉ Sel²` (complete 2-descent, `rank_second.full_selmer2`), or
  `⟨b, c(τ)⟩ = −1` for some `b ∈ Sel²` (`rank_second.DFull.pairing`; `b` over a basis of `Sel²` modulo the images of
  the points found; the killing pair is checked symmetric, the pairing alternating on `c(τ)` and trivial on a point
  image). A pairing-killed tuple is the class of no rational point (the 2-covering of `E_jk` with class `c(τ)` is everywhere locally soluble and has no rational point; the pairing is its Brauer–Manin obstruction from `Sel²`).

## Results (computed; `t10_summary.txt`)

**Answer: sometimes, almost always.** Of 1,766 residual tuples that a good prime removes, **1,762 are also
pairing-killed** and **4 are not** (those 4 are removed by rank 0 of `E_08`, resp. `E_LU`; see below).

- Ratios: the 101 of `g2_997_all.txt` that need a good prime (`t10_997_a.txt`, `t10_997_b.txt`) and 100 drawn with
  seed 1 from the 836 of `g2_15000_need_good.txt` (`t10_15000_s100_a.txt`, `_b.txt`). The residual sets have G2's
  sizes for all 201 (598 + 1,168 = 1,766 tuples). Each tuple dies at a good prime ≤ 173 (first killing prime: 17 for
  420, 29 for 716, 37 for 290, 41 for 166, 53 for 76, 61 for 48, 73 for 8, 89 for 10, 97 for 4, 101 for 20, 109 for 4,
  173 for 4), never above the largest good prime of G2's kill set.
- **Every residual class lies in `Sel²` of all six pair curves** (no "outside `Sel²`"), as it must: the 2-covering of class `c(τ)` is locally soluble at the bad places (a residual tuple has local points there) and
  has good reduction at the other primes. So the good prime never enters at the pair level. Every kill comes from the
  pairing.
- **Where the −1's sit.** The 8,930 kill lines (one per killed tuple and pair curve, one `b` each; 4,976 distinct
  (ratio, curve, class)) carry 18,222 local terms −1, an odd number on each line, all at bad primes of the ratio, none
  at ℝ or at a prime only auxiliary for `rank_second` (the primes of `δ`, `det`). This is forced, not evidence: the
  bad primes of every `E_jk` (1,206 curves) are among those of the ratio, the local terms are 1 outside the bad
  primes of `E_jk`, 2 and ℝ, and which of those places carries a −1 depends on the local points and tangent lines
  chosen (review z11 found other places with other choices); only the product is intrinsic.
- Kills per pair curve: `08` 1,280; `0L`, `0U`, `8L`, `8U` 1,694 each; `LU` 874. (`C_0L ≅ C_0U ≅ C_8L ≅ C_8U` by
  `t ↦ −t, 1/t`, but each of the four tests a different class: `d0 dL`, `d0 dU`, `d8 dL`, `d8 dU`.) 199 of the 201 ratios
  have every residual tuple pairing-killed.
- **Not pairing-killed:** (5623,1017): `(2, 1270798, 1, 2)` and `(635399, 1, 1, 2)`. Its partner (3320,2303):
  `(2, 1, 1, 635399)` and `(2, 1, 1270798, 2)`. All four die at the good prime 29. For (5623,1017):
  - On `E_08`, `Sel²` is the image of the torsion (dim 2), so `E_08` has rank 0, and its torsion has order 8 (8 points; the gcd of `#E(F_p)`, `3 ≤ p ≤ 400` good, is 8: review z11). The 8 rational points of `C_08` are the obvious ones, over `t = 0, 1, −1, ∞`, where the four forms are not squares together. So the ratio has no four-form point: the tuples are removed by `E_08` after all, by rank 0, not by the pairing (their class is a torsion class, so no Selmer or Cassels–Tate condition on `E_08` sees them).
  - On `E_0L, E_0U, E_8L, E_8U` (dim `Sel²` = 6, the pairing trivial on `Sel²`) and on `E_LU` (dim 5), the class is in
    the kernel. It is not the class of any point found, even with PARI's `ellrank` at effort 3, whose bounds `[2, 4]`
    and `[1, 1]` are not relied on (`t10_deep.txt`).
  - The partner gives the same picture with `08 ↔ LU` (`E_LU` of (3320,2303) has rank 0 likewise). Three other ratios of the 201, (9191,7688), (9532,4445), (14885,6076), have `E_LU` of rank 0 and are excluded the same way (their tuples are also pairing-killed).
- **Degenerate controls** (`t10_negative.txt`):
  - The known tuples of (4,3), (7,1), (12,5), (17,7), (627,364), (991,263) survive every pair curve.
  - The extra tuples of (12,5) and (17,7), which die at 29, are pairing-killed.
  - The extra tuples of (627,364) and (991,263), which die at 41, are not. They are classes of rational points on
    `E_08` and `E_LU`, and on the other four curves they are in the kernel without being the class of any point found
    (`t10_deep.txt`; PARI's bounds `[4, 4]` there come from 3 points found and are not relied on).
- **The pairing does kill tuples that no good prime kills,** at H1's phantoms. `--h1control` (`t10_h1control.txt`)
  takes (1109503,597671), whose four extra tuples are locally soluble everywhere (H1, review z8). It kills each of them
  on exactly the pair curve that `h1_pairs.txt` names, and the known tuples survive. At non-degenerate ratios this
  cannot be observed, since every residual tuple dies at a good prime.
- **Forged inputs are refused:** a class times 11 is outside `Sel²`, exchanged roots make the constants disagree, a
  false base point is rejected, and with the pairing left out nothing is killed.

**What it means.** For 199 of the 201 ratios, the Cassels–Tate pairing on `Sel²` of the genus-one quotients `C_jk`
removes every tuple that the good primes remove. For the partner pair (5623,1017) / (3320,2303) the pairing leaves
two tuples each, but there `E_08` (resp. `E_LU`) has rank 0 and its 8 points lie over `t ∈ {0, ±1, ∞}`, so these two
ratios have no four-form point at all (review z11). So on these 201 of the 937 ratios that need a good prime, every
tuple that a good prime removes is also removed by a computation on one pair curve over ℚ (the pairing, or rank 0).
That is a candidate global mechanism, not a proof of one: the pairing does not reproduce the good-prime kills of the
degenerate controls (627,364), (991,263), and the residual sets rest on `exact_local.py`. The "E±" curves of the
question were not tried.

## Labels

Computed (exact), resting on the complete 2-descent and the Cassels–Tate pairing as `rank_second.py` computes them
(Cassels 1998; Skorobogatov §6.2; `reviews/review_t2.md`), on Washington's map and the Kummer map (quoted), and for the
residual sets on G1's `exact_local.py` and the repository's test at 2 (tested, not proved). PARI finds points, base
points and conic points; each is checked. Tier D.

## Rerunning

```
python3 t10_pairct.py --negative                 # controls and forged inputs (t10_negative.txt)
python3 t10_pairct.py --h1control                # (1109503,597671) against h1_pairs.txt (t10_h1control.txt)
python3 t10_pairct.py --list ../roadmap-g/g2/g2_997_all.txt --from 0 --to 79 --out 997_a      # t10_997_a.txt
python3 t10_pairct.py --list ../roadmap-g/g2/g2_997_all.txt --from 79 --to 101 --out 997_b    # t10_997_b.txt
python3 t10_pairct.py --list ../roadmap-g/g2/g2_15000_need_good.txt --sample 100 --seed 1 --from 0 --to 50 --out 15000_s100_a
python3 t10_pairct.py --list ../roadmap-g/g2/g2_15000_need_good.txt --sample 100 --seed 1 --from 50 --to 100 --out 15000_s100_b
python3 t10_pairct.py --summarize t10_997_a.txt,t10_997_b.txt,t10_15000_s100_a.txt,t10_15000_s100_b.txt   # t10_summary.txt
python3 t10_pairct.py --deep '5623,1017;3320,2303;627,364;991,263' --effort 3                   # t10_deep.txt
```
Add `--write` to write the outputs (and the hints `t10_hints_*.json`); without it the script only prints.

Times (one process, PARI 2.17.4): about 25 s per ratio on average; the four list runs about 2 h of CPU in all;
`--negative` about 20 s, `--h1control` 5 s with H1's hints, `--deep` 5 min. Memory: up to about 3.5 GB per process.

**Provenance of the outputs.** The run of the 101 ratios was stopped by the out-of-memory killer after 79 ratios.
`t10_997_a.txt` is that run's output, cut from its log: the header line was rewritten to that of `--from 0 --to 79
--out 997_a`, the timing lines (stderr) were dropped, and there is no SUMMARY block and no hints file for it. A rerun
of its first three ratios with the committed script gives the same lines. The other outputs were written by `--write`.
An earlier, slower version of the script tested both sign choices of every class. On the 662 tuples that it and the
committed version both did, its verdicts (the killing pairs and the killing prime) are the same.
