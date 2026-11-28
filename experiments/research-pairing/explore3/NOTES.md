# Exploration 3: the mixed curves, the modular picture, parity (run 11, 27 November 2026, UTC)

Written by an AI model (Claude); not read by a person. `EXPLORATION_3.md` and its scripts and logs are a separate AI
research session's (Claude, Fable 5.2, with two sub-agents), committed unchanged apart from the banner line; the
scripts name the session's own folders (`/tmp/mix`, `/tmp/t10`) and do not run here as committed. Reviewed by a
separate hostile AI agent (z14; `reviews/review_z14.md`, its checks and path-adapted copies in
`../../referee_z14/`). **Tier D; no count changes.**

## What the review found (z14)

- **No disagreement with the repository's own descent and pairing** (`rank_second.full_selmer2`,
  `t10_pairct.DFullX`) on 134 mixed curves of ratios outside the pilot and its partners (80 non-degenerate of
  opposite parity, 40 both odd, 8 square type, 6 twice-square type):
  - the session's general pairing `cassels.py` agrees on 16,401 of 16,401 pairs (a check of the implementation:
    it implements the same formula as `DFull`);
  - **the Rédei form** (§1: when `x² − d y² = β z²` has a point, the pairing of `(1, d, d)` and `(1, β, β)` is a
    product of Legendre symbols of `x0 + y0√d` at the primes of odd exponent in `W = m⁸ + 14m⁴n⁴ + n⁸`) agrees on
    863 of 863 soluble pairs; its independence of the solution and of the sign of `√d` is forced (the odd primes of
    `W` are `≡ 1 (mod 12)` and `W ≡ 1 (mod 8)`), not evidence; it is **tested, not derived** (see exploration 5);
  - `{d : (1, d, d) ∈ Sel²}` lies inside the norms from `ℚ(√W)` on 134 of 134 curves and equals them on none:
    being a norm is necessary, not sufficient ("exactly when", §1, is wrong).
- **"Only the second reading is defined"** (§1) is wrong: Stevenhagen's Theorem 1.1 (= 7.9, arXiv:1806.06250)
  needs all three Hilbert symbols trivial, a symmetric condition, so where `x² − d y² = β z²` has no point neither
  reading is defined (4,279 such pairs in the reviewer's sample, 1,266 with pairing −1).
- **Counts are lines, not cases**: the four mixed curves are one curve per ratio, printed four times, and partner
  ratios give the same curves; "804 curves" are 160 distinct; the 6,776 mixed kill lines are 1,183 distinct
  (ratio, a, b), 1,009 distinct (curve, a, b); 2,760 / 488 / 3,528 are 400 / 79 / 530 cases; the two samples (seeds
  11 and 12: 57 and 160 curves) overlap, so 14,160 and 16,638 are not independent.
- **§2, the modular picture**, checks symbolically (discriminants, six `I4` fibres of `E_08`, full 4-torsion over
  `ℚ(i)(r)`, `E_LU` the twist by −1, twelve `I2` fibres of the mixed curve, a section of infinite order), with two
  rewordings: the eight zeros of `W` are where `j(E_08) = 0` (on the mixed curve they carry `I2` fibres, `j = ∞`);
  the octahedral group `S4 = PSL2(ℤ/4)` permutes the six cusps (the line's automorphism group is `PGL2`; only a
  dihedral group of order 8 is defined over `ℚ`). "`γ ↦ γ⁻¹`" was not checked.
- **§4, parity**: 37 / 34 / 62 and 36 / 47 / 41 / 26 confirmed with `rank_second` on 150 of 150 ratios (an
  independent draw gives 43 / 47 / 77 and 40 / 37 / 41 / 32); "38 of 80" (square type) does not reproduce (the
  height used was not recorded; 40 at `H = 3000`); "should have rank 0 or a large Ш" is mis-stated: an even bound
  `≥ 2` leaves rank 0 with `Ш ⊇ (ℤ/4)²` or rank `≥ 2`; `rank0.py` would count a failed pairing evaluation as −1
  (none occurred).

## What it means

The pairing on the four mixed curves is, on all tested cases, a Rédei symbol whose third entry is the octahedral
form `W`; derived for classes with first component 1 and `β` a norm in exploration 5 (Theorem B), with the local
factors at 2, ∞ and the primes of `d`, `β` not yet explicit. The tuples' own classes have first component 2, which
none of this covers. Exploration 4 and review z16 put the stage-2 data in context (`../explore4/NOTES.md`).
