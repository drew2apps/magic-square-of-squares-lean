# Roadmap I, parts 1 and 2 (run 11, 28 November 2026, UTC)

Written by an AI model (Claude); not read by a person. The folders `i1` … `i8` and the two hand-offs are a separate,
long-running AI session's (Claude, Fable 5.2), committed unchanged apart from banner lines; part 2's files replace
part 1's where both exist. Superseded partial files and logs were left out (2.3 MB: `i2_residual_partial.jsonl`, the
`.part` files, `i4/*.err`, `i4/done_14000`). Reviewed by a separate hostile AI agent (z17; `reviews/review_z17.md`,
its checks in `../referee_z17/`). **Experiment, not counted; tier D.** Nothing here changes the census. The session
continues (I3's new ratios, I4's chunks from 18000, I8 beyond line 360); these are interim where marked.

**Every committed output reproduces under PARI 2.17.4** (CI's version; the session had 2.19.1 and 2.17.2).

## What holds (z17), with the corrections

- **I1** (`i1/`): three designed hunts, **1,548 ratios** (1,569 rounds, 21 of them empty; 20 in design B, not 6),
  no survivor; 29 distinct ratios keep a class, the bad primes alone empty the four-class set of all but
  (651432,64265), which the good prime 29 empties. The report's killing-set section is out of date: 19 done, 10 over
  the memory cap.
- **I2** (`i2/`, complete): `i2_residual.jsonl` lists, for the 3,876 ratios of H2 and G2 that need a good prime, the
  32,474 tuples the bad primes leave and the first good prime that kills each; every set has even size, every killing
  prime is good, prime and at most 229. 50 seeded ratios recomputed by z17 with `conductor-t10`'s own code (20 of them
  G2's): 50 of 50 agree. Part 1's interim count is 2,930 records (it says both 2,929 and 2,932).
- **I3** (`i3/`): 7 of H1's 12 pairs are true phantoms and 5 are killed at the first octic primes listed by review z8;
  the rescan of H1's 20,009 ratios agrees. "Found independently" overstates it: the session's base already carried
  z8's correction. The scan of 50,014 new ratios is unfinished.
- **I4** (`i4/`, interim): G3's pipeline leaves (14017,2856) with its partner (16873,11161), (14640,10679) and
  (16065,2032) open in 14000–15000 and 16000–17000. **The exchange test is vacuous on these data** ("0 failures over
  1,927 points": every point is "no" (828) or trivial (1,099), so the ratio check never ran); on G3's outputs on main
  it gives 1,064 "no" and 1,501 trivial (the report says 767 and 1,098), still 0 failures. REPORT_I4's "111 ratios of
  chunk 18000 done" is 109 in the bundle.
- **I5** (`i5/`): the exact descent (`rank_indep.rank_bound`) gives rank `≤ 3` on 6 and `≤ 5` on 3 of the 9 curves
  behind the rows conditional on a rank being 1, so none closes (run 11's `../conductor-t11/` lowers two of the 5s to
  3 by Fisher's pairings, none to 1). (12703,2767) follows from its partner (7735,4968) — H7's exact rank and F4's exact
  Demjanenko–Manin chain — through the census's partner step: **computed (exact), resting on H7's and F4's quoted
  statements and the partner step**, not "proved" and "unconditional" as the report says.
- **I8** (`i8/`, checkpoint, lines 1–360 of the 3,876 of `i2_needgood_list.txt`, all from H2's exhaustive list):
  the Cassels–Tate pairing on the six pair curves (t10's code) removes 3,160 of 3,166 residual tuples; 11 ratios are
  also excluded by a pair curve of rank 0; 462 tuples are killed only on the mixed curves. z17 reran t10 from scratch on
  the 3 "neither" ratios and 15 seeded ones: 18 of 18 identical. The G2 ratios are lines 2,949–3,876 (928) plus 9
  labelled `h2-hunt` among lines 2,808–2,948, not "lines 2,940 to 3,876".
- **The three "neither" ratios** (6 tuples, in the kernel or the image on all six pair curves) close by
  Demjanenko–Manin on an `E±` curve (`dm_pari`'s Sig±, the same models):
  - (18313,120) and its partner (18433,18193): **conditional on a rank being 1** (PARI's bound; the exact descent
    gives 3; no curve of the isogeny class has full rational 2-torsion), on PARI's saturation at primes `≤ 1000` and on
    floating-point heights;
  - (20056,19825), through its partner (39881,231): **rank exactly 1** on that curve (review z17: it is 2-isogenous to
    two curves with full 2-torsion where H7's step gives `dim Sel² = 7` with pairing rank 4 and `dim Sel² = 3` with
    pairing rank 0), conditional only on the saturation and the heights.
  So the pairing does not replace every good-prime kill: these 6 tuples fall to a good prime and to a rank-1 curve. All three ratios lie in the 484 residue classes of the census, so they are already excluded for a
  magic square at tier A (release panel P2).

## Rerun

Each script prints `# matches` or `# DIFFERS` and exits 0 either way, so CI runs them through its `hm` helper
(job `research`, about 3 minutes): `i1_laws.py` (and `--negative`), `i1_negative.py`, `i2_residual.py --negative`,
`i3_scan.py --negative`, `i4_exchange.py --negative`, `i5_rank.py --negative`, `i8_pairct.py --negative` and
`--compare50`, and z17's `z17_i8count.py`, `z17_i8compare.py`, `z17_i2check.py`, `z17_neither.py a b c e`. Longer, not
in CI: `i5_rank.py` (3–8 minutes), `i3_scan.py 13` (1.5–5 minutes), `z17_neither.py d` (2–4 minutes).
`i8_checkpoint.py 360 3` needs the session's worker files and does not run from the bundle.
