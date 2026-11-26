# g4k: Demjanenko–Manin on the quotient y² = g4·K (run 5, F4)

**Experiment, not counted; tier D of `LEDGER.md`.** Scripts and outputs written by an AI model (Claude, the "Fable"
session of run 5); reviewed by a separate AI agent (`reviews/review_x3.md`: the DM step holds; the closures are
conditional on rank 1). Not read by a person.

For a G1 ratio, w = t² gives a point of C: y² = g4(w)K(w), and t → 1/t a second point with x·x′ = c3² (x = c3·w), so
Demjanenko–Manin applies to C. `g4k_exact.py` does the DM step in exact arithmetic (torsion ℤ/2; an enclosure of ĥ(G);
saturation at every prime ≤ 199; an index bound from the Cremona–Prickett–Siksek height bound on the minimal model; the
DM bound and every candidate). It shows that **if rank C(ℚ) = 1**, every solution is degenerate for (7735,4968),
(6557,324) and (6881,6233). Rank 1 is not proved: it is PARI's `ellrank` upper bound (an unverified program output; its
class groups are GRH-dependent unless certified; for (6557,324) and (6881,6233) only one curve of the isogeny class gives
the bound 1), and the repository's exact descent (`g4k_rank_exact.py`, through `scripts/g1_full/rank_indep.py`) gives 2,
3 and 5. Tier: conditional on rank C(ℚ) = 1 and on the published Silverman and CPS statements quoted (tier E). The
review re-derived the DM step with its own exact code (`experiments/referee_x3/`) and got the same bounds and candidates.

- `g4k_exact.py` (+ `.txt`, `_negative.txt`): the exact chain; PARI only for hints that are checked.
- `g4k_rank_exact.py` (+ `.txt`, `_negative.txt`): the exact rank bounds 2, 3, 5.
- `dm_g4k.py` (+ `.jsonl`, `_negative.txt`): cross-check with `experiments/census2000/dm_pari.py`'s DM step (PARI
  generators, floating-point heights) on these three and the nine other ratios; PARI's rank of g4·K is not [1,1] for any
  of the seven conditional ratios ((6121,959) has [1,3]).
- The session's original `sympy_pari_shim.py` (a PARI stand-in for sympy) is not committed: the scripts use sympy.
