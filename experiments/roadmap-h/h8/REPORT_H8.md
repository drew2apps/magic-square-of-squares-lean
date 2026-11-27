# H8: the PARI-free L-value for (6557,324)

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted.

## Bottom line

- **Not run.** H8 was "only if cores are idle". The 4 cores were busy for the whole run, first with H1, H2 and H3 (load average 9 to 11), later with H2 and G3.
- **Its purpose is met another way (H7, h7/REPORT_H7.md).** H8 was to replace the numerical L(E,1) behind the exclusion of (6557,324) and its partner (6881,6233). H7 shows rank 1 exactly for y² = g4·K at (6881,6233), by a complete 2-descent on its 2-isogenous curve. With F4's exact Demjanenko–Manin chain, every solution of (6881,6233) is then degenerate. (6557,324) is its partner, so it follows by the partner map the census already uses. *(computed, exact; resting on Vélu, Mazur, Silverman and Cremona–Prickett–Siksek, quoted)*
- **The L-value is still worth certifying:** it would give a second, independent exclusion. For (6557,324) itself, H7's step gives only rank ≤ 3.
- **Why it would not have fitted anyway:** lball_exact53.py takes only curves with conductor N ≤ 10¹⁴ (its selection line). The curve behind (6557,324) has N = 11077451294972790, about 1.1 × 10¹⁶ (PARI, as a finder), so it is outside what the script accepts, and the series behind an L-value grows like √N. Raising that limit is the first step if Drew wants the second exclusion.
