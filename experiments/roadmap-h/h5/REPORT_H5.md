# H5: conic points for Θ₄ on E⁺ of (945,248)

> **Read `../NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review) holds.

Written by an AI model (Claude); not read by a person. Experiment, not counted.

## Bottom line

- **Skipped, as the prompt says.** H5 waits for `scripts/g1_full/rank_fisher_theta4_gen.py` and `rank_fisher_theta4_find.py` to be on `main`. On 27 November 2026 at about 06:45 UTC, `git fetch origin main` showed `main` still at 062c54f (release 0.6.0), and neither file is in `scripts/g1_full/`. The `Θ₄` scripts there are `rank_fisher_theta4.py`, `_control.py`, `_k4.py`, `_k6.py`, `_order.py` and `_recheck.py`, all for (72,65) and its controls.
- Nothing was computed. In particular the δ′ of the U-conic over ℚ(√1356488820560645889749762), which "the finder's notes give", is not available here, so the conics could not be set up from the notes as asked.

## For the next run

When the two scripts are on `main`, H5 is a finder job: PARI with GRH (`bnfinit`, `bnfisnorm` or `rnfisnorm`, or the norm equation in L(√δ′)/L), each point checked exactly on its conic and handed back with its field. The pairing itself is computed elsewhere.
