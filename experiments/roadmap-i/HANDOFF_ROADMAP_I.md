# Roadmap I: hand-off, part 1 (27 Nov 2026, about 20:20 UTC)

> **Read `NOTES.md` first.** This is the session's report, committed unchanged apart from this line; where they differ, `NOTES.md` (after review z17) holds.


Written by an AI model (Claude); not read by a person. Experiments, not counted. Nothing here changes the census. Nothing was committed, pushed, posted or sent; Drew publishes.

**Part 1** covers I1, I2 (partial, marked), I3's pairs and rescan, I4 (interim) and I5. The addendum will bring the rest of I2 and its full residual file, I3's 50,014 new ratios, and I4's chunks 15000 to 19000.

- **Base:** main at 97c8564 (runs 7 and 8). Roadmap H's folders are now read from the repository (experiments/roadmap-h), with the old H*_DIR variables kept only as overrides. Roadmap H's corrections are applied: the complete bound for good p ≡ 1 mod 4 is 67,583 and p ≡ 7 mod 8 is below 1154; H1's 12 pairs are 7.
- **Tools:** real sympy 1.14.0; PARI 2.19.1 with cypari2 here, so every output comes from 2.19.1 (CI has 2.17.4). PARI is a finder; everything found is checked exactly.
- **Layout:** the bundle unpacks to experiments/roadmap-i/i1 to i5. Each folder holds scripts, outputs, `--negative` outputs, times.txt and REPORT_I<n>.md. A `.part` file is a partial output its script resumes from.
- **Not started, as asked:** I6 (not needed) and Θ₄ on I5's curves (the conductor's). I7 waits for idle cores.

## For the conductor

1. **I1: no survivor in 1,569 designed ratios (three designs).** 29 distinct ratios keep a class; bad primes alone empty all but (651432,64265), which good 29 empties. Design C, which works backwards from H3's tables, made the bad primes do all the work: 10 class-keeping ratios, every one emptied with P = 0. *(computed, exact)*
2. **I2: the residual file exists (partial).** `i2/i2_residual_partial.jsonl` gives, for 2,929 of the 3,876 ratios that need a good prime, every residual tuple (d0, d8, dL, dU) and the good prime that kills it. |S| is always even, 69% of the sets are cosets, the killing prime is at most 229, and the separating quotient at it always involves at least two classes. *(computed, exact)*
3. **I3 is complete for H1's ratios and confirms z8.** 7 of H1's 12 pairs are true phantoms, all settled by H1; 5 are not (each killed at an octic prime H1 never tested). The rescan of all 20,009 ratios differs from H1 only by those false pairs. *(computed, exact)*
4. **I4: two open rows so far, (14017,2856) and (14640,10679).** The first's partner (16873,11161) is in the chunk now running. The exchange test found 0 failures over every point reported so far. *(computed)*
5. **I5: no conditional row closes by a rank step.** All 9 curves have isogeny classes of 2 curves with one rational 2-torsion point, so H7's step can't apply. These are the 9 curves the conductor will try Θ₄ on. *(computed, exact; parity remark conditional)*
6. **(12703,2767) from its partner (7735,4968):** now noted in main's h7/NOTES.md, so nothing more is needed.
7. **Two small things in main's notes:** the H3 row of roadmap-h/NOTES.md and the (12703,2767) sentence in h7/NOTES.md each repeat a sentence three times.
8. **A practical tip:** PARI under `ulimit -v` can hang forever when its threads can't get stacks. Setting `nbthreads` to 1 fixes it and changes no output.

## Per task

- **I1 (done):** i1/REPORT_I1.md. Designs A (900), B (519), C (150). Laws that bite at bad and small good primes are tabulated, with the irredundant killing sets of the 19 class-keeping ratios from A and B.
- **I2 (partial, 2,932 of 3,876):** i2/REPORT_I2.md. The shapes of the class-keeping ratios that bad primes empty alone are not started.
- **I3 (pairs and rescan done):** i3/REPORT_I3.md. The new-ratio scan runs on 2 chains.
- **I4 (interim):** i4/REPORT_I4.md, with i4_exchange.py for the exchange test. It runs from the repository's experiments/roadmap-i/i4 so G3's relative paths resolve.
- **I5 (done):** i5/REPORT_I5.md. Next step for those rows would be a 4-descent or a 2-descent over the 2-torsion field, which is the conductor's call.
