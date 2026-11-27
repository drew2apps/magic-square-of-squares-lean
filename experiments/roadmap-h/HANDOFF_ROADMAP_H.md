# Roadmap H: hand-off (27 Nov 2026)

> **Read `NOTES.md` first.** This is the session's hand-off, committed unchanged apart from this line; the reviews corrected H1, H2, H3 and parts of H6 and H7, and `NOTES.md` and each folder's `NOTES.md` hold where they differ.

Written by an AI model (Claude); not read by a person. Experiments, not counted. Nothing here changes the census, which stops at m ≤ 1000. Nothing was committed, pushed, posted or sent; Drew publishes.

- **Base:** release 0.6.0 at 062c54f, and the roadmap-G bundle (experiments/roadmap-g). Nothing in 0.6.0 or in the G bundle is redone, except as a cross-check (H2's control block, below).
- **Tools:** real sympy 1.14.0 (no stand-in). PARI 2.19.1 with cypari2 on this machine; CI has PARI 2.17.4 and cypari2 2.2.9, so every output below comes from 2.19.1. PARI is a finder (factorisations, points, conic points, torsion points), and everything found is checked exactly. Where a line says "rank by PARI", the rank is consistent with PARI and not proved.
- **Layout:** the bundle unpacks to `experiments/roadmap-h/h1` to `h8`, plus `g3/` (roadmap G's G3 folder, continued here by H7 with block 13000; its table files are rewritten). Each folder has its scripts, committed outputs, `--negative` outputs where a script certifies something, `times.txt` and `REPORT_H<n>.md` (bottom line, checks, what failed, LEDGER text, CI lines with times). Scripts find the repository through MSS_REPO, or else `../../..`; G1's, G2's, H3's and H4's folders through G1_DIR, G2_DIR, H3_DIR, H4_DIR, or else the relative paths in each script.
- **Labels:** proved / computed (exact) / computed (numerical) / conditional (on what) / open.

## H1: the phantoms of the four-class descent
- **Bottom line:** a scan of 20,009 random degenerate ratios up to m = 10⁹ finds 13 **phantoms** (ratios whose descent keeps tuples beyond those of t = 0, ±1, ∞). With G1's (1109503,597671) and the partner map, they form 12 partner pairs, and the extra tuples pass every local test (through H3). In **11 of the 12 pairs** every extra tuple is the class of no rational point: 7 by an exactly computed rank-0 quotient y² = w·f0·f8 (rank bound taken over the isogeny class), 4 by the pair curves' complete 2-descent and Cassels–Tate pairing. **1 pair is open:** (11875080,6497849) / (18372929,5377231).
- **Label:** computed (exact rank bounds and descents). The step from one side of a pair to the other is exact too: the partner isomorphism is written down, t′ = (t + 1)/(t − 1) with the forms permuted (0, 8, L, U) → (L, U, 0, 8), and checked to carry known tuples to known tuples for all 12 pairs (h1_iso.py). The extra tuples of (1109503,597671) are settled twice, both times on that ratio itself.
- **Files:** h1/ (h1_pairs.py, h1_scan.py, h1_scan_k13 to k29, h1_partners.py, h1_allphantoms.py, h1_rank0.py, h1_rank0_exact.py, h1_holdouts.py, h1_iso.py, their outputs and negatives).
- **Next step:** the open pair. Its pair curves have 3 ≤ rank ≤ 4 and both extra classes are classes of points found on them, so the pair-curve test cannot decide. A Mordell–Weil sieve on two pair curves at once (same t) is the natural next tool.

## H2: survivor hunt
- **Bottom line:** no survivor. Every coprime non-degenerate ratio with 15000 < m ≤ 30000 (205,172,879; 30,434 keep a class) is emptied by its bad primes and good primes ≤ 229. A biased hunt (m from 10⁴ to 10¹²) and an excess-guided hunt built on H3's tables add 2,421 class-keeping ratios, all emptied, the largest cutoff 137. P grows slowly with m (101 to 229 per block of 1000), and by H3 never needs to pass 1153. A rerun of 14900 < m ≤ 15000 agrees with G2 exactly.
- **Label:** computed (exact).
- **Files:** h2/ (h2_hunt.py, h2_design.py, h2_exhaust.py, h2_summary.py, h2_control.py, h2_retry.py, their outputs).
- **Next step:** the exhaustive run past 30,000 costs about 2.5 core-hours per block of 1000 at this size and grows slowly.

## H3: local images at good primes
- **Bottom line:** good primes above 1153 never remove anything (each twist is a smooth genus-17 curve; Hasse–Weil gives it points), and good primes ≡ 3 mod 4 never do either. So the full local four-class set of any ratio is decided by its bad primes, ℝ and the 96 good primes ≡ 1 mod 4 up to 1153. Tables for those primes, computed two ways, show losses at every one. 15,984 of the 49,764 allowed sets are not affine subspaces of F₂⁴.
- **Label:** proved given Hasse–Weil (quoted); tables computed (exact).
- **Files:** h3/ (h3_tables.py, h3_weil.py, h3_sizes.py, tables, checks, negatives).

## H4: G1's three out-of-memory ratios
- **Bottom line:** all three are emptied by their bad primes alone. All 63 of G1's kept ratios are now done: 60 by bad primes alone, 3 with good primes, the largest 61.
- **Label:** computed (exact).
- **Files:** h4/ (h4_split.py, h4_split.txt, h4_negative.txt).

## H5: conic points for Θ₄
- **Skipped, as the prompt says:** the finder scripts are not on `main` (still 062c54f on 27 Nov at about 06:45 UTC). Nothing computed; Θ₄ itself was not computed.

## H6: Demjanenko–Manin and Chabauty on the quotients of (627,364) and (991,263)
- **Bottom line:** no route opens. The four genus-2 pieces through g4 are one curve A up to isomorphism, so Demjanenko–Manin needs rank J(A) ≤ 3; the Richelot descent gives ≤ 8. Every other genus-2 piece has bound 6 to 9. The two ratios mirror each other (991 = 627 + 364, 263 = 627 − 364).
- **Label:** computed (genus-2 bounds exact; genus-1 ranks by PARI).
- **Files:** h6/.

## H7: G3 toward 30,000, and F4's g4·K route
- **Bottom line:** the g4·K route applies to one of G3's four open ratios, (12201,7360) (PARI rank [1,1]; the other three have PARI rank 2 or more). F4's exact Demjanenko–Manin chain shows every solution there is degenerate if rank 1. A new exact step gives that rank: the 2-isogenous curve has full 2-torsion, and the repository's complete 2-descent with the Cassels–Tate pairing gives rank ≤ 1 on it. The same step gives rank 1 exactly for two of F4's three conditional ratios, (7735,4968) and (6881,6233), so F4's closures there no longer depend on PARI's rank. (6557,324) gets only bound 3 on its own curve, but it is the partner of (6881,6233), so it follows. G3 block 13000 < m ≤ 14000 finished (G3's table now runs to 14000, 916 ratios): of its two new open ratios, (13348,3795) closes the same way, exactly; (13899,1180) stays open (PARI rank 2).
- **Beyond H7:** the same descent gives rank exactly 1 for the curve E04 of 90 of the 91 ratios that tier B closes on E04 (h7_e04_rank.py). Their upper rank bound then needs no PARI, Sage or Magma. Their generators and heights are unchanged. (864,151) is not covered. The D1 and E0L ratios (31) are not covered; no D1 curve has full 2-torsion in its class.
- **Label:** computed (exact), resting on Vélu, Mazur, Silverman and Cremona–Prickett–Siksek (quoted); the dm_g4k table is computed (numerical).
- **Files:** h7/ (run_h7_b13000.sh and its h7_b13000_* outputs, run_h7_g4k.sh, h7_g4k.jsonl, run_h7_g4k_exact.sh, h7_g4k_exact.txt, h7_g4k_rank_exact.txt, h7_g4k_rank1.py, h7_e04_rank.py, their outputs and negatives).
- **Next step:** G3's other four open ratios (now including (13899,1180)) have g4·K of PARI rank 2 or more, so they need another route. (7735,4968) was the one open ratio of CENSUS9000.md, and (6557,324)/(6881,6233) rested on a numerical L-value there. If Drew accepts the rank step, both are excluded exactly (the second pair through the partner map).

## H8: ap_one for (6557,324)
- **Not run:** the cores were never idle. **Its purpose is met another way:** H7 gives rank 1 exactly at (6881,6233). That excludes (6881,6233) and, through the partner map, (6557,324), without the numerical L-value. The L-value would still be a second, independent check.

## Ideas tried, and what came of them
- **Hasse–Weil finiteness (H3):** turns "good primes to some cutoff" into a finite, complete local test. It worked, and it bounds H2's P.
- **The shape of the local sets (H3):** a third of them are not affine, so a proof of the reduction theorem from these tables cannot be pure linear algebra over F₂.
- **The Cassels–Tate kernel test (H1):** excludes a class when it is not the class of any element of the kernel of the pairing, even when the rank bounds differ. It is what kills the extras at 302617 and 550548247.
- **Partner pairs (H1):** every phantom's partner is a phantom, and the explicit isomorphism t′ = (t + 1)/(t − 1) makes a proof on either side settle the pair. Two pairs were settled only on the partner side.
- **Rank-0 quotients with exact bounds over the isogeny class (H1):** at (349960367,118678903) the first model gives rank ≤ 2, an isogenous one gives 0.
- **The excess-guided hunt (H2):** H3's tables as a heuristic for where survivors could hide. It finds class-keeping ratios 3.4 times as often, but they are easier: 169 of the 170 whose descent finished are emptied by bad primes alone.
- **s-line quotients (scratch, not in the bundle):** the involutions t → 1/t and t → −1/t of C give the genus-1 quotients E± (from F_0F_8F_L) and the even quartics Q1, Q2 over s = t ± 1/t. On the open pair none has rank 0 (exact bounds 2 to 4 over the isogeny classes of Q1, Q2; PARI rank 2 for Q1 at 18372929; E± PARI bounds [1,3] and [2,2]).
- **g4 quotients (scratch):** w·f0·g4 and the rest have PARI rank at least 1 at the holdouts tried.
- **Richelot (H6):** negative, as above.
- **Rank 1 through a full-2-torsion isogenous curve (H7):** the complete 2-descent and pairing on the 2-isogenous curve with full 2-torsion is sharper than the descent via 2-isogeny on the curve itself. It turned three conditional rank-1 statements into exact ones, and gives the upper bound for 65 of the 66 E04 curves of tier B. It does not apply to the 12 σ curves of Result 4 (no full 2-torsion in the class), and cannot beat PARI's [1,3] on the seven rank-1-route exclusions.
- **Demjanenko–Manin on H1's open pair (scratch):** g4·K has PARI rank 1 there, but F₄ is not one of the four forms of H1's curve. On y² = w·K, which is a quotient of it, t → 1/t is translation by the 2-torsion point (0,0), so DM gives nothing.
- **Not done:** a Mordell–Weil sieve on two pair curves, which would need exact ranks and saturated generators of the pair curves.

## Not done
- **H5:** waits for the finder scripts.
- **H8:** no free cores.
- **Not started, as asked:** Θ₄, and the (627,364) 4-descent over ℚ(i). No census closure or "first" is claimed.
