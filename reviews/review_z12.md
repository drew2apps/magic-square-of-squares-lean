# Review z12: closed formula for the Cassels–Tate pairing on E_08 and E_LU (`pairct_formula`)

Written by an AI model (reviewer z12), not read by a person. 27 November 2026. Scratch and outputs: `/home/claude/conductor/referee_z12/`.

> Committed as written apart from this note. In this repository the reviewer's scratch folder is `experiments/referee_z12/` (paths made relative, `README.md` there); `pf_run` was a copy of `experiments/research-pairing/` with the pilot's `t10_*` files. The brief's own counts (133 + 111) were the conductor's.

## Bottom line

**Yes, on the sample tested.** `ct08.pairing` agrees with `rank_second` (`full_selmer2` and `DFull.pairing`, run through `t10_pairct.PairData` and `DFullX`) on **every** distinct (ratio, curve, a, b) tested: **0 disagreements** in 18,428 (E_08) + 10,490 (E_LU) = 28,918 triples. The sample is 146 ratios on E_08 and 133 on E_LU, out of 201; (5623,1017) and (3320,2303) are included. a runs over the class of every residual tuple (killed or not, both signs; 1,736 + 1,540 triples) and over every element of Sel² of the form (1, δ, δ) (16,692 + 8,950 triples). b runs over the full Sel² basis.

The run was cut short at about 90 minutes, after the out-of-memory killer had stopped it once; this was not a choice of sample.

In the same run:
- ct08 over **all** of Sel² equals the linear extension of rank_second (69,304 + 47,712 checks).
- Seeds 1 and 7 give identical values everywhere; ct08 never returned None.

Weaknesses:
- **[fix]** ct08 leaves out the real place, which is correct only for A > 0.
- **[fix]** The doc never states gcd(m, n) = 1, and the formula fails without it.
- **[fix]** The doc says the component assignment was "fixed by the data". In fact theory forces it, and a data-free test of local constancy rejects the other assignment. Fitting it to the data did not weaken the test.
- **[fix]** Several counts are inflated by repetition.
- **[fix]** The §6 reading of arXiv:2303.05059 is wrong.
- **[fix]** §4's expectation about the mixed curves is stated as fact in the opening paragraph.

There are no blockers.

## Findings

1. **[note] The test the author did not run, now run.** 0 disagreements in the counts above, with rank_second's values split +1/−1 as 12,146 / 6,282 (E_08) and 7,744 / 2,746 (E_LU). Every tuple class was in Sel². The §1 class formula holds on every residual tuple of the sample, on 2,064 (E_08) + 1,904 (E_LU) (tuple, sign) instances:
   - E_08: (1, ±mn·D);
   - E_LU: (2S, ±2D), which becomes (1, ±(m² − n²)D) after multiplying by the image of the point of order 4.

   Commands:
   - `cd referee_z12; python3 fulltest.py full0.jsonl 1,7 m,n,08 ...` (job lists in `jobs0.txt` and `jobs1.txt`; `fulltest_lowmem.py` for the restart);
   - `python3 -I analyse.py`;
   - `python3 -I classcheck.py`.

2. **[fix] The counts are repetitions.**
   - The 1,280 kill lines of E_08 are 133 distinct (ratio, class, b), over 121 ratios. The 874 of E_LU are 111, over 94 ratios. The 133 + 111 in the brief are confirmed.
   - The 21,838 "odd local terms" are 1,321 + 1,094 = 2,415 distinct (ratio, δ, b, p).
   - The 150 controls are drawn from the repeated lines.

   §3 should report distinct counts. Commands: `python3 -I zparse.py`, plus the inline count in the log.

3. **[fix] The real place.** §2 says "The real place contributes 1"; ct08 never computes that term. The claim is true only for A > 0.
   - Test: the same curve with the other rational point (−A, 0) of the first conic, i.e. `ct08.pairing(-A, C, ...)`. This is wrong on exactly the 323 of 3,912 triples with δ < 0 and b_m < 0.
   - Example: at (5623,1017) it pairs the 2-torsion image (1, −1, −1) with b = (1, −1) to −1.
   - By hand: at ∞ with δ < 0, u small, the product is (−1, b_*)(−1, b_m) = 1 if A > 0, and (−1, b_m)_∞ if A < 0.
   - All 201 pilot ratios have m > n, so the pilot is unaffected. State A > 0 (m > n), or include the term at ∞.

4. **[fix] gcd(m, n) = 1 is never stated.**
   - For coprime inputs, the odd-p table matches the local term (my own local points, with constancy checked) on all 9,412 pairs (a, b) tried. This covers every a = (1, δ, δ) and every b in Im_p, at every odd p | ABC, for E_08 and E_LU.
   - Inputs: 66 coprime (m, n) (60 random with m < 40, all parities, and pilot ratios), plus 15 with m < n. All four cases are covered; p = 3 is included.
   - The two solubility claims (p | A ⇒ p ∤ δ and ±δ square; p | C ⇒ δ square) hold on all 791 cases.
   - Non-coprime inputs ((6,3), (9,3), …, with p | gcd): "p | C − A ⇒ 1" fails on 34 of 156, and the solubility claims fail on 28.

   Command: `python3 -I zlocal.py` (output in `zlocal.txt`).

5. **[note] The table by hand.** I re-derived cases 3 and 4 (δ square: u near ∞; −δ square: u near 0) and the b_r = 1 product formula.
   - Both derivations hold at every odd p, not only away from C ± A.
   - The choice of branch in case 4 changes the term by (−1, b_r)_p. That is 1 on Sel², because (1, −1, −1) is a torsion image and the local images are isotropic.
   - I did not re-derive cases 1 and 2 (p | C ∓ A) by hand. They are confirmed only numerically (finding 4).
   - The b_r = 1 formula matches rank_second on 8,130 + 7,012 triples (data at the time of the check). Command: `python3 -I brone.py`.

6. **[fix] The "free bit" is not free.**
   - Theory: on the trivial covering, the fibre of "forget u_i" is {Q, T_i − Q}, so F_2 − F_3 ~ T_1 and F_1 − F_3 ~ T_2. Writing b = b'_1⊗T_1 + b'_2⊗T_2 gives Kummer components b_1 = b'_2 and b_2 = b'_1. The Brauer class is then ∏_i (b_i, l_i), where l_i is the tangent of the conic that omits u_i.
   - That is exactly rank_second's `(b2, l_v/l_w)(b1, l_u/l_w)`. On these curves it is ct08's conv=0, the assignment chosen.
   - The derivation is mine. I did not fetch Cassels 1998 or FSS 2010 to compare their index conventions.
   - Data-free test: my evaluator `zct.py` tries all 6 assignments, with all 8 sign choices, the real place included, and constancy checked at every place.
     - The chosen assignment agrees with rank_second on all 9,350 triples of the first 100 curves and is locally constant throughout.
     - The other assignment (conv=1) coincides with it when b_r = 1. When b_r ≠ 1 it is **not locally constant** on 795 of 1,367 (E_08) and 765 of 992 (E_LU) triples, and it disagrees on 32 + 34 of the rest.
     - So the other assignment is not a well-defined pairing, and it would not have looked "consistent" on the non-kill set.
   - Also, "fails on 222 lines" is E_08 only: 222 of 1,280 lines, which is 15 of 133 distinct. On E_LU it is 112 of 874 lines (15 of 111 distinct).

   Commands: `python3 -I convtest.py conv_part.txt 1,3 part0.jsonl part1.jsonl`; `python3 -I claims.py`.

7. **[note] The 2-adic term is robust.**
   - With 6 more seeds (2, 3, 4, 5, 11, 12345) on 2,826 triples, the global value and the 2-adic term never change, and nothing returns None.
   - PARI's `hilbert` raises "precision too low" on an under-precise p-adic rather than returning a sign, so there is no silent wrong sign.
   - The `u + p^40` line is dead code: `local_point` never returns ±A or ±C, because `issq(0)` is False.
   - The `um ± us` switch is a legitimate change of point (u_* → −u_*).

   Command: `python3 -I seedtest.py part0.jsonl part1.jsonl`.

8. **[note] The other numbers reproduce:** 1,280; 874; 668 of 1,280 with b_r = 1 (83 of 133 distinct); 12,984 = 2,644 + 2,972 + 4,906 + 2,462; 8,854 = 1,676 + 2,528 + 3,210 + 1,440; 21,838; 2-adic −1 on 156 + 122 lines; controls 150 × 5, all as stated; 8,930 = 2,154 + 6,776.

   Commands: `cd pf_run; python3 check_local.py 08|LU; python3 controls.py`.

9. **[fix] Literature.**
   - **arXiv:2303.05059 (S. Zhang).** The abstract requires "$e_1+e_2+e_3=0$". With e1 = e2 = 1 this forces e3 = −2 and a² + b² = 2c². The second curve then has root differences a², b², 2c². It is not a Frey curve of a Pythagorean triple (all three differences squares), and no twist of it is. So "With e1 = e2 = 1 the second curve is our Frey curve" is false. Its other hypotheses (minimal Sel², twists by n with residue conditions) do not fit either.
   - **Naskręcki, arXiv:1210.6933.** The family is y² = x(x − a²)(x − b²) with **legs** a, b, which is a different family from E_08/E_LU (a leg and the hypotenuse). It is related at most by isogeny, which preserves rank but not Sel² or the pairing.
   - **Fisher, RNT 2022** (arXiv:2208.14977) and **Shukla–Stoll** (2302.01640) are correctly described.
   - One web search found no published version of this formula. That is not conclusive. The result is a direct specialisation of Cassels' formula; what is new is the four local cases.

10. **[fix] The §4 hedge does not reach the opening paragraph.** The "Result in one paragraph" says, as fact, "For the four mixed curves … there is no such formula", and §4 says "therefore of Rédei-symbol type". Only the last sentence of §4 hedges. A conic that depends on the tuple does not rule out a closed form (Fisher 2022 solves no conics). Reword both as expectations.

11. **[note] Status words.** "It reproduces all 2,154 kill lines" is true but means 244 distinct tests. §3's own "Not tested" paragraph is now answered by finding 1, though only on 146 + 133 of the 201 ratios. "Closed form at every odd prime" is supported numerically, not proved (cases 1 and 2 have no hand derivation in the doc or here).
