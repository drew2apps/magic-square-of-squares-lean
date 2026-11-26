# (627,364): rank of the curve over ℚ(i), experiment log (26 Nov 2026)

Each result carries one of four labels:
- **proved**
- **computed (exact)**
- **computed (numerical)**
- **conditional**

Nothing here has been reviewed by a person. The scripts and their outputs sit in this folder.

## Setup
- **Curve:** E₁ : δY² = Q(s) over K = ℚ(i), with Q = L1·L3·M1·M2 and the kclose labels.
- **δ = 1:** Q(−i) is a square in K. *(computed, exact)* (`build_curve.gp`)
- **Not a ℚ-curve:** the conductor is not Galois-stable, so this is not a base change. *(computed, exact)*
  - L(E/K, s) is a degree-4 L-function with conductor ≈ 10⁴¹, so analytic rank is out of reach.
  - Kolyvagin does not apply over K.
- **2-power isogeny class over K:** six curves E₁…E₆, with edges 1-2, 1-3, 1-4, 4-5 and 4-6. *(computed, exact)* (`isog.gp`)
  - E₁ and E₄ have torsion ℤ/4 × ℤ/2.
  - E₂, E₃, E₅ and E₆ have torsion ℤ/4.

## What held
1. **Full 2-descent on E₁ and E₄:** dim Sel² = 4 on each, so rank ≤ 2. *(computed, exact)* (`ksel_full2.py`, using the repo's `kdescent`)
   - Local images are found by sampling.
   - Sampling stops at the exact size |E(K_v)/2E(K_v)|.
2. **2-isogeny descent on every edge, both ways:** the bounds are 2, 4 or 8, so none beats 2. *(computed, exact)* (`phi_descent.py`)
   - The best edge is 1-4. Sel^φ on E₄ has dimension 3 (basis 991i, 145, 263), and the torsion gives the class 263·991i.
   - The open question is therefore whether the excess classes 145, 263 and 145·263 come from points: rank = 0, 1 or 2 of them. *(computed, exact)* (`edge14.py`)
3. **Isogeny invariance of the BSD quotient across the class** *(computed: Tamagawa numbers by PARI's Tate algorithm, periods to 60 digits; the ratios come out as exact powers of 2)* (`sha_ratio.gp`)
   - Tamagawa products are 2²², 2¹², 2¹⁶, 2²⁴, 2¹⁸ and 2¹⁸.
   - All six period areas are equal.
   - If rank 0 and Sha is finite, then |Sha(E_j)| / |Sha(E₁)| = **1, 256, 16, 1/4, 4, 4** for j = 1…6.
4. **Consequence** *(conditional on rank 0 and Sha finite; the deduction is standard)*:
   - Sha(E₁)[2] and Sha(E₄)[2] both have dimension 2.
   - So Sha(E₁)[2^∞] ≅ (ℤ/2^a)² and Sha(E₄)[2^∞] ≅ (ℤ/2^b)², with **b = a − 1**.
   - **E₄ has the smallest Sha in the class.**
   - If the Cassels–Tate pairing on Sel² vanishes on all six curves (run 4's claim), then a ≥ 3 and b ≥ 2.
   - A 4-descent with its pairing on E₄ proves rank 0 exactly when b = 2, that is when |Sha(E₄)[2^∞]| = 16.
   - This agrees with run 4's "wall" result (relayed to me, not read directly) and names the one curve it points to: E₄, the other curve with full 2-torsion, whose a-invariants are in `isog_627_364.gp`.

## Where it got stuck
- **Second φ-descent on edge 1-4** (`second_phi.py`): 4-coverings of the φ-coverings, with exact local solubility tested by residue-disc recursion.
  - For δ = 145, the local data at 2, 5, 29, 37, 181, 263 and 991 is complete.
  - It ran out of room at the split prime 23490241, where the covering has bad reduction and the discs mod p have to be enumerated.
  - I stopped there. By item 4, the excess classes are expected to survive this kind of descent, so it can't decide the rank.
- **Rank 2:** heights of generators are expected to be far too large for a point search (run 4's BSD estimate, not redone here).

## What would decide it
1. **A 4-descent with the Cassels–Tate pairing on E₄ over ℚ(i).** This is Fisher-style work and needs Magma or a large new implementation. If b = 2, rank is 0 and (627,364) closes, modulo the census's other inputs.
2. **Cyclic 4-isogeny descent: tried, does not decide it** (section below).
3. **A conditional decision.** The repo's CVV Mordell–Weil-sieve route already treats (627,364) conditionally.

## Cyclic 4-isogeny descent over ℚ(i) (added later on 26 Nov 2026)
Code: `four/fourdesc.py`, `four/run4.py`, data from `four/gen4.gp`; output `four/main4_seed1.txt`.
- **Method.** μ₄ ⊂ K, so a cyclic 4-isogeny ψ: A → B with a K-rational generator Q of its kernel, and its dual with generator T′, gives Kummer maps into K*/K*⁴. They are g_{T′} on B and g_Q on A, with div g_T = 4(T) − 4(O). Local images are found by sampling until their size meets the exact local product (16 at odd places, 256 at 1+i). At odd places it also has to equal 4·c(B)/c(A). The Selmer groups then come from linear algebra over ℤ/4, and r ≤ log₄(|Sel^ψ|·|Sel^ψ̂|/16).
- **Validation** *(computed, exact; the test curves and their runs are not in this folder)*: on six test curves of known rank over K (ranks 0, 1, 1, 2, 2, 3), the bound equals the true rank every time.
- **On the (627,364) class** *(computed, exact)*: all eight 4-isogenies with rational kernel generators were run, and every local image met its exact size.

| ψ | \|Sel^ψ\| | \|Sel^ψ̂\| | rank bound |
|---|---|---|---|
| E₁ → E₆ and E₁ → E₅ | 8 | 128 | 3 |
| E₃ → E₄ | 512 | 2 | 3 |
| E₂ → E₄ | 8192 | 2 | 5 |

The reverse directions give the same numbers. **None beats the 2-descent bound of 2.** So this descent sees less than the 2-descent here, which fits Sha(E₁)[2] being divisible by 2 (section "What held", item 4).
- **Conditional reading** *(conditional on rank 0)*: |Sha(E₁)[ψ]|·|Sha(E₆)[ψ̂]| = 64.
- **Side check on (609,200)/(809,409).** Both rest on rank 0 for one curve E₀/ℚ, for which only Cassels–Tate pairings give rank 0 (Sha[2] has dimension 4). E₀ has K-rank 1 when the twist's point is included, and I hoped a 4-isogeny descent over K could give the bound 1 as an independent check. It can't: no curve in E₀'s six-curve K-isogeny class has a rational point of order 4, so there is no cyclic 4-isogeny with a rational kernel generator *(computed, exact; `four/g609.gp`)*.

**Still needed for (627,364):** a genuine 4-descent with the Cassels–Tate pairing on E₄ over ℚ(i) (item 1 above). That needs Magma or a large new implementation; PARI's `ellrank` works only over ℚ, and E₄ is not a ℚ-curve, so restriction of scalars does not help.

## Cassels–Tate pairing on Sel² over ℚ(i) (added 26 Nov 2026, afternoon)
Code: `ctp/ctp.gp` (library), `ctp/run627b.gp` (E₄), `ctp/run627c.gp` (E₁); logs `ctp/run627b.log`, `ctp/run627c.log`.

- **Method.** The Cassels pairing on Sel² for curves with full 2-torsion, written as a sum of local Hilbert symbols (Cassels; Shukla–Stoll's formulation). For each pair of Selmer classes the program:
  - finds a K-point on each 2-covering's conics, using `rnfisnorm` plus a unit reduction;
  - samples a local point at every bad or Selmer-relevant place;
  - adds the Hilbert symbols.
- **New shortcut.** Good places outside the bad set are handled by Hilbert reciprocity (the `cassels2` routine). Their total equals Σ over the bad set of (Z_i, β′_i), where Z_i generates the denominator ideal of the covering's coefficients. This removes the need to factor 100-digit norms. *(The argument is the session's own. Its hypotheses: the local term is independent of the local point, and a point avoiding the three conic points modulo P exists when N(P) ≥ 13, which holds over ℚ(i) once 3, 5, 7 are added to the bad set and over ℚ once 11 is. A separate AI review read it and found it correct under these hypotheses; it agrees with the direct method on all 9 test curves over ℚ.)*
- **Validation** *(computed)*: drivers `ctp/validation/val_Q.gp` and `val_Qi.gp`, logs `val_Q.log` and `val_Qi.log`, run with the current code on 26 Nov.
  - **Over ℚ, 9 curves:** the direct method and the reciprocity shortcut give the same matrix on all 9. The pairing rank equals dim Sel² − 2 − rank on all 9, where PARI's `ellrank` pins the rank. Six curves have pairing rank 2 (Sha[2] of dimension 2) and three have rank 0.
  - **Over ℚ(i), 8 curves base-changed from ℚ:** the pairing rank equals dim Sel² − 2 − (rank E(ℚ) + rank E^(−1)(ℚ)) on all 8 (this check holds only when Ш has no element of order 4; [0,−22,−33] over ℚ(i), with analytic Ш of order 16 and pairing 0, is an exception a separate review found). Two curves, [0,195,−376] and [0,333,−238], have pairing rank 2 and six have rank 0.
  - **Did not finish:** [0,240,−283] over ℚ(i), where the conic solver ran out of time, so it is left out.
  - **Consistency check:** each entry was computed twice with different random local points. The per-place values differ, but the totals agree.
- **E₁** *(computed; local points by p-adic sampling, all symbols exact)*: all six entries ⟨k,l⟩ (k<l) on the Sel² basis are 0, each computed twice with no disagreement. So the pairing vanishes on Sel²(E₁), as run 4 found.
- **E₄** *(computed)*:
  - **Selmer basis:** S1…S4 (in `ctp/run627b.log`).
  - **Torsion image:** span{S3·S4, S1·S4}, found by testing squares at 300 split primes (`ctp/torsion_E4.gp`). That test is heuristic in the "is a square" direction.
  - **Deciding value:** the pairing vanishes on torsion, so everything turns on the single value ⟨S2,S4⟩ = ⟨S1,S2⟩ = ⟨S2,S3⟩.
  - **Evidence 1, `ctp/run627.out`:** a partial run of `run627.gp`, stopped by me after 8 of 16 entries, with no summary line and no E₁ part. ⟨1,1⟩ … ⟨2,4⟩ are each computed twice and are all 0, including the deciding ⟨1,2⟩, ⟨2,3⟩ and ⟨2,4⟩. *(I said earlier that this log was lost. That was wrong: the output flushed when the run stopped.)*
  - **Evidence 2, `ctp/run627b.log`:** a separate run of `run627b.gp` with the deciding entries first, unit reduction off and fresh random local points. ⟨1,2⟩, ⟨2,4⟩ and ⟨2,3⟩ are each 0 twice. The remaining check entries are still running; see the log for its current state.
  - `ctp/run627b_attempt1.log` is an earlier start of the same script, with unit reduction on. I stopped it after nearly two hours without a single entry, and it contains no pairing values.
- **Reading** *(conditional on E₄'s rank being 0 and its Sha finite; the deduction is standard)*:
  - A zero pairing on E₄ means Sha(E₄)[2] ⊂ 2·Sha(E₄), so b ≥ 2 and a ≥ 3 in the notation of "What held", item 4.
  - Sel² with its pairing cannot prove rank 0. This independently confirms run 4's wall, with a separate implementation.
- **Still needed:** a pairing on Sel⁴(E₄), i.e. 4-coverings of E₄ over ℚ(i) and their Cassels–Tate pairing. If b = 2, that closes (627,364). It is a large step beyond this code: 4-coverings are intersections of two quadrics, so it needs a Merriman–Siksek–Smart second descent over ℚ(i) and then a pairing on those.
