# Findings to share (26 Nov 2026)

An AI model (Claude) wrote this note, and no person has reviewed it. Nothing has been committed or pushed to the repository. None of it touches the v6 `lean_rank` work (Lemma B, Theorems 1 and 2, Lemma D/E1 over ℚ(√2), Stoll's x−T map). Labels: **proved**, **computed (exact)**, **computed**, **conditional**, **open**.

## 1. Census beyond m = 1000 (tier D, computed; not in Lean)
- **1000 < m ≤ 9000:**
  - every coprime ratio is excluded except **(7735,4968)**, which is **open**, subject to: PARI's rank bounds throughout; Demjanenko–Manin generators saturated only at primes ≤ 1000; for nine route exclusions ((3321,2440), (5073,2336), (6588,805), (7811,3060), (8364,4123), (8687,4416), and through partners (5761,881), (7393,5783), (7409,2737)) the conditions (a), (b) of `CENSUS5000.md`; (6557,324), (6881,6233) on a numerical L-value;
  - seven ratios are excluded **conditionally on a rank being 1** (root number −1, PARI bounds [1,3]): (3540,2581), (4440,2231), (6121,959), (6671,2209), (6480,4961), (6903,5104), (8249,3720).
- **Pipeline:** congruences → the repo's class criterion → second descent → rank-0 quotients (`global.py`, unchanged) → Demjanenko–Manin → the rank-1 route. Write-ups: `census-wide/CENSUS5000.md`, `CENSUS7000.md`, `CENSUS9000.md`.
- **(7735,4968):** every usable curve has rank ≥ 2. A new idea that might close it is Chabauty or Demjanenko–Manin on its rank-1 quotients g4·K and w·g4·K; `global.py` uses only rank-0 quotients.
- **9000 < m ≤ 11000:** running. 6536 ratios pass the congruences and 612 are of square or twice-square type.

## 2. (627,364): Cassels–Tate pairing over ℚ(i) (computed)
Write-up: `open-ratio/WRITEUP.md`, last section. Code: `open-ratio/ctp/`.
- **New PARI implementation** of the Cassels–Tate pairing on Sel² for full 2-torsion over ℚ and ℚ(i).
  - A Hilbert-reciprocity shortcut avoids factoring huge norms. The argument is unreviewed but agrees with the direct method.
  - Validation (`open-ratio/ctp/validation/`, drivers plus logs):
    - over ℚ, 9 curves: the shortcut equals the direct method on all 9, and the pairing rank matches dim Sel² − 2 − rank on all 9 (six with pairing rank 2);
    - over ℚ(i), 8 curves: the pairing rank matches on all 8 (two with pairing rank 2).
- **E₁:** the pairing on Sel² is identically 0, with every entry computed twice. This independently confirms run 4.
- **E₄**, the curve with the smallest Sha in the class:
  - the torsion image in Sel² was identified, so everything turns on a single value;
  - it is 0 in two separate runs (`ctp/run627.out`, 8 of 16 entries, and `ctp/run627b.log`), each entry computed twice. A separate AI review also recomputed it with fresh local points: 0.
- **Reading** *(conditional on rank 0 and Sha finite)*: Sha(E₄)[2^∞] ⊇ (ℤ/4)². Sel² with its pairing cannot prove rank 0.
- **What would close it:** the pairing on Sel⁴(E₄). If |Sha(E₄)[2^∞]| = 16, that proves rank 0. It needs 4-coverings over ℚ(i).
- **Earlier in the day** *(computed, exact)*: cyclic 4-isogeny descent over ℚ(i) on all eight rational-kernel 4-isogenies gives bounds 3 to 5, none better than 2 (`open-ratio/four/`).
