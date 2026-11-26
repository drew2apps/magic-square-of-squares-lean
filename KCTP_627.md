# The corner K-curves of (627,364): the Cassels–Tate pairing over ℚ(i) (experiment)

**Status.** Experiment, 25 November 2026. Not counted in the census, and it changes no count. Written by an AI model (Claude); not read by a person, not reviewed. Script: `scripts/second_descent/kctp.py` (exact arithmetic; PARI and Sage only as search tools or as guides, as marked).

## The question

`kclose.py` attaches to each of the four corner quartics `L1·L3·M1·M2`, `L1·L3·M3·M4`, `L2·L4·M1·M2`, `L2·L4·M3·M4` of (627,364) (labels of `kclose.py`, `s = −it`) a curve `E_J: y² = (x − e1)(x − e2)(x − e3)` over `K = ℚ(i)`. A solution of the ratio gives a point of `E_J(K)`. An earlier AI slice (25 November 2026) found that the 2-Selmer bound of these curves over `K` is 2, that `j ∉ ℚ`, and that PARI's root number is +1, so the rank over `K` is at most 2, and 0 or 2 if Ш is finite (parity). Rank 0 for one of them would close (627,364) by the finite torsion check of `kclose.py`. The question was which of the two it is.

## What was computed (`kctp.py`, output `kctp.txt`)

- **The curves.** The four curves are rebuilt in pure Python from `m, n`. The script checks that `Q_J(±i)` are squares in `K`, so the class of `kclose.py` is 1, and that the Weierstrass model is exactly `(x − e1)(x − e2)(x − e3)`. It also checks, exactly, that `j ∉ ℚ`. Two exact isomorphisms over `K` hold: `L2L4M3M4 ≅ L1L3M1M2` and `L2L4M1M2 ≅ L1L3M3M4`. In addition, `L1L3M3M4` is isomorphic to the complex conjugate of `L1L3M1M2`. So the four curves have one rank between them. The labels are checked too: `q₊ = M1·M4` and `q₋ = M2·M3`, so `L1·L3·M1·M2` is not `q₊·L1·L3`. This supports the remark of that slice on `SECOND_DESCENT.md` §6.
- **Sel²(E/K)**, by the complete 2-descent. It reuses the local images and class labels of `kdescent.py`, now with the Selmer elements written out explicitly. Each curve has `dim Sel² = 4`. The torsion image has dimension 2: it is spanned by `T1`, `T2` and the points of order 4 over the degenerate solutions `s = ∓i`. The points over `s = ∓i` have order 4, and twice each of them is `T3`; this is checked exactly with the group law. Sage gives `E(K)_tors = ℤ/2 × ℤ/4`, as a guide.
- **The Cassels–Tate pairing on Sel².** It uses Cassels' conic formula `⟨b, D_a⟩ = Σ_v (b₂, l_v l_w)_v + (b₁, l_u l_w)_v`. This is the formula of `scripts/g1_full/rank_second.py`, ported from ℚ to `K`:
  - points on the conics over `ℤ[i]` come from a lattice search and are checked exactly;
  - the tangent lines are primitive;
  - local points on `D_a` come from Hensel lifting in `ℤ_p`, `ℤ_p[i]` or `ℤ₂[i]`, with the precision asserted;
  - Hilbert symbols at odd places use the tame formula, and at `1 + i` a 4×4 form obtained from the product formula (checked against PARI's `nfhilbert` and the product formula on random pairs).

  The places that can contribute are those above 2 and above the primes of `a`, `b` and the `e_i − e_j`. A purity argument for the other places is written in the script header. It is not proved by a program, but the script checks it at the 15 places above 3, 7, …, 61: there the local term is 0, i.e. `v(l_v l_w(P))` is even.

  The matrix is also recomputed with the tangent lines taken at second rational points of all twelve conics. As a second computation that uses no Hilbert symbols, the script solves an affine system over `F₂` for each basis class `a` of Sel². It looks for `λ ∈ K(S,2)²` that agree, at every place of `S`, with the classes of `(l_v l_w, l_u l_w)` at local points of `D_a`; such a `λ` exists exactly when `D_a` lifts to an everywhere locally soluble 4-covering. On the controls, a class lifts exactly when its row of the pairing is zero (checked on every curve).
- **Result.** On all four curves the pairing matrix on Sel² is identically zero, with either choice of conic points, and every basis class lifts to an everywhere locally soluble 4-covering. It is symmetric and alternating, and the torsion images lie in its kernel. The same holds on the 2-isogenous curve `E/⟨T3⟩`. It is the only curve 2-isogenous to `E` with all of its 2-torsion over `K`, and the curves 2-isogenous to `E/⟨T3⟩` other than `E` lack full 2-torsion over `K` (checked). So this computation gives **rank E_J(K) ≤ 2 again.** The rank is 0 or 2 if `Sha(E_J/K)[2^∞]` is finite (the pairing is alternating), which agrees with PARI's root number +1. **It does not decide between 0 and 2.** If the rank is 0, then `Sha(E_J/K)[2] = (ℤ/2)²` lies in `2·Sha[4]`. **(627,364) is not closed.**
- **Route (b), points.** For each of the three non-torsion cosets of Sel², the script searches the 2-coverings in the models `L_j(σ, τ) = β_j y_j²` over `ℤ[i]`. Every admissible `β₄` is tried, with `|Re y|, |Im y| ≤ 60` (numpy sieve, then an exact test; `kctp_search.txt`). No point of infinite order was found. The same models find all eight torsion points at `|y| ≤ 4`. Sage's `simon_two_descent` over `K` (a guide) gave no answer in 7 CPU-minutes. PARI's conductor norm (a guide) is about 5.7·10³⁹. Point search is therefore not expected to reach generators, if they exist. Analytic methods are out of reach.

## Checks (outputs `kctp_controls.txt`, `kctp_negative.txt`)

- **Base changes of `y² = x(x + a)(x + b)` over ℚ.** The rank over `K` is the rank of `E` plus the rank of `E⁻¹`, taken from PARI (a guide). The bound from Sel² and the pairing is never below that rank. It closes the gap on 16 of 17 curves with excess 2. The pairing is zero on all 12 curves whose Selmer bound equals the rank.
- **Curves over `K` with `j ∉ ℚ`, `dim Sel² = 4`.** Lower bounds come from points found and checked exactly. On every curve the bound is sound, and the images of the points pair to zero with every `D_a`. The pairing has rank 2 (bound 0) on 6 of 10. It is zero on 4 of 10; on 3 of those, two independent points exist, so the rank is exactly 2.
- **`--negative`.** 29 forged inputs, all refused:
  - a point off a conic;
  - a forged conic hint;
  - a non-tangent line (the local term stops being constant);
  - each class `(π, 1)` outside Sel² (the local term is not constant), 15 of them;
  - classes outside Sel² used as 2-coverings (no local points);
  - an incomplete factorisation hint;
  - a forged Selmer element;
  - the 4-covering lift with a non-tangent line;
  - a lift claimed for a control class whose pairing row is non-zero;
  - a flipped pairing entry.

## What it rests on, and what is not done

- **Theory quoted, not proved here:**
  - the Cassels–Tate pairing and its Brauer–Manin description;
  - constancy of the local term for Selmer classes (checked only at 6 points × 8 signs per place);
  - the purity argument for places outside `S` (checked only at 15 places);
  - the parity statement behind "0 or 2", which also assumes that `Sha[2^∞]` is finite. Unconditionally the result is only rank ≤ 2.

  The arithmetic is exact.
- **Not done:**
  - a descent deeper than Sel² with its pairing: a 4-descent, or a descent by the cyclic 4-isogeny with kernel `⟨P4⟩`. Since `μ₄ ⊂ K`, that descent is Kummer theory over `K`, and it is the natural next step;
  - Sel² of the curves `E/⟨T1⟩` and `E/⟨T2⟩`, which need arithmetic over a quadratic extension of `K`;
  - the torsion check of `kclose.py`, which is not needed, since rank 0 is not certified.

## Reproduce

    cd scripts/second_descent
    python3 kctp.py                # kctp.txt, about 3 minutes
    python3 kctp.py --controls     # kctp_controls.txt (PARI for the ranks over Q, a guide)
    python3 kctp.py --negative     # kctp_negative.txt
    python3 kctp.py --search=60    # kctp_search.txt, about 4 minutes
    (add --write to rewrite the outputs and kctp_hints.json; --nopari uses only the hints: kctp.txt is then reproduced
    except the line of the PARI nfhilbert check)
