# Hostile review of sess_t5 (scripts/open_routes/routes.py, reports/t5.md), focused on the SECOND_DESCENT §6 label claim

Reviewer: slice v5, run 3. Written by an AI model (Claude); not read by a person. Script:
`v5/scripts/second_descent/review_t5_labels.py` (sympy only; output `review_t5_labels.txt`), written from
`kclose.factors` and SECOND_DESCENT.md §2 without using routes.py.

## Verdict: t5's claim is CORRECT. SECOND_DESCENT.md §6 mislabels the four bound-2 quartics.

What was checked (exact, sympy):
1. With `kclose.factors` (`M1 = ps − p̄, M2 = p̄s + p, M3 = ps + p̄, M4 = p̄s − p`, `p = m + in`) and §2's
   `q₊ = Ss² − 2Ds + S`, `q₋ = Ss² + 2Ds + S` (`D = m² − n²`): `M1·M4 = q₊` and `M2·M3 = q₋` identically in `m, n, s`;
   `M1·M2 = Ss² + 4imn·s − S`, which equals `−GL/2` at `s = −it`; no other product of two M's is `q₊` or `q₋`.
2. `kdist_627_364.txt` (the table §6 quotes): the four quartics with bound 2 are `L1L3M1M2, L1L3M3M4, L2L4M1M2,
   L2L4M3M4`. The quartics §6 names, `q₊L1L3 = L1L3M1M4`, `q₋L1L3 = L1L3M2M3`, `q₊L2L4 = L2L4M1M4`,
   `q₋L2L4 = L2L4M2M3`, all have bound **11** in the same table.
3. Of the 210 quartics, 46 are rational (up to a constant) in `s` and 10 in `t`; the four bound-2 quartics are in
   neither set. Stronger than t5 says: their j-invariants (from the quartic invariants I, J, exact in ℚ(i)) are
   **not rational** (`Re j ≈ 963.59`, `Im j ≈ ±928.21`), so their Jacobians are not base changes of curves over ℚ in
   any coordinate (t5 says "not base changes in the s- or t-structure" and, in its log, "j not rational"; confirmed).
   Negative control: the same test gives j ∈ ℚ for `q₊L1L3` etc.
4. New detail: the four bound-2 curves are one curve up to isomorphism and conjugation: `s ↦ −s` sends
   `L1L3M1M2` to `L2L4M3M4` (`L1 ↔ L2, L3 ↔ L4, M1 → −M3, M2 → −M4`), and complex conjugation of the coefficients sends
   `L1L3M1M2` to `L1L3M3M4` (`M1 ↔ M4, M2 ↔ M3`); the j-invariants pair up accordingly (j, j̄, j̄, j). So "rank 0 for
   one of them" is the same as for all four (over ℚ(i) conjugate curves have equal rank).

Consequence (t5 says this too, correctly): §6's conclusion "as far as 2-descent sees, every route through a single
curve of genus 1 is used up, over ℚ and over ℚ(i)" does not follow. Under the mislabel the bound-2 curves were base
changes of `E₁` (rank 1 over ℚ), so their ℚ(i)-rank was ≥ 1; in fact they are a different curve over ℚ(i), with
2-Selmer bound 2, so rank 0 or 2 over ℚ(i) (PARI root number +1 over ℚ(i), a guide, not rechecked here; even rank also
follows if Sha is finite), and rank 0 would close (627,364) outright by the torsion check of kclose.py. Theorem A of
LAST_PAIR.md is not affected: it defines `Q₁ = q₊(m − ns)(n + ms)`, `Q₂ = q₊(m + ns)(n − ms)` by formula, over ℚ.

## Points against t5's proposed doc text
- "So these four quartics are `L1·L3·GL`-type products": mixes the variables (L1L3 is in s, GL in t). Say instead
  "`L1·L3·M1·M2` with `M1·M2 = −GL/2` at `s = −it`".
- "(the 2-descent over ℚ(i) is weak on base changes)": an explanation, not checked. What is checked is only that
  `q₊L1L3` has bound 11 over ℚ(i) while §6 reports ranks 1 and 2 for the curve and its twist by −1 over ℚ. A large gap
  between the ℚ(i)-Selmer group of a base change and the two ℚ-Selmer groups is possible (norm-index terms at the
  many split bad primes, Kramer's formula), but nobody has computed it here; drop the parenthesis or mark it a guess.
- "with `j ∉ ℚ` and root number `+1` (PARI)": j ∉ ℚ is now exact (this review); the root number is PARI only.
- Add that the four are one curve up to `s ↦ −s` and conjugation (item 4).

## Proposed replacement for the SECOND_DESCENT.md §6 sentence
"The four with bound 2 are `L1·L3·M1·M2`, `L1·L3·M3·M4`, `L2·L4·M1·M2`, `L2·L4·M3·M4` (labels of `kclose.py`, where
`q₊ = M1·M4` and `q₋ = M2·M3`; an earlier version of this paragraph called them `q±·L1·L3`, `q±·L2·L4`, which is
wrong: those have bound 11). They are one curve up to `s ↦ −s` and complex conjugation; its j-invariant is not
rational (`review_t5_labels.py`), so it is not defined over ℚ in any coordinate, and its bound 2 means rank 0 or 2 over
`ℚ(i)` (root number +1 by PARI). Rank 0 would close (627,364) by the torsion check of `kclose.py`; 2-descent cannot
decide it. Theorem A of `LAST_PAIR.md` uses `q₊·L1·L3` and `q₊·L2·L4` over ℚ and is not affected." (Found by t5;
checked by an independent AI review, 25 Nov 2026; not read by a person.)

## Not reviewed
The rest of routes.py (rank ranges, root numbers, the u-structure quotients, the K-quartic distribution for the other
ratios) was not re-run; spot checks only: its rationality criterion for s (line 432: the M's must come as `M1M4` and/or
`M2M3`) gives exactly the 46 quartics found here; its parity argument (root number +1 and rank ≥ 1 ⇒ 2^∞-Selmer corank
≥ 2, Monsky; Dokchitser–Dokchitser) is the standard p-parity theorem and is used correctly; "no corner form splits over
ℚ(√2)" is right (as forms in t, F0 and F8 have discriminant −4m²n², GL and GU −16(m² − n²)²: all four split over
ℚ(i) and over no other quadratic field).
