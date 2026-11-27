# H6: quotients of the (ℤ/2)⁵-cover for (627,364) and (991,263) (Roadmap H, run 7)

Written by an AI model (Claude) in a separate session; reviewed by a separate AI agent (z9, `reviews/review_z9.md`, its
own checks in `../../referee_z9/`); not read by a person. **Experiment, not counted; tier D of `LEDGER.md`.**
`REPORT_H6.md` is the session's report, unchanged; these notes take precedence.

- **The two ratios are partners, and their covers are isomorphic:** with `(m′, n′) = (m + n, m − n)` and
  `t′ = (t + 1)/(t − 1)`, `(t − 1)²·F′(t′)` is `L, U, 4F4, 16F0, 16F8` for `F′ = F0, F8, F4, L, U` (checked with sympy
  by review z9). So the equal bounds for the two ratios check that the code does not depend on the model; they are not
  new data. The report's "That the two covers are isomorphic is not proved here" is superseded.
- **The four pieces F0·F4·L, F0·F4·U, F8·F4·L, F8·F4·U are one curve A** (identities for `σ: t → −t` and `τ: t → 1/t`
  checked with sympy), so `J(cover) ~ A⁴ × B` by the quoted decomposition of the Jacobian of a `(ℤ/2)⁵`-cover into
  those of its hyperelliptic quotients (a written argument, not "proved"). Demjanenko–Manin through those four maps
  would need `rank J(A)(ℚ) ≤ 3`; a descent via a Richelot isogeny, with local images certified by the product formula
  of Tate local duality, gives 8 (review z9 re-derived the logic and the certification and reran it: identical). The
  twist of the dual curve is fixed by point counts over `F_p` at 8 to 12 good primes.
- **Genus-1 ranks are PARI's**, and some of PARI's lower bounds rest on parity (Ш finite), not on points: for
  (627,364), `y² = L·U` (the same curve as `y² = F0·F8` of (991,263)) has no known point of infinite order, so its
  rank is 1 only if Ш is finite (rank 0 would exclude the ratio outright, but by the 2-parity theorem for elliptic curves over ℚ,
  quoted, it would need Ш[2^∞] of odd corank, so Ш infinite); the class of `y² = F0·L`, with four maps from the cover, has three known independent
  points and rank 4 only if Ш is finite (Demjanenko–Manin through it would need rank 3, which by the same theorem
  needs Ш infinite; this remark is the conductor's; the root numbers +1 and −1 were checked with PARI's `ellrootno` by the pre-publication referee). With PARI's ranks, no Demjanenko–Manin or Chabauty route opens.
- Computed: genus-2 bounds in exact arithmetic on quoted theory (the Richelot construction, the Kummer map, local duality,
  unramified images at good odd primes; PARI's `nfislocalpower` decides local squares for quadratic points); genus-1
  ranks by PARI 2.19.1, and `h6_ranks.py 627 364` gives the same table under 2.17.4.

```
python3 h6_split.py 627 364         # 1 s
python3 h6_ranks.py --negative      # 1 s
python3 h6_richelot.py --negative   # 1 s
python3 h6_richelot.py --controls   # 14 s
```
