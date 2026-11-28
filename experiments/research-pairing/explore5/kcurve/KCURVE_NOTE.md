# Corner K-curve for (627,364) over K = Q(i): point search, 27 Nov 2026

> **Read `../NOTES.md` first.** This is the session's note, committed unchanged apart from this line; where they differ, `../NOTES.md` (after review) holds.


Written by an AI model (Claude); unreviewed. Time box about 45 minutes. **Result: the rank is not decided.**
No non-torsion point was found, so rank 2 is not proved; nothing here is evidence for rank 0 either.

## Model (checked by exact arithmetic in `selmer.py`, assertion on the polynomial identity)

Quartic: `y^2 = Q(s) = (m - n s)(n + m s)(p s - pb)(pb s + p)`, `m=627, n=364, p=m+in, pb=m-in`
(`Q(i)` and `Q(-i)` are squares, so no twist is needed). With `r = -p/pb` (root of M2), `A = Q'(r)`,
`s = r + A/x`, `x = (c^2 X + shift)/den^2`, `c = den = 525625`,
`shift = -2783595387469704535000 + 28543239158509955880000 i`:

    E : Y^2 = X (X - e2)(X - e3),   e2 = 273990441250,   e3 = -102929464375 + 118967496648 i

`j(E)` equals the `j` of `ellfromeqn(y^2 - Q)` (PARI); conductor norm 5683468881658146768135645923544141966400;
`E(K)_tors = Z/4 x Z/2`.

2-isogenous curve by `T3 = (e3,0)`:

    E' : y'^2 = x'(x'^2 + A1 x' + e2^2),  A1 = 2(e2 - 2 e3) = 959698740000 - 475869986592 i
    dual isogeny E' -> E:  X = e3 + (y'/(2x'))^2

`E'` also has full 2-torsion over K (other roots `-63145210128 - 36055447954 i`, `-896553529872 + 511925434546 i`),
torsion `Z/4 x Z/2`.

## 2-Selmer group (code: r16 `kdescent.py`, experiment code with sampled local images; not re-derived here)

`dim Sel^2(E/K) = 4` (so rank <= 2, as recorded before); the 16 elements `(b1,b2,b3)` are in `selmer.log`.
The two points of order 4 found by the search lie in classes 11 and 2, so the torsion classes are {0,2,9,11}.
Observation from the table: each of the three non-torsion cosets has representatives of the form `(d,d,1)`:
`d = 991i, 263` (coset A), `d = 145, 37791785i` (coset B), `d = 38135, 143695i` (coset C).
For such a class a point of E(K) is the image of a point of E'(K) with `x' = d * square` (shown by hand:
`X - e1 = x' * square` under the dual isogeny). So if the rank is 2, generators modulo 2E(K)+torsion come from E'(K),
and searching on E' is cheaper (heights halve). `dim Sel^2(E'/K) = 4` as well (`selmer_Eprime.log`).

## Searches (C sieve `search.c` / `search2.c`, 80 degree-one places of Z[i] above p = 13..., p = 1 mod 4; exact check in PARI)

Notation: `U = a + b i`, `W = c + d i`, box H means `0 <= a <= H, |b| <= H, 1 <= c <= H, 0 <= d <= H`
(all `U/W` up to sign with a representative in the box).

1. On E, all 16 Selmer elements and all three coordinates: `X = e_j + b_j U^2/W^2` for the 34 distinct pairs `(j, b_j)`
   (`tasks.txt`), **box H = 100**. Hits (`hits_H100.txt`, checked by `verify.py`): only the two x-coordinates of the
   points of order 4, `X = 105422615561 - 18027723977 i` and `X = -311281544311 + 255962717273 i`
   (canonical height 0). No non-torsion point.
2. On E', the six non-torsion classes `x' = d U^2/W^2`, `d in {991i, 145, 263, 38135, 143695i, 37791785i}`,
   **box H = 280** (`isohits_H280.txt`): no hit at all. (The two torsion classes d = 1, 260633i were run only to H = 60: no hit; that log was not kept.)

So: no point P' of E'(K) with `x'(P') = d (U/W)^2`, d as listed, `U/W` of box height <= 280; and no point of E(K) with
`X - e_j = b_j (U/W)^2`, box height <= 100. These are statements about naive coordinates only; no bound on the
canonical height of a generator was derived (the height-estimate step was not done).

Test of the sieve: it finds the points of order 4 (both binaries, including the generalised `search2` on the same
task). It was not tested on a planted non-torsion point of E'.

## Not done / where it stops

- Simon's `ell.gp` (shipped with passagemath) loads in cypari2 but `bnfellrank` stops with a PARI 2.17 type error
  ("forbidden comparison t_INT, t_POL"); not debugged. No `gp`, no Magma.
- Full 2-coverings as reduced quartics (needs conics over K solved and minimisation/reduction over K), 4-coverings,
  and a third isogeny step from E' (possible: E' has full 2-torsion; its Selmer table is in `selmer_Eprime.log`) were not done.
  These are the next steps; each roughly halves the height to be searched.
- Structure remark (by hand, not machine-checked): in `t = i s`, `M1*M2` is, up to sign, the rational form L(t), so
  `Q = L1*L3*L(t)`; the fibre product of E and its conjugate over the t-line is a genus-3 curve over Q whose Jacobian is
  isogenous to `Res_{K/Q}(E)` times the curve `y^2 = F0*F8`. Since E is not a Q-curve, `Res E` is not a product of curves
  over Q, so no elliptic curve over Q supplies points; this route gave nothing.
- A proof of rank 0 would need a second (in fact, by the earlier note, a third) descent over K; out of reach here.

## Files

`selmer.py`, `selmer.log`, `selmer.json`, `selmer_Eprime.py`, `selmer_Eprime.log`, `selmer_Eprime.json`, `tasks.py`, `tasks.txt`,
`tasks_iso.txt`, `tasks_iso6.txt`, `search.c`, `search2.c`, `run.sh`, `run_iso.sh`, `verify.py`, `verify_iso.py`,
`hits_H100.txt`, `isohits_H280.txt`, `ell.gp`/`ellcommon.gp`/`ellQ.gp`/`qfsolve.gp` (copies of Simon's scripts, unused).
