# The pairing on `E_08` and `E_LU`: the local terms derived (run 11, 27 November 2026, UTC)

Written by an AI model (Claude); not read by a person. `PAIRING_FORMULA_2.md` and its scripts are a separate AI
research session's second note (Claude, Fable 5.2), committed unchanged apart from its banner line, a regenerated `shape.txt`
and one `sys.path` line removed from two scripts (below). Reviewed by a separate hostile AI agent (z13; `reviews/review_z13.md`, its checks in
`../../referee_z13/`), not by a person. No census count changes.

## What the note claims

For `E: y² = x(x − A²)(x − C²)`, `(A, B, C)` a primitive Pythagorean triple with `C` odd, and classes
`a = (1, δ, δ)`, `b` in `Sel²`: a closed formula (2) for Cassels' pairing with no local point and no conic, under
**H2** (if `A` is odd, 32 divides `C − A` or `C + A`; true for `E_08`, `E_LU` at every ratio with a 2-adic point),
with a proof of every local term given Cassels' theorem; the subgroup `{b ∈ Sel² : b_r = 1}` described from
factorisations alone (`S′`, the image of the Selmer group of the 2-isogeny with kernel `(0, 0)`), and the pairing on
it as a matrix of symbols (`Θ`, (3)). §5 and §7 are observations on tuples with trivial pairing and on the radical of
`Θ`.

## What the review found (z13)

- **The proofs hold, given Cassels' theorem**, with three gaps in the write-up, each closed by the reviewer:
  1. the good primes `p ≡ 1 (mod 4)` with `δ` a non-square unit need an argument ("all entries are units for a
     suitable point"): the Hasse bound gives a point for `p ≥ 13`, and `p = 5` never occurs because 5 divides `ABC`
     for every Pythagorean triple;
  2. **`δ` must be squarefree** in (2): the terms depend on the representative (at `E_08`, (535,23), `δ = 1` written
     as `7²` gives −1);
  3. in (iv) the threshold is `v(u + A) ≤ v(C + A) − 1`, resp. `− 3` at 2.
- **H2 is a real hypothesis**: outside it (2) fails on half the pairs of every curve tested (16 control curves, e.g.
  `E_LU` at (725,497)); `S′` and `Θ` still agreed there, unexplained.
- **Independent test, outside the data the session used** (the first review's tables): 362 curves from 181 ratios
  not in the pilot, within H2: (2) equals `rank_second`'s pairing on **194,504 of 194,504** pairs (16,296 of them
  `rank_second`'s own values on a basis, the rest by multiplicativity in `b`); `S′ = {b_r = 1}` on 362 of 362;
  `Θ` equals the pairing on 96,632 of 96,632. `A` odd on 168 distinct curves, `v₂(C ∓ A) = 5` exactly on 86;
  `dim Sel²` 2 to 8; square and twice-square ratios included.
- **Every number of the note reproduces.** Corrections: the 720 and 1,434 lines of §4 are pilot kill lines, 90 and
  154 distinct cases; "agree with the committed values" means the first review's stored `ct08` values; `shape.txt`
  lacked one line.
- **§5** follows in two lines from `S′`, (3) and the definition of the radical, so "not proved" understates it; but
  "candidate answer" overstates it, since it decides nothing where the radical is larger than `±1` (about two thirds
  of the ratios). A torsion class is the class of a rational point and is not removed by rank 0. (3320,2303) is the
  partner of (5623,1017), so the pilot has one example, seen twice, not two.
- **§7**: "the reason is exact" holds for the square-type statement once one adds why `2mn` is not `±` a square
  (`a⁴ + 4b⁴ = □` has no solution); the non-degenerate part is an explanation, not a proof, and it needs
  `(1, 2mn, 2mn) ∈ S′` (true on 32 of 93 of the reviewer's opposite-parity `E_08` curves, in the radical on 16 of
  them). `(2√β, S)_p` is `(2/p)` times the quartic symbol.
- **References**: Fisher, *Higher descents on an elliptic curve with a rational 2-torsion point*, Math. Comp. 86
  (2017), 2493–2518, Theorem 2.1 (the pairing on the isogeny Selmer group, whose kernel is the next descent); the
  Monsky analogy is fair, but the matrices have 2-adic entries too. A separate prior-art search (run 11) found no
  source for this formula; the nearest work is on the sibling family `a² + b² = 2c²` (Z. Wang, arXiv:1703.05909;
  Z. Wang and S. Zhang, arXiv:2303.05058), and Y. Mokrani (M.Sc., Montréal 2020) gives Monsky matrices for any curve
  in Legendre form.

## Revision 3 (`PAIRING_FORMULA_3_2026-11-27.md`) and its check (z13)

The session revised the note after the review. Review z13 checked the revision (`reviews/review_z13.md`, section
"Revision 3"; `../../referee_z13/h2_rerun_z13.py`):

- Applied correctly: (1) the good-prime step, now with no bound on `p` (Hasse–Weil on the smooth genus-1 reduction and
  four exclusions at every `𝔽_p`-point; better than the reviewer's version); (2) `δ` and the components squarefree;
  (3) the threshold in (iv); (4) H2 necessary, with the counterexample (32 of 64 pairs wrong at `E_LU` of (725,497);
  recomputed); (6) §5 rewritten; (7) §7's Fermat line and `(2/p)` times the quartic symbol; (8) Fisher 2017 in full.
- (5) The citation of Fisher–Schaefer–Stoll 2010 (§3, §5, Theorem 5.2) is right and so is the conclusion, but one
  sentence is false as written: under FSS's normalisation the function `f_j` is not the tangent form of the conic
  that omits `u_j`; FSS's triple is `(l_u, l_v, l_w)` times the common factor `l_u l_v l_w`, which drops out since
  `b_1 b_2 b_3` is a square, so `∏ (l_j, b_j)` is FSS's pairing and formula (1) holds. "The only one of the three
  defined over `K(θ_j)`" fails for cyclic cubics; the intrinsic reason is that the conic omitting `u_j` is the
  quotient by `R ↦ T_j − R`.
- New error: the Prior art paragraph's "correction" is wrong. arXiv:2303.05059 is S. Zhang, *On a comparison of
  Cassels pairings of different elliptic curves* (as the first note said); arXiv:2303.05058 is Z. Wang and S. Zhang,
  *On the quadratic twist of elliptic curves with full 2-torsion*. Both exist; the sentence calling 2303.05058 "the
  correct one" should be dropped. The prior-art search was the conductor's, not the reviewer's.
- Not applied from the first review: §4 still counts kill lines (720 lines are 90 distinct cases, 1,434 are 154) and
  says "committed values" for the first review's stored `ct08` values; at non-square ratios §7 should say that
  `(1, 2mn, 2mn)` must lie in `S′` before (2) applies to it; "matrices of Legendre symbols" omits the 2-adic entries.
- `h2_check.py` does not run as committed (it imports `pipeline` from the session's folders).

## Tier

Formula (2) and the description of `S′` and `Θ` are **tier C** (a written proof here, reviewed by a separate AI
agent, read by no person), resting on: Cassels' formula (1998) with the assignment of the components of `b` to the
three tangent forms (derived by review z12, matched with `rank_second`, and matched with Fisher–Schaefer–Stoll 2010,
§3 and Theorem 5.2, by review z13 on revision 3; Cassels 1998 itself was not opened by anyone); Cassels' local constancy; Fisher–Schaefer–Stoll 2010 (Cassels' pairing
is the Cassels–Tate pairing); the standard sizes of the local images. With the three fixes above. §5 and §7 are tier D.

## Rerun

The scripts read the first review's tables or the pilot's files from the current folder and import `parse.py`,
`ct08.py` of `../`:

    cd experiments/referee_z12 && export PYTHONPATH=../research-pairing/note2:../research-pairing
    python3 ../research-pairing/note2/closed_formula.py | diff - ../research-pairing/note2/closed_formula.txt   # 67 s
    python3 ../research-pairing/note2/isogeny_part.py   | diff - ../research-pairing/note2/isogeny_part.txt     # 25 s
    python3 ../research-pairing/note2/shape.py          | diff - ../research-pairing/note2/shape.txt            # 4 s
    python3 ../research-pairing/note2/radical.py        | diff - ../research-pairing/note2/radical.txt          # 16 s
    cd ../conductor-t10
    python3 ../research-pairing/note2/two.py  | diff - ../research-pairing/note2/two.txt                         # 3 s
    python3 ../research-pairing/note2/two2.py | diff - ../research-pairing/note2/two2.txt                        # 2 s

(the line `sys.path.insert(0,'/tmp/t10')` of `closed_formula.py` and `shape.py` was removed, so that no folder outside the repository can supply their imports; that is the only change to the session's scripts). The review's checks: `../../referee_z13/`.
