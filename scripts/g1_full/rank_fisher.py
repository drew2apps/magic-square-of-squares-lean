#!/usr/bin/env python3
"""Fisher's higher descents (arXiv:1509.03234) for the five curves whose rank bound 1 rests on Magma alone: how far
the descent by 2-isogenies gets in exact arithmetic, and where Magma's 1 must come from.

    python3 rank_fisher.py              the five curves (E+ of (760,39), (945,248), (72,65), (920,129); D1 of (631,71))
    python3 rank_fisher.py --negative   altered inputs that must be refused
    (add --write to rewrite rank_fisher.txt / rank_fisher_negative.txt; without it the script only prints)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.  It imports rank_indep.py (exact
first and second descent) and rank_second.py (the curves D_lam, their tangent lines and the pairing <., D_lam>).

NOTATION (Fisher, Section 2).  E: y^2 = x(x^2 + a x + b), T = (0, 0); E': Y^2 = X(X^2 + a'X + b'), a' = -2a,
b' = a^2 - 4b, T' = (0, 0); phi: E -> E' with kernel <T>, phihat the dual.  S_1 = S^(phi)(E/Q) (rank_indep's
"E' side": classes d' | b', the quartics C'_d', torsors under E); S'_1 = S^(phihat)(E'/Q) (rank_indep's "E side":
classes d | b, the quartics C_d, torsors under E').  S_m (resp. S'_m) is the image in S_1 (resp. S'_1) of
Sel(m): Sel(2n) = S^(2^n)(E'/Q), Sel(2n+1) = S^(2^n phi)(E/Q) (and with E, E' swapped for the primed groups).
    [Thm 2.1]  for m odd, alternating pairings Theta_m on S_m x S_m and Theta'_m on S'_m x S'_m, with kernels
               S_{m+1}, S'_{m+1}; for m even, a pairing Theta_m: S_m x S'_m with left kernel S_{m+1} and right kernel
               S'_{m+1}.  The images of the rational points lie in every S_m, S'_m and in every kernel.
    [Lemma 2.2] the 2^n-descent bound: rank E(Q) <= dim S_{2n-1} + dim S'_{2n} - 2 (and the same with primes swapped;
               the -2 is dim E(Q)[phi] + dim E'(Q)[phihat]).
    [Sect. 3]  for x in Sel(m) with covering curve C and y in the first Selmer group of the other side,
               Theta_m = <x, y> = sum_v (f(P_v), y)_v, f in Q(C) a pushout: div f = 2 b, b a rational divisor in the
               class of the rational 2-torsion point of Pic^0(C); the Selmer set
               S(C, f) = { xi : xi = f(P_v) mod Im(delta_v) for all v }, i.e. the xi for which the double cover
               D_xi: xi w^2 = f of C (the next level) has points everywhere, is empty or a coset of the first Selmer
               group of C's own side.
These are the statements of the paper used here (read through summaries of arXiv:1509.03234v1; the PDF could not be
downloaded); they rest on the Cassels-Tate pairing and are not proved here.

WHAT IS COMPUTED, EXACTLY, FOR EACH CURVE.
 (1) rank_indep: the first descent (|S_1|, |S'_1| as upper bounds) and the second descent: S_2 <= <known, y0>,
     S'_2 <= <known, y1> (y0, y1 the open classes; known = the images of T, T' and the known point G).
 (2) The lifts of y0 and y1 to level 2 (curves D_lam: X^2 = lam U, Y^2 = lam V) with a checked local point at every
     place of S (all bad primes of D_lam, since the v4 fix of rank_second.py) and at R: so y0 is in S_2 and y1 in
     S'_2, and dim S_2 = dim <known, y0>, dim S'_2 = dim <known, y1> exactly.
 (3) Theta_2(y0, y1) = <D_lam, y0> = +1 for a checked lift D_lam of y1 (rank_second's formula, recomputed here).
     By Thm 2.1, S_3 = S_2 and S'_3 = S'_2.
 (4) NEW, the level-3 lifts, constructively (Fisher Section 3 with the pushout f = l_X l_Y of Lemma 7.1: the
     tangent planes at rational points of the two rational singular quadrics X^2 = lam U, Y^2 = lam V of the pencil
     of D_lam).  For each checked lift D_lam: S3 = R, 2, the primes of S above, the bad primes of E and E', and
     3, 5, 7, 11, 13.  At each place of S3 one sampled value f(P_v) and the local image Im(delta_v) give the coset
     f(P_v) Im(delta_v) (Fisher Sect. 3); every other sampled value is checked to lie in it.  Im(delta_v) is
     rank_indep's local image of this side's first descent where every possible class is certified, the unit
     classes at an odd prime of good reduction, and, where rank_indep leaves classes undecided (D1 at 2), the
     unique pair of subgroups with |Im delta_phi| |Im delta_phihat| = 8 (p = 2; 4 at odd p) that contains the
     certified classes and the differences of values of f (pin_images).  The F_2-linear system "class_v(xi) in the
     coset at every v" for xi in <-1, S3> gives the Selmer set S(D_lam, f).  For an xi of it, a sampled P_v with
     xi f(P_v) a non-zero square of Q_v is found at every place of S3 (exact square test of an integer known modulo
     p^40, the precision asserted; at R an exact sign; sampling by random points and by a breadth-first walk
     through the small discs where the forms vanish).  For p not in S3 (so p >= 17, p does not divide 2 lam, the
     discriminants, delta and det of the tangent decompositions, or b b'): D_lam has good reduction, l_X l_Y
     reduces to a form with divisor 2(F_X + F_Y) on it, xi is a unit, and xi w^2 = l_X l_Y is an etale double cover
     of the reduction, connected because F_X + F_Y - H reduces to the reduction of the 2-torsion point (non-zero,
     p odd of good reduction); a smooth genus-one curve over F_p has p + 1 - 2 sqrt(p) > 8 points (p >= 17), so one
     lies off the 4 zeros of f and Hensel lifts it.  So C3 = {X^2 = lam U, Y^2 = lam V, xi w^2 = l_X l_Y} is
     everywhere locally soluble: the class lifts to Sel(3).  An independent, constructive check of (3).
       Cross-check: by Poitou-Tate the Selmer set of D_lam is non-empty iff <D_lam, eta> = 0 for every eta of the
     dual first Selmer group (the other side's S_1).  These pairings are evaluated (rank_second's Hilbert-symbol
     formula) for a basis and compared; a control on E+ of (945,248) (lifts of the trivial class) must show both
     outcomes.
 (5) The 4-descent bound, exactly: Theta_3 is alternating on S_3 and Theta'_3 on S'_3, each with the images of the
     rational points (codimension 1) in its kernel; an alternating form whose kernel has codimension <= 1 is zero.
     So S_4 = S_3, S'_4 = S'_3, and Lemma 2.2 with n = 2 gives exactly 2 + 3 - 2 = 3 (dimensions of the two sides
     in either order).  The 4-descent cannot give less than 3 on these curves.  A bound 1 (Magma) needs
     Theta_4(y0, y1) != 0: Fisher's fifth descent (Sect. 7).  (Theta_4 is the fifth descent, not the 8-descent,
     which is the sixth; with one rational 2-torsion point Magma's method stops at S_5.)
 (6) Theta_4 as one pairing of explicit curves.  On the side whose first Selmer group S_1 equals S_2 (the E' side
     of (760,39), (72,65), (920,129); the E side of (945,248) and of D1), the obstruction to lifting a level-3 lift C3 of y to level 4 is its
     pairing with S_1 = <images of rational points, y>, which is 0 (Thm 2.1; Theta_3(y, y) = 0, alternating).  With
     Theta_4(alpha_1, beta_1) = <alpha_4, beta_1> and <psi x, z> = <x, psihat z> this gives
     Theta_4(y0, y1) = <[C3], [D]>, the Cassels-Tate pairing in Sha(E) (resp. Sha(E')) of the checked C3 and a
     checked level-2 lift D = D_lam of the other side.  Both are written to rank_fisher_curves.json (--write).
NOT COMPUTED HERE: Theta_4.  It is not blocked.  Fisher Sect. 7: the level-3 curve C3 satisfies hypothesis (19)
(Im rho in <(12),(34)>: the singular quadrics of its pencil are the pairs over Q(sqrt disc U) and Q(sqrt disc V));
Lemma 7.1 then gives a pushout form on C3 from tangent planes at points of the two conjugate cones (a conic over a
quadratic field), Thm 7.2 the level-4 curve C4 as a quadric intersection, Thm 7.3-7.4 a pushout on C4 from the
companion C4^dagger, and Theta_4(y0, y1) = sum_v (f(P_v), y)_v (Sect. 3).  Hypothesis (19) fails only one step
later, at C4 (Sect. 9).  This is presumably how Magma's bound 1 is reached.  Not implemented here (it needs conics
over quadratic fields); see rank_fisher5.py for the first steps.
"""
import os, sys, json, random
from fractions import Fraction as Fr
from math import gcd, isqrt

WRITE = '--write' in sys.argv
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rank_indep as RI
import rank_second as R2

HINTS_FILE = os.path.join(here, 'rank_fisher_hints.json')
_BASE = {}

def load_hints():
    """the hints of rank_indep.py and rank_second.py (PARI's conic points and factorisations, each checked on use),
    then this script's own extra ones, so that a run needs no PARI."""
    if os.path.exists(RI.HINTS_FILE):
        for k, v in json.load(open(RI.HINTS_FILE)).items(): RI.HINTS[k] = v
    for f in (R2.HINTS_FILE, HINTS_FILE):
        if os.path.exists(f):
            H = json.load(open(f))
            R2.HINTS['conic'].update(H.get('conic', {}))
            for k, v in H.get('rank_indep', {}).items():
                for kk, vv in v.items(): RI.HINTS.setdefault(k, {}).setdefault(kk, vv)
        if f == R2.HINTS_FILE:
            _BASE['conic'] = set(R2.HINTS['conic'])
            _BASE['ri'] = {k: set(v) for k, v in RI.HINTS.items() if isinstance(v, dict)}

def save_hints():
    if not WRITE: return
    extra_c = {k: v for k, v in R2.HINTS['conic'].items() if k not in _BASE.get('conic', set())}
    extra_r = {k: {kk: vv for kk, vv in RI.HINTS.get(k, {}).items() if kk not in _BASE.get('ri', {}).get(k, set())}
               for k in ('factor', 'conic', 'points')}
    json.dump({'conic': extra_c, 'rank_indep': extra_r}, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)

OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

NPREC = R2.NPREC
# the equations of the curves written to rank_fisher_curves.json, in the conventions of rank_second.py
EQ2 = ('D_lam: X^2 = lam*U(r,s), Y^2 = lam*V(r,s) in P^3(r:s:X:Y); U = [u0, u1, u2] means u0*s^2 + u1*r*s + u2*r^2 '
       '(rank_second.qv), the same for V; l_X = [c0, c1, c2] means c0*r + c1*s + c2*X, the tangent line at the point '
       'p_X = (r, s, X) of X^2 = lam*U(r,s); l_Y, p_Y the same with Y, V')
EQ3 = ('C3: X^2 = lam*U(r,s), Y^2 = lam*V(r,s), xi*w^2 = l_X(r,s,X)*l_Y(r,s,Y) (conventions as for D_lam); a local '
       'point [r, s, sx, sy] at v is (r : s : sx*X : sy*Y) with X, Y the square roots of lam*U(r,s), lam*V(r,s) in Q_v '
       '(positive at v = inf; at a prime p, rank_second.padic_sqrt)')
SMALL = [3, 5, 7, 11, 13]       # primes below 17 always in S3 (the Hasse bound argument needs p >= 17)
NSAMP = 6                       # sampled points (r : s) per place, each with the four sign choices

# ---------------------------------------------------------------- classes of values of the pushout
def val_mod(c, r, s, W, p):
    """c0 r + c1 s + c2 W modulo p^NPREC, W a p-adic square root known modulo p^NPREC."""
    return (c[0] * r + c[1] * s + c[2] * W) % p ** NPREC

def cls_of_mod(y, p):
    """class in Q_p*/Q_p*^2 (rank_second.sq_cls bit vector) of an element of Z_p known modulo p^NPREC."""
    assert y % p ** NPREC != 0, 'p-adic precision (zero modulo p^N)'
    v = RI.vp(y % p ** NPREC, p)
    assert v <= NPREC - 6, 'p-adic precision'
    u = ((y % p ** NPREC) // p ** v) % (8 if p == 2 else p)
    return R2.sq_cls(p ** v * u, p)

def point_data(D, p, r, s_):
    X = R2.padic_sqrt(D.lam * R2.qv(D.U, r, s_), p, NPREC)
    Y = R2.padic_sqrt(D.lam * R2.qv(D.V, r, s_), p, NPREC)
    return X, Y

def f_class(D, p, pt):
    """class of f(P) = l_X(P) l_Y(P) at the point P = (r : s : sx X : sy Y) of D_lam(Q_p) (p prime) or R."""
    r, s_, sx, sy = pt
    if p == 'inf':
        X2, Y2 = D.lam * R2.qv(D.U, r, s_), D.lam * R2.qv(D.V, r, s_)
        sg = R2._sign_lin(D.Y['l'], r, s_, Y2, sy) * R2._sign_lin(D.X['l'], r, s_, X2, sx)
        return 1 if sg < 0 else 0
    X, Y = point_data(D, p, r, s_)
    return cls_of_mod(val_mod(D.X['l'], r, s_, sx * X, p), p) ^ cls_of_mod(val_mod(D.Y['l'], r, s_, sy * Y, p), p)

def sample(D, p, rng, want=NSAMP):
    """[(class of f(P), P)] for sample points P of D_lam(Q_v); every P is checked to be on D_lam(Q_v)."""
    out = []
    if p == 'inf':
        for r in R2._real_samples(D.U, D.V):
            X2, Y2 = D.lam * R2.qv(D.U, r, 1), D.lam * R2.qv(D.V, r, 1)
            if X2 > 0 and Y2 > 0:
                for sx in (1, -1):
                    for sy in (1, -1):
                        try: out.append((f_class(D, 'inf', (r, 1, sx, sy)), (r, 1, sx, sy)))
                        except AssertionError as ex:
                            if 'real zero' not in str(ex): raise
        return out
    pts = R2.local_points(D.lam, D.U, D.V, p, rng, want)
    for (r, s_) in pts:
        assert R2._is_padic_square(D.lam * R2.qv(D.U, r, s_), p) and R2._is_padic_square(D.lam * R2.qv(D.V, r, s_), p)
        for sx in (1, -1):
            for sy in (1, -1):
                try: out.append((f_class(D, p, (r, s_, sx, sy)), (r, s_, sx, sy)))
                except AssertionError as ex:
                    if 'precision' not in str(ex): raise
    return out

def tree_points(D, p, rng, maxpts=80, maxnodes=3000):
    """points (r, s) of D_lam(Q_p) from the discs t0 + p^k Z_p of the two charts that rank_second.local_points
    explores (a disc is cut further where a form without its content vanishes mod p), breadth first, at most 6
    points per depth and chart, without the random phase: these reach the small discs where the rarer classes of f live.  Candidates only;
    sample_tree tests every point exactly."""
    out = []
    lam, U, V = D.lam, D.U, D.V
    def good(r, s_):
        if r == 0 and s_ == 0: return False
        A, B = lam * R2.qv(U, r, s_), lam * R2.qv(V, r, s_)
        return A != 0 and B != 0 and R2._is_padic_square(A, p) and R2._is_padic_square(B, p)
    from collections import deque
    for chart in (0, 1):
        polys = [[lam * F[0], lam * F[1], lam * F[2]] if chart == 0 else [lam * F[2], lam * F[1], lam * F[0]] for F in (U, V)]
        dq = deque([(0, 1 if chart == 0 else p, 0)])
        nodes = 0; per = {}
        while dq and nodes < maxnodes and len(out) < maxpts:
            nodes += 1
            t0, sc, depth = dq.popleft()
            hs = []
            for g in polys:
                h = RI.ptrim(RI.psub_lin(g, t0, sc))
                k = min(RI.vp(c_, p) for c_ in h if c_)
                hs.append([c_ // p ** k for c_ in h])
            if p < 3000: cands = range(p)
            else:
                cands = set()
                for h in hs:
                    hh = RI.ptrim(h, p)
                    if len(hh) > 1: cands.update(RI.roots_mod(hh, p))
                cands.update(rng.randrange(p) for _ in range(10))
                cands = sorted(cands)
            for c_ in cands:
                t = t0 + sc * c_
                r, s_ = (t, 1) if chart == 0 else (1, t)
                if per.get(depth, 0) < 6 and gcd(r, s_) == 1 and (r, s_) not in out and good(r, s_):
                    out.append((r, s_)); per[depth] = per.get(depth, 0) + 1
                if depth < 10 and ((p == 2 and depth < 5) or any(RI.peval(h, c_) % p == 0 for h in hs)):
                    dq.append((t, sc * p, depth + 1))
    return out

def sample_tree(D, p, rng):
    out = []
    for (r, s_) in tree_points(D, p, rng):
        assert R2._is_padic_square(D.lam * R2.qv(D.U, r, s_), p) and R2._is_padic_square(D.lam * R2.qv(D.V, r, s_), p)
        for sx in (1, -1):
            for sy in (1, -1):
                try: out.append((f_class(D, p, (r, s_, sx, sy)), (r, s_, sx, sy)))
                except AssertionError as ex:
                    if 'precision' not in str(ex): raise
    return out

# ---------------------------------------------------------------- F_2 linear algebra
def rref(vecs):
    """reduced echelon basis {pivot bit: vector} of the span of the bit vectors."""
    B = {}
    for v in vecs:
        for pb, b in B.items():
            if v >> pb & 1: v ^= b
        if v:
            pv = v.bit_length() - 1
            for pb in list(B):
                if B[pb] >> pv & 1: B[pb] ^= v
            B[pv] = v
    return B

def red(B, x):
    for pb, b in B.items():
        if x >> pb & 1: x ^= b
    return x

def solve_f2(eqs, n):
    """eqs: [(mask over n unknowns, rhs bit)].  Returns (particular solution, nullspace basis) or None."""
    rows = []
    for m, r in eqs:
        for (pm, pr, pv) in rows:
            if m >> pv & 1: m ^= pm; r ^= pr
        if m == 0:
            if r: return None
            continue
        pv = m.bit_length() - 1
        rows = [((pm ^ m, pr ^ r, q) if pm >> pv & 1 else (pm, pr, q)) for (pm, pr, q) in rows]
        rows.append((m, r, pv))
    piv = {q: (pm, pr) for (pm, pr, q) in rows}
    x0 = 0
    for q, (pm, pr) in piv.items():
        if pr: x0 |= 1 << q
    free = [i for i in range(n) if i not in piv]
    null = []
    for f in free:
        v = 1 << f
        for q, (pm, pr) in piv.items():
            if pm >> f & 1: v |= 1 << q
        null.append(v)
    for m, r in eqs:
        assert bin(m & x0).count('1') % 2 == r
        for v in null: assert bin(m & v).count('1') % 2 == 0
    return x0, null

# ---------------------------------------------------------------- the level-3 Selmer set of one D_lam
def places_S3(D, R):
    ps = set(D.S(1)) | set(R['bad']) | set(SMALL)
    return sorted(ps)

def local_image_bits(v, R, side):
    """Im delta_v (the local image of the first descent of this side, which is the twisting group of the double
    covers C3 -> D_lam) as a set of rank_second.sq_cls bit vectors: from rank_indep's local images where every
    possible class is certified (else None); the unit classes at an odd prime of good reduction of E and E'; at R
    the d = +-1 with C_d real (rank_indep.real_ok)."""
    A, Bc, loc = (R['a'], R['b'], R['locS']) if side == 'E' else (R['ap'], R['bp'], R['locSp'])
    if v == 'inf': return {R2.sq_cls(d, 'inf') for d in (1, -1) if RI.real_ok(d, A, Bc)}
    if (side, v) in R.get('exact', {}): return R['exact'][(side, v)]
    if v in loc:
        poss, cert = loc[v]
        if poss == cert:
            return {R2.sq_cls(dl, v) for dl in RI.class_reps(v) if RI.qp_class(dl, v) in cert}
        # not decided by rank_indep at its depth: return the bounds cert <= Im <= poss (as bit vectors); the
        # caller adds the differences of sampled values of f (which lie in Im by Fisher Section 3) to the lower one
        return ('bounds', {R2.sq_cls(dl, v) for dl in RI.class_reps(v) if RI.qp_class(dl, v) in cert},
                {R2.sq_cls(dl, v) for dl in RI.class_reps(v) if RI.qp_class(dl, v) in poss})
    assert v != 2 and v not in R['bad']
    return {0, 2}

def pin_images(R, Dsides, rng):
    """At a prime p where rank_indep leaves a local image undecided on a side: Im(delta_phi,p) and
    Im(delta_phihat,p) contain the certified classes and the differences of values of f on a D_lam of the side
    (Fisher Section 3), and |Im(delta_phi,p)| |Im(delta_phihat,p)| = |E(Q_p)/2E(Q_p)| |E'(Q_p)[phihat]/phi(E(Q_p)[2])|
    = 8 at p = 2 and 4 at odd p here (E(Q_p)[2] = Z/2 iff -1 is not a square in Q_p, since a^2 - 4b = -square;
    |E(Q_p)/2E(Q_p)| = |E(Q_p)[2]| |2|_p^-1).  When the two lower bounds multiply to that number, both images are
    known exactly.  Dsides: {side: a checked D_lam of that side}."""
    R.setdefault('exact', {})
    for p in R['bad']:
        und = [s for s, loc in (('E', R['locS']), ("E'", R['locSp'])) if loc[p][0] != loc[p][1]]
        if not und: continue
        lows = {}
        for side, D in Dsides.items():
            loc = R['locS'] if side == 'E' else R['locSp']
            cert = {R2.sq_cls(dl, p) for dl in RI.class_reps(p) if RI.qp_class(dl, p) in loc[p][1]}
            poss = {R2.sq_cls(dl, p) for dl in RI.class_reps(p) if RI.qp_class(dl, p) in loc[p][0]}
            sm = sample(D, p, rng) + sample_tree(D, p, rng)
            a0 = sm[0][0]
            B = rref(list(cert) + [c ^ a0 for c, _ in sm])
            span = set()
            for k in range(1 << len(B)):
                x = 0
                for i, bb in enumerate(B.values()):
                    if k >> i & 1: x ^= bb
                span.add(x)
            assert span <= poss, ('a certified class outside the possible ones', p, side)
            lows[side] = span
        if len(lows) == 2:
            mod2 = 0 if p == 2 else 1
            prod = 8 if p == 2 else 4
            if len(lows['E']) * len(lows["E'"]) == prod:
                for side in ('E', "E'"): R['exact'][(side, p)] = lows[side]
                R.setdefault('pinned', []).append((p, len(lows['E']), len(lows["E'"])))

def selmer_set_level3(D, R, side, seed=0, max_rounds=2):
    """Fisher Section 3: S(D_lam, f) = { xi : class_v(xi) in f(P_v) Im(delta_v) for every v } with f = l_X l_Y,
    computed for xi in <-1, S3> from ONE sampled value f(P_v) per place and the local image Im(delta_v) (exact,
    from rank_indep); every further sampled value is checked to lie in that coset (a check of the formula).  If it
    is non-empty, an xi of it whose class is attained by a sampled point at every place is sought (more sampling if
    needed), giving checked points of C3.  Returns nonempty True (xi and points), False (empty), None (an image not
    known exactly), and 'theory' (non-empty, but no sampled point for some place)."""
    rng = random.Random(seed + abs(D.lam) % 997)
    S3 = places_S3(D, R)
    places = ['inf'] + S3
    samples, coset = {}, {}
    for v in places:
        samples[v] = sample(D, v, rng)
        assert samples[v], ('no local point of D_lam found', v, D.lam)
        I = local_image_bits(v, R, side)
        if isinstance(I, tuple):
            # Im(delta_v) lies between the certified classes and the possible ones; the differences of the values
            # of f lie in it (Fisher Section 3); it is known exactly when the group they generate is as large as
            # the set of possible classes
            lowc, high = I[1], I[2]
            for rnd in range(2):
                a0 = samples[v][0][0]
                low = rref(list(lowc) + [c ^ a0 for c, _ in samples[v]])
                if 1 << len(low) == len(high): break
                samples[v] += sample_tree(D, v, rng)
            if 1 << len(low) != len(high): return dict(lam=D.lam, places=places, nonempty=None, dim=None)
            I = set()
            for k in range(1 << len(low)):
                x = 0
                for i, bb in enumerate(low.values()):
                    if k >> i & 1: x ^= bb
                I.add(x)
            assert I == high, 'the local image is not the set of possible classes'
        a0 = samples[v][0][0]
        coset[v] = {a0 ^ i for i in I}
        for (c, _) in samples[v]:
            assert c in coset[v], ('a value of f outside f(P) Im(delta): f is not a pushout', v, D.lam)
    gens = [-1] + S3
    n = len(gens)
    eqs = []
    for v in places:
        a0 = min(coset[v])
        B = rref([c ^ a0 for c in coset[v]])
        assert len(coset[v]) == 1 << len(B)
        gv = [red(B, R2.sq_cls(g, v)) for g in gens]
        t = red(B, a0)
        nb = 3 if v == 2 else (1 if v == 'inf' else 2)
        for j in range(nb):
            m = 0
            for i in range(n):
                if gv[i] >> j & 1: m |= 1 << i
            eqs.append((m, t >> j & 1))
    sol = solve_f2(eqs, n)
    if sol is None:
        return dict(lam=D.lam, places=places, nonempty=False, dim=None)
    x0, null = sol
    for rnd in range(max_rounds + 1):
        for k in range(1 << len(null)):
            e = x0
            for i, v_ in enumerate(null):
                if k >> i & 1: e ^= v_
            xi = 1
            for i, g in enumerate(gens):
                if e >> i & 1: xi *= g
            cert = {}
            for v in places:
                cv = R2.sq_cls(xi, v)
                assert cv in coset[v]
                hit = [P for (c, P) in samples[v] if c == cv]
                if not hit: break
                cert[v] = hit[0]
            else:
                return dict(lam=D.lam, places=places, nonempty=True, dim=len(null), xi=xi, cert=cert, D=D)
        if rnd < max_rounds:
            for v in places:
                if v != 'inf' and len(set(c for c, _ in samples[v])) < len(coset[v]):
                    samples[v] += sample_tree(D, v, rng) if rnd == 0 else sample(D, v, rng, 3 * NSAMP)
                    for (c, _) in samples[v]:
                        assert c in coset[v], ('a value of f outside f(P) Im(delta): f is not a pushout', v, D.lam)
    return dict(lam=D.lam, places=places, nonempty='theory', dim=len(null))

def check_C3_point(D, xi, v, pt):
    """the certificate: xi f(P) is a non-zero square of Q_v at the point P of D_lam(Q_v), so (P, w) is a point of
    C3: X^2 = lam U, Y^2 = lam V, xi w^2 = l_X l_Y over Q_v.  Exact (p-adic precision asserted)."""
    r, s_, sx, sy = pt
    if v == 'inf':
        X2, Y2 = D.lam * R2.qv(D.U, r, s_), D.lam * R2.qv(D.V, r, s_)
        assert X2 > 0 and Y2 > 0, 'not a real point of D_lam'
        sg = R2._sign_lin(D.Y['l'], r, s_, Y2, sy) * R2._sign_lin(D.X['l'], r, s_, X2, sx)
        assert sg * (1 if xi > 0 else -1) > 0, 'xi f(P) is not a positive real number'
        return True
    p = v
    assert R2._is_padic_square(D.lam * R2.qv(D.U, r, s_), p) and R2._is_padic_square(D.lam * R2.qv(D.V, r, s_), p), \
        'not a point of D_lam(Q_p)'
    X, Y = point_data(D, p, r, s_)
    assert (X * X - D.lam * R2.qv(D.U, r, s_)) % p ** NPREC == 0 and (Y * Y - D.lam * R2.qv(D.V, r, s_)) % p ** NPREC == 0
    # l_X(P), l_Y(P) are known modulo p^NPREC; each valuation is asserted <= NPREC - 6, so its unit part is known
    # modulo p^6 (>= 8 at p = 2); xi is an exact integer.  Then xi l_X(P) l_Y(P) is a square of Q_p iff the total
    # valuation is even and the product of the unit parts is a square modulo 8 (p = 2) or modulo p.
    w, u = RI.vp(xi, p), xi // p ** RI.vp(xi, p)
    for c, W in ((D.X['l'], sx * X), (D.Y['l'], sy * Y)):
        y = val_mod(c, r, s_, W, p)
        assert y != 0, 'p-adic precision (zero modulo p^N)'
        vy = RI.vp(y, p)
        assert vy <= NPREC - 6, 'p-adic precision'
        w += vy; u *= y // p ** vy
    assert w % 2 == 0 and (u % 8 == 1 if p == 2 else RI.legendre(u % p, p) == 1), 'xi f(P) is not a square in Q_p'
    return True

def certify_C3(res):
    D, xi = res['D'], res['xi']
    # the argument for the primes outside S3 needs xi to be a unit there
    assert xi != 0 and all(p in res['places'] for p in RI.factor_int(abs(xi))), 'xi has a prime outside S3'
    for v in res['places']:
        assert check_C3_point(D, xi, v, res['cert'][v])
    return True

# ---------------------------------------------------------------- one curve
def log2(n):
    k = 0
    while (1 << (k + 1)) <= n: k += 1
    assert 1 << k == n, 'not a power of 2'
    return k

def analyse(name, a, b, G, lam_cap=12, cross=True):
    say('%s:  y^2 = x(x^2 + a x + b),  a = %d,  b = %d' % (name, a, b))
    ap, bp = -2 * a, a * a - 4 * b
    assert bp < 0 and isqrt(-bp) ** 2 == -bp, 'a^2 - 4b is not minus a square'
    r, rec = RI.rank_bound(a, b, [G], RI.find_points([0, ap, 0, bp, 0]), say=lambda s: None)
    R = RI.first_descent(a, b)
    assert len(rec['E']['open']) == 1 and len(rec["E'"]['open']) == 1
    y1 = rec['E']['open'][0][0]; y0 = rec["E'"]['open'][0][0]
    for (dd, AA, BB, loc) in ((y1, a, b, R['locS']), (y0, ap, bp, R['locSp'])):
        for p in R['bad']: assert R2.certified(dd, AA, BB, p, loc), ('class not certified locally soluble', dd, p)
        assert RI.real_ok(dd, AA, BB)
    kS, kSp = log2(len(R['S'])), log2(len(R['Sp']))
    kH, kHp = log2(len(rec['E']['known'])), log2(len(rec["E'"]['known']))
    say('  (1) first descent: dim S\'_1 <= %d (E side, classes d | b), dim S_1 <= %d (E\' side, classes d\' | b\'); '
        'second descent: dim S\'_2 <= %d, dim S_2 <= %d; known: %d, %d; open classes y1 = %d (E side), y0 = %d (E\' side)'
        % (kS, kSp, rec['E']['log2'], rec["E'"]['log2'], kH, kHp, y1, y0))
    assert rec['E']['log2'] == kH + 1 and rec["E'"]['log2'] == kHp + 1
    out = dict(name=name)
    level3 = {}
    if any(R['locS'][p][0] != R['locS'][p][1] or R['locSp'][p][0] != R['locSp'][p][1] for p in R['bad']):
        Ds = {}
        for side, (dd, AA, BB) in (('E', (y1, a, b)), ("E'", (y0, ap, bp))):
            P0 = RI.find_conic_point(dd, AA, BB // dd)
            surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
            for lam in sorted(surv, key=abs):
                try:
                    D = R2.DLam(lam, info['U'], info['V'], bad=info['primes']); D.pairing(1); Ds[side] = D; break
                except AssertionError as ex:
                    if 'no conic point' not in str(ex) and 'no local point' not in str(ex): raise
        pin_images(R, Ds, random.Random(7))
        say('  local images left undecided by rank_indep, now known exactly from |Im phi| |Im phihat| = 8 (p = 2) or 4: '
            + (', '.join('p = %d (sizes %d, %d)' % t for t in R.get('pinned', [])) or 'none'))
    for side, (dd, AA, BB, other, kown, kother, grp_other) in (('E', (y1, a, b, y0, kS, kSp, R['Sp'])),
                                                              ("E'", (y0, ap, bp, y1, kSp, kS, R['S']))):
        c = BB // dd
        P0 = RI.find_conic_point(dd, AA, c)
        surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
        U, V = info['U'], info['V']
        R2.check_conic_forms(dd, AA, c, U, V, P0[2])
        # basis of the other side's first Selmer group (for the cross-check), with the class of its 2-torsion point first
        tors_other = RI.sqfree_part(bp if side == 'E' else b)
        Bo = rref_ints([tors_other, other] + sorted(grp_other, key=abs))
        say('  %s side: class %d; lifts D_lam not excluded by rank_indep\'s second descent: %d; used here: up to %d'
            % (side, dd, len(surv), lam_cap))
        rows = []
        for lam in sorted(surv, key=abs):
            if len(rows) >= lam_cap: break
            if len(rows) >= 4 and all(x[2]['nonempty'] is None for x in rows): break
            try:
                D = R2.DLam(lam, U, V, bad=info['primes'])
            except AssertionError as ex:
                if 'no conic point' not in str(ex): raise
                continue
            try:
                th2, _ = D.pairing(other)                  # level 2: local points at every place of S, and Theta_2
            except AssertionError as ex:
                if 'no local point found' not in str(ex): raise
                continue
            assert not D.special_rational
            res = selmer_set_level3(D, R, side)
            if res['nonempty'] is True:
                certify_C3(res)
            if len(rows) < 2:
                # the outcome must not depend on the sampled points: again with other random points
                res2 = selmer_set_level3(D, R, side, seed=101)
                assert (res2['nonempty'] is False) == (res['nonempty'] is False) and res2['dim'] == res['dim'], \
                    ('the level-3 Selmer set depends on the sample', lam)
                if res2['nonempty'] is True: certify_C3(res2)
            pv = {}
            if cross:
                for eta in Bo:
                    try: pv[eta] = D.pairing(eta)[0]
                    except AssertionError as ex:
                        if 'not constant' not in str(ex) and 'no local point' not in str(ex): raise
                        pv[eta] = None                      # eta not shown to be in the Selmer group
                pred = all(x == 1 for x in pv.values() if x is not None)
                ne = res['nonempty']
                agree = ((ne is not False) == pred) if ne is not None and None not in pv.values() else None
            else: agree = None
            rows.append((lam, th2, res, pv, agree))
            say('    lam = %d: level 2 checked (Theta_2 = <D_lam, %d> = %+d); level 3: %s%s'
                % (lam, other, th2,
                   ('Selmer set non-empty (a coset of dimension %d; dim of this side\'s S_1 bound: %d); xi = %d, '
                    'C3 checked locally soluble at R and %d primes' % (res['dim'], kown, res['xi'], len(res['places']) - 1))
                   if res['nonempty'] is True else {False: 'Selmer set empty', None: 'not decided (a local image not known exactly)',
                                                     'theory': 'Selmer set non-empty (dimension %s), but no sampled local point for every place of one xi' % res['dim']}[res['nonempty']],
                   ('; pairings with a basis of the other side\'s S_1: %s -> %s' %
                    (' '.join('%d:%s' % (k_, '+1' if v_ == 1 else ('-1' if v_ == -1 else '?')) for k_, v_ in pv.items()),
                     'agrees' if agree else ('DISAGREES' if agree is False else 'not compared'))) if cross else ''))
        assert rows, 'no checked lift'
        assert all(th == 1 for (_, th, _, _, _) in rows), 'Theta_2 is not trivial'
        assert all(ag is not False for (*_, ag) in rows), 'the level-3 Selmer set disagrees with the pairing'
        ok3 = [x for x in rows if x[2]['nonempty'] is True]
        undecided = [x for x in rows if x[2]['nonempty'] is None]
        if not ok3:
            # allowed only when every lift tried is undecided because a local image is not known exactly (D1 at 2)
            assert undecided and len(undecided) == len(rows), 'no lift to level 3 found'
            say('    (no constructive level-3 lift on this side: a local image of the first descent at some place is not'
                ' decided by rank_indep; the theory (Thm 2.1, from (3)) still gives the lift)')
        level3[side] = dict(dd=dd, lifts=len(rows), lift3=len(ok3), lam=ok3[0][0] if ok3 else None,
                            xi=ok3[0][2]['xi'] if ok3 else None, dims=sorted(set(x[2]['dim'] for x in ok3)))
        def dl_json(D):
            return dict(equations=EQ2, lam=str(D.lam), U=[str(c_) for c_ in D.U], V=[str(c_) for c_ in D.V],
                        l_X=[str(c_) for c_ in D.X['l']], l_Y=[str(c_) for c_ in D.Y['l']],
                        p_X=[str(c_) for c_ in D.X['P']], p_Y=[str(c_) for c_ in D.Y['P']])
        level3[side]['level2_curve'] = dl_json(rows[0][2]['D'] if 'D' in rows[0][2] else R2.DLam(rows[0][0], U, V, bad=info['primes']))
        if ok3:
            r3 = ok3[0][2]
            level3[side]['level3_curve'] = dict(dl_json(r3['D']), xi=str(r3['xi']), equations=EQ3,
                                                local_points={str(v): [str(x_) for x_ in P] for v, P in r3['cert'].items()})
    # (5) the bookkeeping
    dS2, dS2p = rec["E'"]['log2'], rec['E']['log2']        # S_2 (E' side), S'_2 (E side): exact after (2)
    dS3, dS3p = dS2, dS2p                                    # Theta_2 = 0 (Thm 2.1 (ii)); level-3 lifts checked in (4)
    assert dS3 - kHp <= 1 and dS3p - kH <= 1                 # codim of the known images in S_3, S'_3
    dS4, dS4p = dS3, dS3p                                    # alternating Theta_3, Theta'_3 with kernel of codim <= 1
    b4 = dS3 + dS4p - 2; b4p = dS3p + dS4 - 2
    say('  (5) dim S_2 = dim S_3 = dim S_4 = %d,  dim S\'_2 = dim S\'_3 = dim S\'_4 = %d  =>  4-descent bound '
        '(Lemma 2.2, n = 2): %d + %d - 2 = %d (the other order: %d)' % (dS3, dS3p, dS3, dS4p, b4, b4p))
    assert b4 == b4p == 3
    say('  => %s: the descent by 2-isogenies through level 4 gives exactly 3; Magma\'s 1 needs Theta_4(y0, y1) != 0 '
        '(Fisher\'s fifth descent, Sect. 7)' % name)
    # Theta_4 as one Cassels-Tate pairing of explicit curves: on the side whose first Selmer group equals S_2
    # (every level-3 lift of it lifts to level 4: the obstruction pairs it with S_1 = S_2 = <images of rational points, y>, and
    # Theta_3(y, y) = 0), Theta_4(y0, y1) = <[C3 of that side], [D_lam of the other side]> (Fisher Sect. 2-3 and
    # <psi x, y> = <x, psihat y>).
    if kSp == dS2: red4 = ("E'", 'E')            # alpha side = E' side: C3 is a torsor under E, D_lam (E side) too
    elif kS == dS2p: red4 = ('E', "E'")
    else: red4 = None
    out.update(level3={k: {kk: vv for kk, vv in v.items()} for k, v in level3.items()}, bound4=b4, red4=red4)
    if red4 and level3[red4[0]]['lift3']:
        say('  Theta_4(y0, y1) = <[C3], [D]> in Sha(%s): C3 = {X^2 = lam U, Y^2 = lam V, xi w^2 = l_X l_Y} the checked level-3 '
            'lift of the %s side (lam = %d, xi = %d), D = the checked level-2 lift of the %s side (lam = %d); '
            'both curves are in rank_fisher_curves.json' % ('E' if red4[0] == "E'" else "E'", red4[0], level3[red4[0]]['lam'],
                                                          level3[red4[0]]['xi'], red4[1], int(level3[red4[1]]['level2_curve']['lam'])))
    return out

def rref_ints(vals):
    """squarefree integers forming an F_2-basis of the span of vals (in order of first appearance)."""
    keyp = sorted(set(p for v in vals if abs(v) > 1 for p in RI.factor_int(abs(v))))
    gens = [-1] + keyp
    def vec(v):
        x = 0
        if v < 0: x |= 1
        for i, p in enumerate(keyp):
            if RI.vp(abs(v), p) % 2: x |= 1 << (i + 1)
        return x
    basis, B = [], {}
    for v in vals:
        x = vec(v)
        y = red(B, x)
        if y:
            basis.append(v)
            B = rref(list(B.values()) + [x])
    return basis

def curves():
    out = []
    for (m, n) in ((760, 39), (945, 248), (72, 65), (920, 129)):
        a, b, G = R2.sigma_curve(m, n, 'E+')
        out.append(('%d,%d E+' % (m, n), a, b, G))
    a, b, G = R2.d1_curve(631, 71)
    out.append(('631,71 D1', a, b, G))
    return out

def control(say=say):
    """A control with an answer known beforehand: on E+ of (945,248) the lifts D_lam of the trivial class C_1 of
    the E side pair to -1 with some classes of S_1 that rank_indep's second descent excludes (rank_second.control).
    By Poitou-Tate the level-3 Selmer set of such a D_lam must then be empty, and non-empty when every pairing with
    a basis of S_1 is +1.  The Selmer sets below are computed only from the classes of f at local points (complete
    cosets), the pairings only from Hilbert symbols: two different computations."""
    a, b, G = R2.sigma_curve(945, 248, 'E+')
    ap, bp = -2 * a, a * a - 4 * b
    R = RI.first_descent(a, b)
    P0 = RI.find_conic_point(1, a, b)
    surv, info = RI.second_descent(1, a, b, P0, stop_at_first=False)
    R2.check_conic_forms(1, a, b, info['U'], info['V'], P0[2])
    basis = rref_ints([RI.sqfree_part(bp)] + sorted(R['Sp'], key=abs))
    certS1 = [e for e in basis if all(R2.certified(e, ap, bp, p, R['locSp']) for p in R['bad']) and RI.real_ok(e, ap, bp)]
    assert certS1 == basis, 'a basis element of S_1 is not certified locally soluble'
    say('control 945,248 E+: lifts D_lam of the trivial class of the E side (%d); basis of S_1 (E\' side): %s' % (len(surv), basis))
    n_empty = n_full = 0
    for lam in sorted(surv, key=abs):
        D = R2.DLam(lam, info['U'], info['V'], bad=info['primes'])
        res = selmer_set_level3(D, R, 'E')
        if res['nonempty'] is True: certify_C3(res)
        pv = [D.pairing(e)[0] for e in basis]
        pred = all(x == 1 for x in pv)
        st = {True: 'non-empty (xi = %s, C3 checked)' % res.get('xi'), False: 'empty', None: 'not decided',
              'theory': 'non-empty (no sampled point at every place)'}[res['nonempty']]
        agree = res['nonempty'] is not None and (res['nonempty'] is not False) == pred
        say('  lam = %d: pairings %s; level-3 Selmer set %s: %s' % (lam, pv, st, 'agrees' if agree else 'DISAGREES'))
        assert agree, 'the level-3 Selmer set disagrees with the pairing'
        n_empty += res['nonempty'] is False; n_full += res['nonempty'] is not False
    assert n_empty >= 1 and n_full >= 1, 'the control did not test both outcomes'
    say('  control: %d empty, %d non-empty, all agree with the pairings' % (n_empty, n_full))

def negative(say=say):
    """altered inputs that must be refused."""
    a, b, G = R2.sigma_curve(72, 65, 'E+')
    R = RI.first_descent(a, b)
    P0 = RI.find_conic_point(2, a, b // 2)
    surv, info = RI.second_descent(2, a, b, P0, stop_at_first=False)
    U, V = info['U'], info['V']
    D = R2.DLam(sorted(surv, key=abs)[0], U, V, bad=info['primes'])
    res = selmer_set_level3(D, R, 'E')
    assert res['nonempty'] is True and certify_C3(res)
    n = ok = 0
    def refuse(label, fn):
        nonlocal n, ok
        n += 1
        try: fn()
        except AssertionError as ex:
            ok += 1; say('  refused: %s (%s)' % (label, str(ex)[:110])); return
        say('  NOT REFUSED: %s' % label)
    xi = res['xi']
    p0 = next(p for p in res['places'] if p != 'inf' and RI.vp(xi, p) == 0)
    refuse('a forged xi (times a prime %d of S3) with the same local points' % p0,
           lambda: certify_C3(dict(res, xi=xi * p0)))
    refuse('a forged xi (times -1): the real point must refuse it',
           lambda: check_C3_point(D, -xi, 'inf', res['cert']['inf']))
    q = next(p for p in res['places'] if p != 'inf' and p > 2)
    r_, s_, sx, sy = res['cert'][q]
    bad_pt = next((r_ + k, s_) for k in range(1, 400)
                  if not (R2._is_padic_square(D.lam * R2.qv(U, r_ + k, s_), q) and R2._is_padic_square(D.lam * R2.qv(V, r_ + k, s_), q)))
    refuse('a forged local point at %d (not on D_lam(Q_%d))' % (q, q),
           lambda: check_C3_point(D, xi, q, (bad_pt[0], bad_pt[1], sx, sy)))
    Dbad = R2.DLam(D.lam, U, V, bad=info['primes'])
    Dbad.Y = dict(Dbad.Y, l=[Dbad.Y['l'][0] + 1, Dbad.Y['l'][1], Dbad.Y['l'][2] + 1])
    refuse('a non-tangent line in place of l_Y (f is not a pushout: its classes are not a coset)',
           lambda: (lambda r: (_ for _ in ()).throw(AssertionError('no xi')) if r['nonempty'] is not True else certify_C3(r))(selmer_set_level3(Dbad, R, 'E')))
    refuse('a lam outside the survivors of the second descent (no local point somewhere)',
           lambda: R2.DLam(sorted(surv, key=abs)[0] * 3, U, V, bad=info['primes']).pairing(1))
    ell = next(q_ for q_ in range(17, 10 ** 4) if RI.is_prime_td(q_) and q_ not in res['places'])
    refuse('a forged xi (times the prime %d outside S3, where no point is checked)' % ell,
           lambda: certify_C3(dict(res, xi=xi * ell)))
    # a forged local image: at a prime of S3 where Im(delta) has more than one class, keep only the trivial class
    qf = next(p for p in res['places'] if p != 'inf' and p in R['bad'] and len(local_image_bits(p, R, 'E')) > 1)
    Rf = dict(R, exact={('E', qf): {0}})
    refuse('a forged local image at %d (only the trivial class): the sampled values of f leave the coset' % qf,
           lambda: (lambda r: (_ for _ in ()).throw(AssertionError('no xi')) if r['nonempty'] is not True else certify_C3(r))(selmer_set_level3(D, Rf, 'E')))
    U2 = [U[0] + 1, U[1], U[2]]
    refuse('a mutated coefficient of U', lambda: R2.check_conic_forms(2, a, b // 2, U2, V, P0[2]))
    say('negative tests: %d of %d refused' % (ok, n))
    assert ok == n
    return ok, n

def assert_eq(x, y):
    assert x == y, 'record inconsistent: %s != %s' % (x, y)

def main():
    say('rank_fisher.py: Fisher\'s higher descents (arXiv:1509.03234) on the five Magma-alone curves, levels 1-4, exact')
    control()
    res = []
    for (name, a, b, G) in curves():
        res.append(analyse(name, a, b, G))
    say('summary: ' + '; '.join('%s: 4-descent bound %d (level-3 lifts checked: E side %d of %d, E\' side %d of %d)'
                                % (r['name'], r['bound4'], r['level3']['E']['lift3'], r['level3']['E']['lifts'],
                                   r['level3']["E'"]['lift3'], r['level3']["E'"]['lifts']) for r in res))
    say('Theta_4 (Fisher\'s fifth descent, Sect. 7) is not computed here: Magma\'s bound 1 still rests on Magma alone.')
    if WRITE:
        open(os.path.join(here, 'rank_fisher.txt'), 'w').write('\n'.join(OUT) + '\n')
        js = [dict(name=r['name'], theta4_reduction=r['red4'],
                   sides={s: {k: v for k, v in r['level3'][s].items() if k in ('dd', 'level2_curve', 'level3_curve')}
                          for s in r['level3']}) for r in res]
        json.dump(js, open(os.path.join(here, 'rank_fisher_curves.json'), 'w'), indent=1)

if __name__ == '__main__':
    load_hints()
    if '--negative' in sys.argv:
        say('rank_fisher.py --negative')
        negative()
        if WRITE: open(os.path.join(here, 'rank_fisher_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
        save_hints()
    else:
        main()
        save_hints()
