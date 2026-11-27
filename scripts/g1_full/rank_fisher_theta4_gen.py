#!/usr/bin/env python3
"""Fisher's fifth-descent pairing Theta_4 (arXiv:1509.03234 Sects. 2, 3, 7), the computation of
rank_fisher_theta4.py made GENERAL: it takes any curve of rank_fisher_curves.json (the five Magma-alone curves)
together with points of its two cone conics (hints, re-checked exactly), and computes
    Theta_4(lifted class, paired class) = sum_v (g0(P_v), y)_v
on the level-4 lift C4 of the side named first in the curve's 'theta4_reduction' (rank_fisher.py: the side whose
first Selmer group equals its S_2), with y the open class of the other side.

    python3 rank_fisher_theta4_gen.py <name> [--side=E|E'] [--write]   name: 72,65 (positive control) | 920,129 |
                                                                       760,39 | 631,71 --side=E' | 945,248 --side=E'
    python3 rank_fisher_theta4_gen.py <name> [--side=...] --negative [--write]
    (outputs only with --write: rank_fisher_theta4_gen_<m>_<n>.txt, ..._negative.txt)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.  Reviewed by a separate AI agent on
27 November 2026 (reviews/review_y5.md): holds with gaps.  Every committed output and CI line reduces on the E' side (the lifted class on E', as for
(72,65)): for (631,71) and (945,248) by --side=E', for the other curves because that is rank_fisher.py's side.  The
branch for a reduction on the E side (the else of prepare()) is run by no output and no CI line: untested,
and its K2 list holds only the class of T'.

RESULTS (rank_fisher_theta4_gen_<m>_<n>.txt): Theta_4 != 0 (product of the local symbols -1) on E+ of (72,65)
(positive control: rank_fisher_theta4.py's value reproduced), of (920,129) (Theta_4(-23, 156252226)), of (760,39)
(Theta_4(-799, 761)) and on the D1 curve of (631,71) (Theta_4(-627214, 31836878641), with --side=E': output
rank_fisher_theta4_gen_631_71_Eprime.txt); for each, K2 (the classes of the rational points pair to 0), K3 (two
more level-4 twists) and K5 (other k_1, k_3) agree.  So rank <= (dim S_4 - 1) + (dim S'_4 - 1) - 2 = 1 for these
curves without Magma, on the written arguments of rank_fisher_theta4.py (tier C, as for (72,65)).
THE SIDE.  rank_fisher.py names first the side whose first Selmer group equals its S_2, where every level-3 lift
lifts to level 4.  --side= takes the other side; then the level-3 lift of rank_fisher_curves.json is used only if
its level-4 Selmer set is found non-empty, and the C4 used is checked locally soluble at every place of S (points
found) and at R, the good-reduction argument covering the rest; Fisher's Thm 3.1(i) (independence of the lift,
needing y in S'_4) is then what makes <[C4], y> equal Theta_4(lifted class, y).  For the D1 curve of (631,71) the
E side's cone fields Q(sqrt 49449047761), Q(sqrt 31836878641) were out of reach of the finder, the E' side's
Q(sqrt +-627214) were not.
(945,248), added 27 November 2026 (run 7, t9): with --side=E' (output rank_fisher_theta4_gen_945_248_Eprime.txt)
Theta_4(-5820647, 2) != 0, K2, K3, K5 agree; rank <= (2 - 1) + (3 - 1) - 2 = 1.  Two differences from the other
curves: the conic points over Q(sqrt 5820647), Q(sqrt -5820647) come from rank_fisher_theta4_norm.py (a Legendre-type
descent over L, rnfisnorm, a quaternion reduction of the point), and the json's E'-side C3 (xi = 1271) has an empty
level-4 Selmer set, so its twist xi = 17794 = 1271 * 14 (14 in S_1 of that side; the hint key '945,248 E+|E'|xi')
is used; prepare() asserts the ratio is in S_1 and nothing rests on the choice (C4 is checked locally soluble).
The E side (Q(sqrt 1356488820560645889749762), class group (2,2,2,2), a regulator too large for generators) was
not attempted.  For (945,248) choose_Xi_lin also retries more seeds at a place where the first seed finds no point (p =
97 here); this fires on no other curve, whose outputs are unchanged.
THE CONIC POINTS (rank_fisher_theta4_gen_hints.json) were found by rank_fisher_theta4_find.py with
conic_lattice_find.py (finders: PARI, lattice reduction and enumeration), those of (945,248) by
rank_fisher_theta4_norm.py; they are inputs, re-checked exactly by rank_fisher5.analyse (on the cone, tangency
along the ruling, f and f_dag rational, q and c by exact linear algebra).
NEW EXACT CHECK (not in rank_fisher_theta4.py): the primes p where f0 mod p lies in span(QA, QB) mod p (where the
good-reduction argument for the double cover would not apply) are computed (gcd of the 3 x 3 minors) and put into S;
on the three curves they are already in B0 u {p < 47}.

WHAT IS TAKEN FROM rank_fisher_theta4.py UNCHANGED (imported): the exact checks (I1) (identities modulo I(C4) by
exact linear algebra), (I2) (resultants), (I3') (the ruling discriminant), (I4) (two classes of g0 on C4(Q_p)); the
pushout g0 = R3 / (W N1) (build_g0, generic in the fields L = Q(sqrt D1), L' = Q(sqrt D3)); the local points
(PARI p-adics, precision asserted) and the real points (interval enclosures); the Hilbert symbols; the
degenerate-prime test; the place set S = {R} u {p < 47} u B0 u primes(Xi) u {degenerate primes} u primes(y).
WHAT IS NEW HERE: (a) setup() without the (72,65)-specific assertion (D1, D3) = (-2, 2); (b) the curve, the
side, the paired class y and the images of the rational points are read from rank_fisher_curves.json and
rank_second / rank_indep instead of being hard-coded; for a reduction on the E side (untested, see above) the roles of E and E' are
exchanged (the 'working curve' is E': y^2 = x(x^2 - 2a x + a^2 - 4b), whose X-classes are the paired group);
(c) B0 = {2} u bad primes of the working curve u primes of the pencil discriminant (factored completely);
(d) the level-4 twists Xi by F_2-linear algebra on the sampled classes (closed up to cosets, as rank_fisher_theta4's
choose_Xi does by enumeration); (e) the labels of the printed fields.
THE WRITTEN ARGUMENTS of rank_fisher_theta4.py (evenness of div g0 from (I1)-(I3'); the good-reduction argument
outside S; Fisher's Thm 2.1, Lemma 2.2 and the Sect. 3 formula, read through summaries) are used unchanged; they
do not depend on the curve beyond what is checked here.
"""
import os, sys, json, time
from fractions import Fraction as Fr
from math import isqrt, gcd
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_fisher5 as R5
import rank_second as R2
import rank_indep as RI
import rank_fisher_theta4 as T
from rank_fisher_theta4 import pari

WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

HINTS_GEN = os.path.join(here, 'rank_fisher_theta4_gen_hints.json')
c, d, e, g, W = T.c, T.d, T.e, T.g, T.W
VS = T.VS

# dim S_4 (E' side), dim S'_4 (E side), from rank_fisher.txt, line '(5) dim S_2 = dim S_3 = dim S_4 = ...'
DIMS = {'760,39': (2, 3), '945,248': (2, 3), '72,65': (2, 3), '920,129': (2, 3), '631,71': (3, 2)}

def curve_index(name):
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    for k, cv in enumerate(js):
        if cv['name'].split(' ')[0] == name: return k, cv
    raise KeyError(name)

def base_curve(cvname):
    rat, typ = cvname.split(' ')
    m, n = [int(x) for x in rat.split(',')]
    if typ == 'D1': return R2.d1_curve(m, n)
    return R2.sigma_curve(m, n, typ)

def load_hints(cv, side, mutate=None):
    """the conic points: rank_fisher5_hints.json (already in R5.HINTS: (72,65)) or rank_fisher_theta4_gen_hints.json;
    and, under the key '<key>|xi', another level-3 twist xi of the json's C3 (only (945,248), E' side; see prepare)."""
    key = '%s|%s' % (cv['name'], side)
    H = json.load(open(HINTS_GEN)) if os.path.exists(HINTS_GEN) else {}
    for k_ in (key, key + '|V'):
        if k_ in H: R5.HINTS[k_] = H[k_]
    assert key in R5.HINTS and key + '|V' in R5.HINTS, 'no conic points (hints) for %s' % key
    if mutate == 'hintZ':
        h = dict(R5.HINTS[key]); z0 = Fr(h['Z'][0]) + 1; h['Z'] = [str(z0), h['Z'][1]]; R5.HINTS[key] = h
    xi2 = H.get(key + '|xi')
    return key, (int(xi2) if xi2 is not None else None)

def setup_gen(R, mut=None):
    """rank_fisher_theta4.setup without the assertion (D1, D3) = (-2, 2)."""
    QA, QB = sp.expand(R['QA']), sp.expand(R['QB'])
    for Q_ in (QA, QB):
        assert all(sp.Rational(k).q == 1 for k in sp.Poly(Q_, *VS).coeffs())
    dU, dV = int(R['dU']), int(R['dV'])
    D1, D3 = T.sqfree(dU), T.sqfree(dV)
    m1, m3 = isqrt(dU // D1), isqrt(dV // D3)
    assert m1 * m1 * D1 == dU and m3 * m3 * D3 == dV and D1 != D3
    al1, B1 = T.split_t(R['pushout']['l1'], D1, sp.Rational(1, m1))
    al3, B3 = T.split_t(R['pushout3']['l3'], D3, sp.Rational(1, m3))
    B1 = sp.expand(B1 * m1 * m1); B3 = sp.expand(B3 * m3 * m3)
    f = sp.expand(R['pushout']['f']); fd = sp.expand(R['pushout3']['f'])
    assert sp.expand(al1 ** 2 - D1 * B1 ** 2 - f) == 0, 'f != N(l_1)'
    assert sp.expand(al3 ** 2 - D3 * B3 ** 2 - fd) == 0, 'f_dag != N(l_3)'
    q = sp.expand(R['comp']['q']); cc = Fr(str(R['comp']['c']))
    if mut == 'q': q = sp.expand(q + e * g)
    assert T.in_ideal_C3(sp.expand(f * fd - sp.Rational(cc.numerator, cc.denominator) * q * q), QA, QB), \
        'f f_dag - c q^2 is not in (QA, QB)'
    s0, f0 = T.poly_int(f, VS)
    return dict(QA=QA, QB=QB, D1=D1, D3=D3, al1=al1, B1=B1, al3=al3, B3=B3, f=f, fd=fd, q=q, c=cc,
                s0=Fr(int(s0.p), int(s0.q)), f0=f0)

def pencil_primes(D, known):
    """the primes of the pencil discriminant, factored completely (the cofactor after the known primes is factored
    by PARI and every factor is checked prime with isprime)."""
    A2 = sp.hessian(D['QA'], VS); B2 = sp.hessian(D['QB'], VS)
    x_, y_ = sp.symbols('x_ y_')
    F = sp.Poly(sp.expand((x_ * A2 + y_ * B2).det()), x_, y_)
    co = [int(F.coeff_monomial(x_ ** (4 - i) * y_ ** i)) for i in range(5)]
    disc = int(pari.poldisc(pari('Pol([' + ','.join(str(k) for k in co) + '])')))
    assert disc != 0 and co[0] != 0
    n = abs(disc); ps = []
    for p in known:
        if n % p == 0: ps.append(p)
        while n % p == 0: n //= p
    if n > 1:
        Fa = pari.factor(n)
        for i in range(len(Fa[0])):
            q_ = int(Fa[0][i]); assert pari.isprime(q_), 'incomplete factorisation of the pencil discriminant'
            ps.append(q_)
    return len(str(abs(disc))), sorted(set(ps))

def gf2_solve(rows, rhs, n):
    """all solutions e in F_2^n of rows . e = rhs (rows: lists of 0/1)."""
    M = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for col in range(n):
        k = next((i for i in range(r, len(M)) if M[i][col]), None)
        if k is None: continue
        M[r], M[k] = M[k], M[r]
        for i in range(len(M)):
            if i != r and M[i][col]: M[i] = [(x + y) % 2 for x, y in zip(M[i], M[r])]
        piv.append(col); r += 1
    if any(row[n] and not any(row[:n]) for row in M): return []
    free = [j for j in range(n) if j not in piv]
    sols = []
    for mask in range(1 << len(free)):
        e_ = [0] * n
        for i, j in enumerate(free): e_[j] = mask >> i & 1
        for i, col in enumerate(piv):
            e_[col] = (M[i][n] + sum(M[i][j] * e_[j] for j in free)) % 2
        sols.append(e_)
    return sols

def choose_Xi_lin(M, D, C, places, say=say):
    """Xi in <-1, places> with Xi f0 (P_v) a square for sampled P_v at every place of places and at R: the sampled
    classes of f0 at p are closed up to the affine F_2-span they generate (T.coset_closure), and the condition
    'class(Xi) in that coset' is solved by linear algebra (rank_fisher_theta4.choose_Xi does the same by
    enumerating all 2^(1 + #places) candidates)."""
    K = {p: T.f0_classes(M, D, C, p) for p in places}
    for p in places:
        # a place without a sampled point: more seeds (added for (945,248), where the default seed finds no point
        # at p = 97; it changes nothing where the first seed finds points; a sampler only)
        sd = 2
        while not K[p] and sd <= 8:
            K[p] = T.f0_classes(M, D, C, p, seed=sd); sd += 1
    Rs = T.real_signs(M, D, C=C)
    if not all(K[p] for p in places): return [], K, Rs     # a place without a sampled point: no twist found
    gens = [-1] + list(places)
    n = len(gens)
    rows, rhs = [], []
    for p in places:
        ks = [T.cls_vec(k, p) for k in K[p]]
        k0 = ks[0]
        span = []
        for k in ks[1:]:
            dv = [(x + y) % 2 for x, y in zip(k, k0)]
            # reduce against span
            for b in span:
                lead = next(i for i, x in enumerate(b) if x)
                if dv[lead]: dv = [(x + y) % 2 for x, y in zip(dv, b)]
            if any(dv): span.append(dv)
        dim = len(k0)
        # annihilator of span: vectors a with a . s = 0 for all s in span
        ann = gf2_solve(span, [0] * len(span), dim) if span else [[1 if i == j else 0 for i in range(dim)] for j in range(dim)]
        # a basis of ann (rank)
        basis = []
        for a in ann:
            v = list(a)
            for b in basis:
                lead = next(i for i, x in enumerate(b) if x)
                if v[lead]: v = [(x + y) % 2 for x, y in zip(v, b)]
            if any(v): basis.append(v)
        gv = [T.cls_vec(T.int_class(q_, p), p) for q_ in gens]
        for a in basis:
            rows.append([sum(ai * gi for ai, gi in zip(a, gv[j])) % 2 for j in range(n)])
            rhs.append(sum(ai * ki for ai, ki in zip(a, k0)) % 2)
    rs = set(Rs)
    if rs == {1}: rows.append([1] + [0] * (n - 1)); rhs.append(0)
    elif rs == {-1}: rows.append([1] + [0] * (n - 1)); rhs.append(1)
    elif not rs: return [], K, Rs          # no real point sampled: no twist found (as rank_fisher_theta4.choose_Xi)
    sols = []
    for e_ in gf2_solve(rows, rhs, n):
        Xi = 1
        for j, q_ in enumerate(gens):
            if e_[j]: Xi *= q_
        sols.append(Xi)
    return sorted(sols, key=lambda z: (abs(z), z)), K, Rs

def relabel(say_, D):
    """rank_fisher_theta4.check_g0 prints the (72,65) fields (-2, 2) in its labels; the identities it checks use
    D1, D3 (D['D1'], D['D3']).  This wrapper prints the right labels."""
    def s_(txt):
        txt = txt.replace("N1 N1' = -2 delta1^2", "N1 N1' = (%d) delta1^2" % D['D1'])
        txt = txt.replace("R3 R3' = 2 W^2 delta3^2", "R3 R3' = (%d) W^2 delta3^2" % D['D3'])
        return say_(txt)
    return s_

def f0_degenerate_primes(D):
    """the primes p for which f0 mod p lies in span(QA, QB) mod p (then f0 vanishes on the reduced curve and the
    good-reduction argument does not apply): the primes of the last elementary divisor of the 3 x 10 integer
    matrix of the coefficients of QA, QB, f0: the primes of the gcd of its 3 x 3 minors (exact)."""
    mons = [VS[i] * VS[j] for i in range(4) for j in range(i, 4)]
    rows = [[int(sp.Poly(Q_, *VS).coeff_monomial(m_)) for m_ in mons] for Q_ in (D['QA'], D['QB'], D['f0'].as_expr())]
    from itertools import combinations
    gg = 0
    for cols in combinations(range(10), 3):
        m3 = [[rows[i][j] for j in cols] for i in range(3)]
        det = (m3[0][0] * (m3[1][1] * m3[2][2] - m3[1][2] * m3[2][1]) - m3[0][1] * (m3[1][0] * m3[2][2] - m3[1][2] * m3[2][0])
               + m3[0][2] * (m3[1][0] * m3[2][1] - m3[1][1] * m3[2][0]))
        gg = gcd(gg, det)
    assert gg != 0, 'QA, QB, f0 are linearly dependent'
    return T.primes_of(gg)

def theta4_value(M, D, C, Xi, ys, B0, label='', say=say, kmul=False, kap=None):
    """rank_fisher_theta4.theta4_for with generic labels: returns ({y: product}, S, data)."""
    G = T.build_g0(D, Xi, kap1=kap[0] if kap else None, kap3=kap[1] if kap else None)
    if kmul:
        mul = lambda k, Dq: (k[0] * (1 + Dq) + Dq * k[1] * 2, k[0] * 2 + k[1] * (1 + Dq))
        G = T.build_g0(D, Xi, kap1=mul(G['kap1'], D['D1']), kap3=mul(G['kap3'], D['D3']))
    say('%sC4: Xi W^2 = f0 with Xi = %d; k_1 = %s + %s sqrt(%d) (N = %s mod squares: xi4), k_3 = %s + %s sqrt(%d) '
        '(N = c xi4 mod squares); sigma_1, tau_3 rational (%d/%d, %d/%d digits)'
        % (label, Xi, G['kap1'][0], G['kap1'][1], D['D1'], T.sqfree(G['Nk1']), G['kap3'][0], G['kap3'][1], D['D3'],
           len(str(G['sig1'].numerator)), len(str(G['sig1'].denominator)), len(str(G['tau3'].numerator)),
           len(str(G['tau3'].denominator))))
    T.check_g0(D, G, say=relabel(say, D))
    Nz = T.normalise(G)
    ps, bad = T.degenerate_primes(D, Nz, say=say)
    small = [p for p in range(2, 47) if pari.isprime(p)]
    extra = sorted(set(q_ for y in ys for q_ in T.primes_of(y)))
    fbad = f0_degenerate_primes(D)
    say('  f0 mod p in span(QA, QB) mod p only for p in %s (exact: gcd of the 3 x 3 minors); these are put into S' % fbad)
    S = sorted(set(small + list(B0) + T.primes_of(Xi) + bad + extra + fbad))
    say('  S = {R} u %s' % S)
    terms, classes = T.local_terms(M, D, C, Nz, Xi, S, ys)
    rt, sg = T.real_term(M, D, Nz, Xi, ys, C=C)
    say('  R: a point of C4(R) (interval enclosure), g0^0 there %s'
        % ('> 0' if sg == 1 else ('< 0' if sg == -1 else 'of uncertified sign (all y > 0: term +1)')))
    p4 = next((p for p in S if len(classes[p]) > 1), None)
    assert p4 is not None, '(I4) not shown: g0 might be a constant times a square'
    say('  (I4) at p = %d g0^0 takes the classes %s on C4(Q_p): [B] = T\'' % (p4, classes[p4]))
    out = {}
    for y in ys:
        nontriv = [p for p in S if terms[p][y] == -1] + (['R'] if rt[y] == -1 else [])
        prod = rt[y]
        for p in S: prod *= terms[p][y]
        out[y] = prod
        say('  sum_v (g0(P_v), %d)_v: the term is -1 exactly at %s; product = %+d' % (y, nontriv if nontriv else 'no place', prod))
    for p in S:
        say('    p = %d: classes of g0^0 on C4(Q_p) %s; terms %s' % (p, classes[p], ' '.join('%d:%+d' % (y, terms[p][y]) for y in ys)))
    return out, S, dict(G=G, Nz=Nz)

SIDE = None          # --side=E or --side=E' : reduce on that side instead of theta4_reduction[0]

def prepare(name, mutate=None, quiet=False):
    sy = (lambda s: None) if quiet else say
    k, cv = curve_index(name)
    side, other = cv['theta4_reduction']
    if SIDE is not None and SIDE != side: side, other = other, side
    S = cv['sides'][side]; C = S['level3_curve']
    lifted = int(S['dd']); y = int(cv['sides'][other]['dd'])
    _, xi2 = load_hints(cv, side, mutate=mutate)
    if xi2 is not None:
        # (945,248), E' side: the json's C3 (xi = 1271) has an empty level-4 Selmer set, its twist xi2 = 1271 * 14
        # (14 in S_1 of this side: same level-2 lift D_lam, another member of its level-3 coset) has not
        # (rank_fisher_theta4_norm.py).  Nothing rests on that choice: C4 is checked everywhere locally soluble below.
        C = dict(C, xi=str(xi2))
    R = R5.analyse(cv['name'], side, C, dd=lifted)
    assert R['pushout'] is not None and R['pushout3'] is not None and R['comp'] is not None, \
        'rank_fisher5 steps (6)-(8) not completed (conic points missing or rejected)'
    D = setup_gen(R)
    a, b, Gpt = base_curve(cv['name'])
    if side == "E'":
        aw, bw = a, b
        imgs = [T.sqfree(b), T.sqfree(Gpt[0])]
    else:
        aw, bw = -2 * a, a * a - 4 * b
        imgs = [T.sqfree(bw)]
    RD = RI.first_descent(aw, bw)
    if xi2 is not None:
        assert side == "E'" and T.sqfree(xi2 * int(S['level3_curve']['xi'])) in RD['Sp'], \
            'the level-3 twist xi2 is not in the coset xi S_1 of this side'
    return dict(k=k, cv=cv, side=side, other=other, C=C, lifted=lifted, y=y, R=R, D=D, a=a, b=b, G=Gpt,
                aw=aw, bw=bw, imgs=imgs, RD=RD, xi2=xi2, xi_json=int(S['level3_curve']['xi']))

def main(name):
    t0 = time.time()
    P = prepare(name)
    cv, side, other, C, D, RD = P['cv'], P['side'], P['other'], P['C'], P['D'], P['RD']
    y, imgs = P['y'], P['imgs']
    say('rank_fisher_theta4_gen.py: Fisher\'s Theta_4 on %s (arXiv:1509.03234 Sects. 2, 3, 7), exact' % cv['name'])
    say('E: y^2 = x(x^2 + a x + b), a = %d, b = %d; known point with x = %s' % (P['a'], P['b'], P['G'][0]))
    say('the lifted class %d (%s side) -> C3\' (rank_fisher5 steps (1)-(8) re-run, conic points re-checked); '
        'L = Q(sqrt %d), L\' = Q(sqrt %d)' % (P['lifted'], side, D['D1'], D['D3']))
    if P['xi2'] is not None:
        say('  level-3 lift: the twist xi = %d of the json\'s C3 (xi = %d there; the ratio %d is in S_1 of this side, so '
            'the same D_lam and another member of its level-3 coset); the json\'s C3 has no level-4 twist '
            '(rank_fisher_theta4_norm.py); C4 below is checked everywhere locally soluble'
            % (P['xi2'], P['xi_json'], T.sqfree(P['xi2'] * P['xi_json'])))
    say('working curve (whose x-classes are the paired group): y^2 = x(x^2 + %d x + %d) (%s); paired class y = %d'
        % (P['aw'], P['bw'], 'E' if side == "E'" else "E'", y))
    for p in RD['bad']:
        assert R2.certified(y, P['aw'], P['bw'], p, RD['locS']), 'y is not certified in the local image at %s' % p
    assert RI.real_ok(y, P['aw'], P['bw'])
    say('  y is in the local image at R and at the bad primes %s (rank_indep certified images); y in S\'_1, and '
        'S\'_4 = S\'_2 contains y (rank_fisher.py, levels 2-3)' % RD['bad'])
    nd, ps = pencil_primes(D, known=sorted(set([2] + RD['bad'])))
    B0 = sorted(set([2] + RD['bad'] + ps))
    say('  pencil discriminant (%d digits): primes %s; B0 = %s' % (nd, ps, B0))
    M = T.Model(D['QA'], D['QB'])
    places = sorted(set([2, 3, 5, 7, 11, 13] + B0))
    sols, K, Rs = choose_Xi_lin(M, D, C, places)
    say('  level-4 twists Xi (sampled classes closed up to cosets, F_2 linear algebra) over R and %s: %d: %s'
        % (places, len(sols), sols[:16]))
    assert sols, 'no level-4 twist found'
    ys = [y] + [z for z in imgs if z != 1]
    res, S, dat = theta4_value(M, D, C, sols[0], ys, B0)
    for z in imgs:
        if z != 1: assert res[z] == 1, 'CONTROL (K2) FAILED: the class %d of a rational point pairs non-trivially' % z
    say('  (K2) the classes %s of the rational points pair to 0, as they must' % [z for z in imgs if z != 1])
    th = res[y]
    for Xi2 in sols[1:3]:
        r2, _, _ = theta4_value(M, D, C, Xi2, ys, B0, label='(K3) ')
        assert all(r2[z] == 1 for z in imgs if z != 1), 'CONTROL (K2) FAILED for the twist %d' % Xi2
        assert r2[y] == th, 'CONTROL (K3) FAILED: the value depends on the level-4 lift'
    say('  (K3) the same value for the level-4 lifts %s' % sols[1:3])
    r5, _, _ = theta4_value(M, D, C, sols[0], ys, B0, label='(K5) other k_1, k_3: ', kmul=True)
    assert r5[y] == th and all(r5[z] == 1 for z in imgs if z != 1), 'CONTROL (K5) FAILED'
    say('  (K5) the same value with k_1 (1 + sqrt D1)^2, k_3 (1 + sqrt D3)^2 (another g0)')
    say('')
    say('RESULT: Theta_4(%d, %d) on %s = %s (the product of the local Hilbert symbols is %+d)'
        % (P['lifted'], y, cv['name'], '1/2 (non-zero)' if th == -1 else '0', th))
    if th == -1:
        dS, dSp = DIMS[cv['name'].split(' ')[0]]
        say('  => the left and right kernels S_5, S\'_5 drop by one each: dim %d and %d (rank_fisher.py: dim S_4 = %d, '
            'dim S\'_4 = %d, each the known classes and one open class), and Lemma 2.2 (n = 3) gives rank <= %d + %d - 2 '
            '= 1; with the known point of infinite order, rank = 1 for %s, without Magma (written arguments as listed '
            'in the docstring)' % (dS - 1, dSp - 1, dS, dSp, dS - 1, dSp - 1, cv['name']))
    else:
        say('  => Theta_4 = 0: the bound stays 3 at this level')
    return th

def negative(name):
    say('rank_fisher_theta4_gen.py --negative on %s: altered inputs that must be refused' % name)
    n = ok = 0
    def refuse(label, fn):
        nonlocal n, ok
        n += 1
        try: fn()
        except AssertionError as ex:
            msg = ex.args[0] if ex.args else type(ex).__name__
            msg = msg[0] if isinstance(msg, tuple) else msg
            msg = str(msg).replace("N1 N1' = -2 delta1^2", "N1 N1' = D1 delta1^2").replace("R3 R3' = 2 W^2 delta3^2", "R3 R3' = D3 W^2 delta3^2")
            ok += 1; say('  refused: %s (%s)' % (label, msg[:90])); return
        say('  NOT REFUSED: %s' % label)
    quiet = lambda s: None
    refuse('a conic point moved off its cone (Z -> Z + 1 in the hint)', lambda: prepare(name, mutate='hintZ', quiet=True))
    R5.HINTS.clear(); R5.HINTS.update(json.load(open(R5.HINTS_FILE)) if os.path.exists(R5.HINTS_FILE) else {})
    P = prepare(name, quiet=True)
    D, C, y = P['D'], P['C'], P['y']
    refuse('q replaced by q + e g (f f_dag - c q^2 leaves the ideal)', lambda: setup_gen(P['R'], mut='q'))
    M = T.Model(D['QA'], D['QB'])
    RD = P['RD']
    nd, ps = pencil_primes(D, known=sorted(set([2] + RD['bad'])))
    B0 = sorted(set([2] + RD['bad'] + ps))
    places = sorted(set([2, 3, 5, 7, 11, 13] + B0))
    sols, K, Rs = choose_Xi_lin(M, D, C, places, say=quiet)
    Xi = sols[0]
    for m_, what in (('sigma', 'sigma_1 doubled in N1'), ('tau', 'tau_3 tripled in R3')):
        refuse(what, lambda m_=m_: T.check_g0(D, T.build_g0(D, Xi, mut=m_), say=quiet))
    refuse('k_1 = 1 + sqrt(D1) (norm 1 - D1, not xi4 modulo squares in general): sigma_1 not rational',
           lambda: T.build_g0(D, Xi, kap1=(Fr(1), Fr(1))))
    # a forged twist: Xi * q for the first prime q of places with Xi q outside the computed set
    def forged():
        bad_q = next(q_ for q_ in places if Xi * q_ not in sols)
        Kp = {p: T.f0_classes(M, D, C, p) for p in places}
        for p in places:
            Kc = T.coset_closure(Kp[p], p)
            assert T.int_class(Xi * bad_q, p) in Kc, 'no point of %d W^2 = f0 over Q_%d (sampled classes)' % (Xi * bad_q, p)
    refuse('a forged twist Xi q outside the Selmer coset', forged)
    def wrong_y():
        z = next(q_ for q_ in (3, 5, 7, 11, 13, 17, 19, 23)
                 if not all(R2.certified(q_, P['aw'], P['bw'], p, RD['locS']) for p in RD['bad']))
        for p in RD['bad']:
            assert R2.certified(z, P['aw'], P['bw'], p, RD['locS']), '%d is not in the local image at %s' % (z, p)
    refuse('a wrong paired class (not Selmer)', wrong_y)
    G = T.build_g0(D, Xi); Nz = T.normalise(G)
    Nbad = dict(Nz, R30=sp.Poly(sp.expand(Nz['R30'].as_expr() * c), *(VS + (W,))),
                N10=sp.Poly(sp.expand(Nz['N10'].as_expr() * e), *(VS + (W,))))
    # places where y is not a local square (elsewhere every symbol (., y)_p is +1 and nothing can be seen)
    cand = sorted(set([q_ for q_ in range(2, 200) if pari.isprime(q_)] + places))
    pl = [q_ for q_ in cand if T.int_class(y, q_) != (0, 1 if q_ == 2 else 0)][:6]
    refuse('g0 replaced by g0 c / e (not a pushout form): the local term depends on the point (at %s)' % pl,
           lambda: T.local_terms(M, D, C, Nbad, Xi, pl, [y], say=quiet))
    def low_prec():
        old = T.prec
        T.prec = lambda p: 4
        try: T.local_terms(M, D, C, Nz, Xi, [2], [y], say=quiet)
        finally: T.prec = old
    refuse('the p-adic precision lowered to 4 digits at p = 2', low_prec)
    say('negative tests: %d of %d refused' % (ok, n))
    assert ok == n
    return ok, n

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    name = args[0] if args else '72,65'
    tag = name.replace(',', '_')
    sides = [a[7:] for a in sys.argv if a.startswith('--side=')]
    if sides:
        SIDE = sides[0]; tag += '_' + SIDE.replace("'", 'prime')
    if '--negative' in sys.argv:
        negative(name)
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_%s_negative.txt' % tag), 'w').write('\n'.join(OUT) + '\n')
    else:
        main(name)
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_%s.txt' % tag), 'w').write('\n'.join(OUT) + '\n')
