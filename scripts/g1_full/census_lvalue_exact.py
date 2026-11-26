#!/usr/bin/env python3
"""L(E,1) > 0 without PARI for the closing curves of the 19 census closures of G1_FULL_SYSTEM.md Result 3 that "use
L(E,1) != 0 (Kolyvagin) where 2-descent gives only [0, 2]".  Written by an AI model; not read by a person.

The pairs and curves: census_lvalue_identify.py (PARI-guided) recovers the 19 pairs from the history of
open_pairs_1000.txt and finds, for each, the curve of global.py's twelve on which ellrank gives [0, 2], w = +1,
N <= 10^14 and PARI's lfun(E,1) is clearly nonzero, with no torsion point giving a solution.  The table CLOSURES below
copies its result (one curve per pair: the smallest conductor).  The model is global.py's [0, c2, 0, c1 c3, c0 c3^2]
for y^2 = P(w), rebuilt here from (m, n) in pure Python (checked against global.py with --pari).

The certificate is the machinery of ../second_descent/lball_exact53.py and lball_exact.py, imported, not changed:
  N by Tate's algorithm (lball_exact.local_data, exact); a_p at bad p from Tate, at good p < 1000 by counting points,
  above by the exact order certificate of lball_exact53 (Hasse + Lagrange, unique multiple), counting as fallback;
  a_n by lball_exact.an_from_ap; w = +1 from F = sum a_n exp(-2 pi n/sqrt N) != 0 (w = -1 forces F = 0);
  L(E,1) by lball_exact.evaluate, integer interval arithmetic with directed rounding, tail bound |a_n| <= 2n.
  The guard (the enclosures at the cut-offs y = 1 and y = 6/5 meet, and the Fricke residual at 6/5 contains 0) is
  asserted: it is what refuses a wrong N, w or a_p (Tate's algorithm itself is trusted).
Certified line: L(E,1) >= lower > 0 with N and the a_n computed here and w = +1 from F != 0.  With modularity, Carayol,
Atkin-Lehner and Kolyvagin (or Kato), rank E(Q) = 0.
The torsion step (--torsion, since 26 November 2026; before, global.py's elltors): E(Q)_tors of each of the 19 curves
by torsion_exact.py (#E(Q)_tors divides B = gcd #E(F_p) over the good odd p < 1000, reduction being injective on
torsion for odd p of good reduction (quoted: Silverman VII.3.1, IV.6.4); the torsion points exhibited exactly by
halving; when they number fewer than B, the rational 2-torsion is counted exactly and each point of T with no half in
T is shown not to be a double modulo some good p), then global.py's own classify_t and analyse (compiled from
global.py's source, without its PARI function torsion_w) on every torsion point, with rank 0 from the certificate
above.  So the 19 closures use PARI nowhere.  PARI's elltors agrees on all 19 (--torsion --pari, a check, not used).

usage: census_lvalue_exact.py --list                      the 19 curves, N (Tate), w (local), terms needed
       census_lvalue_exact.py [--pari] [--only i,j | --pairs 'm,n;m,n' | --maxN X] [--jsonl out.jsonl]   certify (index in the N-sorted list)
       census_lvalue_exact.py --merge a.jsonl ... [--write] census_lvalue_exact.txt from the JSON lines
       census_lvalue_exact.py --negative [--write]          forged N, w, a_p routed through lball_exact.evaluate: refused
       census_lvalue_exact.py --torsion [--pari] [--write]  the torsion step, census_lvalue_exact_torsion.txt (--pari: elltors check)
       census_lvalue_exact.py --torsion-negative [--write]  forged torsion certificates refused, larger torsion recognised,
                                                            a planted solution flagged: census_lvalue_exact_torsion_negative.txt
(700,529), N = 8.7e13, is done in chunks by census_lvalue_big.py (same certificate, a_n by lball_exact53.an_array); its
record is merged here like the others.
Outputs: census_lvalue_exact.txt, census_lvalue_exact.jsonl, census_lvalue_exact_timing.txt, census_lvalue_exact_negative.txt,
census_lvalue_exact_torsion.txt, census_lvalue_exact_torsion_negative.txt."""
import sys, os, json, math, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'second_descent'))
import lball_exact as LX
import lball_exact53 as L53

# pair -> closing curve (census_lvalue_identify.txt); index order is fixed by sorting on N below
CLOSURES = [((112, 15), 'w*K'), ((127, 97), 'w*f0*f8'), ((240, 161), 'w*K'), ((260, 69), 'w*K'), ((288, 175), 'w*K'),
            ((329, 191), 'w*f0*f8'), ((352, 135), 'w*K'), ((401, 79), 'w*f0*f8'), ((463, 113), 'w*f0*f8'), ((487, 217), 'w*f0*f8'),
            ((527, 336), 'w*K'), ((544, 33), 'w*K'), ((575, 48), 'w*f0*f8'), ((577, 511), 'w*f0*f8'), ((623, 527), 'w*K'),
            ((697, 696), 'w*f0*f8'), ((700, 529), 'w*f0*f8'), ((863, 191), 'w*f0*f8'), ((943, 576), 'w*f0*f8')]
NMAX = 10**14

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r

def poly(m, n, name):
    """global.curves(m, n)[name], for the products used here"""
    W = [0, 1]; f0 = [m*m, n*n]; f8 = [n*n, m*m]; S = m*m + n*n; g4 = [2*S, 2*S]
    K = [S*S, 2*S*S - 16*m*m*n*n, S*S]
    fac = {'w': W, 'f0': f0, 'f8': f8, 'g4': g4, 'K': K}
    P = [1]
    for x in name.split('*'): P = pmul(P, fac[x])
    return P

def model(P):
    """global.torsion_w's model of y^2 = P(w): cubic c0..c3, or 0, c1..c4 (w = 1/u)"""
    if len(P) == 5: assert P[0] == 0; c = [P[4], P[3], P[2], P[1]]
    else: c = P
    c0, c1, c2, c3 = c
    return (0, c2, 0, c1 * c3, c0 * c3 * c3)

# ---------------------------------------------------------------- exact 2-isogeny chains (partner pairs)
def c46(a, b):
    """c4, c6 of y^2 = x^3 + a x^2 + b x"""
    return 16*a*a - 48*b, -64*a**3 + 288*a*b

def isomorphic(E1, E2):
    """y^2 = x^3 + a x^2 + b x, Q-isomorphic iff c4' = u^4 c4, c6' = u^6 c6 for some u in Q* (c4 c6 != 0 here)"""
    from fractions import Fraction as Fr
    (c4, c6), (d4, d6) = c46(*E1), c46(*E2)
    assert c4 and c6 and d4 and d6
    s = Fr(d6 * c4, c6 * d4)                                  # = u^2
    if s <= 0: return False
    r = (math.isqrt(s.numerator), math.isqrt(s.denominator))
    if r[0]**2 != s.numerator or r[1]**2 != s.denominator: return False
    return d4 == s*s*c4 and d6 == s**3*c6

def two_isogenous(E):
    """every curve 2-isogenous to y^2 = x(x^2 + a x + b) over Q: for each rational root r of x(x^2 + a x + b), move
    (r, 0) to (0, 0) (a' = 3r + a, b' = 3r^2 + 2ar + b) and apply Velu: y^2 = x(x^2 - 2a' x + a'^2 - 4b')"""
    a, b = E; roots = [0]; D = a*a - 4*b
    if D >= 0 and math.isqrt(D)**2 == D:
        d = math.isqrt(D)
        for num in (-a + d, -a - d):
            assert num % 2 == 0; roots.append(num // 2)
    out = []
    for r in sorted(set(roots)):
        a1 = 3*r + a; b1 = 3*r*r + 2*a*r + b
        assert b1 != 0 and a1*a1 - 4*b1 != 0
        out.append((-2*a1, a1*a1 - 4*b1))
    return out

def isogeny_chain(src, dst, depth=4):
    """a chain of rational 2-isogenies from src to a curve Q-isomorphic to dst (both y^2 = x^3 + a x^2 + b x), by
    breadth-first search; None if there is none of length <= depth"""
    frontier = [(src, [src])]
    for k in range(depth + 1):
        nxt = []
        for E, path in frontier:
            if isomorphic(E, dst): return path
            for F in two_isogenous(E): nxt.append((F, path + [F]))
        frontier = nxt
    return None

# ---------------------------------------------------------------- the torsion step without PARI (--torsion)
import torsion_exact as TE

def global_py():
    """global.py's own pure-Python functions (pmul, issq_int, issq, sqrt_fr, curves, classify_t, analyse), compiled from
    its source text without running its module code (which imports cypari2).  torsion_w, its one PARI function, is
    not loaded: torsion_step supplies it."""
    import ast
    from fractions import Fraction
    tree = ast.parse(open(os.path.join(HERE, 'global.py')).read())
    keep = ('pmul', 'issq_int', 'issq', 'sqrt_fr', 'curves', 'classify_t', 'analyse')
    ns = {'Fr': Fraction, 'math': math, 'json': json, 'sys': sys}; got = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in keep:
            exec(compile(ast.Module([node], []), 'global.py', 'exec'), ns); got.append(node.name)
    assert sorted(got) == sorted(keep), got
    return ns

def torsion_ws(P, pts):
    """the w-values of the points pts of the model of y^2 = P(w): the lines of global.torsion_w that map a point
    (X, Y) of y^2 = X^3 + c2 X^2 + c1 c3 X + c0 c3^2 to w (O -> w = infinity (None) or 0 for the quartic case)"""
    from fractions import Fraction as Fr
    if len(P) == 5: assert P[0] == 0; c = [P[4], P[3], P[2], P[1]]; inv = True
    else: c = P; inv = False
    c3 = c[3]; ws = []
    for Q in pts:
        if Q is None: ws.append(Fr(0) if inv else None)
        else:
            x = Fr(Q[0]) / c3
            ws.append((None if x == 0 else 1 / x) if inv else x)
    return ws

def torsion_step(pair, name, a, G, T=None, primes=None):
    """E(Q)_tors of the closing curve, certified by torsion_exact (gcd of #E(F_p) over the good odd p < 1000, the
    points exhibited, the 2-part completed by the halving test modulo p), then global.py's own analyse() with this
    torsion in place of PARI's elltors and rank 0 (from the L-value certificate) in place of ellrank.  T, primes:
    forged inputs for the negative tests."""
    m, n = pair; P = poly(m, n, name)
    assert G['curves'](m, n)[name] == P, 'poly() differs from global.curves'
    E = TE.Curve(a[1], a[3], a[4]); assert tuple(a) == E.model()
    if primes is None: primes = TE.good_odd_primes(E, 1000)
    if T is None: T = TE.find_torsion(E, TE.bound(E, primes))
    cert = TE.certify(E, T, primes)
    rec = dict(pair=pair, curve=name, cert=cert, T=T)
    if not cert['ok']: return rec
    ws = torsion_ws(P, T)
    def torsion_w(poly_):
        if poly_ == P: return 0, 0, ws, len(T)       # rank 0 by the L-value certificate; ws as global.torsion_w
        return 0, 1, [], 0                           # the other eleven curves of global.curves: not used here
    G['torsion_w'] = torsion_w
    res = [x for x in G['analyse'](m, n) if x['curve'] == name]; assert len(res) == 1
    verdicts = {}
    for w in ws:
        if w is None: verdicts['w = inf'] = 'degenerate'; continue
        if w < 0 or not G['issq'](w): continue
        t = G['sqrt_fr'](w)
        for tt in {t, -t}: verdicts['t = %s' % tt] = G['classify_t'](m, n, tt)
    rec.update(ws=ws, closes=res[0]['closes'], alert=res[0]['alert'], verdicts=verdicts)
    return rec

def torsion_line(rec, nprimes):
    c = rec['cert']; pair = rec['pair']
    if not c['ok']: return '%-8s %-8s torsion NOT CERTIFIED: %s' % ('%d %d' % pair, rec['curve'], '; '.join(c['reasons']))
    xs = {}
    for Q in rec['T']:
        if Q is not None: xs[Q[0]] = xs.get(Q[0], 0) + 1
    l1 = '%-8s %-8s B = %d (gcd of #E(F_p), %d good odd p < 1000); E(Q)_tors = %s, order %d: %s' % (
        '%d %d' % pair, rec['curve'], c['B'], nprimes, TE.fmt_structure(c['structure']), c['order'], '; '.join(c['reasons']))
    l2 = '                  x of the torsion points (number of points): %s' % ', '.join('%s (%d)' % (x, k) for x, k in sorted(xs.items()))
    l3 = '                  torsion t (w = t^2) and global.py classify_t: %s; global.py analyse: %s' % (
        ', '.join('%s %s' % kv for kv in sorted(rec['verdicts'].items())) or 'none',
        'closes' if rec['closes'] and not rec['alert'] else ('SOLUTION' if rec['alert'] else 'does NOT close'))
    return '\n'.join((l1, l2, l3))

def torsion_all(write, pari):
    G = global_py(); rows = the_curves()
    rank0 = {}                                       # rank 0 of each curve: its line of census_lvalue_exact.txt
    for l in open(os.path.join(HERE, 'census_lvalue_exact.txt')):
        w = l.split()
        if len(w) > 3 and w[0].isdigit() and w[1].isdigit():
            how = 'by the series' if 'certified positive' in l and 'guard passes' in l else ('by isogeny' if l.rstrip().endswith('certified by isogeny') else None)
            rank0[((int(w[0]), int(w[1])), w[2])] = how
    lines = ['census_lvalue_exact.py --torsion: E(Q)_tors of the 19 closing curves without PARI (torsion_exact.py), and', 
             'global.py\'s classify_t/analyse (loaded from global.py\'s source) on every torsion point', '']
    done = 0
    for N, pair, name, a, LD in rows:
        rec = torsion_step(pair, name, a, G)
        E = TE.Curve(a[1], a[3], a[4]); nprimes = len(TE.good_odd_primes(E, 1000))
        assert rank0.get((pair, name)), ('rank 0 of this curve is not certified in census_lvalue_exact.txt', pair, name)
        lines.append(torsion_line(rec, nprimes) + '; rank 0: census_lvalue_exact.txt, %s' % rank0[(pair, name)])
        assert rec['cert']['ok'], ('torsion not certified', pair)
        assert not rec['alert'], ('a torsion point gives a SOLUTION', pair)
        assert rec['closes'], ('global.py analyse does not close', pair)
        done += 1
        if pari:                                     # a check, not used: global.py's torsion_w (PARI elltors)
            import importlib
            sys.path.insert(0, HERE); GP = importlib.import_module('global')
            lo, hi, wsp, Tp = GP.torsion_w(poly(pair[0], pair[1], name))
            same = Tp == len(rec['T']) and sorted(wsp, key=str) == sorted(rec['ws'], key=str)
            print('  PARI elltors (check, not used): #E(Q)_tors = %d, the same torsion w: %s' % (Tp, same), flush=True)
            assert same
    lines.append('torsion step without PARI: %d of %d curves certified; no torsion point gives a solution; with rank 0 '
                 '(census_lvalue_exact.txt) global.py\'s analyse closes each pair' % (done, len(rows)))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    if write: open(os.path.join(HERE, 'census_lvalue_exact_torsion.txt'), 'w').write(out)

def torsion_negative(write):
    """forged torsion certificates and inputs, each must be refused; curves of larger torsion recognised; a planted
    solution flagged"""
    G = global_py(); rows = {tuple(r[1]): r for r in the_curves()}
    lines = ['census_lvalue_exact.py --torsion-negative: forged torsion inputs, each must be refused; controls', '']
    nref = ntot = 0
    def refused(label, rec):
        nonlocal nref, ntot
        ok = rec['cert']['ok']; ntot += 1; nref += not ok
        lines.append('forged %-58s %s' % (label + ':', 'NOT REFUSED' if ok else 'refused (%s)' % rec['cert']['reasons'][-1]))
    # a Z/4 curve and a Z/2 x Z/8 curve of the census
    for pair in ((112, 15), (575, 48)):
        N, _, name, a, LD = rows[pair]; E = TE.Curve(a[1], a[3], a[4]); primes = TE.good_odd_primes(E, 1000)
        T = TE.find_torsion(E, TE.bound(E, primes)); P4 = [Q for Q in T if Q is not None and E.order(Q, 16) >= 4][0]
        tag = '(%d,%d) %s' % (pair[0], pair[1], name)
        rec = torsion_step(pair, name, a, G); lines.append('true certificate %s: E(Q)_tors = %s, certified %s' % (tag, TE.fmt_structure(rec['cert']['structure']), rec['cert']['ok']))
        assert rec['cert']['ok']
        refused('%s, an extra point (x, y + 1) off the curve' % tag, torsion_step(pair, name, a, G, T=T + [(P4[0], P4[1] + 1)]))
        refused('%s, a torsion point left out' % tag, torsion_step(pair, name, a, G, T=[Q for Q in T if Q != P4]))
        refused('%s, E(Q)_tors claimed to be E(Q)[2]' % tag, torsion_step(pair, name, a, G, T=[Q for Q in T if E.order(Q, 2) is not None]))
        refused('%s, the bad prime %d among the primes' % (tag, min(p for p in LD['rows'] if p > 2)),
                torsion_step(pair, name, a, G, T=T, primes=sorted(primes[:20] + [min(p for p in LD['rows'] if p > 2)])))
        refused('%s, p = 2 among the primes' % tag, torsion_step(pair, name, a, G, T=T, primes=[2] + primes[:20]))
    # a point of infinite order listed as torsion: y^2 = x^3 - 2x, (2, 2) (2P has x = 9/4: not integral)
    E = TE.Curve(0, -2, 0); pr = TE.good_odd_primes(E, 400)
    for lab, T in (('y^2 = x^3 - 2x: (2, 2) (infinite order) with O, (0, 0)', [None, (0, 0), (2, 2), (2, -2)]),
                   ('y^2 = x^3 - 2x: the multiples kP, |k| <= 3, of P = (2, 2)', [E.mul(k, (2, 2)) for k in range(-3, 4)])):
        c = TE.certify(E, T, pr); ntot += 1; nref += not c['ok']
        lines.append('forged %-58s %s' % (lab + ':', 'NOT REFUSED' if c['ok'] else 'refused (%s)' % c['reasons'][-1]))
    # curves with larger torsion (Cremona's tables): recognised
    for lab, (a2, a4, a6), want in TE.CONTROLS:
        E = TE.Curve(a2, a4, a6); T, c = TE.torsion(E, TE.good_odd_primes(E, 400))
        ok = c['ok'] and c['structure'] == want
        lines.append('control %-6s %-30s E(Q)_tors = %-12s (tables: %s): %s' % (lab, E.model(), TE.fmt_structure(c['structure']) if c['ok'] else '-',
                     TE.fmt_structure(want), 'recognised' if ok else 'NOT RECOGNISED'))
        assert ok
    E = TE.Curve(-4, -125120, -16869104); T, c = TE.torsion(E, TE.good_odd_primes(E, 400)); ntot += 1; nref += not c['ok']
    lines.append('control 11a2 (trivial torsion, B = %d from its 5-isogeny): %s' % (c['B'], 'refused, no claim (%s)' % c['reasons'][-1] if not c['ok'] else 'CERTIFIED (wrong)'))
    # a planted solution: no solution is known, so one is simulated by a forged square test in global.py's namespace
    # (on (575,48), whose torsion has w-values other than 0, 1, infinity)
    N, pair, name, a, LD = rows[(575, 48)]
    Gf = global_py(); Gf['issq'] = lambda q: True
    rec = torsion_step(pair, name, a, Gf); ntot += 1; nref += rec['alert'] and not rec['closes']
    lines.append('planted solution (575,48) w*f0*f8, global.py\'s issq forged to accept every value: classify_t %s; analyse: %s' % (
        ', '.join('%s %s' % kv for kv in sorted(rec['verdicts'].items())), 'flags SOLUTION, does not close' if rec['alert'] and not rec['closes'] else 'NOT FLAGGED'))
    try: assert not rec['alert'], ('a torsion point gives a SOLUTION', pair); fired = False     # the assertion of --torsion
    except AssertionError as ex: fired = True; lines.append('  the assertion of --torsion fires: %s' % (ex.args[0],))
    ntot += 1; nref += fired
    lines.append('negative tests: %d forged inputs, %d refused' % (ntot, nref))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert nref == ntot
    if write: open(os.path.join(HERE, 'census_lvalue_exact_torsion_negative.txt'), 'w').write(out)

def the_curves():
    out = []
    for pair, name in CLOSURES:
        a = model(poly(pair[0], pair[1], name)); LD = LX.local_data(a)
        out.append((LD['N'], pair, name, a, LD))
    out.sort(key=lambda t: (t[0], t[1]))
    return out

def fN(LD):
    return ' * '.join('%d^%d' % (p, T['f']) if T['f'] > 1 else str(p) for p, T in LD['rows'].items() if T['f'])

def certify(idx, N, pair, name, a, LD, pari):
    t = {}; t0 = time.time()
    Mmax = LX.max_terms(N)
    ap, st = L53.ap_table_fast(a, LD, Mmax); t['ap'] = time.time() - t0
    an = LX.an_from_ap(ap, LD, Mmax); t['an'] = time.time() - t0 - t['ap']
    PI = LX.pi_interval(); t1 = time.time()
    ev = LX.evaluate(an, N, 1, PI, Mmax); t['series'] = time.time() - t1
    F = ev['F']; wF = 1 if F[0] > 0 or F[1] < 0 else None
    L1 = ev['L_y1']; L2 = ev['L_y2']
    r = dict(idx=idx, pair=list(pair), curve=name, model=list(a), N=N, fN=fN(LD), w_local=LD['w_local'], w_fricke=wF,
             F=[LX.lower_str(F[0], 3), LX.upper_str(F[1], 3)], terms=list(ev['M']), Mmax=Mmax,
             L_y1=[L53.dec(L1[0]), L53.dec(L1[1], up=True)], L_y2=[L53.dec(L2[0]), L53.dec(L2[1], up=True)],
             guard=ev['guard'], positive=(L1[0] > 0 and wF == 1), ap_stats=st)
    if pari:
        import importlib
        t2 = time.time(); P = L53.get_pari()
        sys.path.insert(0, HERE); G = importlib.import_module('global')
        gm = G.curves(pair[0], pair[1])[name]
        E = P.ellminimalmodel(P.ellinit(list(a)))
        Np = int(P.ellglobalred(E)[0]); wp = int(P.ellrootno(E))
        pa = [int(x) for x in P.ellan(E, Mmax)]
        Lp = float(P.lfun(P.lfuncreate(E), 1))
        from fractions import Fraction as Fr
        r['pari'] = dict(poly=gm == poly(pair[0], pair[1], name), N=Np == N, w=wp == wF, an=pa == an[1:Mmax + 1], n=Mmax,
                         lfun=Lp, lfun_inside=Fr(r['L_y1'][0]) <= Fr(Lp) <= Fr(r['L_y1'][1]) or abs(Lp - float(r['L_y1'][0])) < 1e-8)
        t['pari'] = time.time() - t2
    t['total'] = time.time() - t0
    r['time'] = {k: round(v, 1) for k, v in t.items()}
    return r

def line(r):
    return ('%-8s %-8s N = %d = %s  w = %s (F in [%s, %s]; local %s)  terms %d  L(E,1) in [%s, %s]  guard %s  %s' % (
        '%d %d' % tuple(r['pair']), r['curve'], r['N'], r['fN'], '+1' if r['w_fricke'] == 1 else 'NOT DETERMINED', r['F'][0], r['F'][1],
        {1: '+1', -1: '-1', None: 'not implemented'}[r['w_local']], r['terms'][0], r['L_y1'][0], r['L_y1'][1],
        'passes' if r['guard'] else 'FAILS', 'certified positive' if r['positive'] and r['guard'] else 'NOT CERTIFIED'))

def run(sel, pari, jsonl):
    rows = the_curves(); fh = open(jsonl, 'a') if jsonl else None
    for i in sel:
        N, pair, name, a, LD = rows[i]
        r = certify(i, N, pair, name, a, LD, pari)
        print(line(r), flush=True)
        if fh: fh.write(json.dumps(r) + '\n'); fh.flush()
        assert r['positive'] and r['guard'], ('NOT CERTIFIED', pair)
        if 'pari' in r: assert all(r['pari'][k] for k in ('poly', 'N', 'w', 'an', 'lfun_inside')), ('PARI cross-check disagrees', pair, r['pari'])

def merge(paths, write):
    recs = {}; rows = the_curves(); index = {tuple(r[1]): i for i, r in enumerate(rows)}
    for p in paths:
        for l in open(p):
            r = json.loads(l); i = index[tuple(r['pair'])]; r['idx'] = i
            assert r['curve'] == rows[i][2] and r['N'] == rows[i][0] and tuple(r['model']) == tuple(rows[i][3]) and r['positive'] and r['guard'], ('record does not match', r['pair'])
            if 'pari' in r: assert all(r['pari'][k] for k in ('N', 'w', 'an')), ('PARI cross-check disagrees', r['pair'])
            recs[i] = r
    # check, not used: PARI's lfun(E,1) from census_lvalue_identify.json (a float) must lie in the enclosure (to 1e-9)
    from fractions import Fraction as Fr
    ident = {}
    for r in json.load(open(os.path.join(HERE, 'census_lvalue_identify.json'))):
        for c in r['curves']:
            if c.get('L') is not None: ident[(tuple(r['pair']), c['curve'])] = (c['L'], c['N'], c['w'])
    lines = []
    for i in range(len(rows)):
        if i in recs:
            r = recs[i]; L, Np, wp = ident[(tuple(r['pair']), r['curve'])]
            inside = Fr(r['L_y1'][0]) - Fr(1, 10**9) <= Fr(L) <= Fr(r['L_y1'][1]) + Fr(1, 10**9)
            assert inside and Np == r['N'] and wp == r['w_fricke'], ('PARI identification disagrees', r['pair'])
            px = ''
            if 'pari' in r and 'primes' in r['pari']: px = '; PARI ellap at all %d primes p <= %d and ellan to 10^6 agree' % (r['pari']['primes'], r['pari']['n'])
            elif 'pari' in r: px = '; PARI ellan: all %d a_n agree' % r['pari']['n']
            lines.append(line(r) + '  (PARI lfun %.9f, N, w agree%s)' % (L, px))
    nd = [i for i in range(len(rows)) if i not in recs]
    viaiso = []
    for i in nd:                 # a curve linked by an exact chain of 2-isogenies to a curve certified above
        N, pair, name, a, LD = rows[i]
        for j in sorted(recs):
            b = recs[j]['model']
            if recs[j]['N'] != N: continue
            ch = isogeny_chain((b[1], b[3]), (a[1], a[3]))
            if ch is not None:
                lines.append('%-8s %-8s N = %d  rank 0: a chain of %d rational 2-isogenies (Velu, exact) from the curve of (%d,%d) certified above '
                             '(L(E,1) the same)  certified by isogeny' % ('%d %d' % pair, name, N, len(ch) - 1, recs[j]['pair'][0], recs[j]['pair'][1]))
                viaiso.append(i); break
    nd = [i for i in nd if i not in viaiso]
    lines.append('certified %d of %d by the series (L(E,1) > 0, w = +1, guard passes), %d more by an exact isogeny chain to one of them; not done: %s' % (
        sum(1 for r in recs.values() if r['positive'] and r['guard']), len(rows), len(viaiso),
        ', '.join('(%d,%d) %s N=%d' % (rows[i][1] + (rows[i][2], rows[i][0])) for i in nd) or 'none'))
    tl = ['pair  curve  N  terms(max)  a_p: bad/counted/certified/fallback/points  seconds: a_p, a_n, series, PARI, total']
    for i in sorted(recs):
        r = recs[i]; s = r['ap_stats']; tm = r['time']
        tl.append('%d %d  %s  %d  %d(%d)  %d/%d/%d/%d/%d  %s, %s, %s, %s, %s' % (r['pair'][0], r['pair'][1], r['curve'], r['N'], r['terms'][0], r['Mmax'],
                  s['bad'], s['counted'], s['certified'], s['fallback'], s['points'], tm.get('ap'), tm.get('an'), tm.get('series'), tm.get('pari', '-'), tm.get('total')))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out); sys.stdout.write('\n'.join(tl) + '\n')
    if write:
        open(os.path.join(HERE, 'census_lvalue_exact.txt'), 'w').write(out)
        open(os.path.join(HERE, 'census_lvalue_exact_timing.txt'), 'w').write('\n'.join(tl) + '\n')
        with open(os.path.join(HERE, 'census_lvalue_exact.jsonl'), 'w') as fh:
            for i in sorted(recs): fh.write(json.dumps(recs[i]) + '\n')

def negative(write):
    """forged N, w and a_p, routed through lball_exact.evaluate itself (guard and Fricke residual): each must be refused"""
    rows = the_curves(); lines = ['census_lvalue_exact.py --negative: forged inputs through lball_exact.evaluate, each must be refused', '']
    PI = LX.pi_interval(); nref = ntot = 0
    for i in (0,):
        N, pair, name, a, LD = rows[i]
        bad = sorted(p for p, T in LD['rows'].items() if T['f'])
        lN = bad[-1]; nN = next(q for q in (3, 5, 7, 11, 13) if N % q)
        Mmax = max(LX.max_terms(Nf) for Nf in (N, N // lN, N * nN, 2*N))
        ap, st = L53.ap_table_fast(a, LD, Mmax); an0 = LX.an_from_ap(ap, LD, Mmax)
        good = [p for p in sorted(ap) if p >= 5 and p not in LD['rows']]
        pmid = next(p for p in good if p > math.isqrt(N) // 4)
        ev = LX.evaluate(an0, N, 1, PI, Mmax)
        lines.append('curve (%d,%d) %s, N = %d:' % (pair[0], pair[1], name, N))
        lines.append('  true inputs: L(y=1) %s, L(y=6/5) %s, Fricke residual %s: guard %s' % (LX.iv(ev['L_y1'], 10), LX.iv(ev['L_y2'], 10), LX.iv(ev['R'], 6), 'passes' if ev['guard'] else 'FAILS'))
        assert ev['guard'] and ev['L_y1'][0] > 0
        bp = bad[0]; bT = LD['rows'][bp]; q5 = 5 if 5 not in LD['rows'] else 7
        forged = [('N/%d (%d left out)' % (lN, lN), N // lN, 1, None), ('%dN (%d added)' % (nN, nN), N * nN, 1, None),
                  ('2N', 2*N, 1, None), ('w = -1', N, -1, None),
                  ('a_%d one more' % q5, N, 1, {q5: ap[q5] + 1}),
                  ('a_%d one more (p near sqrt(N)/4)' % pmid, N, 1, {pmid: ap[pmid] + 1}),
                  ('a_%d = %d (bad prime, %s)' % (bp, 1 if bT['ap'] != 1 else -1, bT['red']), N, 1, {bp: 1 if bT['ap'] != 1 else -1})]
        for nm, Nf, wf, ov in forged:
            an = an0 if ov is None else LX.an_from_ap({**ap, **ov}, LD, Mmax)
            e = LX.evaluate(an, Nf, wf, PI, Mmax)
            ok = e['guard']
            lines.append('  forged %-36s L(y=1) %s, L(y=6/5) %s, residual %s: %s' % (nm + ':', LX.iv(e['L_y1'], 10), LX.iv(e['L_y2'], 10),
                         LX.iv(e['R'], 6), 'guard passes (NOT REFUSED)' if ok else 'refused (%s)' % ', '.join(
                             x for x, y in (('enclosures disjoint', not e['overlap']), ('residual excludes 0', not e['R_has_0'])) if y)))
            ntot += 1; nref += not ok
    # the isogeny step: a non-partner curve of another conductor and a mutated partner (b + 1) must not be linked
    cur = {tuple(r[1]): r[3] for r in rows}
    A = cur[(112, 15)]; B = cur[(127, 97)]; C = cur[(260, 69)]
    ok = isogeny_chain((A[1], A[3]), (B[1], B[3])) is not None
    lines.append('isogeny chain (112,15) -> (127,97) (true partners): %s' % ('found' if ok else 'NOT FOUND')); assert ok
    for nm, dst in (('(112,15) -> (260,69) (not isogenous)', (C[1], C[3])), ('(112,15) -> (127,97) with b + 1', (B[1], B[3] + 1)),
                    ('(112,15) -> (127,97) with a + 2', (B[1] + 2, B[3])), ('(112,15) -> twist by -1 of (127,97)', (-B[1], B[3]))):
        r = isogeny_chain((A[1], A[3]), dst)
        lines.append('forged isogeny %-40s %s' % (nm + ':', 'refused (no chain)' if r is None else 'NOT REFUSED'))
        ntot += 1; nref += r is None
    lines.append('negative tests: %d forged inputs, %d refused' % (ntot, nref))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert nref == ntot
    if write: open(os.path.join(HERE, 'census_lvalue_exact_negative.txt'), 'w').write(out)

if __name__ == '__main__':
    args = sys.argv[1:]
    def opt(k, d=None): return args[args.index(k) + 1] if k in args else d
    if '--list' in args:
        for i, (N, pair, name, a, LD) in enumerate(the_curves()):
            print(i, pair, name, N, fN(LD), 'w_local', LD['w_local'], 'terms', LX.max_terms(N), 'model', a)
    elif '--merge' in args: merge([x for x in args[args.index('--merge') + 1:] if not x.startswith('--')], '--write' in args)
    elif '--torsion-negative' in args: torsion_negative('--write' in args)
    elif '--torsion' in args: torsion_all('--write' in args, '--pari' in args)
    elif '--negative' in args: negative('--write' in args)
    else:
        rows = the_curves(); sel = list(range(len(rows)))
        if opt('--only'): sel = [int(x) for x in opt('--only').split(',')]
        elif opt('--maxN'): sel = [i for i in sel if rows[i][0] <= float(opt('--maxN'))]
        elif opt('--pairs'):
            want = [tuple(int(x) for x in q.split(',')) for q in opt('--pairs').split(';')]
            sel = [i for i in sel if tuple(rows[i][1]) in want]; assert len(sel) == len(want)
        run(sel, '--pari' in args, opt('--jsonl'))
