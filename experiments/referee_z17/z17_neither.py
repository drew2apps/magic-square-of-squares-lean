"""Referee z17 (written by an AI model; not read by a person): the three "neither" closures of roadmap I's I8.
(a) the curve DM ran on (census2000/dm_pari.py's Sig+ of (18313,120), Sig- of (20056,19825), E of (39881,231)) against
    I8's "E+"/"E-" (scripts/g1_full/sigma_exact.model(m, n, +-2), as neither/epm.py builds them): isomorphic over Q?
    (c4, c6 scaled by u^4, u^6 with u^2 rational and a square);
(b) PARI (this container's version) on each closing curve: ellrank effort 1 [lo, hi], ellrootno, elltors, the isogeny
    class (ellisomat) and the number of rational 2-torsion points of each member (exact: integer roots of the cubic);
(c) dm_pari.dm_curve rerun on the three closing curves and compared, field by field (without 'sec'), with
    roadmap-i/i8/neither/dm_*.jsonl;
(d) exact rank bounds: for every member of the class with full rational 2-torsion, H7's step (roadmap-h/h1
    h1_rank0_exact.bound_full: complete 2-descent and the Cassels-Tate pairing of rank_second.py); for the others,
    rank_indep.rank_bound (descent via 2-isogeny and second descent), as neither/ri.py;
(e) PARI's point on each closing curve checked exactly: on the curve and of infinite order (k G != O, k <= 12).
usage: python3 -I z17_neither.py PART [--write]    PART in a b c d e  (z17_neither_<PART>.txt)"""
import sys, os, io, json, contextlib, time
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
os.environ['MSS_REPO'] = REPO
argv = sys.argv; sys.argv = [argv[0]]
for d in (('scripts', 'g1_full'), ('experiments', 'census2000'), ('experiments', 'roadmap-h', 'h1')):
    sys.path.insert(0, os.path.join(REPO, *d))
with contextlib.redirect_stdout(io.StringIO()):
    import sigma_exact as SX
    import dm_pari as DM
    import rank_indep as RI
    import h1_rank0_exact as R
sys.argv = argv
pari = DM.pari; pari.default('nbthreads', 1)
PART = argv[1]
OUT = []
def say(s=''): OUT.append(s); print(s, flush=True)
say('z17 (part %s): the "neither" closures of I8. Written by an AI model (reviewer z17); not read by a person. PARI %s.' % (PART, list(pari.version())[:3]))
CL = [((18313, 120), 'Sig+', 2), ((20056, 19825), 'Sig-', -2), ((39881, 231), 'E', None)]
def dm_entry(m, n, nm):
    return [c for c in DM.all_curves(m, n) if c[4] == nm][0]
def c46(ai):
    E = pari.ellinit(ai); return Fr(str(E[9])), Fr(str(E[10]))
def iso(ai, bi):
    (c4, c6), (d4, d6) = c46(ai), c46(bi)
    if c4 == 0 or c6 == 0: return None
    u2 = (d6 / c6) / (d4 / c4)
    from math import isqrt
    sq = u2 > 0 and isqrt(u2.numerator) ** 2 == u2.numerator and isqrt(u2.denominator) ** 2 == u2.denominator
    return u2 * u2 == d4 / c4 and u2 ** 3 == d6 / c6 and sq      # u2 not a square: a quadratic twist, not isomorphic
def roots2(ai):
    a1, a2, a3, a4, a6 = [int(v) for v in ai]
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    import sympy
    X = sympy.Symbol('X')
    return [b2, 8 * b4, 16 * b6], sorted(int(r) for r in sympy.Poly(X ** 3 + b2 * X ** 2 + 8 * b4 * X + 16 * b6, X).ground_roots() if r.is_Integer)

if PART == 'a':
    for (m, n), nm, s in CL:
        ai = dm_entry(m, n, nm)[0]
        if s is None: say('(%d,%d) %s: a1..a6 = %s (y^2 = x(x + 2Sm^2)(x + 2Sn^2): full rational 2-torsion by its form)' % (m, n, nm, ai)); continue
        M = SX.model(m, n, s); bi = [0, M['a2'], 0, M['a4'], M['a6']]
        say('(%d,%d): dm_pari %s %s; sigma_exact.model(s = %d) %s; isomorphic over Q: %s; j equal: %s'
            % (m, n, nm, ai, s, bi, iso(ai, bi), pari.ellinit(ai)[12] == pari.ellinit(bi)[12]))
        o = SX.model(m, n, -s); ci = [0, o['a2'], 0, o['a4'], o['a6']]
        say('  and with sigma_exact.model(s = %d): isomorphic %s' % (-s, iso(ai, ci)))
elif PART == 'b':
    for (m, n), nm, s in CL:
        ai = dm_entry(m, n, nm)[0]; E = pari.ellinit(ai)
        r = pari.ellrank(E, 1)
        say('(%d,%d) %s: ellrank effort 1 [%d,%d] s=%d points %d; rootno %d; torsion %s; conductor digits %d'
            % (m, n, nm, int(r[0]), int(r[1]), int(r[2]), len(r[3]), int(pari.ellrootno(E)), [int(v) for v in pari.elltors(E)[1]], len(str(abs(int(pari.ellglobalred(E)[0]))))))
        lo, hi = DM.rank_bounds(E)
        say('  dm_pari.rank_bounds (efforts 0, 1 over the isogeny class): [%d,%d]' % (lo, hi))
        cls = [[int(v) for v in pari.ellminimalmodel(pari.ellinit(c))[:5]] for c in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]]
        say('  isogeny class (PARI): %d curves; rational 2-torsion points per curve %s' % (len(cls), [len(roots2(c)[1]) for c in cls]))
elif PART == 'c':
    ref = {}
    for f in ('dm_18313_120.jsonl', 'dm_20056_19825.jsonl'):
        for l in open(os.path.join(REPO, 'experiments', 'roadmap-i', 'i8', 'neither', f)):
            o = json.loads(l); ref[(tuple(o['pair']), o['curve'].replace('-partner', ''))] = o
    for (m, n), nm, s in CL:
        ai, rel, C, tt, _ = dm_entry(m, n, nm)
        o = DM.dm_curve(m, n, ai, rel, C, tt, nm); o = json.loads(json.dumps(o, default=str))
        r = dict(ref[((m, n), nm)]); r.pop('sec'); r.pop('of'); r['curve'] = nm
        diff = sorted(k for k in set(o) | set(r) if o.get(k) != r.get(k))
        say('(%d,%d) %s: %s' % (m, n, nm, json.dumps(o, sort_keys=True)))
        say('  fields that differ from the committed record: %s' % (diff or 'none'))
        for k in diff:
            if isinstance(o.get(k), float): say('    %s: here %r, committed %r' % (k, o.get(k), r.get(k)))
elif PART == 'd':
    for (m, n), nm, s in CL:
        ai = dm_entry(m, n, nm)[0]; E = pari.ellinit(ai)
        cls = [[int(v) for v in pari.ellminimalmodel(pari.ellinit(c))[:5]] for c in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]]
        best = None
        for c in cls:
            sc, rs = roots2(c); t0 = time.time()
            if len(rs) == 3:
                e = rs[0]; a, b = 3 * e + sc[0], 3 * e * e + 2 * sc[0] * e + sc[1]
                es = R.roots_of(a, b); rb, info = R.bound_full(es)
                say('(%d,%d) %s: class member %s, full 2-torsion: H7 step (bound_full): %s, rank <= %d' % (m, n, nm, c, info, rb)); sys.stderr.write('%.0f s\n' % (time.time() - t0))
            elif len(rs) == 1:
                full = [0, sc[0], 0, sc[1], sc[2]]
                e, a, b = RI.shift_to_two_torsion(sc[0], sc[1], sc[2]); ap, bp = -2 * a, a * a - 4 * b
                pE = [(P[0] - e, P[1]) for P in RI.find_points(full)]
                buf = []
                rb, rec = RI.rank_bound(a, b, pE, RI.find_points([0, ap, 0, bp, 0]), say=buf.append)
                say('(%d,%d) %s: class member %s, one 2-torsion point: rank_indep.rank_bound rank <= %d' % (m, n, nm, c, rb)); sys.stderr.write('%.0f s\n' % (time.time() - t0))
                for x in buf: say('    ' + str(x).strip())
            else: continue
            best = rb if best is None else min(best, rb)
        say('  exact bound over the class: rank <= %s' % best)
elif PART == 'e':
    # PARI's point (ellrank effort 1, a finder) checked exactly: on the curve, and k G != O for k <= 12 (Mazur), by the
    # group law in exact rationals on y^2 = x^3 + a2 x^2 + a4 x + a6
    def add(A, P, Q):
        a2 = A[1]
        if P is None: return Q
        if Q is None: return P
        if P[0] == Q[0] and P[1] == -Q[1]: return None
        lam = (3 * P[0] ** 2 + 2 * a2 * P[0] + A[3]) / (2 * P[1]) if P == Q else (Q[1] - P[1]) / (Q[0] - P[0])
        x3 = lam * lam - a2 - P[0] - Q[0]
        return (x3, lam * (P[0] - x3) - P[1])
    for (m, n), nm, s in CL:
        ai = dm_entry(m, n, nm)[0]; E = pari.ellinit(ai); r = pari.ellrank(E, 1)
        G = (Fr(str(r[3][0][0])), Fr(str(r[3][0][1])))
        on = G[1] ** 2 == G[0] ** 3 + ai[1] * G[0] ** 2 + ai[3] * G[0] + ai[4]
        R = None; inf = True
        for k in range(1, 13):
            R = add(ai, R, G)
            if R is None: inf = False; break
        say('(%d,%d) %s: PARI point on the curve (exact): %s; k G != O for k <= 12 (exact group law): %s, so rank >= 1' % (m, n, nm, on, inf))
s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'z17_neither_%s.txt' % PART)
if '--write' in argv: open(f, 'w').write(s)
elif os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
