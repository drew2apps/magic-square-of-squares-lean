#!/usr/bin/env python3
"""F4: Demjanenko-Manin on the quotient C: Y^2 = g4(s) K(s) (s = t^2) of the G1 system, through the relation s -> 1/s.
Written by an AI model (Claude); not read by a person.  Experiment, not counted.

Why there is a relation.  With s = t^2, g4 = 2S(1 + s) = 4 F4 and K = S^2 + (2S^2 - 16 m^2 n^2) s + S^2 s^2 = 4 L U, so a
solution t gives a point of C at s = t^2.  The system is invariant under t -> 1/t (F0 <-> F8, the others times 1/t^2),
so it also gives a point at s = 1/t^2.  (Equivalently: w*g4*K is the same curve through (s, Y) -> (1/s, Y/s^2), so the
two rank-1 quotients g4*K and w*g4*K of census-wide CENSUS9000.md are one curve, and a solution gives two points on it.)
They are the two images of a point of the genus-2 curve D: y^2 = g4(u^2) K(u^2) under its two elliptic quotient maps
(u, y) -> (u^2, y) and (u, y) -> (1/u^2, y/u^3); J(D) is isogenous to C x C, so the two maps are independent and
Demjanenko-Manin applies when rank C(Q) = 1.
Weierstrass model of the cubic c3 s^3 + c2 s^2 + c1 s + c0: x = c3 s, y = c3 Y, y^2 = x^3 + c2 x^2 + c1 c3 x + c0 c3^2,
so the relation is x x' = c3^2 (kappa = c3), and t = +-sqrt(x / c3).

The Demjanenko-Manin step itself is census2000/dm_pari.py's dm_curve, unchanged (PARI ranks and generator, saturation
at primes <= 1000, so 'conditional on PARI generators'; height bounds in floats with a margin).

    python3 dm_g4k.py                 run the 12 ratios, check the output equals dm_g4k.jsonl
    python3 dm_g4k.py --write         write dm_g4k.jsonl
    python3 dm_g4k.py --negative      controls (with --write: dm_g4k_negative.txt)
    python3 dm_g4k.py m,n [m,n ...]   one JSON record per ratio
"""
import sys, os, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'census2000'))
import dm_pari as DM
from fractions import Fraction as Fr
x = DM.x

def g4k_curve(m, n):
    S = m * m + n * n
    g4 = [2 * S, 2 * S]; K = [S * S, 2 * S * S - 16 * m * m * n * n, S * S]
    P = [0] * 4
    for i, u in enumerate(g4):
        for j, v in enumerate(K): P[i + j] += u * v
    c0, c1, c2, c3 = P
    ainv = [0, c2, 0, c1 * c3, c0 * c3 * c3]
    def tC(X):
        w = Fr(X) / c3
        if w < 0 or not DM.issq(w): return []
        r = DM.sqrt_fr(w); return [r, -r] if r else [r]
    # relation x x' = c3^2:  -c3^2 + x x' = 0
    return ainv, (-(c3 * c3) + 0 * x, x), 2 * math.log(abs(c3)), tC, 'g4*K'

RATIOS = ['7735,4968', '12703,2767', '6557,324', '6881,6233',
          '3540,2581', '6121,959', '4440,2231', '6671,2209', '6480,4961', '6903,5104', '12007,1799', '8249,3720']

def run(arg, rel_scale=1, fake=None, rel_forge=None):
    m, n = (int(v) for v in arg.split(','))
    ainv, rel, C, tt, nm = g4k_curve(m, n)
    if rel_scale != 1: rel = (rel[0] * rel_scale, rel[1])          # forged relation x x' = rel_scale c3^2
    if rel_forge is not None: rel = rel_forge
    keep = DM.classify_t
    if fake is not None: DM.classify_t = lambda mm, nn, t: 'SOLUTION' if t == fake else keep(mm, nn, t)
    try: r = DM.dm_curve(m, n, ainv, rel, C, tt, nm)
    finally: DM.classify_t = keep
    r['rootno'] = int(DM.pari.ellrootno(DM.pari.ellinit(ainv)))
    for k in ('hG', 'B'):                                          # floats only guide amax; print them rounded
        if k in r: r[k] = round(r[k], 3)
    r.pop('sec', None)
    return r

def negative():
    out = []
    r = run('7735,4968')
    out.append('positive control: (7735,4968) unchanged: status %s, t found %s (must be closed, with t = 1 and -1 degenerate)' % (r['status'], sorted(r['t'])))
    ok1 = r['status'] == 'closed' and sorted(r['t']) == ['-1', '1']
    r = run('7735,4968', rel_forge=(-x, 1 + 0 * x))
    out.append("relation forged to x' = x, which every point satisfies: status %s (must be identity: no finite list exists)" % r['status'])
    ok2 = r['status'] == 'identity'
    r = run('7735,4968', fake=Fr(1))
    out.append('classifier forged to call t = 1 a solution: status %s (must be ALERT)' % r['status'])
    ok3 = r['status'] == 'ALERT'
    r = run('3540,2581')
    out.append('a ratio whose g4*K has PARI rank [2,2]: status %s (must be "rank not 1")' % r['status'])
    ok4 = r['status'] == 'rank not 1'
    out.append('negative checks: %s' % ('all as required (4 of 4)' if ok1 and ok2 and ok3 and ok4 else 'FAILED'))
    hb = DM.height_bounds; DM.height_bounds = lambda E: tuple(3 * v for v in hb(E))
    try:
        for a in ('7735,4968', '6557,324', '6881,6233'):
            r = run(a)
            out.append('robustness, not a control: %s with the height bound tripled: amax %d, %d candidates, t %s, status %s'
                       % (a, r['amax'], r['ncand'], sorted(r['t']), r['status']))
    finally: DM.height_bounds = hb
    return out

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--negative' in sys.argv:
        lines = negative(); txt = '\n'.join(lines) + '\n'; print(txt, end='')
        if '--write' in sys.argv: open(os.path.join(HERE, 'dm_g4k_negative.txt'), 'w').write(txt)
        sys.exit(0)
    txt = ''.join(json.dumps(run(a), default=str) + '\n' for a in (args or RATIOS))
    print(txt, end='')
    if '--write' in sys.argv: open(os.path.join(HERE, 'dm_g4k.jsonl'), 'w').write(txt)
    elif not args:
        assert txt == open(os.path.join(HERE, 'dm_g4k.jsonl')).read(), 'output differs from dm_g4k.jsonl'
        print('# matches dm_g4k.jsonl')
