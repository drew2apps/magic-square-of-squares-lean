"""The five curves over Q behind the seven census-wide ratios that are 'conditional on rank 1' (1000 < m <= 9000).

Written by an AI model (Claude); not read by a person.  Experiment, not counted.
    route E1 curves (route_pari.py): family(m, n)['+13'] twisted by d, for (3540,2581,-11740078), (4440,2231,-13342),
        (6480,4961,709342), (8249,3720,-108415202)  [(6121,959), (6671,2209) are the partners of the first two]
    Demjanenko-Manin D1 curve of (12007,1799) (dm_pari.py; rank_second.d1_curve gives the same model), for (6903,5104)
Each is returned as integers (a, b) of y^2 = x(x^2 + a x + b) with a rational point G of infinite order (a hint from
PARI's ellrank; on_curve is checked exactly by the caller)."""
import sys, os, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.abspath(os.path.join(HERE, '..', '..'))
HINTS = os.path.join(HERE, 'f2_curves_hints.json')
CURVES = [('E1 of (3540,2581), d = -11740078', 3540, 2581, -11740078),
          ('E1 of (4440,2231), d = -13342', 4440, 2231, -13342),
          ('E1 of (6480,4961), d = 709342', 6480, 4961, 709342),
          ('E1 of (8249,3720), d = -108415202', 8249, 3720, -108415202),
          ('D1 of (12007,1799)', 12007, 1799, None)]

def weier(m, n, d):
    """(a2, a4, a6) with a1 = a3 = 0."""
    if d is None:
        sys.path.insert(0, os.path.join(EXP, '..', 'scripts', 'g1_full'))
        S = m * m + n * n; D2 = (m * m - n * n) ** 2
        c3 = 4 * S * S * m * m * n * n; c2 = -16 * S * m ** 3 * n ** 3; c1 = 4 * S * S * D2; c0 = -16 * S * m * n * D2
        return c1, c0 * c2, c0 * c0 * c3
    sys.path.insert(0, os.path.join(EXP, 'census-wide'))
    from sd_pari import family, weierstrass, pari
    Q, r = family(m, n)['+13']; ain, A = weierstrass(d * Q, r)
    Em = pari.ellminimalmodel(pari.ellinit([pari(x.numerator) / x.denominator for x in ain]))
    a = [int(Em[i]) for i in range(5)]
    assert a[0] == 0 and a[2] == 0
    return a[1], a[3], a[4]

def two_torsion_form(a2, a4, a6):
    import cypari2; pari = cypari2.Pari()
    rts = [int(r) for r in pari.nfroots(None, pari('x^3 + (%d)*x^2 + (%d)*x + (%d)' % (a2, a4, a6)))]
    assert len(rts) == 1; e = rts[0]; assert e ** 3 + a2 * e * e + a4 * e + a6 == 0
    return e, 3 * e + a2, 3 * e * e + 2 * a2 * e + a4

def curves(find=False):
    hints = json.load(open(HINTS)) if os.path.exists(HINTS) else {}
    out = []
    for name, m, n, d in CURVES:
        a2, a4, a6 = weier(m, n, d)
        e, a, b = two_torsion_form(a2, a4, a6)
        if name not in hints:
            assert find, 'no hint for ' + name
            import cypari2; pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9)
            R = pari.ellrank(pari.ellinit([0, a, 0, b, 0]), 4)
            P = R[3][0]
            hints[name] = [str(P[0]), str(P[1])]
        G = (Fr(hints[name][0]), Fr(hints[name][1]))
        assert G[1] ** 2 == G[0] * (G[0] ** 2 + a * G[0] + b), 'G is not on the curve'
        out.append((name, a, b, G))
    if find: json.dump(hints, open(HINTS, 'w'), indent=1)
    return out

if __name__ == '__main__':
    for name, a, b, G in curves(find='--find' in sys.argv):
        print(name, 'a =', a, 'b =', b, 'G =', G[0])
