#!/usr/bin/env python3
"""covers.py (conductor-t13): the 2-coverings of the corner curve E = L1L3M1M2 of (627,364) over K = Q(i), and of
E3 = E/<T3>, as binary quartics y^2 = G(t1, t2) over Z[i], one per non-trivial coset of Sel^2 modulo the torsion
image (and one per torsion class, as a control).  EXPERIMENT; written by an AI model (build agent t13), not read by a
person.

    python3 covers.py [--write]      -> covers.txt, covers.json (the quartics and the maps back to E)

Construction (exact, Sage over K).  For a class a = (d1, d2) of Sel^2 (d3 = the class of d1 d2, c^2 = d1 d2 d3), the
2-covering is x = e1 + d1 u^2/z^2 = e2 + d2 v^2/z^2 = e3 + d3 w^2/z^2.  A K-point P0 of the conic
d1 u^2 - d2 v^2 + (e1 - e2) z^2 = 0 (kctp_recheck.conic_point: lattice search, checked exactly) parametrises the
conic, (u : v : z) = X(t1, t2) quadratic, and the covering becomes W^2 = G(t) := d3 (d1 u(t)^2 + (e1 - e3) z(t)^2),
W = d3 w.  G is made integral and primitive up to squares, minimised at every prime dividing 2 (e_i - e_j) or the
d_k (moves G -> pi^-2 G(pi t1 + r t2, t2) at the multiple roots r of G mod pi, and G -> pi^-2 G; a move is kept when it
lowers the valuation of the invariants I, J), and reduced (translations t1 -> t1 + c t2, c in Z[i], and t -> (t2, -t1),
greedy in the size sum |g_k|^2 / binom(4, k)).  The total change of variable M (2x2 over K) is stored, so that a point
t of the final quartic gives x = e1 + d1 (u/z)^2 with (u : v : z) = X(M t); the search (search.py) re-checks every
point exactly on E (f(x) a square in K).
Nothing here certifies anything about E(K): it only builds search models.
"""
from sage.all__sagemath_schemes import *
import sys, os, json, random, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'scripts', 'second_descent'))
import kctp_recheck as R          # existing, AI-written; imported, not modified

K = R.K; I = R.I; OK = K.ring_of_integers()
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

S2 = PolynomialRing(K, ['t1', 't2']); t1, t2 = S2.gens()

def invariants(G):
    """I, J of the binary quartic G = a t1^4 + b t1^3 t2 + c t1^2 t2^2 + d t1 t2^3 + e t2^4."""
    a, b, c, d, e = [G.monomial_coefficient(t1 ** (4 - k) * t2 ** k) for k in range(5)]
    Iv = 12 * a * e - 3 * b * d + c ** 2
    Jv = 72 * a * c * e + 9 * b * c * d - 27 * a * d ** 2 - 27 * e * b ** 2 - 2 * c ** 3
    return Iv, Jv

def coeffs(G): return [G.monomial_coefficient(t1 ** (4 - k) * t2 ** k) for k in range(5)]
def from_coeffs(cs): return sum(c * t1 ** (4 - k) * t2 ** k for k, c in enumerate(cs))
def subst(G, M):
    """G(M (t1, t2)^T) with M = [[a, b], [c, d]]: t1 -> a t1 + b t2, t2 -> c t1 + d t2."""
    (a, b), (c, d) = M
    return G(a * t1 + b * t2, c * t1 + d * t2)
def mmul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]

def cval(G, P):
    cs = [c for c in coeffs(G) if c != 0]
    return min(c.valuation(P) for c in cs)

def make_integral(G):
    """G * s^2 with s in K so that G is integral, then divide by the square part of the content (s recorded)."""
    cs = [c for c in coeffs(G) if c != 0]
    den = lcm([c.denominator() for c in cs]); s = K(den)
    G = G * s ** 2
    g = K.ideal(coeffs(G))
    fac = g.factor()
    q = prod([P.gens_reduced()[0] ** (e // 2) for P, e in fac]) if len(fac) else K(1)
    G = G / q ** 2; s = s / q
    assert all(c in OK for c in coeffs(G))
    return G, s

def multiple_roots(G, P):
    """candidate centres of a move at P: with c = min(v_P(coefficients)) in {0, 1}, the roots of G / pi^c mod P in
    P^1 of the residue field of multiplicity >= 2 - c, as lifts r (None for infinity)."""
    F = K.residue_field(P); pi = P.gens_reduced()[0]
    c = cval(G, P)
    if c >= 2: return []
    cs = [x / pi ** c for x in coeffs(G)]
    Rf = PolynomialRing(F, 'z'); z = Rf.gen()
    g = sum(F(x) * z ** (4 - k) for k, x in enumerate(cs))
    need = 2 - c
    out = [F.lift(r) for r, mlt in g.roots() if mlt >= need]
    if 4 - g.degree() >= need:      # root at infinity of multiplicity 4 - deg
        out.append(None)
    return out

def move(G, P, r):
    pi = P.gens_reduced()[0]
    M = [[pi, r], [0, 1]] if r is not None else [[1, 0], [0, pi]]
    return subst(G, M), M

def minimise_at(G, P, Mtot, s, depth=4):
    """greedy minimisation at P (see the module docstring); returns G, Mtot, s and the number of reductions."""
    pi = P.gens_reduced()[0]
    red = 0
    while True:
        c = cval(G, P)
        if c >= 2:
            k = c // 2; G = G / pi ** (2 * k); s = s / pi ** k; red += k; continue
        Iv, Jv = invariants(G)
        vI = Iv.valuation(P) if Iv != 0 else 99; vJ = Jv.valuation(P) if Jv != 0 else 99
        if vI < 4 or vJ < 6: break
        # breadth-first over move sequences: keep models with content >= 2 after the move (divide by pi^2)
        frontier = [(G, Mtot, s)]; found = None
        for dep in range(depth):
            new = []
            for H, MH, sH in frontier:
                rts = multiple_roots(H, P)
                for r in rts:
                    H2, M = move(H, P, r)
                    if cval(H2, P) >= 2:
                        H3 = H2 / pi ** 2; s3 = sH / pi; M3 = mmul(MH, M)
                        if cval(H3, P) >= 2:
                            found = (H3 / pi ** 2, M3, s3 / pi); break
                        new.append((H3, M3, s3))
                if found: break
            if found: break
            frontier = new[:16]
            if not frontier: break
        if not found: break
        G, Mtot, s = found; red += 1
    return G, Mtot, s, red

def size(G):
    cs = coeffs(G)
    return sum(float(abs(CC(c))) ** 2 / binomial(4, k) for k, c in enumerate(cs))

def covariant(G):
    """the point (z, t) of hyperbolic 3-space minimising sum_j log((|z - a_j|^2 + t^2) / t) over the roots a_j of
    G(z, 1) (a reduction covariant of the binary form; floating point, used only to choose the change of variable)."""
    import numpy as np
    from scipy.optimize import minimize
    cs = coeffs(G)
    Rz = PolynomialRing(ComplexField(400), 'z'); z = Rz.gen()
    g = sum(ComplexField(400)(CC(c)) * z ** (4 - k) for k, c in enumerate(cs))
    rts = [complex(r) for r in g.roots(multiplicities=False)]
    assert len(rts) == 4
    def f(v):
        x, y, lt = v; t = np.exp(lt)
        return sum(np.log(((x - r.real) ** 2 + (y - r.imag) ** 2 + t * t) / t) for r in rts)
    c0 = sum(rts) / 4; sc = max(1e-30, sum(abs(r - c0) for r in rts) / 4)
    best = None
    for st in ([c0.real, c0.imag, np.log(sc)], [rts[0].real, rts[0].imag, np.log(sc)]):
        res = minimize(f, st, method='Nelder-Mead', options={'xatol': 1e-10, 'fatol': 1e-12, 'maxiter': 20000})
        if best is None or res.fun < best.fun: best = res
    x, y, lt = best.x
    return complex(x, y), float(np.exp(lt))

def reduce_quartic(G, Mtot):
    """move the covariant point into the fundamental domain of PSL2(Z[i]) (|re z|, |im z| <= 1/2, |z|^2 + t^2 >= 1)
    by translations t1 -> t1 + c t2 and the inversion (t1, t2) -> (-t2, t1); exact on G, floating point only in the
    choice of the moves."""
    for it in range(500):
        if coeffs(G)[0] == 0:
            G = subst(G, [[1, 1], [0, 1]]); Mtot = mmul(Mtot, [[1, 1], [0, 1]]); continue
        zc, tc = covariant(G)
        if max(abs(zc.real), abs(zc.imag)) > 0.5 + 1e-9:
            c = K(round(zc.real)) + K(round(zc.imag)) * I
            G = subst(G, [[1, c], [0, 1]]); Mtot = mmul(Mtot, [[1, c], [0, 1]]); continue
        if abs(zc) ** 2 + tc ** 2 < 1 - 1e-9:
            M = [[0, -1], [1, 0]]
            G = subst(G, M); Mtot = mmul(Mtot, M); continue
        break
    return G, Mtot

def parametrise(co, P0):
    """X(t) = Q(D) P0 - 2 B(P0, D) D, D = t1 e_a + t2 e_b, for the diagonal conic sum co_k X_k^2."""
    for a, b in ((0, 1), (0, 2), (1, 2)):
        Mx = matrix(K, [P0, [1 if j == a else 0 for j in range(3)], [1 if j == b else 0 for j in range(3)]])
        if Mx.det() != 0: break
    D = [t1 * (1 if j == a else 0) + t2 * (1 if j == b else 0) for j in range(3)]
    QD = sum(co[j] * D[j] ** 2 for j in range(3)); BPD = sum(co[j] * P0[j] * D[j] for j in range(3))
    X = [QD * P0[j] - 2 * BPD * D[j] for j in range(3)]
    assert sum(co[j] * X[j] ** 2 for j in range(3)) == 0
    return X

def build(es, d1, d2, d3, rng, label):
    e1, e2, e3 = es
    c2 = d1 * d2 * d3; assert c2.is_square()
    co = [d1, -d2, e1 - e2]
    P0 = R.conic_point(co, rng)
    X = parametrise(co, P0)
    G0 = d3 * (d1 * X[0] ** 2 + (e1 - e3) * X[2] ** 2)
    G, s = make_integral(G0)
    Mtot = [[K(1), K(0)], [K(0), K(1)]]
    D = 2 * (e1 - e2) * (e1 - e3) * (e2 - e3) * d1 * d2 * d3
    Iv0, Jv0 = invariants(G)
    gIJ = K.ideal([Iv0, Jv0])                      # primes where the model can be non-minimal: v(I), v(J) > 0
    primes = sorted(set(P for P, _ in K.ideal(D).factor()) | set(P for P, _ in gIJ.factor()), key=lambda P: P.norm())
    reds = {}
    for P in primes:
        G, Mtot, s, red = minimise_at(G, P, Mtot, s)
        if red: reds[str(P.gens_reduced()[0])] = red
    G, Mtot = reduce_quartic(G, Mtot)
    # exact consistency: G(t) = s^2 G0(M t)
    assert subst(G0, Mtot) * s ** 2 == G
    assert all(c in OK for c in coeffs(G))
    Iv, Jv = invariants(G)
    return dict(label=label, d=[d1, d2, d3], P0=P0, X=X, G=G, M=Mtot, s=s, reds=reds, I=Iv, J=Jv)

def kstr(x):
    x = K(x); return [str(QQ(x[0])), str(QQ(x[1]))]

def curve_data(es, rng, name):
    sel = R.Selmer(es, rng)
    E = R.curve_from_e(es)
    tor = []
    T = [E([ee, 0]) for ee in es]
    for t in T:
        for H in t.division_points(2):
            tor.append(H)
    timg = []
    for x in [es[0], es[1]] + [H.xy()[0] for H in tor]:
        dd1, dd2 = R.kummer(x, es); timg.append(sel.vec_of(dd1) + sel.vec_of(dd2))
    rt = R.f2_rank(timg)
    say('%s: e = %s' % (name, [str(x) for x in es]))
    say('  dim Sel^2 = %d, torsion image dim %d (T1, T2 and %d points of order 4)' % (sel.dim, rt, len(tor)))
    # a complement of the torsion image in Sel^2
    comp = []
    for b in sel.basis:
        if R.f2_rank(timg + comp + [b]) > R.f2_rank(timg + comp): comp.append(b)
    assert len(comp) == sel.dim - rt
    tbasis = []
    for v in timg:
        if R.f2_rank(tbasis + [v]) > R.f2_rank(tbasis): tbasis.append(v)
    return sel, E, comp, tbasis, rt

def main():
    rng = random.Random(20261127)
    t0 = time.time()
    say('covers.py (conductor-t13): 2-coverings of E = L1L3M1M2 of (627,364) and of E3 = E/<T3> over K = Q(i).')
    say('Written by an AI model (build agent t13); not read by a person.  Experiment; builds search models only.')
    Q, fs = R.quartic(627, 364, 'L1L3M1M2'); es = R.full2_model(Q, fs)
    E = R.curve_from_e(es)
    E3 = E.isogeny(E([es[2], 0])).codomain()
    import importlib
    KI = importlib.import_module('kctp_isog')
    es3 = KI.short_full2(E3)
    out = {'curves': {}}
    for name, ee in (('E', es), ('E3', es3)):
        sel, Em, comp, tbasis, rt = curve_data(ee, rng, name)
        cl = []
        # non-trivial cosets modulo the torsion image: comp[0], comp[1], comp[0] + comp[1]
        vecs = [('coset A', comp[0]), ('coset B', comp[1]), ('coset A+B', [(x + y) % 2 for x, y in zip(comp[0], comp[1])])]
        # torsion classes as controls (their coverings carry torsion points)
        vecs += [('torsion t1', tbasis[0]), ('torsion t2', tbasis[1])]
        for lab, v in vecs:
            dd1, dd2, dd3 = sel.pair(v)
            assert sel.contains(dd1, dd2)
            c = build(ee, dd1, dd2, dd3, rng, lab)
            nG = [float(abs(CC(x))) for x in coeffs(c['G'])]
            say('  %s: d = %s' % (lab, [str(x) for x in c['d']]))
            say('    conic point %s; minimisation steps %s' % ([str(x) for x in c['P0']], c['reds']))
            zc, tc = covariant(c['G'])
            say('    reduced G: max |g_k| = %.3e, covariant (z, t) = (%.3f%+.3fi, %.3e); N(4I^3 - J^2) = %s' % (max(nG), zc.real, zc.imag, tc, factor(ZZ((4 * c['I'] ** 3 - c['J'] ** 2).norm()))))
            cl.append({'label': lab, 'd': [kstr(x) for x in c['d']], 'G': [kstr(x) for x in coeffs(c['G'])],
                       'M': [[kstr(x) for x in row] for row in c['M']], 's': kstr(c['s']),
                       'covariant': [float('%.6f' % float(zc.real)), float('%.6f' % float(zc.imag)), float('%.6e' % float(tc))],
                       'X': [[kstr(Xk.monomial_coefficient(m)) for m in (t1 ** 2, t1 * t2, t2 ** 2)] for Xk in c['X']]})
        out['curves'][name] = {'e': [kstr(x) for x in ee], 'classes': cl}
    print('time %.0f s' % (time.time() - t0), file=sys.stderr)
    if WRITE:
        json.dump(out, open(os.path.join(HERE, 'covers.json'), 'w'), indent=1, sort_keys=True)
        open(os.path.join(HERE, 'covers.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
