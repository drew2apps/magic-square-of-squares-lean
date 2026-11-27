#!/usr/bin/env python3
"""A machine test (not a proof) of the one written step of rank_fisher_theta4_gen.py that is about g0 itself: that
div g0 is even on C4 (div g0 = 2B).  Weil reciprocity on the reduction of C4 modulo a prime p of good reduction:
for h = l1 / l2 (l1, l2 linear forms involving W), g0(div h) = h(div g0); if div g0 = 2B with B defined over F_p,
the right side is h(B)^2, a square in F_p^*.  So for every pair of hyperplane sections Z1 = C4 n {l1 = 0},
Z2 = C4 n {l2 = 0} (8 points each over the algebraic closure), N(Z1) / N(Z2) must be a square, where
N(Z) = prod_{P in Z} g0(P) is computed exactly in F_p[e] / (h(e)) from a lex Groebner basis in shape position
(sympy, modulus p).  Control: g0 c / e (an odd divisor) must fail for some pairs.  (reviews/review_w4.md, FIX 4,
did this by hand for (72,65); here for the general module's curves.)

    python3 rank_fisher_theta4_gen_weil.py [--write]        (rank_fisher_theta4_gen_weil.txt)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.  A test modulo two primes is
evidence for the written argument, not a proof over Q.
"""
import os, sys, random
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import sympy as sp
import rank_fisher_theta4 as T
import rank_fisher_theta4_gen as TG
from rank_fisher_theta4 import pari, c, d, e, g, W

WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

def norm_mod(Gexpr, phis, h, p):
    """prod over the roots of h of G(phi_c(t), phi_d(t), t): the determinant of multiplication by G in F_p[t]/(h)."""
    tt = sp.Symbol('tt')
    H = sp.Poly(h, tt, modulus=p)
    n = H.degree()
    Gt = sp.Poly(Gexpr, tt, modulus=p).rem(H)
    rows = []
    cur = sp.Poly(1, tt, modulus=p)
    for i in range(n):
        v = (Gt * cur).rem(H)
        co = [int(v.coeff_monomial(tt ** j)) % p for j in range(n)]
        rows.append(co)
        cur = (cur * sp.Poly(tt, tt, modulus=p)).rem(H)
    M = pari.matrix(n, n, [x for r in rows for x in r])
    return int(pari.lift(pari.matdet(pari.Mod(1, p) * M)))

def section_norms(D, Nz, Xi, p, lin, num, den):
    """N(num) and N(den) over C4 n {W = lin(c, d, e, g)} with g = 1 (points with g = 0 are avoided: the ideal must have
    degree 8 in shape position)."""
    tt = sp.Symbol('tt')
    L = lin
    Q3 = sp.expand(Xi * L ** 2 - D['f0'].as_expr())
    eqs = [sp.expand(D['QA'].subs(g, 1)), sp.expand(D['QB'].subs(g, 1)), sp.expand(Q3.subs(g, 1))]
    GB = sp.groebner(eqs, c, d, e, modulus=p, order='lex')
    polys = list(GB.exprs)
    if len(polys) != 3: return None
    he = [q for q in polys if q.free_symbols <= {e}]
    if len(he) != 1: return None
    h = sp.Poly(he[0], e, modulus=p)
    if h.degree() != 8 or not h.is_sqf: return None
    phis = {}
    for q in polys:
        P_ = sp.Poly(q, c, d, e, modulus=p)
        if P_.degree(c) == 1 and P_.degree(d) == 0 and P_.coeff_monomial(c) == 1 and len(q.free_symbols) <= 2:
            phis[c] = sp.expand(c - q)
        if P_.degree(c) == 0 and P_.degree(d) == 1 and P_.coeff_monomial(d) == 1 and len(q.free_symbols) <= 2:
            phis[d] = sp.expand(d - q)
    if c not in phis or d not in phis: return None
    sub = {c: phis[c], d: phis[d], g: 1, W: L.subs({c: phis[c], d: phis[d], g: 1})}
    out = []
    for F_ in (num, den):
        Ge = sp.expand(F_.subs(sub, simultaneous=True)).subs(e, tt)
        out.append(norm_mod(Ge, None, sp.expand(he[0]).subs(e, tt), p))
    return out

def test(name, p, npairs, seed, side=None):
    TG.SIDE = side
    P = TG.prepare(name, quiet=True)
    D, C, RD = P['D'], P['C'], P['RD']
    nd, ps = TG.pencil_primes(D, known=sorted(set([2] + RD['bad'])))
    B0 = sorted(set([2] + RD['bad'] + ps))
    M = T.Model(D['QA'], D['QB'])
    places = sorted(set([2, 3, 5, 7, 11, 13] + B0))
    sols, K, Rs = TG.choose_Xi_lin(M, D, C, places, say=lambda s: None)
    Xi = sols[0]
    G = T.build_g0(D, Xi); Nz = T.normalise(G)
    num, den = Nz['R30'].as_expr(), sp.expand(W * Nz['N10'].as_expr())
    numb, denb = sp.expand(num * c), sp.expand(den * e)          # g0 c / e: an odd divisor
    assert Xi % p and all(x % p for x in (2,)) and p not in B0
    rng = random.Random(seed * 1000 + p)
    vals, valsb = [], []
    tries = 0
    while len(vals) < npairs and tries < 6 * npairs:
        tries += 1
        lin = sum(rng.randrange(p) * v_ for v_ in (c, d, e, g))
        r = section_norms(D, Nz, Xi, p, lin, num, den)
        rb = section_norms(D, Nz, Xi, p, lin, numb, denb) if r else None
        if not r or not rb or 0 in r or 0 in rb: continue
        vals.append(r[0] * pow(r[1], -1, p) % p); valsb.append(rb[0] * pow(rb[1], -1, p) % p)
    assert len(vals) >= 2, 'not enough usable sections'
    sq = lambda x: pow(x, (p - 1) // 2, p) == 1
    ok = sum(1 for i in range(1, len(vals)) if sq(vals[i] * pow(vals[0], -1, p) % p))
    okb = sum(1 for i in range(1, len(valsb)) if sq(valsb[i] * pow(valsb[0], -1, p) % p))
    n = len(vals) - 1
    say('  %s, Xi = %d, p = %d: %d pairs of hyperplane sections: g0(Z1) / g0(Z2) a square for %d of %d; control '
        'g0 c / e: %d of %d' % (P['cv']['name'], Xi, p, n, ok, n, okb, n))
    return ok == n, okb < n

def main():
    say('rank_fisher_theta4_gen_weil.py: Weil reciprocity modulo p, a test of div g0 = 2B on C4 (not a proof)')
    allok = True
    for (name, side) in (('920,129', None), ('760,39', None), ('631,71', "E'"), ('72,65', None)):
        ctl = False
        for p in (10007, 100003):
            a_, b_ = test(name, p, 12, 1, side)
            allok = allok and a_; ctl = ctl or b_
        assert ctl, 'THE ODD-DIVISOR CONTROL NEVER FAILED on %s' % name
    assert allok, 'A WEIL TEST FAILED'
    say('all tests passed: every ratio for g0 a square; the odd-divisor control g0 c / e fails on every curve')
    if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_weil.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    main()
