#!/usr/bin/env python3
"""FINDER, not certifying: points of the cone conics Gamma_1 (over L = Q(sqrt disc U)) and Gamma_3 (over
L' = Q(sqrt disc V)) of the level-3 curves C3' of rank_fisher5.py, for any curve of rank_fisher_curves.json.
Written as hints to rank_fisher_theta4_gen_hints.json (with --write), in the format of rank_fisher5_hints.json;
rank_fisher5.analyse (called by rank_fisher_theta4_gen.py) re-checks every hint exactly on its cone.

    python3 rank_fisher_theta4_find.py <curve index> [U|V] [--write]

Method (PARI is a finder only): the conic h0 Z^2 + b20 e^2 + b11 e g + b02 g^2 = 0 over L becomes
E^2 - A1 Z^2 = B1 g^2 (E = 2 b20 e + b11 g, A1 = -4 b20 h0, B1 = b11^2 - 4 b20 b02); A1, B1 are reduced modulo
squares (A1 = A2 sa^2: ideal factorisation from the factored norm, idealprimedec, idealval, bnfisprincipal); then
conic_lattice_find.find_point (a Holzer-type lattice in O_L^3, short vectors by qfminim, and a zero on the O_L-span
of two short vectors when B(v1, v2)^2 - Q(v1) Q(v2) is a square in L).  (PARI rnfisnorm was tried first: its
bnfinit of the quartic L(sqrt B2), discriminant about 1e34, did not finish in 15 minutes.)
Written by an AI model (Claude); not read by a person."""
import os, sys, json, time
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import sympy as sp
import rank_fisher5 as F
import cypari2
pari = cypari2.Pari(); pari.allocatemem(3 * 10 ** 9)
pari.default('parisizemax', 6 * 10 ** 9)
HF = os.path.join(here, 'rank_fisher_theta4_gen_hints.json')
t = sp.Symbol('t')

def sqfree_split(n):
    fa = pari.factor(abs(n)); D0 = -1 if n < 0 else 1; m = 1
    for p_, ex in zip(fa[0], fa[1]):
        D0 *= int(p_) ** (int(ex) % 2); m *= int(p_) ** (int(ex) // 2)
    assert m * m * D0 == n
    return D0, m

class QF:
    def __init__(self, D0):
        self.D0 = D0
        self.pol = pari('y^2 - (%d)' % D0)
        self.bnf = pari.bnfinit(self.pol, 1)
    def el(self, a, b):
        return pari('Mod((%s) + (%s)*y, y^2 - (%d))' % (a, b, self.D0))
    def sqred(self, z):
        """z = zz * s^2 with zz integral and (zz) = (odd part of (z)) * (product of class-group generators)^2:
        I = prod P^floor(e/2) = (s) prod g_i^c_i (bnfisprincipal), zz = z / s^2."""
        bnf = self.bnf
        zl = z.lift()
        den = pari.denominator(pari.content(zl))
        zz = pari('Mod(%s, y^2 - (%d))' % (zl * den * den, self.D0)); s = pari('Mod(1, y^2 - (%d))' % self.D0) / den
        n = abs(int(pari.norm(zz)))
        fa = pari.factor(n)
        I = pari.idealhnf(bnf, 1)
        for i in range(len(fa[0])):
            p = int(fa[0][i])
            for Pr in pari.idealprimedec(bnf, p):
                ex = int(pari.idealval(bnf, zz.lift(), Pr))
                if ex >= 2: I = pari.idealmul(bnf, I, pari.idealpow(bnf, Pr, ex // 2))
        ip = pari.bnfisprincipal(bnf, I, 0)
        if any(int(c_) for c_ in ip):
            # a small prime ideal J (norm < 2000, prime to z) in the class of I: I / J = (al), (z / al^2) = odd part J^2
            J = None
            for q_ in range(2, 2000):
                if not pari.isprime(q_) or n % q_ == 0: continue
                for Pq in pari.idealprimedec(bnf, q_):
                    if pari.bnfisprincipal(bnf, Pq, 0) == ip: J = Pq; break
                if J is not None: break
            if J is not None: I = pari.idealdiv(bnf, I, J)
        ip = pari.bnfisprincipal(bnf, I)
        al = pari.nfbasistoalg(bnf, ip[1])
        zz = zz / al ** 2; s = s * al
        assert zz * s ** 2 == z
        return zz, s

def _sqred_rational(self, z):
    """z = zz * s^2 with s rational: rational square factors removed (inert p: p^(2 floor(e/2)); split p:
    p^(2 min); ramified p: p^(2 floor(e/4)))."""
    nf = self.bnf
    zl = z.lift()
    den = pari.denominator(pari.content(zl))
    zz = pari('Mod(%s, y^2 - (%d))' % (zl * den * den, self.D0)); s = pari('Mod(1, y^2 - (%d))' % self.D0) / den
    n = abs(int(pari.norm(zz)))
    fa = pari.factor(n)
    for i in range(len(fa[0])):
        p = int(fa[0][i])
        dec = pari.idealprimedec(nf, p)
        ex = [int(pari.idealval(nf, zz.lift(), Pr)) for Pr in dec]
        if len(dec) == 2: k = min(ex) // 2
        elif int(dec[0][2]) == 2: k = ex[0] // 4
        else: k = ex[0] // 2
        if k: zz = zz / p ** (2 * k); s = s * p ** k
    assert zz * s ** 2 == z
    return zz, s
QF.sqred_rational = _sqred_rational

def find(k, PAIR, verbose=True):
    js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
    cv = js[k]; side = os.environ.get('FIND_SIDE') or cv['theta4_reduction'][0]
    S = cv['sides'][side]; C = S['level3_curve']
    R = F.analyse(cv['name'], side, C, dd=int(S['dd']))
    if PAIR == 'U': h0, Bq, v1, v2, dd = R['h0'], R['Bq'], F.e, F.g, R['dU']
    else: h0, Bq, v1, v2, dd = R['coneV'][3], R['coneV'][4], F.c, F.d, R['dV']
    D0, m = sqfree_split(dd)
    K = QF(D0)
    t0 = time.time()
    def topari(expr):
        P_ = sp.Poly(sp.expand(expr), t)
        return K.el(P_.coeff_monomial(1), '(%s)*%d' % (P_.coeff_monomial(t), m))
    P = sp.Poly(Bq, v1, v2)
    b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (v1 ** 2, v1 * v2, v2 ** 2)]
    h0p = topari(h0)
    assert b20 != 0
    A1 = -4 * b20 * h0p; B1 = b11 ** 2 - 4 * b20 * b02
    try:
        A2, sa = K.sqred(A1); B2, sb = K.sqred(B1)
    except cypari2.PariError:
        # generators of the square parts not computable (large regulator): rational squares only
        A2, sa = K.sqred_rational(A1); B2, sb = K.sqred_rational(B1)
        if verbose: print('  (square parts reduced by rational squares only: generators not available)', flush=True)
    nA, nB = abs(int(pari.norm(A2))), abs(int(pari.norm(B2)))
    if verbose: print('%s %s: L = Q(sqrt %d), class group %s; N(A2) %d digits, N(B2) %d digits (%.1fs)'
                      % (cv['name'], PAIR, D0, K.bnf.bnf_get_cyc(), len(str(nA)), len(str(nB)), time.time() - t0), flush=True)
    if os.environ.get('FIND_SIZES_ONLY'): return None, dict(nA=len(str(nA)), nB=len(str(nB)))
    import conic_lattice_find as CL
    big = os.environ.get('FIND_BIG')
    if big: pt = CL.find_point(D0, A2, B2, tries=int(big), nvec=20000, Rmul=40, npair=2500, verbose=verbose)
    else: pt = CL.find_point(D0, A2, B2, tries=20, nvec=3000, Rmul=8, npair=400, verbose=verbose)
    assert pt is not None, 'the lattice finder found no point'
    E_, Zr, gr = pt          # E^2 - A2 Z^2 - B2 g^2 = 0
    E_ = pari.Mod(pari.lift(E_), K.pol); Zr = pari.Mod(pari.lift(Zr), K.pol); gr = pari.Mod(pari.lift(gr), K.pol)
    assert E_ ** 2 - A2 * Zr ** 2 == B2 * gr ** 2
    Z = Zr / sa; gg = gr / sb
    assert E_ ** 2 - A1 * Z ** 2 == B1 * gg ** 2
    ee = (E_ - b11 * gg) / (2 * b20)
    assert h0p * Z ** 2 + b20 * ee ** 2 + b11 * ee * gg + b02 * gg ** 2 == 0
    def out(z):
        z = z.lift(); return [str(pari.polcoef(z, 0, 'y')), str(pari.polcoef(z, 1, 'y'))]
    if PAIR == 'U': hint = dict(D0=D0, m=m, Z=out(Z), e=out(ee), g=out(gg))
    else: hint = dict(D0=D0, m=m, Z=out(Z), c=out(ee), d=out(gg))
    if verbose: print('  point found (%.1fs), %d digits' % (time.time() - t0, max(len(s_) for v in hint.values() if isinstance(v, list) for s_ in v)), flush=True)
    return '%s|%s' % (cv['name'], side) + ('|V' if PAIR == 'V' else ''), hint

if __name__ == '__main__':
    k = int(sys.argv[1]); PAIR = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] in ('U', 'V') else 'U'
    key, hint = find(k, PAIR)
    outs = [a[6:] for a in sys.argv if a.startswith('--out=')]
    if outs: json.dump({key: hint}, open(outs[0], 'w'), indent=1, sort_keys=True)
    if '--write' in sys.argv:
        H = json.load(open(HF)) if os.path.exists(HF) else {}
        H[key] = hint
        json.dump(H, open(HF, 'w'), indent=1, sort_keys=True)
        print('written', key)
