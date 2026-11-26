"""FINDER, not certifying: a point of Gamma_1 over L = Q(sqrt -2) for one level-3 curve, written as a hint to
rank_fisher5_hints.json (with --write); rank_fisher5.py re-checks it exactly.  PARI (cypari2) factors the
discriminant and reduces the conic coefficients modulo squares; rank_fisher5_legendre.py solves the reduced conic
by Legendre descent in Z[sqrt -2].  Usage: python3 rank_fisher5_find.py 2 "E'" [--pair=V] [--write]
Written by an AI model (Claude); not read by a person."""
import sys, json, time
import os
here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
import sympy as sp
import rank_fisher5 as F
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
k = int(sys.argv[1]); side = sys.argv[2]
PAIR = 'V' if '--pair=V' in sys.argv else 'U'      # U: Gamma_1 over L = Q(sqrt disc U); V: Gamma_3 over L' = Q(sqrt disc V)
js = json.load(open(os.path.join(here, 'rank_fisher_curves.json')))
C = js[k]['sides'][side]['level3_curve']
R = F.analyse(js[k]['name'], side, C)
t = sp.Symbol('t'); e, g = F.e, F.g
dU = R['dU']; t0 = time.time()
if PAIR == 'V':
    red_, t, Mw_, h0V, BqV = R['coneV']
    dU = R['dV']; e, g = F.c, F.d           # the conic h0 Z^2 + Bq(c, d) = 0 over L'
    R = dict(R, h0=h0V, Bq=BqV)
fa = pari.factor(abs(dU))
D0 = -1 if dU < 0 else 1; m = 1
for p_, ex in zip(fa[0], fa[1]):
    D0 *= int(p_) ** (int(ex) % 2); m *= int(p_) ** (int(ex) // 2)
assert m * m * D0 == dU
print('L = Q(sqrt(%d)), t = %d sqrt(D0)' % (D0, m), flush=True)
def topari(expr):
    P = sp.Poly(sp.expand(expr), t)
    a = P.coeff_monomial(1); b = P.coeff_monomial(t)
    return pari('Mod((%s) + (%s)*%d*y, y^2 - (%d))' % (a, b, m, D0))
P = sp.Poly(R['Bq'], e, g)
b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (e**2, e*g, g**2)]
h0 = topari(R['h0'])
A1 = -4 * b20 * h0; B1 = b11**2 - 4 * b20 * b02          # E^2 - A1 Z^2 = B1 g^2, E = 2 b20 e + b11 g
bnf = pari.bnfinit(pari('y^2 - (%d)' % D0), 1)
print('class number', bnf.bnf_get_no(), time.time() - t0, flush=True)
def sqred(z):
    """z = zr * s^2 with zr small: returns (zr, s)"""
    z = z.lift(); s = pari('Mod(1, y^2 - (%d))' % D0)
    dd = pari.denominator(pari.content(z)); zz = pari('Mod(%s, y^2 - (%d))' % (z * dd * dd, D0)); s = s / dd
    fz = pari.idealfactor(bnf, zz.lift())
    for i in range(int(pari.matsize(fz)[0])):
        Pr, ex = fz[i, 0], int(fz[i, 1])
        if ex >= 2:
            ip = pari.bnfisprincipal(bnf, Pr)
            assert pari.length(pari.vecsort(ip[0])) == 0 or all(int(c_) == 0 for c_ in ip[0])
            pi = pari.nfbasistoalg(bnf, ip[1])
            zz = zz / pi ** (2 * (ex // 2)); s = s * pi ** (ex // 2)
    return zz, s
A2, sa = sqred(A1); B2, sb = sqred(B1)        # A1 = A2 sa^2 ... (A2 = A1 / (dd^-1 pi^k)^2)
print('A2', A2, 'B2', B2, time.time() - t0, flush=True)
# E^2 - A1 Z^2 = B1 g^2 with A1 = A2 sa^-2?  we have A2 = A1 * (1/sa)^2 ... check:
assert A2 == A1 / sa ** 2 and B2 == B1 / sb ** 2
import subprocess
def ab(z):
    z = z.lift(); return int(pari.polcoef(z, 0, 'y')), int(pari.polcoef(z, 1, 'y'))
assert D0 in (-2, 2)
out_ = subprocess.run([sys.executable, '-B', os.path.join(here, 'rank_fisher5_legendre.py')] + [str(v) for v in ab(A2) + ab(B2)],
                      capture_output=True, text=True, timeout=600, env=dict(os.environ, LEG_D=str(-D0)))
sol = json.loads(out_.stdout)
el = lambda v: pari('Mod(%d + %d*y, y^2 - (%d))' % (v[0], v[1], D0))
xx, yy, zz = el(sol['x']), el(sol['y']), el(sol['z'])
assert xx ** 2 == A2 * yy ** 2 + B2 * zz ** 2
E = xx; Z = yy / sa; gg = zz / sb
assert E ** 2 - A1 * Z ** 2 == B1 * gg ** 2
ee = (E - b11 * gg) / (2 * b20)
assert h0 * Z ** 2 + b20 * ee ** 2 + b11 * ee * gg + b02 * gg ** 2 == 0
def out(z):
    z = z.lift(); return [str(pari.polcoef(z, 0, 'y')), str(pari.polcoef(z, 1, 'y'))]
hint = dict(D0=D0, m=m, Z=out(Z), e=out(ee), g=out(gg)) if PAIR == 'U' else dict(D0=D0, m=m, Z=out(Z), c=out(ee), d=out(gg))
print('HINT', json.dumps(hint), flush=True)
if '--write' in sys.argv:
    hf = os.path.join(here, 'rank_fisher5_hints.json')
    H = json.load(open(hf)) if os.path.exists(hf) else {}
    H['%s|%s' % (js[k]['name'], side) + ('|V' if PAIR == 'V' else '')] = hint
    json.dump(H, open(hf, 'w'), indent=1, sort_keys=True)
