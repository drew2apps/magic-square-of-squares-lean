import sys, json, time; sys.path.insert(0, '.')
import t11_lib as L
import sympy as sp
import rank_fisher5 as R5
import rank_fisher_theta4_find as TF
from cysignals.alarm import alarm, AlarmInterrupt, cancel_alarm
fn, name, side, dd, PAIR = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
C3 = json.load(open(fn)); R = L.analyse3(name, side, C3, dd)
pari_ = TF.pari; t = TF.t
if PAIR == 'U': h0, Bq, v1, v2, ddd = R['h0'], R['Bq'], R5.e, R5.g, R['dU']
else: h0, Bq, v1, v2, ddd = R['coneV'][3], R['coneV'][4], R5.c, R5.d, R['dV']
t0 = time.time()
D0, m = TF.sqfree_split(ddd); print('sqfree_split', time.time() - t0, D0, flush=True)
K = TF.QF(D0); print('bnfinit', time.time() - t0, flush=True)
def topari(expr):
    P_ = sp.Poly(sp.expand(expr), t)
    return K.el(P_.coeff_monomial(1), '(%s)*%d' % (P_.coeff_monomial(t), m))
P = sp.Poly(Bq, v1, v2)
b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (v1 ** 2, v1 * v2, v2 ** 2)]
h0p = topari(h0)
A1 = -4 * b20 * h0p; B1 = b11 ** 2 - 4 * b20 * b02
for z in (A1, B1):
    zl = z.lift(); den = pari_.denominator(pari_.content(zl)); n = abs(int(pari_.norm(z * den * den)))
    fa = pari_.factor(n, 10**6)
    print('norm digits', len(str(n)), 'small-factor part leaves cofactor', [len(str(int(p))) for p in fa[0]][-3:], flush=True)
C = L.curve('Sig+(12204,9853)'); Rd = L.RI.first_descent(C['a'], C['b'])
import rank_fisher_theta4_gen as G
D = dict(QA=sp.expand(R['QA']), QB=sp.expand(R['QB']))
nd, ps = G.pencil_primes(D, known=sorted(set([2] + Rd['bad'])))
known = sorted(set(Rd['bad'] + ps + list(L.RI.factor_int(abs(int(C3['xi'])))) + (list(L.RI.factor_int(abs(int(C3['lam'])))) if abs(int(C3['lam'])) > 1 else [])))
print('known primes', known)
for z in (A1, B1):
    zl = z.lift(); den = pari_.denominator(pari_.content(zl)); n = abs(int(pari_.norm(z * den * den)))
    for p in known + [q for q in range(2, 10**5) if pari_.isprime(q)]:
        while n % p == 0: n //= p
    print('cofactor digits after known primes and p < 1e5:', len(str(n)), flush=True)
