import sys, json; sys.path.insert(0, '.')
import t11_lib as L
import sympy as sp
import rank_fisher5 as R5
import rank_fisher_theta4_norm as NM
fn, name, side, dd, PAIR = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5]
C3 = json.load(open(fn))
R = L.analyse3(name, side, C3, dd)
pari_ = NM.pari; t_ = NM.t_
if PAIR == 'U': h0, Bq, v1, v2, ddd = R['h0'], R['Bq'], R5.e, R5.g, R['dU']
else: h0, Bq, v1, v2, ddd = R['coneV'][3], R['coneV'][4], R5.c, R5.d, R['dV']
D0, m = NM.sqfree_split(ddd)
def topari(expr):
    P_ = sp.Poly(sp.expand(expr), t_)
    return pari_('Mod((%s) + (%s)*%d*y, y^2 - (%d))' % (P_.coeff_monomial(1), P_.coeff_monomial(t_), m, D0))
P = sp.Poly(Bq, v1, v2)
b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (v1 ** 2, v1 * v2, v2 ** 2)]
h0p = topari(h0)
A1 = -4 * b20 * h0p; B1 = b11 ** 2 - 4 * b20 * b02
print('D0', D0, 'm digits', len(str(m)))
print('A1', len(str(A1)), 'B1', len(str(B1)))
K = NM.Fld(D0)
for z in (A1, B1):
    zz, s, J, Jn = NM.sqred(K, z)
    print('sqred: N(J) digits', len(str(Jn)), 'zz digits', len(str(pari_.lift(zz))), 'norm zz', len(str(abs(int(pari_.norm(zz))))))
    print('   factor N', pari_.factor(abs(int(pari_.norm(zz)))))
