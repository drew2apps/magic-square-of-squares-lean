import sys, os, io, contextlib
from fractions import Fraction as Fr
from math import isqrt
REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'census-wide')); sys.path.insert(0, os.path.join(REPO, 'scripts', 'g1_full'))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'census2000'))
_argv = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import sd_pari as SD
    import rank_indep as RI
    import dm_pari as DM
sys.argv = _argv
pari = SD.pari
CURVES = [('E1(3540,2581)', 'E1', 3540, 2581, -11740078), ('E1(4440,2231)', 'E1', 4440, 2231, -13342),
          ('E1(6480,4961)', 'E1', 6480, 4961, 709342), ('E1(8249,3720)', 'E1', 8249, 3720, -108415202),
          ('E1(10512,145)+', 'E1', 10512, 145, 220962238), ('E1(10512,145)-', 'E1', 10512, 145, -220962238),
          ('D1(12007,1799)', 'D1', 12007, 1799, None), ('Sig+(12204,9853)', 'Sig+', 12204, 9853, None),
          ('Sig+(22057,2351)', 'Sig+', 22057, 2351, None)]
for lab, src, m, n, d in CURVES:
    if src == 'E1':
        Q, r = SD.family(m, n)['+13']; ain, A = SD.weierstrass(d * Q, r)
    else:
        ain = [Fr(int(v)) for v in [c[0] for c in DM.all_curves(m, n) if c[4] == src][0]]
    print(lab, 'denoms', [x.denominator for x in ain], 'a1,a3', ain[0], ain[2])
    a2, a4, a6 = ain[1], ain[3], ain[4]
    E0 = pari.ellinit([pari(x.numerator)/x.denominator for x in ain]); Em = pari.ellminimalmodel(E0)
    mi = [int(v) for v in Em[:5]]
    print('  minimal model', mi)
    a1_, a2_, a3_, a4_, a6_ = mi
    b2, b4, b6 = a1_*a1_ + 4*a2_, 2*a4_ + a1_*a3_, a3_*a3_ + 4*a6_
    if a1_ == 0 and a3_ == 0: A2, A4, A6 = a2_, a4_, a6_
    else: A2, A4, A6 = b2, 8*b4, 16*b6
    e, a, b = RI.shift_to_two_torsion(A2, A4, A6)
    bp = a*a - 4*b
    print('  a digits', len(str(abs(a))), 'b digits', len(str(abs(b))), 'bp<0?', bp < 0, '-bp square?', bp < 0 and isqrt(-bp)**2 == -bp, 'sqfree(bp)', RI.sqfree_part(bp) if len(str(abs(bp)))<0 else '')
