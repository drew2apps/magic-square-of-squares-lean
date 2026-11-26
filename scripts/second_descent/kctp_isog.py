#!/usr/bin/env python3
"""kctp_isog.py: the isogeny class over K = Q(i) of the corner curve L1L3M1M2 of (627,364), why no 2-descent on any
curve of it can prove rank 0, and the BSD scale of the curve.  EXPERIMENT; written by an AI model (Claude), not read by
a person.

    python3 kctp_isog.py [--write]          -> kctp_isog.txt
    python3 kctp_isog.py --simon [--write]  -> also Sage's simon_two_descent bounds for the four other curves (slow)

EXACT PARTS (certifying, given the quoted theory):
  * E3 = E/<T3> (Velu, exact) has all of its 2-torsion over K; kctp_recheck's new code gives dim Sel^2(E3/K) = 4 and
    every class lifts to an everywhere locally soluble 4-covering (pairing zero), as kctp.py found.
  * The curves 2-isogenous to E other than E3 (kernels T1, T2) and those 2-isogenous to E3 other than E (kernels the
    two other points of order 2 of E3) have exactly one point of order 2 over K (exact: the 2-division polynomial).
LEMMA (written proof; quoted theory: the Cassels-Tate pairing is alternating and non-degenerate on a finite Sha, and
  compatible with isogenies).  Let E be the corner curve and suppose rank E(K) = 0 and Sha(E/K)[2^inf] finite.  Then
  Sha(E'/K)[2] != 0 for each of the six curves E' of the isogeny class (E, E/<T1>, E/<T2>, E3 = E/<T3>, and the two
  curves E3/<S1>, E3/<S2> 2-isogenous to E3 other than E).  Hence dim Sel^2(E'/K) >= dim E'(K)[2] + 2 for each: no
  2-descent on any curve of the class gives rank 0.
  Proof.  Sel^2(E/K) = E(K)/2E(K) + Sha(E)[2] with dim E(K)/2E(K) = 2 (torsion Z/2 x Z/4, rank 0) and dim Sel^2 = 4,
  so Sha(E)[2] = (Z/2)^2.  The pairing is zero on Sel^2, so Sha(E)[2] is orthogonal to itself; for finite Sha[2^inf]
  the orthogonal of Sha[2] is 2 Sha[2^inf], so Sha(E)[2] lies in 2 Sha(E): [2] != 0 on Sha(E)[2^inf].  For phi: E -> E'
  of degree 2, [2] = phi^_* phi_* on Sha(E), so phi_* != 0 on Sha(E)[2^inf]: Sha(E')[2^inf] != 0, so Sha(E')[2] != 0.
  Finiteness of Sha[2^inf] is an isogeny invariant, so Sha(E')[2^inf] = A x A (alternating pairing) and dim Sha(E')[2]
  is even, >= 2.  The same argument with E3 in place of E (E3: dim Sel^2 = 4, torsion image of dimension 2, pairing
  zero -- checked below) gives E3/<S1>, E3/<S2>.  With rank 0, dim Sel^2(E') = dim E'(K)[2] + dim Sha(E')[2].
  WHAT THE LEMMA DOES NOT EXCLUDE.  For the four curves with one point of order 2 over K (E/<T1>, E/<T2>, E3/<S1>,
  E3/<S2>) the level-2 pairing on Sel^2(E'/K) has not been computed (it needs the quadratic extension K(E'[2]) of K).
  Under rank 0 it is non-zero on all of Sha(E')[2] exactly when Sha(E')[2^inf] is killed by 2; then it gives rank 0.
  REFINEMENT (section 5 of the output; uses the ratio |Sha(E3)| / |Sha(E)| = B(E) / B(E3) = 2^-2 from Cassels'
  isogeny invariance, the ratio computed in floating point): under rank 0, Sha(E)[2^inf] = (Z/2^n)^2 with n >= 3 and
  Sha(E3)[2^inf] = (Z/2^(n-1))^2.  This excludes the level-2 pairing on E/<T1>, E/<T2> and a 4-descent on E as ways
  to prove rank 0; a 4-descent on E3, or the level-2 pairing on E3/<S1>, E3/<S2>, could only succeed if n = 3.
GUIDES (not certifying): Sage's isogeny class over K (6 curves, degrees 1, 2, 4, 8: all powers of 2; no curve with
  j = conj(j(E)), so E is not a Q-curve if the class is complete); conductor, Tamagawa numbers and the period lattice
  (numerical), giving the BSD scale; the 3-division field.
"""
from sage.all__sagemath_schemes import *
import sys, os, random, time, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kctp_recheck as R

K = R.K; I = R.I
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

def two_torsion_x(F):
    """the x-coordinates over K of the points of order 2 of F (from the 2-division polynomial, exact)."""
    a1, a2, a3, a4, a6 = F.a_invariants()
    Rx = PolynomialRing(K, 'x'); x = Rx.gen()
    psi2 = 4 * x ** 3 + (a1 ** 2 + 4 * a2) * x ** 2 + (2 * a4 + a1 * a3) * 2 * x + (a3 ** 2 + 4 * a6)
    return [r for r, mult in psi2.roots()]

def short_full2(F):
    """F with all of its 2-torsion over K, as y^2 = (x - e1)(x - e2)(x - e3) with e_k in Z[i]."""
    Fs = F.short_weierstrass_model()
    a4, a6 = Fs.a4(), Fs.a6()
    Rx = PolynomialRing(K, 'x'); x = Rx.gen()
    rts = [r for r, mlt in (x ** 3 + a4 * x + a6).roots()]
    assert len(rts) == 3
    den = lcm([r.denominator() for r in rts])
    es = [r * den ** 2 for r in rts]
    assert R.curve_from_e(es).is_isomorphic(F)
    return es

def main():
    rng = random.Random(20261126)
    say('kctp_isog.py: the isogeny class over K = Q(i) of L1L3M1M2 of (627,364).  Written by an AI model; not read by a')
    say('person.  Experiment, not counted.')
    Q, fs = R.quartic(627, 364, 'L1L3M1M2'); es = R.full2_model(Q, fs); E = R.curve_from_e(es)
    T = [E([ee, 0]) for ee in es]
    # exact Velu isogenies
    phis = [E.isogeny(t) for t in T]
    E1, E2, E3 = [p.codomain() for p in phis]
    say('1. The 2-isogenous curves (Velu, exact): points of order 2 over K: E/<T1>: %d, E/<T2>: %d, E/<T3>: %d' % (
        len(two_torsion_x(E1)), len(two_torsion_x(E2)), len(two_torsion_x(E3))))
    assert len(two_torsion_x(E1)) == 1 and len(two_torsion_x(E2)) == 1 and len(two_torsion_x(E3)) == 3
    T3s = [E3.lift_x(xx) for xx in two_torsion_x(E3)]
    # the dual kernel: the image of T1 (= image of T2)
    dual = phis[2](T[0]); assert dual != 0 and 2 * dual == 0 and phis[2](T[1]) == dual
    others = [t for t in T3s if t != dual]
    E3b = [E3.isogeny(t).codomain() for t in others]
    say('   the curves 2-isogenous to E3 other than E (kernels the two points of order 2 of E3 other than the image of T1):')
    say('   points of order 2 over K: %s' % [len(two_torsion_x(F)) for F in E3b])
    assert all(len(two_torsion_x(F)) == 1 for F in E3b)
    assert E3.isogeny(dual).codomain().is_isomorphic(E)
    curves = {'E': E, 'E/<T1>': E1, 'E/<T2>': E2, 'E3 = E/<T3>': E3, 'E3/<S1>': E3b[0], 'E3/<S2>': E3b[1]}
    # Sage's isogeny class (guide)
    t0 = time.time(); C = E.isogeny_class()
    say('2. Sage isogeny_class over K (guide): %d curves, degree matrix rows %s (%.0f s)' % (
        len(C), [list(r) for r in C.matrix().rows()], time.time() - t0))
    assert len(C) == 6 and all(set(ZZ(d).prime_factors()) <= {2} for r in C.matrix().rows() for d in r if d != 0)
    match = all(any(F.is_isomorphic(G) for G in C.curves) for F in curves.values())
    say('   the six curves above are, up to isomorphism, the six curves of the class: %s' % match)
    assert match
    jbar = E.j_invariant().conjugate() if hasattr(E.j_invariant(), 'conjugate') else None
    jb = K(list(E.j_invariant())[0]) - K(list(E.j_invariant())[1]) * I
    say('   a curve of the class with j = conj(j(E)): %s (so E is not a Q-curve, if the class is complete)' % any(G.j_invariant() == jb for G in C.curves))
    # E3: Sel^2 and the lift criterion with the new code (exact)
    es3 = short_full2(E3)
    t0 = time.time()
    sel3 = R.Selmer(es3, rng)
    timg = []
    halves = []
    E3m = R.curve_from_e(es3)
    for t in [E3m([ee, 0]) for ee in es3]:
        try: halves += [H for H in t.division_points(2)]
        except Exception: pass
    for xx in [es3[0], es3[1]] + [H.xy()[0] for H in halves]:
        d1, d2 = R.kummer(xx, es3); timg.append(sel3.vec_of(d1) + sel3.vec_of(d2))
    tor3 = E3.torsion_subgroup().invariants()
    say('3. E3 = E/<T3> (model with all 2-torsion over K, e in Z[i]): dim Sel^2(E3/K) = %d, torsion image dim %d' % (
        sel3.dim, R.f2_rank(timg)))
    say('   (Sage torsion of E3, guide: %s)' % (tor3,))
    assert sel3.dim == 4 and R.f2_rank(timg) == 2
    lifts = []
    for b in sel3.basis:
        ok, dims, _, _, _ = R.lift_data(sel3, sel3.pair(b), rng)
        expect = {str(pl.P.gens_reduced()[0]): (4 if pl.two else 2) for pl in sel3.places}
        assert dims == expect
        lifts.append(ok)
    say('   lift criterion for the 4 basis classes: %s -> the pairing on Sel^2(E3/K) is zero  (%.0f s)' % (lifts, time.time() - t0))
    assert all(lifts)
    say('4. LEMMA (proof in the header): if rank E(K) = 0 and Sha(E/K)[2^inf] is finite, then Sha(E\'/K)[2] != 0 for')
    say('   all six curves E\' of the class, so dim Sel^2(E\'/K) >= dim E\'(K)[2] + 2 for each: no 2-descent on any curve')
    say('   of the isogeny class proves rank 0.  The level-2 pairing is zero on Sel^2 of E and of E3 (both computed).  On')
    say('   the four curves with one point of order 2 over K it is not computed; under rank 0 it would prove rank 0')
    say('   exactly when Sha(E\'/K)[2^inf] is killed by 2, which the lemma does not exclude (open route).')
    # BSD scale (guide)
    Em = E.global_minimal_model()
    cond = Em.conductor(); Nc = cond.norm()
    tam = [(ZZ(P.norm()), str(Em.local_data(P).kodaira_symbol()), Em.local_data(P).tamagawa_number()) for P in cond.prime_factors()]
    ctam = prod(t for _, _, t in tam)
    emb = K.embeddings(ComplexField(80))[0]
    w1, w2 = Em.period_lattice(emb).basis()
    area = abs((w1.conjugate() * w2).imag())
    Om = 2 * area
    B = Om * ctam / (64 * 2)
    say('5. BSD scale (guide, floating point).  Normalisation check on base changes of y^2 = x(x + a)(x + b) of rank 0')
    say('   over Q for E and E^-1: L(E/K,1) = L(E,1) L(E^-1,1) (PARI lfun) against B = 2 area prod c_v / (|tors|^2 * 2):')
    import cypari2
    pari = cypari2.Pari()
    for a_, b_ in [(1, 2), (2, 3), (1, 4)]:
        LQ = 1.0
        for sg in (1, -1):
            Eq = pari.ellinit([0, sg * (a_ + b_), 0, a_ * b_, 0]); assert int(pari.ellrank(Eq)[1]) == 0
            LQ *= float(pari.lfun(Eq, 1))
        Ec = EllipticCurve(K, [0, a_ + b_, 0, a_ * b_, 0]).global_minimal_model()
        cc = prod(Ec.local_data(P).tamagawa_number() for P in Ec.conductor().prime_factors())
        u1, u2 = Ec.period_lattice(K.embeddings(ComplexField(80))[0]).basis()
        Bc = 2 * abs((u1.conjugate() * u2).imag()) * cc / (Ec.torsion_order() ** 2 * 2)
        say('     a=%d b=%d: L(E/K,1) = %.6f, B = %.6f, ratio %.4f (= |Sha(E/K)|)' % (a_, b_, LQ, float(Bc), LQ / float(Bc)))
        assert abs(LQ / float(Bc) - round(LQ / float(Bc))) < 1e-6
    say('   The target: global minimal model %s.' % (Em.ainvs(),))
    nf = pari.nfinit('y^2 + 1')
    ai = [pari('Mod(%s, y^2 + 1)' % str(c).replace('i', 'y')) for c in Em.ainvs()]
    w = int(pari.ellrootno(pari.ellinit(ai, nf)))
    say('   global root number of E/K (PARI ellrootno over the number field; a guide): %+d, so the analytic rank is even' % w)
    say('   conductor norm %s = %s' % (Nc, factor(Nc)))
    say('   Kodaira symbols and Tamagawa numbers %s; product %s = 2^%d' % (tam, ctam, ZZ(ctam).valuation(2)))
    say('   area of the period lattice %.4e; Omega = 2 area = %.4e; |E(K)_tors|^2 = 64; sqrt|d_K| = 2' % (area, Om))
    say('   B = Omega prod c_v / (64 * 2) = %.3e.  BSD over K: L(E/K,1) = B |Sha| if the rank is 0, and' % B)
    say('   L\'\'(E/K,1)/2 = B Reg |Sha| if it is 2.  So rank 0 needs |Sha| = %.2e L(E/K,1), and rank 2 needs' % (1 / B))
    say('   Reg |Sha| = %.2e L\'\'(E/K,1)/2.  For L-values of moderate size this means |Sha| of order 10^5 to 10^6 (rank 0)' % (1 / B))
    say('   or a regulator of that order (rank 2, Sha = 1): generators of canonical height of order 10^3, out of reach of')
    say('   any point search on 2- or 4-coverings (these reach canonical height about 2 h and 4 h for points of naive')
    say('   height h on the covering).  The order of magnitude only; L(E/K, s) itself is out of reach')
    say('   (N_{K/Q}(conductor) d_K^2 = %.1e, about 10^20 Dirichlet coefficients).' % (float(Nc) * 16))
    say('   B for the six curves of the class (same L-function; Cassels: the BSD quotient is an isogeny invariant, so with')
    say('   rank 0 and finite Sha, |Sha(E\')| / |Sha(E)| = B(E) / B(E\') exactly; floating point here, a guide):')
    def Bof(F):
        F = F.global_minimal_model()
        cF = prod(F.local_data(P).tamagawa_number() for P in F.conductor().prime_factors())
        v1, v2 = F.period_lattice(K.embeddings(ComplexField(80))[0]).basis()
        return 2 * abs((v1.conjugate() * v2).imag()) * cF / (F.torsion_order() ** 2 * 2), cF, F.torsion_order()
    BE = Bof(E)[0]; ks = []; kd = {}
    for name, F in curves.items():
        Bf, cF, tF = Bof(F)
        ratio = float(BE / Bf); k = int(round(math.log2(ratio)))
        say('     %-12s prod c_v = 2^%d, |tors| = %d, B = %.4e, B(E)/B = %.6f = 2^%d' % (name, ZZ(cF).valuation(2), tF, float(Bf), ratio, k))
        assert abs(ratio / 2 ** k - 1) < 1e-8
        ks.append(k); kd[name] = k
    say('   exponents %s: %s' % (ks, 'all even, as needed for rank 0 with finite Sha (|Sha| a square for each curve)' if all(k % 2 == 0 for k in ks)
        else 'SOME ODD: rank 0 with finite Sha would be impossible'))
    say('   Exact form of the ratios for the three curves 2-isogenous to E: Sage\'s Velu isogeny from the minimal model is')
    say('   normalised (phi^* omega\' = omega), so the period lattice of the codomain contains that of E with index 2 (area')
    say('   ratio 1/2), and the minimal model of the codomain changes the area by N(u)^(+-1) (u from the K-isomorphism; the')
    say('   sign of the exponent, fixed by the convention, is read off the floating-point area ratio).  With the Tamagawa')
    say('   numbers (Tate\'s algorithm) and the torsion orders, B(E\')/B(E) is then a rational number:')
    Emin = E.global_minimal_model()
    aE = abs((lambda v: (v[0].conjugate() * v[1]).imag())(Emin.period_lattice(K.embeddings(ComplexField(100))[0]).basis()))
    tamq = lambda F: prod(F.local_data(P).tamagawa_number() for P in F.conductor().prime_factors())
    exact = {}
    for r in Emin.two_division_polynomial().roots(multiplicities=False):
        Cc = Emin.isogeny(Emin.lift_x(r)).codomain(); Cm = Cc.global_minimal_model(); u = Cc.isomorphism_to(Cm).u
        v = Cm.period_lattice(K.embeddings(ComplexField(100))[0]).basis(); fl = abs((v[0].conjugate() * v[1]).imag()) / aE
        Nu = QQ(K(u).norm())
        conv = [c for c in (Nu / 2, 1 / (2 * Nu)) if abs(float(c) / float(fl) - 1) < 1e-12]
        assert len(conv) == 1
        q = QQ(Emin.torsion_order() ** 2 * tamq(Cm)) / (Cm.torsion_order() ** 2 * tamq(Emin)) * conv[0]
        nm = [n_ for n_, F in curves.items() if F.is_isomorphic(Cm)][0]
        exact[nm] = q
        say('     %-12s N(u) = %s, area ratio %s, B(E\')/B(E) = %s' % (nm, Nu, conv[0], q))
    assert exact['E3 = E/<T3>'] == 4
    assert kd['E3 = E/<T3>'] == -2
    say('   CONSEQUENCE (from the ratio for E3, which is 2^-2, and the zero pairings on Sel^2 of E and E3).  Suppose rank 0')
    say('   and Sha finite.  dim Sha(E)[2] = dim Sha(E3)[2] = 2, so Sha(E)[2^inf] = (Z/2^n)^2 and Sha(E3)[2^inf] =')
    say('   (Z/2^n3)^2 (a finite group with an alternating non-degenerate pairing is C x C); the zero pairings give n, n3 >= 2;')
    say('   the ratio gives 2 n3 = 2 n - 2.  Hence n3 = n - 1 >= 2, so n >= 3: Sha(E/K)[2^inf] contains (Z/8)^2.  Then every')
    say('   class of Sha(E)[2] = 2^(n-1) Sha(E)[2^inf] lifts to Sha[8], so a 4-descent on E with its pairing cannot prove rank')
    say('   0 either; on E3 a 4-descent with its pairing proves rank 0 exactly when n3 = 2.  For E/<T1>, E/<T2> (2-isogenous')
    say('   to E): Sha[2^inf] killed by 2 would give 4 Sha(E)[2^inf] = 0 (2 = phi^ phi on E, phi phi^ on E\'), i.e. n = 2:')
    say('   excluded, so their level-2 pairing cannot prove rank 0.  For E3/<S1>, E3/<S2> it could only if n3 = 2.')
    say('   (The ratio 2^-2 for E3 is the exact value B(E3)/B(E) = 4 above; floating point only fixes the convention.')
    say('   Quoted theory: Cassels\' isogeny invariance of the BSD quotient, and the properties of the pairing.)')
    # the 3-division field
    psi3 = Em.division_polynomial(3)
    fac3 = psi3.factor()
    say('6. 3-descent: the 3-division polynomial over K factors as degrees %s (exact).' % [g.degree() for g, _ in fac3])
    assert len(fac3) == 1 and fac3[0][0].degree() == 4
    Lr = K.extension(fac3[0][0].monic(), 'a'); a = Lr.gen()
    fx = a ** 3 + Em.a4() * a + Em.a6() if Em.a1() == 0 and Em.a2() == 0 and Em.a3() == 0 else None
    sq = fx.is_square() if fx is not None else None
    say('   f(x) at a root x of psi_3 is a square in K(x): %s, so a point of order 3 generates a field of degree %d' % (
        sq, 4 if sq else 8))
    say('   over K (%d over Q).  A 3-descent (Schaefer-Stoll) needs class groups and units in that field; if rank 0 and' % (8 if sq else 16))
    say('   Sha[3] = 0, Sel^3(E/K) = 0 would prove rank 0.  Not attempted.')
    if '--simon' in sys.argv:
        say('7. Sage simon_two_descent upper bounds (guide; GRH in PARI\'s class groups):')
        for name, F in curves.items():
            t0 = time.time()
            try:
                r = F.simon_two_descent(lim1=2, lim3=4, limtriv=2, maxprob=20, limbigprime=30)
                say('   %s: lower %s, upper %s (%.0f s)' % (name, r[0], r[1], time.time() - t0))
            except Exception as ex:
                say('   %s: failed (%s)' % (name, str(ex)[:60]))

if __name__ == '__main__':
    main()
    if WRITE:
        open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kctp_isog.txt'), 'w').write('\n'.join(OUT) + '\n')
