#!/usr/bin/env python3
"""kctp_isog_exact.py: the period (area) ratios B(E')/B(E) of kctp_isog.py section 5 for the three curves 2-isogenous to
the corner curve E = L1L3M1M2 of (627,364) over K = Q(i), with the floating-point choice removed.  EXPERIMENT; written
by an AI model (Claude), not read by a person.

    python3 kctp_isog_exact.py [--write]             -> kctp_isog_exact.txt
    python3 kctp_isog_exact.py --negative [--write]  -> kctp_isog_exact_negative.txt

WHAT kctp_isog.py DOES.  For phi: Emin -> C (Velu, Sage), C' the global minimal model of C, u the parameter of the
K-isomorphism C -> C', it takes area(C')/area(Emin) to be N(u)/2 or 1/(2 N(u)), "the sign of the exponent, fixed by the
convention, read off the floating-point area ratio" (the period lattices at 100 bits, tolerance 1e-12).
THE EXACT ARGUMENT USED HERE (quoted facts: Silverman AEC VI.3-VI.5 and VI.4 / III.5).
  (i) For a Weierstrass model W over C with invariant differential w = dx/(2y + a1 x + a3), the period lattice L(W) of
      w has g2(L) = c4/12, g3(L) = c6/216, so Delta(L) := g2^3 - 27 g3^2 = Delta(W).  Delta is homogeneous of weight 12:
      Delta(mu L) = mu^-12 Delta(L).
 (ii) Two models W, W' of the same curve over C have L(W') = mu L(W) for some mu in C* (their differentials differ by a
      scalar).  So area L(W') / area L(W) = |mu|^2 = |Delta(W) / Delta(W')|^(1/6).  For z in K and the embedding
      K -> C, |z|^2 = N_{K/Q}(z).  Hence (area C' / area C)^12 = N(Delta(C)) / N(Delta(C')), a rational number, and
      the ratio is its positive 12th root -- whatever convention Sage uses for u.  (Checked: it is a 12th power, and
      Delta(C') = u^-12 Delta(C) holds exactly for Sage's u, so the ratio is N(u).)
(iii) An isogeny phi: C/L -> C/L' with phi^* w' = c w is z -> c z, and c L is a sublattice of L' of index deg phi.  Velu's
      isogeny is normalised, c = 1 (checked exactly: Sage's phi.is_normalized(), which compares phi^* w' with w by the
      formal group); so area L(C) = area L(Emin) / 2.
  So area(C') / area(Emin) = (N(Delta(C)) / N(Delta(C')))^(1/12) / 2, exactly.  With Omega = 2 area for the complex
  place of K (the same factor for all curves), Tamagawa products (Sage's Tate's algorithm) and torsion orders (Sage),
  B(C')/B(E) = area ratio * prod c(C') / prod c(Emin) * |Emin(K)_tors|^2 / |C'(K)_tors|^2, a rational number.  It
  is compared with kctp_isog.py's printed values (1/16, 1/256 and, for E3 = E/<T3>, 4) and with Sage's floating-point
  area ratio (a check, not used).
TRUSTED.  Sage for exact arithmetic in K, Velu's formulas, global minimal models, isomorphisms, Tate's algorithm and the
torsion subgroup over K; the quoted facts (i)-(iii); Cassels' isogeny invariance of the BSD quotient (used by
kctp_isog.py, not here).  No floating point enters a value printed as exact.
--negative: the opposite convention (1/(2 N(u))), a non-normalised isogeny (phi followed by a rescaling by 1 + i), a
quadratic twist of C' by 1 + i in place of C', and a forged parameter u = 1 + 2i must each be refused.
"""
from sage.all__sagemath_schemes import *
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kctp_recheck as R

K = R.K; I = R.I
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

PRINTED = {'E/<T1>': QQ(1) / 256, 'E/<T2>': QQ(1) / 16, 'E3 = E/<T3>': QQ(4)}      # kctp_isog.txt, section 5

def twelfth_root(q):
    q = QQ(q); assert q > 0
    n, d = ZZ(q.numerator()), ZZ(q.denominator())
    rn, rd = n.nth_root(12, truncate_mode=True), d.nth_root(12, truncate_mode=True)
    if not (rn[1] and rd[1]): return None
    return QQ(rn[0]) / QQ(rd[0])

def exact_ratio(Emin, phi, Cm, iso):
    """area(Cm)/area(Emin) by (i)-(iii); returns (ratio, details) or (None, reason)"""
    C = phi.codomain()
    if not phi.is_normalized(): return None, 'the isogeny is not normalised (phi^* w\' != w)'
    if iso.domain() != C or iso.codomain() != Cm: return None, 'the isomorphism does not go from the codomain to the model'
    if C.is_isomorphic(Cm) is False: return None, 'C\' is not isomorphic to C over K'
    u = K(iso.u)
    if Cm.discriminant() != C.discriminant() / u ** 12: return None, 'Delta(C\') != u^-12 Delta(C)'
    rho = twelfth_root(QQ(K(C.discriminant()).norm()) / QQ(K(Cm.discriminant()).norm()))
    if rho is None: return None, 'N(Delta(C)) / N(Delta(C\')) is not a 12th power'
    return rho / 2, dict(Nu=QQ(u.norm()), rho=rho)

def tam(F): return prod(F.local_data(P).tamagawa_number() for P in F.conductor().prime_factors())

def setup():
    Q, fs = R.quartic(627, 364, 'L1L3M1M2'); es = R.full2_model(Q, fs); E = R.curve_from_e(es)
    Emin = E.global_minimal_model()
    T = [E([ee, 0]) for ee in es]
    names = {}
    for k, t in enumerate(T):
        names[k] = ['E/<T1>', 'E/<T2>', 'E3 = E/<T3>'][k]
    cods = {names[k]: E.isogeny(t).codomain() for k, t in enumerate(T)}
    return E, Emin, cods

def float_ratio(Emin, Cm):
    emb = K.embeddings(ComplexField(100))[0]
    ar = lambda F: abs((lambda v: (v[0].conjugate() * v[1]).imag())(F.period_lattice(emb).basis()))
    return ar(Cm) / ar(Emin)

def main(write):
    say('kctp_isog_exact.py: B(E\')/B(E) for the three curves 2-isogenous to L1L3M1M2 of (627,364) over Q(i), exact.')
    say('Written by an AI model; not read by a person.  Experiment, not counted.')
    E, Emin, cods = setup()
    say('Emin = %s (global minimal model of E); prod c_v = %s; |E(K)_tors| = %d' % (Emin.ainvs(), tam(Emin), Emin.torsion_order()))
    res = {}
    for r in Emin.two_division_polynomial().roots(multiplicities=False):
        phi = Emin.isogeny(Emin.lift_x(r)); C = phi.codomain(); Cm = C.global_minimal_model(); iso = C.isomorphism_to(Cm)
        nm = [n_ for n_, F in cods.items() if F.is_isomorphic(Cm)]; assert len(nm) == 1; nm = nm[0]
        ratio, d = exact_ratio(Emin, phi, Cm, iso)
        assert ratio is not None, d
        q = QQ(Emin.torsion_order() ** 2 * tam(Cm)) / (Cm.torsion_order() ** 2 * tam(Emin)) * ratio
        fl = float_ratio(Emin, Cm)
        alt = 1 / (2 * d['Nu'])
        say('%-12s Velu normalised: True; N(u) = %s; Delta(C\') = u^-12 Delta(C) exactly; (N Delta(C) / N Delta(C\'))^(1/12) = %s;'
            % (nm, d['Nu'], d['rho']))
        say('             area(C\')/area(Emin) = %s (exact; = N(u)/2: %s); prod c_v = %s, |tors| = %d; B(E\')/B(E) = %s (kctp_isog.txt: %s)'
            % (ratio, ratio == d['Nu'] / 2, tam(Cm), Cm.torsion_order(), q, PRINTED[nm]))
        say('             check, not used: Sage\'s floating-point area ratio %.15f; the rejected alternative 1/(2 N(u)) = %s' % (float(fl), alt))
        assert ratio == d['Nu'] / 2 and q == PRINTED[nm] and abs(float(fl) / float(ratio) - 1) < 1e-12
        res[nm] = q
    assert res['E3 = E/<T3>'] == 4
    say('')
    say('RESULT: the three ratios of kctp_isog.py section 5 hold exactly (B(E3)/B(E) = 4, i.e. |Sha(E3)|/|Sha(E)| = 2^-2 under')
    say('rank 0 and finite Sha, by Cassels\' isogeny invariance); the convention N(u)/2 that kctp_isog.py read off the')
    say('floating-point area ratio is the one forced by the discriminants.  No floating point is used for this.')
    if write: open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kctp_isog_exact.txt'), 'w').write('\n'.join(OUT) + '\n')

def negative(write):
    say('kctp_isog_exact.py --negative: inputs that must be refused')
    E, Emin, cods = setup()
    ntot = nref = 0
    r = [rr for rr in Emin.two_division_polynomial().roots(multiplicities=False)
         if Emin.isogeny(Emin.lift_x(rr)).codomain().global_minimal_model().is_isomorphic(cods['E3 = E/<T3>'])][0]
    phi = Emin.isogeny(Emin.lift_x(r)); C = phi.codomain(); Cm = C.global_minimal_model(); iso = C.isomorphism_to(Cm)
    ratio, d = exact_ratio(Emin, phi, Cm, iso); assert ratio == 1
    # (a) the opposite convention
    alt = 1 / (2 * d['Nu']); ntot += 1; nref += alt != ratio
    say('forged convention 1/(2 N(u)) = %s for E3: %s' % (alt, 'refused (the discriminants force %s)' % ratio if alt != ratio else 'NOT REFUSED'))
    # (b) a non-normalised isogeny: phi followed by the K-isomorphism C -> Cs, Cs = C rescaled by u = 1 + i
    Cs = C.change_weierstrass_model(1 + I, 0, 0, 0)
    psi = C.isomorphism_to(Cs) * phi
    Csm = Cs.global_minimal_model()
    rr, why = exact_ratio(Emin, psi, Csm, Cs.isomorphism_to(Csm))
    ntot += 1; nref += rr is None
    say('forged non-normalised isogeny (Velu followed by a rescaling by 1 + i; scaling factor %s): %s' % (
        psi.scaling_factor(), 'refused (%s)' % why if rr is None else 'NOT REFUSED'))
    # (c) a quadratic twist of C' by 1 + i in place of C'
    Ct = Cm.quadratic_twist(1 + I)
    try:
        isot = C.isomorphism_to(Ct); rt, why = exact_ratio(Emin, phi, Ct, isot)
    except Exception as ex:
        rt, why = None, 'no K-isomorphism C -> C\' (%s)' % str(ex).split('\n')[0][:60]
    ntot += 1; nref += rt is None
    say('forged C\' = twist of the minimal model by 1 + i: %s' % ('refused (%s)' % why if rt is None else 'NOT REFUSED'))
    # (d) the isomorphism to a different model of C (a non-minimal one) with the minimal model's Delta claimed
    Cn = C.change_weierstrass_model(1 + 2 * I, 0, 0, 0)
    class Fake:                                   # an isomorphism C -> Cm that claims Sage's u of C -> Cn
        def __init__(s_, f): s_.u = f.u; s_._d, s_._c = C, Cm
        def domain(s_): return s_._d
        def codomain(s_): return s_._c
    rf, why = exact_ratio(Emin, phi, Cm, Fake(C.isomorphism_to(Cn)))
    ntot += 1; nref += rf is None
    say('forged u = 1 + 2i (the parameter of C -> another model, claimed for C -> C\'): %s' % ('refused (%s)' % why if rf is None else 'NOT REFUSED'))
    say('positive control, the true data for E3: area ratio %s, holds' % ratio)
    say('negative tests: %d forged inputs, %d refused' % (ntot, nref))
    assert ntot == nref
    if write: open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'kctp_isog_exact_negative.txt'), 'w').write('\n'.join(OUT) + '\n')

if __name__ == '__main__':
    if '--negative' in sys.argv: negative('--write' in sys.argv)
    else: main('--write' in sys.argv)
