#!/usr/bin/env python3
"""kctp_isodescent.py: 2-isogeny descents over K = Q(i) in the isogeny class of the corner curve E = L1L3M1M2 of
(627,364), and what they say about the routes to rank 0 left open by kctp_isog.py.  EXPERIMENT; written by an AI model
(Claude), not read by a person.

    python3 kctp_isodescent.py [--write]      -> kctp_isodescent.txt

For A: y^2 = x (x^2 + a x + b) over K with the kernel point (0, 0), phi: A -> B: Y^2 = X (X^2 - 2 a X + a^2 - 4 b), the
Selmer groups Sel^phi(A) (local images of B(K_v) under X -> X mod squares, (0,0) -> a^2 - 4b) and Sel^phihat(B)
(local images of A(K_v) under x -> x, (0,0) -> b) are subgroups of K(S,2), computed with kctp_recheck's local square
classes.  The local images are found from local points; their sizes are the known ones (quoted, Schaefer's lemma):
    |B(K_v)/phi A(K_v)| = |A(K_v)[phi]| * c_v(B)/c_v(A) * |alpha|_v^-1,
with alpha the leading coefficient of phi on the formal groups of the minimal models (|alpha|_v = 1 at odd v; at
1 + i, alpha = u from Sage's normalised Velu isogeny and the K-isomorphism of its codomain to the minimal model,
omega_min = u omega), and alpha alpha^ = 2 for the dual.  Tamagawa numbers from Sage (Tate's algorithm).  The sampled
image is a subgroup of the true one; reaching the known size proves equality (given the quoted size); 1500 further
local points per place are then sampled and must not enlarge it (a check of the quoted size from below).  Check: Cassels' formula
    dim Sel^phi(A) - dim Sel^phihat(B) = sum_{v in S} (log2 |B(K_v)/phi A(K_v)| - 1) - 1     (the -1: the complex place).
CONSEQUENCE.  For B with [2] = phi phihat (B -> A -> B):  dim Sel^2(B/K) <= dim Sel^phihat(B) + dim Sel^phi(A) - 1
(the -1: A(K)[phi] / phihat(B(K)[2]) has order 2 when B(K)[2] = ker phihat).  kctp_isog.py: under rank 0 (Sha finite)
the level-2 pairing on B = E3/<S_j> can prove rank 0 only if Sha(B)[2^inf] = (Z/2)^8, i.e. dim Sel^2(B/K) = 9.
"""
from sage.all__sagemath_schemes import *
import sys, os, random, math, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kctp_recheck as R

K = R.K; I = R.I; OK_ = R.OK_
WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(str(s))

def ainv_curve(a, b): return EllipticCurve(K, [0, a, 0, b, 0])

def sample_points(pl, a, b, rng, target, extra, maxit=60000, more=1500):
    """classes of x in K_v*/K_v*^2 for points (x, y) of y^2 = x (x^2 + a x + b) over K_v, until dimension target,
    then `more` further local points, which must not enlarge the span (a check that the quoted size is not too small)."""
    rows = [pl.cls(c) for c in extra if c != 0]
    def rk(rs): return R.f2_rank(rs) if rs else 0
    k = 0; after = 0
    while k < maxit:
        k += 1
        t = rng.randint(0, 3); z = K(rng.randint(-300, 300) + rng.randint(-300, 300) * I); j = rng.randint(0, 10)
        x = [z, pl.pi ** j * z, z / pl.pi ** (2 * j), pl.pi ** (2 * j) * z][t]
        if x == 0: continue
        fx = x * (x * x + a * x + b)
        if fx == 0 or not pl.is_square(fx): continue
        v = pl.cls(x)
        if rk(rows + [v]) > rk(rows): rows.append(v)
        if rk(rows) >= target:
            after += 1
            if after > more: break
    assert rk(rows) <= target, 'local image larger than the quoted size at %s' % pl.P
    return [list(r) for r in matrix(GF(2), rows).row_space().basis()] if rows else []

def selmer_isog(places, gens, imgs):
    """kernel of K(S,2) -> prod_v (K_v*/K_v*^2) / image_v."""
    rows = []
    for pl in places:
        W = imgs[pl.P]
        Wm = matrix(GF(2), W) if W else matrix(GF(2), 0, pl.w)
        ann = Wm.right_kernel().basis() if W else matrix.identity(GF(2), pl.w).rows()
        res = [pl.cls(g) for g in gens]
        for h in ann:
            rows.append([sum(int(h[t]) * res[j][t] for t in range(pl.w)) % 2 for j in range(len(gens))])
    return matrix(GF(2), rows).right_kernel().dimension()

def tam(F): return {P: F.local_data(P).tamagawa_number() for P in F.conductor().prime_factors()}

def descent(a, b, rng, label, forge=None):
    """phi: A -> B with A: y^2 = x(x^2 + a x + b); returns (dim Sel^phi(A), dim Sel^phihat(B))."""
    A = ainv_curve(a, b); c2 = a * a - 4 * b; B = ainv_curve(-2 * a, c2)
    phi = A.isogeny(A([0, 0]))
    assert phi.codomain().is_isomorphic(B)
    Am = A.global_minimal_model(); Bm = B.global_minimal_model()
    # alpha: from the normalised Velu isogeny on the minimal model of A and the isomorphism to the minimal model of B
    T0 = [P for P in Am(0).division_points(2) if P != 0 and Am.isogeny(P).codomain().is_isomorphic(B)]
    assert len(T0) == 1
    ps = Am.isogeny(T0[0]); C = ps.codomain(); u = C.isomorphism_to(Bm).u
    P2 = K.ideal(1 + I)
    va = K(u).valuation(P2); vah = 2 - va          # alpha = u, alpha alpha^ = 2
    D = 2 * b * c2
    Ps = sorted(set(P for P, _ in K.ideal(D).factor()), key=lambda P: (P.norm(), str(P.gens_reduced())))
    places = [R.Place(P) for P in Ps]
    gens = [K(I)] + [pl.pi for pl in places]
    tA = tam(Am); tB = tam(Bm)
    imgs_phi = {}; imgs_hat = {}; tsum = 0; tgt = []
    for pl in places:
        cA = next((c for P, c in tA.items() if P == pl.P), 1); cB = next((c for P, c in tB.items() if P == pl.P), 1)
        t_phi = 1 + int(round(math.log2(cB / cA))) + (va if pl.two else 0)        # log2 |B(K_v)/phi A(K_v)|
        t_hat = 1 + int(round(math.log2(cA / cB))) + (vah if pl.two else 0)       # log2 |A(K_v)/phihat B(K_v)|
        assert 2 ** t_phi * cA == 2 * cB * 2 ** (va if pl.two else 0)
        if forge == 'size' and not pl.two and not any(f[0] != 2 for f in tgt): t_phi += 1   # a forged local size (first odd place)
        # images: phi-descent map on B(K_v): X -> X, (0,0) -> c2;  phihat on A(K_v): x -> x, (0,0) -> b
        imgs_phi[pl.P] = sample_points(pl, -2 * a, c2, rng, t_phi, [c2])
        imgs_hat[pl.P] = sample_points(pl, a, b, rng, t_hat, [b])
        assert len(imgs_phi[pl.P]) == t_phi, '%s: phi image %d, expected %d at %s' % (label, len(imgs_phi[pl.P]), t_phi, pl.P)
        assert len(imgs_hat[pl.P]) == t_hat, '%s: phihat image %d, expected %d at %s' % (label, len(imgs_hat[pl.P]), t_hat, pl.P)
        tsum += t_phi - 1; tgt.append((ZZ(pl.P.norm()), t_phi, t_hat))
    s_phi = selmer_isog(places, gens, imgs_phi)
    s_hat = selmer_isog(places, gens, imgs_hat)
    if forge == 'selmer': s_hat += 1                                                          # a forged Selmer dimension
    cas = (s_phi - s_hat == tsum - 1)
    # exact B(B)/B(A): area ratio N(u)/2 (normalised Velu isogeny: lattice index 2; omega_min = u omega), checked in
    # floating point; Tamagawa products and torsion orders from Sage
    emb = K.embeddings(ComplexField(100))[0]
    ar = lambda F: abs((lambda v: (v[0].conjugate() * v[1]).imag())(F.period_lattice(emb).basis()))
    area_q = QQ(K(u).norm()) / 2
    assert abs(float(ar(Bm) / ar(Am)) / float(area_q) - 1) < 1e-12
    Bratio = area_q * QQ(prod(tB.values())) / QQ(prod(tA.values())) * QQ(Am.torsion_order() ** 2) / QQ(Bm.torsion_order() ** 2)
    say('  %s: v_(1+i)(alpha) = %d; local sizes (N(v), log2 |B/phi A|, log2 |A/phihat B|): %s' % (label, va, tgt))
    say('     dim Sel^phi(A) = %d, dim Sel^phihat(B) = %d; Cassels\' formula (%d - %d = %d): %s' % (
        s_phi, s_hat, s_phi, s_hat, tsum - 1, cas))
    say('     B(B)/B(A) = %s (exact; area ratio N(u)/2 = %s checked in floating point)' % (Bratio, area_q))
    assert cas
    return s_phi, s_hat, Bratio

def main():
    rng = random.Random(20261128)
    say('kctp_isodescent.py: 2-isogeny descents over K = Q(i) in the isogeny class of L1L3M1M2 of (627,364).')
    say('Written by an AI model; not read by a person.  Experiment, not counted.')
    Q, fs = R.quartic(627, 364, 'L1L3M1M2'); es = R.full2_model(Q, fs); E = R.curve_from_e(es)
    say('Controls: three base changes y^2 = x(x + p)(x + q) (kernel (0,0): a = p + q, b = p q); the local sizes must be')
    say('reached and not exceeded, and Cassels\' formula must hold.  (For (14,74): rank_K 0, Sha(A)[2] = (Z/2)^2 with')
    say('non-zero pairing (kctp_recheck --controls), B(B)/B(A) = 4, so Sha(B)[2^inf] = 0 and dim Sel^2(B) = 1, inside')
    say('the interval [dim Sel^phihat - 1, dim Sel^phihat + dim Sel^phi - 1] printed for it.)')
    for p, q in [(14, 74), (41, 80), (1, 7)]:
        sp, sh, Br = descent(K(p + q), K(p * q), rng, 'control y^2 = x(x+%d)(x+%d)' % (p, q))
        say('     dim Sel^2(B) in [%d, %d]' % (sh - 1, sh + sp - 1))
    say()
    say('The target.  E: y^2 = (x - e1)(x - e2)(x - e3); for the kernel (e_k, 0), x -> x + e_k gives a = 2 e_k - e_i - e_j,')
    say('b = (e_k - e_i)(e_k - e_j).')
    res = {}
    for k in range(3):
        i_, j_ = [t for t in range(3) if t != k]
        a = 2 * es[k] - es[i_] - es[j_]; b = (es[k] - es[i_]) * (es[k] - es[j_])
        res['E/<T%d>' % (k + 1)] = descent(a, b, rng, 'E -> E/<T%d>' % (k + 1))
    E3 = E.isogeny(E([es[2], 0])).codomain()
    Es = E3.short_weierstrass_model(); a4, a6 = Es.a4(), Es.a6()
    Rx = PolynomialRing(K, 'x'); x = Rx.gen()
    rts = [r for r, mlt in (x ** 3 + a4 * x + a6).roots()]; assert len(rts) == 3
    den = lcm([r.denominator() for r in rts]); e3s = [r * den ** 2 for r in rts]
    E3m = R.curve_from_e(e3s); assert E3m.is_isomorphic(E3)
    nS = 0
    for k in range(3):
        i_, j_ = [t for t in range(3) if t != k]
        a = 2 * e3s[k] - e3s[i_] - e3s[j_]; b = (e3s[k] - e3s[i_]) * (e3s[k] - e3s[j_])
        Bc = ainv_curve(-2 * a, a * a - 4 * b)
        if Bc.is_isomorphic(E): continue                      # the dual of E -> E3
        nS += 1
        res['E3/<S%d>' % nS] = descent(a, b, rng, 'E3 -> E3/<S%d>' % nS)
    say()
    rE3 = res['E/<T3>'][2]; assert rE3 == 4
    # cross-check with kctp_recheck's direct Sel^2 of E and E3 (both with all 2-torsion over K):
    # E -phi-> E3 -phihat-> E:  Sel^phi(E) -> Sel^2(E) -> Sel^phihat(E3), kernel E3(K)[phihat]/phi(E(K)[2]) = 0;
    # E3 -phihat-> E -phi-> E3: Sel^phihat(E3) -> Sel^2(E3) -> Sel^phi(E), kernel E(K)[phi]/phihat(E3(K)[2]) = 0.
    sp, sh, _ = res['E/<T3>']
    dE = R.Selmer(es, rng).dim; dE3 = R.Selmer(e3s, rng).dim
    okx = (sp <= dE <= sp + sh) and (sh <= dE3 <= sp + sh)
    say('Cross-check with kctp_recheck (direct 2-descent): dim Sel^2(E) = %d in [%d, %d], dim Sel^2(E3) = %d in [%d, %d]: %s' % (
        dE, sp, sp + sh, dE3, sh, sp + sh, okx))
    assert okx
    say('Exact dim Sel^2 of the four curves with one point of order 2 over K: Sel^phihat(B) -> Sel^2(B) -> Sel^phi(A) has')
    say('kernel of dimension 1, so dim Sel^phihat(B) - 1 <= dim Sel^2(B) <= dim Sel^phihat(B) + dim Sel^phi(A) - 1.')
    summary = {}
    for name, (sp, sh, Br) in res.items():
        if name == 'E/<T3>': continue
        lo, hi = sh - 1, sh + sp - 1
        ratio_to_E = Br if name.startswith('E/') else Br * rE3          # B(E')/B(E)
        k = -ZZ(QQ(ratio_to_E).numerator()).valuation(2) + ZZ(QQ(ratio_to_E).denominator()).valuation(2)
        summary[name] = (lo, hi, k)
        say('  %-9s dim Sel^2 in [%d, %d]%s;  |Sha(B)| / |Sha(E)| = B(E)/B(B) = 2^%d (exact)' % (
            name, lo, hi, ' (= %d)' % lo if lo == hi else '', k))
    say()
    say('CONSEQUENCES, if rank E(K) = 0 and Sha is finite (kctp_isog.py: Sha(E)[2^inf] = (Z/2^n)^2 with n >= 3).')
    say('For each B: dim Sha(B)[2] = dim Sel^2(B) - 1 = 2 m (m pairs of cyclic factors, Sha(B)[2^inf] = (prod_j Z/2^a_j)^2,')
    say('a_j >= 1, sum a_j = n + k/2).  The level-2 pairing on Sel^2(B) proves rank 0 iff all a_j = 1; a 4-descent with its')
    say('pairing iff all a_j <= 2.')
    for name, (lo, hi, k) in summary.items():
        assert lo == hi and (lo - 1) % 2 == 0
        m = (lo - 1) // 2
        # level 2 needs sum a_j = m, i.e. n = m - k/2; level 4 needs n + k/2 <= 2 m
        n2 = m - k // 2; n4 = 2 * m - k // 2
        say('  %-9s m = %d, sum a_j = n + %d:  level 2 needs n = %d (%s);  level 4 needs n <= %d (%s)' % (
            name, m, k // 2, n2, 'excluded' if n2 < 3 else 'possible', n4, 'excluded' if n4 < 3 else 'possible for n <= %d' % n4))

def negative():
    rng = random.Random(5)
    say('kctp_isodescent.py --negative: forged inputs that must be refused (on the control y^2 = x(x + 1)(x + 7)).')
    n = 0
    for forge, what in [('size', 'a local image size one larger than the true one'),
                        ('selmer', 'a Selmer dimension one too large (Cassels\' formula)')]:
        try:
            descent(K(8), K(7), rng, 'forged (%s)' % forge, forge=forge); say('  NOT REFUSED: %s' % what)
        except AssertionError as ex:
            n += 1; say('  refused: %s (%s)' % (what, str(ex)[:60]))
    say('%d of 2 refused.' % n); assert n == 2

if __name__ == '__main__':
    if '--negative' in sys.argv:
        negative(); fn = 'kctp_isodescent_negative.txt'
    else:
        main(); fn = 'kctp_isodescent.txt'
    if WRITE:
        open(os.path.join(os.path.dirname(os.path.abspath(__file__)), fn), 'w').write('\n'.join(OUT) + '\n')
