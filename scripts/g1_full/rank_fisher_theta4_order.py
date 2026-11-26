#!/usr/bin/env python3
"""rank_fisher_theta4_order.py: the known point G of E+ of (72,65) has infinite order, in exact arithmetic, with no
PARI, no Sage, no floating point and no height bound.  Written by an AI model; not read by a person.

    python3 rank_fisher_theta4_order.py [--write]              -> rank_fisher_theta4_order.txt
    python3 rank_fisher_theta4_order.py --negative [--write]   -> rank_fisher_theta4_order_negative.txt

rank_fisher_theta4.py concludes "rank E(Q) <= 1 + 2 - 2 = 1; with G of infinite order, rank E(Q) = 1".  Until now the
clause "G of infinite order" rested on sigma_exact.py (item 2: a positive lower bound for the canonical height, in a
script whose index bound uses a recorded floating-point lambda).  This script replaces that clause by three exact
certificates, each sufficient on its own:
  (1) reduction bound.  For odd p not dividing Delta, reduction is injective on E(Q)_tors (Silverman, AEC VII.3.1 and
      IV.6.4; quoted), so #E(Q)_tors divides B = gcd #E~(F_p) over the good odd p < 1000.  If G had finite order d,
      then d | B and B G = O.  Here k G != O for every k <= B (exact group law), in particular B G != O.
  (2) orders modulo p.  If G had finite order d, its reduction would have order exactly d modulo every good odd p (an
      injective homomorphism on the finite group <G>).  Here the order of G~ in E~(F_p) differs between two good odd
      primes (orders computed exactly in F_p).
  (3) Nagell-Lutz.  On y^2 = x^3 + a x^2 + b x with a, b integers, a point of finite order has integer coordinates
      (Nagell-Lutz in the form of Silverman-Tate, Rational Points on Elliptic Curves, Ch. II, for y^2 = x^3 + ax^2 +
      bx + c; quoted).  Here G and 2G are integral but 3G has a non-integral x-coordinate, so 3G, hence G, has infinite
      order.
The curve is rebuilt from (m, n) = (72, 65) by sigma_exact.model (pure Python), moved to y^2 = x(x^2 + a x + b) by its
unique rational root (counted exactly: the quadratic cofactor has non-square discriminant), and compared with the
a, b printed in rank_fisher_theta4.txt; G is the generator recorded in sigma_evidence.txt (the same file
rank_fisher_theta4.py reads through rank_indep.evidence), moved by the same shift, and checked to lie on the curve.
E(Q)_tors itself is certified by torsion_exact.py (printed; not needed for (1)-(3)).
--negative: G replaced by each torsion point, by a point off the curve, by points of order 4, 8 and 12 on curves of
known torsion (15a8, 15a4, 90c3), and the prime list forged with a bad prime and with p = 2: each of the three
certificates must refuse each.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sigma_exact as SE
import torsion_exact as TE

M, N = 72, 65

def curve_and_G():
    Mo = SE.model(M, N, 2)                                   # E+: s = 2
    E0 = TE.Curve(Mo['a2'], Mo['a4'], Mo['a6'])
    T2 = TE.two_torsion(E0)
    ok, why = TE.two_torsion_exact(E0, T2, TE.good_odd_primes(E0, 1000))
    assert ok and len(T2) == 2, ('expected exactly one rational root of the cubic', why)
    e = T2[1][0]; assert E0.f(e) == 0
    a, b = 3*e + Mo['a2'], 3*e*e + 2*Mo['a2']*e + Mo['a4']
    assert a.denominator == 1 and b.denominator == 1
    E = TE.Curve(int(a), int(b), 0)
    ev = None
    for line in open(os.path.join(HERE, 'sigma_evidence.txt')):
        w = line.split()
        if len(w) > 2 and w[0] == '%d,%d' % (M, N) and w[1] == 'E+':
            f = dict(x.split('=', 1) for x in w[2:]); ev = tuple(Fr(v) for v in f['G'].split(':'))
    assert ev is not None and E0.on(ev), 'G of sigma_evidence.txt is not on the model'
    G = (ev[0] - e, ev[1]); assert E.on(G)
    head = open(os.path.join(HERE, 'rank_fisher_theta4.txt')).read().split('\n')[1]
    assert head.startswith('E: y^2 = x(x^2 + a x + b), a = %d, b = %d; known point G with x(G) = %s' % (E.a2, E.a4, G[0])), head
    return E, G, e, why

def order_mod(E, P, p, cap):
    Pp = TE.reduce_pt(P, p)
    if Pp is None: return 1
    Q = Pp
    for k in range(1, cap + 1):
        if Q is None: return k
        Q = TE.add_mod(E, Q, Pp, p)
    return None

def certificates(E, G, primes):
    """the three certificates for 'G has infinite order'; each returns (holds, text).  A forged prime list is refused."""
    out = {}
    if not E.on(G): return {k: (False, 'G is not on the curve') for k in ('bound', 'orders', 'nagell-lutz')}
    try: B = TE.bound(E, primes)
    except AssertionError as ex: return {k: (False, 'bound refused: %s' % (ex.args[0],)) for k in ('bound', 'orders', 'nagell-lutz')}
    Q = G; first = None
    for k in range(1, B + 1):
        if Q is None: first = k; break
        Q = E.add(Q, G)
    out['bound'] = (first is None, 'B = gcd #E~(F_p) over %d good odd primes p < 1000 = %d; k G %s' % (
        len(primes), B, '!= O for every k <= B' if first is None else '= O for k = %d' % first))
    ords = {}
    for p in primes[:12]: ords[p] = order_mod(E, G, p, p + 1 + 2 * int(p ** 0.5) + 2)
    diff = len(set(ords.values())) > 1
    out['orders'] = (diff, 'order of G~ in E~(F_p) for p = %s: %s (%s)' % (
        ', '.join(str(p) for p in ords), ', '.join(str(v) for v in ords.values()), 'not all equal' if diff else 'ALL EQUAL'))
    found = None; Q = G
    for k in range(1, 7):                                    # any multiple of a torsion point is torsion, so integral
        if Q is None: break
        if Fr(Q[0]).denominator != 1 or Fr(Q[1]).denominator != 1: found = (k, Q); break
        Q = E.add(Q, G)
    out['nagell-lutz'] = (found is not None, ('%d G has x = %s/%s (denominator %d != 1: not integral)' % (found[0], Fr(found[1][0]).numerator,
                          Fr(found[1][0]).denominator, Fr(found[1][0]).denominator)) if found else 'k G integral (or O) for every k <= 6')
    return out

def main(write):
    E, G, e, why = curve_and_G(); primes = TE.good_odd_primes(E, 1000)
    L = ['rank_fisher_theta4_order.py: G of infinite order on E+ of (72,65), exact (no PARI, no floating point, no heights)', '']
    L.append('E: y^2 = x(x^2 + a x + b), a = %d, b = %d (sigma_exact.model(72, 65, +2) moved by its rational root e = %s; %s;' % (E.a2, E.a4, e, why))
    L.append('   the same a, b as rank_fisher_theta4.txt)')
    L.append('G = (%s, %s) (sigma_evidence.txt, moved by e), on the curve' % G)
    T, c = TE.torsion(E, primes)
    assert c['ok']
    L.append('E(Q)_tors = %s = {%s} (torsion_exact.py: %s)' % (TE.fmt_structure(c['structure']), ', '.join(TE.show(P) for P in T), '; '.join(c['reasons'])))
    cs = certificates(E, G, primes)
    for k, lab in (('bound', '(1) reduction bound'), ('orders', '(2) orders modulo p'), ('nagell-lutz', '(3) Nagell-Lutz')):
        L.append('%-22s %s: %s' % (lab, cs[k][1], 'G of infinite order' if cs[k][0] else 'NOT SHOWN'))
        assert cs[k][0]
    L.append('')
    L.append('RESULT: G has infinite order (each of (1), (2), (3) suffices).  With rank E(Q) <= 1 from rank_fisher_theta4.txt')
    L.append('(Theta_4(y0, y1) != 0, Fisher Lemma 2.2), rank E+(72,65)(Q) = 1.')
    out = '\n'.join(L) + '\n'; sys.stdout.write(out)
    if write: open(os.path.join(HERE, 'rank_fisher_theta4_order.txt'), 'w').write(out)

def negative(write):
    E, G, e, why = curve_and_G(); primes = TE.good_odd_primes(E, 1000)
    T, c = TE.torsion(E, primes); assert c['ok']
    L = ['rank_fisher_theta4_order.py --negative: inputs that must be refused', '']
    ntot = nref = 0
    cases = [('G replaced by the torsion point %s' % TE.show(P), P, primes) for P in T if P is not None]
    cases += [('G replaced by (x(G), y(G) + 1), off the curve', (G[0], G[1] + 1), primes),
              ('primes with the bad prime %d' % next(p for p in range(3, 1000, 2) if E.D % p == 0),
               G, sorted(primes[:30] + [next(p for p in range(3, 1000, 2) if E.D % p == 0)])),
              ('primes with p = 2', G, [2] + primes[:30])]
    for lab, P, pr in cases:
        cs = certificates(E, P, pr)
        for k in ('bound', 'orders', 'nagell-lutz'):
            ntot += 1; nref += not cs[k][0]
            L.append('%-50s %-12s %s: %s' % (lab + ':', k, cs[k][1], 'NOT REFUSED' if cs[k][0] else 'refused'))
    for lab, (a2, a4, a6), want in TE.CONTROLS:              # points of order 4, 8 and 12 on curves of known torsion
        if lab not in ('15a8', '15a4', '90c3'): continue
        Ec = TE.Curve(a2, a4, a6); pc = TE.good_odd_primes(Ec, 1000); Tc, cc = TE.torsion(Ec, pc); assert cc['ok']
        n = max(Ec.order(P, 16) for P in Tc); P = [P for P in Tc if Ec.order(P, 16) == n][0]
        cs = certificates(Ec, P, pc)
        for k in ('bound', 'orders', 'nagell-lutz'):
            ntot += 1; nref += not cs[k][0]
            L.append('%-50s %-12s %s: %s' % ('a point of order %d on %s %s:' % (n, lab, Ec.model()), k, cs[k][1], 'NOT REFUSED' if cs[k][0] else 'refused'))
    cs = certificates(E, G, primes)
    L.append('positive control, the true G: %s' % ', '.join('%s %s' % (k, 'holds' if v[0] else 'FAILS') for k, v in cs.items()))
    assert all(v[0] for v in cs.values())
    L.append('negative tests: %d forged inputs, %d refused' % (ntot, nref))
    out = '\n'.join(L) + '\n'; sys.stdout.write(out)
    assert ntot == nref
    if write: open(os.path.join(HERE, 'rank_fisher_theta4_order_negative.txt'), 'w').write(out)

if __name__ == '__main__':
    if '--negative' in sys.argv: negative('--write' in sys.argv)
    else: main('--write' in sys.argv)
