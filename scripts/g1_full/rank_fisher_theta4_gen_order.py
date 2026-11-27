#!/usr/bin/env python3
"""G has infinite order on the curves of rank_fisher_theta4_gen.py, by exact integer arithmetic only: if G were a
torsion point of order k, its reduction modulo every prime p of good reduction with p > 2 would have order exactly k
(reduction is injective on torsion there), so two such primes at which G mod p has different orders prove that G has
infinite order.  The orders are computed by repeated addition on y^2 = x(x^2 + a x + b) over F_p (no PARI).

    python3 rank_fisher_theta4_gen_order.py [--write] [--negative]

--negative: a torsion point (T = (0, 0), order 2) must give the same order at every prime (so it is not certified).
Written by an AI model (Claude); not read by a person.  Experiment, not counted.
"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from fractions import Fraction as Fr
import rank_second as R2

WRITE = '--write' in sys.argv
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

def add(P, Q, a, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    if P == Q: lam = (3 * x1 * x1 + 2 * a * x1 + B_[0]) * pow(2 * y1, -1, p) % p
    else: lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - a - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

B_ = [0]
def order_mod(G, a, b, p):
    B_[0] = b % p
    P0 = (int(G[0].numerator * pow(G[0].denominator, -1, p) % p), int(G[1].numerator * pow(G[1].denominator, -1, p) % p))
    assert (P0[1] ** 2 - P0[0] * (P0[0] ** 2 + a * P0[0] + b)) % p == 0
    P, k = P0, 1
    while P is not None:
        P = add(P, P0, a % p, p); k += 1
        assert k <= p + 1 + 2 * int(p ** 0.5) + 2
    return k

def good(p, a, b, G):
    return p > 2 and b % p and (a * a - 4 * b) % p and G[0].denominator % p and G[1].denominator % p

def certify(label, a, b, G, need=2):
    orders = {}
    for p in [q for q in range(3, 400) if all(q % r for r in range(2, int(q ** 0.5) + 1))]:
        if not good(p, a, b, G): continue
        orders[p] = order_mod(G, a, b, p)
        if len(set(orders.values())) >= need: break
    ok = len(set(orders.values())) >= 2
    say('%s: orders of G modulo the good primes tried: %s -> %s' % (label, ', '.join('%d: %d' % kv for kv in sorted(orders.items())),
        'two different orders: G has infinite order' if ok else 'no two different orders'))
    return ok

def main():
    say('rank_fisher_theta4_gen_order.py: G has infinite order (reductions modulo good primes, integer arithmetic)')
    curves = [('920,129 E+', R2.sigma_curve(920, 129, 'E+')), ('760,39 E+', R2.sigma_curve(760, 39, 'E+')),
              ('631,71 D1', R2.d1_curve(631, 71))]
    for label, (a, b, G) in curves:
        G = (Fr(G[0]), Fr(G[1]))
        assert G[1] ** 2 == G[0] * (G[0] ** 2 + a * G[0] + b), 'G is not on the curve'
        assert certify(label, a, b, G)

def negative():
    say('rank_fisher_theta4_gen_order.py --negative: a torsion point must not be certified')
    a, b, G = R2.d1_curve(631, 71)
    T = (Fr(0), Fr(0))
    ok = certify('631,71 D1, the point T = (0, 0) of order 2', a, b, T, need=2)
    assert not ok, 'a torsion point was certified'
    say('negative test: refused (all orders equal)')

if __name__ == '__main__':
    if '--negative' in sys.argv:
        negative()
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_order_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
    else:
        main()
        if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_order.txt'), 'w').write('\n'.join(OUT) + '\n')
