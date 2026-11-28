#!/usr/bin/env python3
"""Claim (D) of the separate session's note (27 November 2026, section 3): the corner curve L1 L3 M1 M2 of (627,364)
over K = Q(i) is not a Q-curve.  Here the argument is made to rest on the j-invariant alone, in exact arithmetic.

Written by an AI model (build agent t12); not read by a person.  Experiment, tier D (the argument in NOTES.md is
tier C).  Nothing here changes a count.

    python3 -I kcurve_exact.py [--write]        kcurve_exact.txt
    python3 -I kcurve_exact.py --negative       forged inputs that must be refused

The curve: y^2 = Q(s) = L1 L3 M1 M2 = (m - n s)(n + m s)(p s - pb)(pb s + p), p = m + i n, pb = m - i n (labels of
scripts/second_descent/kclose.py).  j of its Jacobian from the invariants of the binary quartic,
j = 6912 I^3 / (4 I^3 - J^2) (I = 12ae - 3bd + c^2, J = 72ace + 9bcd - 27ad^2 - 27b^2e - 2c^3), in Q(i) with Fractions.
For a prime p = 1 mod 4 and pi = a + b i with a^2 + b^2 = p: v_pi(j) by exact division in Z[i].

The criterion (proof in NOTES.md): if some prime pi of K has v_pi(j) < 0 and v_pibar(j) >= 0, then E is not
isogenous over Kbar to its Galois conjugate E^c, so E is not a Q-curve.  (E has potentially multiplicative reduction
at pi; j(E^c) = conj(j(E)), so v_pi(j(E^c)) = v_pibar(j(E)) >= 0 and E^c has potentially good reduction at pi;
potentially good reduction is invariant under isogeny over any finite extension (Neron-Ogg-Shafarevich, quoted).)
It also gives non-CM (CM j-invariants are algebraic integers, quoted).  No conductor, no Tate's algorithm and no PARI
is used for the criterion.  PARI (cypari2) is used only for a cross-check: its j agrees with this one, and its
conductor (ellglobalred over nf, Tate's algorithm) has exponent 1 at the prime where v(j) < 0 and 0 at its conjugate.
"""
import sys, os
from fractions import Fraction as Fr
import sympy as sp

WRITE = '--write' in sys.argv; NEG = '--negative' in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

# Gaussian rationals as pairs (re, im) of Fractions
def gadd(u, v): return (u[0] + v[0], u[1] + v[1])
def gsub(u, v): return (u[0] - v[0], u[1] - v[1])
def gmul(u, v): return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])
def gconj(u): return (u[0], -u[1])
def ginv(u):
    nr = u[0] * u[0] + u[1] * u[1]; assert nr != 0
    return (u[0] / nr, -u[1] / nr)
def gdiv(u, v): return gmul(u, ginv(v))
def G(a, b=0): return (Fr(a), Fr(b))
def gsc(c, u): return (c * u[0], c * u[1])

def polymul(P, Q):
    R = [G(0)] * (len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        for j, b in enumerate(Q): R[i + j] = gadd(R[i + j], gmul(a, b))
    return R

def quartic(m, n):
    """coefficients [e, d, c, b, a] (low to high) of Q(s) = (m - n s)(n + m s)(p s - pb)(pb s + p)"""
    p = G(m, n); pb = G(m, -n)
    L1 = [G(m), G(-n)]; L3 = [G(n), G(m)]
    M1 = [gsc(-1, pb), p]; M2 = [p, pb]
    return polymul(polymul(L1, L3), polymul(M1, M2))

def jinv(Q):
    e, d, c, b, a = Q
    I = gadd(gsub(gsc(12, gmul(a, e)), gsc(3, gmul(b, d))), gmul(c, c))
    J = gsub(gsub(gsub(gadd(gsc(72, gmul(gmul(a, c), e)), gsc(9, gmul(gmul(b, c), d))), gsc(27, gmul(a, gmul(d, d)))),
                  gsc(27, gmul(e, gmul(b, b)))), gsc(2, gmul(c, gmul(c, c))))
    I3 = gmul(I, gmul(I, I))
    den = gsub(gsc(4, I3), gmul(J, J))
    assert den != G(0), 'singular quartic'
    return gdiv(gsc(6912, I3), den), I, J

def gauss_prime(p):
    """pi = a + b i with a^2 + b^2 = p (p prime, p = 1 mod 4); a > b > 0"""
    assert p % 4 == 1 and sp.isprime(p)
    g = 2
    while pow(g, (p - 1) // 2, p) != p - 1: g += 1
    x = pow(g, (p - 1) // 4, p)                        # x^2 = -1 mod p
    # Euclid in Z[i] on (p, x + i)
    def gdivmod(u, v):
        q = gdiv(u, v); q = (Fr(round(q[0])), Fr(round(q[1])))
        return q, gsub(u, gmul(q, v))
    u, v = G(p), G(x, 1)
    while v != G(0):
        _, r = gdivmod(u, v); u, v = v, r
    a, b = abs(int(u[0])), abs(int(u[1]))
    assert a * a + b * b == p
    return (a, b) if a > b else (b, a)

def vpi(alpha, pi, p):
    """valuation at pi of a nonzero Gaussian integer alpha (exact): alpha / pi = alpha * conj(pi) / p"""
    assert alpha != (0, 0)
    k = 0
    u = (int(alpha[0]), int(alpha[1])); a, b = pi
    while True:
        w = (u[0] * a + u[1] * b, u[1] * a - u[0] * b)     # u * conj(pi)
        if w[0] % p or w[1] % p: return k
        u = (w[0] // p, w[1] // p); k += 1

def vj(j, pi, p):
    from math import lcm
    w = lcm(j[0].denominator, j[1].denominator)
    alpha = (j[0] * w, j[1] * w)
    assert alpha[0].denominator == 1 and alpha[1].denominator == 1
    vw = 0
    while w % p == 0: w //= p; vw += 1
    return vpi(alpha, pi, p) - vw

def criterion(j, primes):
    rows = []; hit = []
    for p in primes:
        a, b = gauss_prime(p)
        v1 = vj(j, (a, b), p); v2 = vj(j, (a, -b), p)
        rows.append((p, (a, b), v1, v2))
        if (v1 < 0) != (v2 < 0) and min(v1, v2) < 0: hit.append(p)
    return rows, hit

def run():
    m, n = 627, 364
    say('kcurve_exact.py: the corner curve y^2 = L1 L3 M1 M2 of (%d,%d) over K = Q(i)' % (m, n))
    Q = quartic(m, n)
    say('  Q(s) = %s' % ' + '.join('(%s, %s)s^%d' % (c[0], c[1], k) for k, c in enumerate(Q)))
    say('         (coefficients as (real part, imaginary part), s^0 first)')
    j, I, J = jinv(Q)
    say('  j = 6912 I^3/(4I^3 - J^2) in Q(i): j is not rational (Im j != 0): %s' % (j[1] != 0))
    assert j[1] != 0
    jc = gconj(j)
    jQc, _, _ = jinv([gconj(c) for c in Q])
    say('  j of the conjugate quartic L1 L3 M3 M4 (= conj of Q) is conj(j): %s' % (jQc == jc))
    PR = [37, 181, 23490241, 395249137273]
    rows, hit = criterion(j, PR)
    for p, (a, b), v1, v2 in rows:
        say('  p = %d = %d^2 + %d^2: v_pi(j) = %d at pi = %d + %di, %d at its conjugate' % (p, a, b, v1, a, b, v2))
    say('  primes where exactly one of pi, conj(pi) has v(j) < 0: %s' % hit)
    assert hit, 'no prime satisfies the criterion'
    # every prime p = 1 mod 4 dividing the norm of the denominator of j: which ones are asymmetric
    from math import lcm
    w = lcm(j[0].denominator, j[1].denominator)
    alpha = (int(j[0] * w), int(j[1] * w))
    nden = sp.factorint(w)
    say('  the common denominator of Re j, Im j has the primes %s' % sorted(nden))
    allp = sorted(set(q for q in nden if q % 4 == 1))
    rows2, hit2 = criterion(j, allp)
    for p, (a, b), v1, v2 in rows2:
        say('    p = %d: v(j) = %d, %d at the two primes above p' % (p, v1, v2))
    say('  asymmetric split primes among them: %s' % hit2)
    say('  => E is non-CM and not isogenous over Kbar to its conjugate: not a Q-curve (NOTES.md, section D)')
    return j, Q

def pari_crosscheck(j, Q):
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(2 * 10**9)
    pari('nf = nfinit(w^2 + 1)')
    pol = ' + '.join('(%s + (%s)*w)*x^%d' % (c[0], c[1], k) for k, c in enumerate(Q))
    pari('E = ellinit(ellfromeqn(y^2 - (%s)), nf)' % pol)
    jp = pari('lift(E.j)')
    jm = pari('(%s) + (%s)*w' % (j[0], j[1]))
    say('PARI cross-check (a guide): PARI\'s j of ellfromeqn equals the j above: %s' % bool(pari('%s == %s' % (jp, jm))))
    G = pari('G = ellglobalred(E); G[1]')
    for p in (37, 181, 23490241, 395249137273):
        dec = pari('idealprimedec(nf, %d)' % p)
        ex = [int(pari('idealval(nf, G[1], idealprimedec(nf, %d)[%d])' % (p, k + 1))) for k in range(len(dec))]
        vjp = [int(pari('nfeltval(nf, lift(E.j), idealprimedec(nf, %d)[%d])' % (p, k + 1))) for k in range(len(dec))]
        say('  p = %d: conductor exponents at the two primes above p %s, v(j) there %s' % (p, ex, vjp))
    nm = pari('idealnorm(nf, G[1])')
    say('  conductor norm %s = %s' % (nm, pari('factor(%s)' % nm)))

def negative():
    res = []
    def t(name, fn):
        try: fn(); res.append(False); say('  NOT refused: ' + name)
        except AssertionError as ex: res.append(True); say('  refused: %s (%s)' % (name, ex.args[0] if ex.args else ''))
    m, n = 627, 364
    Q = quartic(m, n); j, _, _ = jinv(Q)
    # a base change from Q: y^2 = (m - n s)(n + m s)(s - 2)(s + 3) has rational j: symmetric valuations everywhere
    def bc():
        Q2 = polymul(polymul([G(m), G(-n)], [G(n), G(m)]), polymul([G(-2), G(1)], [G(3), G(1)]))
        j2, _, _ = jinv(Q2)
        w = j2[0].denominator
        ps = sorted(q for q in sp.factorint(w) if q % 4 == 1)
        _, hit = criterion(j2, ps)
        assert hit, 'a curve over Q (base change) shows no asymmetric prime'
    t('a base change from Q (rational j) must show no asymmetric prime', bc)
    def conjsym():
        # j * conj(j) is rational: symmetric at every prime
        jj = gmul(j, gconj(j))
        _, hit = criterion(jj, [37, 181, 23490241, 395249137273])
        assert hit, 'N(j) (rational) shows no asymmetric prime'
    t('the norm j conj(j) (rational) must show no asymmetric prime', conjsym)
    t('a composite p (37 * 41 = 1517) refused by gauss_prime', lambda: gauss_prime(1517))
    t('p = 3 mod 4 refused by gauss_prime', lambda: gauss_prime(43))
    def singular():
        Q3 = polymul(polymul([G(m), G(-n)], [G(m), G(-n)]), polymul([G(1), G(1)], [G(2), G(1)]))
        jinv(Q3)
    t('a quartic with a double root (singular)', singular)
    say('negative tests: %d of %d refused' % (sum(res), len(res)))
    return all(res)

if __name__ == '__main__':
    if NEG:
        say('kcurve_exact.py --negative')
        ok = negative()
        if WRITE: open(os.path.join(HERE, 'kcurve_exact_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
        sys.exit(0 if ok else 1)
    j, Q = run()
    if '--nopari' not in sys.argv: pari_crosscheck(j, Q)
    if WRITE: open(os.path.join(HERE, 'kcurve_exact.txt'), 'w').write('\n'.join(OUT) + '\n')
