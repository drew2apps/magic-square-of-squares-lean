#!/usr/bin/env python3
"""Claim (C) of the separate session's note (27 November 2026, section 2): the elliptic quotients of (627,364) fall into
15 isogeny classes (isogeny_classes.py groups them by equal Frobenius traces at 40 primes).  Here every grouping is
made exact: inside a group by an explicit isomorphism or 2-isogeny (sympy, for general m, n); across groups by one
prime where a_p differs (pure Python counts; isogenous curves over Q have the same a_p).  Then the multiplicity of the
classes in the Jacobian of the curve of the five forms: see genus3_mult.py.

Written by an AI model (build agent t12); not read by a person.  Experiment, tier D.  Nothing here changes a count.

    python3 -I genus3_isog.py [m n] [--write]      genus3_isog_627_364.txt (default ratio 627 364)

Facts quoted, not proved here: a nonconstant morphism of genus-1 curves induces an isogeny of their Jacobians of the
same degree; isogenous curves over Q have equal a_p at good p.  All counting is exact.
"""
import sys, os, itertools
import sympy as sp

WRITE = '--write' in sys.argv
args = [a for a in sys.argv[1:] if not a.startswith('--')]
MM, NN = (int(args[0]), int(args[1])) if len(args) >= 2 else (627, 364)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

m, n, x = sp.symbols('m n x')
S = m**2 + n**2
FORMS = {'F0': m**2 + n**2 * x**2, 'F8': n**2 + m**2 * x**2, 'F4': 2 * S * (x**2 + 1),
         'L': (m + n)**2 * (x - 1)**2 + (m - n)**2 * (x + 1)**2, 'U': (m - n)**2 * (x - 1)**2 + (m + n)**2 * (x + 1)**2}
A, B, C, D = m**2 * n**2, m**4 + n**4, S**2, 2 * S**2 - 16 * m**2 * n**2
LUx = 4 * S**2 * (x + 1)**2 - 64 * m**2 * n**2 * x
CURVES = {}
for a, b in itertools.combinations(FORMS, 2): CURVES[a + '*' + b] = FORMS[a] * FORMS[b]
CURVES['E1'] = (x**2 - 4) * (A * x + B) * (C * x + D)
CURVES['E2'] = (x + 2) * (A * x + B) * (C * x + D)
CURVES['E3'] = (x - 2) * (A * x + B) * (C * x + D)
def Epm(nn):
    return (x + 2 * S**3) * (x**2 + 16 * m**2 * nn**2 * (m + nn)**2 * x + 16 * m**2 * nn**2 * (m + nn)**4 * S**2)
CURVES['E+'] = Epm(n); CURVES['E-'] = Epm(-n)
CURVES['[F4F0LU|x=t^2]'] = 2 * S * (x + 1) * (m**2 + n**2 * x) * LUx
CURVES['[F4F8LU|x=t^2]'] = 2 * S * (x + 1) * (n**2 + m**2 * x) * LUx
CURVES['[F0F8LU|x=t^2]'] = (m**2 + n**2 * x) * (n**2 + m**2 * x) * LUx
CURVES['[F4F0F8|x(..)]'] = 2 * S * x * (x + 1) * (m**2 + n**2 * x) * (n**2 + m**2 * x)
CURVES['[F4LU|x(..)]'] = 2 * S * x * (x + 1) * LUx
CURVES['[F4LU|x=t^2]q'] = 2 * S * (x + 1) * LUx
CURVES['[F0LU|x(..)]'] = x * (m**2 + n**2 * x) * LUx
CURVES['[F8LU|x(..)]'] = x * (n**2 + m**2 * x) * LUx
GROUPS = [['F0*F8'], ['F0*F4', 'F8*F4'], ['F0*L', 'F0*U', 'F8*L', 'F8*U'], ['F4*L', 'F4*U'], ['L*U'],
          ['E1', '[F0F8LU|x=t^2]'], ['E2'], ['E3'], ['E+'], ['E-'], ['[F4F0LU|x=t^2]', '[F4F8LU|x=t^2]'],
          ['[F4F0F8|x(..)]'], ['[F4LU|x(..)]', '[F4LU|x=t^2]q'], ['[F0LU|x(..)]'], ['[F8LU|x(..)]']]
assert sorted(sum(GROUPS, [])) == sorted(CURVES)

say('genus3_isog.py: the elliptic quotients of the note\'s isogeny_classes.py at (%d,%d), exactly' % (MM, NN))
say()
# --- inside the groups: explicit maps, generic in m, n --------------------------------------------------------
def mobius_iso(q1, q2, mob, lam=1):
    """q2(x) == lam^2 (c x + d)^4 q1((a x + b)/(c x + d)): (x, y) -> ((ax+b)/(cx+d), y lam^-1 (cx+d)^-2)... an
    isomorphism of the genus-1 curves y^2 = q1 and y^2 = q2 (binary quartic forms of degree <= 4)"""
    a, b, c, d = mob
    e = sp.expand(lam**2 * (c * x + d)**4 * q1.subs(x, (a * x + b) / (c * x + d)))
    return sp.cancel(e - sp.expand(q2)) == 0
ISO = [('F0*F4', 'F8*F4', (0, 1, 1, 0), 'x -> 1/x'), ('F0*L', 'F0*U', (-1, 0, 0, 1), 'x -> -x'),
       ('F0*L', 'F8*L', (0, 1, 1, 0), 'x -> 1/x'), ('F0*L', 'F8*U', (0, -1, 1, 0), 'x -> -1/x'),
       ('F4*L', 'F4*U', (-1, 0, 0, 1), 'x -> -x'), ('[F4F0LU|x=t^2]', '[F4F8LU|x=t^2]', (0, 1, 1, 0), 'x -> 1/x'),
       ('[F4LU|x(..)]', '[F4LU|x=t^2]q', (0, 1, 1, 0), 'x -> 1/x')]
say('Inside the groups (sympy, for general m, n):')
for a, b, mob, txt in ISO:
    ok = mobius_iso(CURVES[a], CURVES[b], mob)
    say('  %-16s ~ %-16s isomorphic by %s, (x, y) -> (x\', y/(cx+d)^2): %s' % (a, b, txt, ok)); assert ok
# the 2-isogeny E_sigma = [F0F8LU|x=t^2] -> E1: (x, y) -> (z, w) = (x + 1/x, y (x - 1/x) / (2x))
ysym, zsym = sp.symbols('y z')
h = CURVES['[F0F8LU|x=t^2]']
zz = x + 1 / x; ww = ysym * (x - 1 / x) / (2 * x)
e = sp.expand(ww**2 - CURVES['E1'].subs(x, zz)).subs(ysym**2, h)
ok = sp.cancel(sp.together(e)) == 0
say('  [F0F8LU|x=t^2] -> E1, (x, y) -> (x + 1/x, y (x - 1/x)/(2x)), a morphism of degree 2 (z = x + 1/x has degree 2'
    ' in x, and on E1 degree 2; the map is the quotient by (x, y) -> (1/x, -y/x^2)): %s' % ok); assert ok
say('  so the two curves of each group have isogenous Jacobians (isomorphic, or 2-isogenous for E1).')
say()

# --- across the groups: a_p ------------------------------------------------------------------------------------
num = lambda e: sp.Poly(sp.expand(e.subs({m: MM, n: NN})), x)
P = {k: [int(c) for c in num(v).all_coeffs()][::-1] for k, v in CURVES.items()}
def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
def ev(cs, t, p):
    r = 0
    for c in reversed(cs): r = (r * t + c) % p
    return r
def good_for(cs, p):
    if cs[-1] % p == 0: return False
    q = sp.Poly(list(reversed([c % p for c in cs])), x, modulus=p)
    return sp.degree(sp.gcd(q, q.diff(x))) == 0
def ap(cs, p):
    N = sum(1 + chi(ev(cs, t, p), p) for t in range(p))
    N += 1 if (len(cs) - 1) % 2 else 1 + chi(cs[-1], p)
    return p + 1 - N
primes = [int(p) for p in sp.primerange(5, 1000)]
APS = {}
def ap_c(k, p):
    if (k, p) not in APS: APS[(k, p)] = ap(P[k], p) if good_for(P[k], p) else None
    return APS[(k, p)]
say('Across the groups: for each pair of groups the first prime p >= 5 at which both are good and a_p differ')
worst = 0; lines = []
for i, j in itertools.combinations(range(len(GROUPS)), 2):
    a, b = GROUPS[i][0], GROUPS[j][0]
    p = next(p for p in primes if ap_c(a, p) is not None and ap_c(b, p) is not None and ap_c(a, p) != ap_c(b, p))
    worst = max(worst, p); lines.append((a, b, p))
for a, b, p in lines:
    say('  %-16s %-16s p = %3d: a_p = %3d, %3d' % (a, b, p, ap_c(a, p), ap_c(b, p)))
say('  all %d pairs of the 15 groups are non-isogenous; the largest prime needed is %d' % (len(lines), worst))
# consistency with the groups themselves (equal a_p inside a group is a consequence of the maps; checked anyway)
for g in GROUPS:
    for k in g[1:]:
        for p in primes[:60]:
            if ap_c(g[0], p) is not None and ap_c(k, p) is not None: assert ap_c(g[0], p) == ap_c(k, p)
say('  (inside each group a_p agree at the good primes < 300, as they must)')
say()

# --- twists ---------------------------------------------------------------------------------------------------
def IJ(cs):
    cs = cs + [0] * (5 - len(cs))
    e_, d_, c_, b_, a_ = [sp.Integer(c) for c in cs]
    I = 12 * a_ * e_ - 3 * b_ * d_ + c_**2
    J = 72 * a_ * c_ * e_ + 9 * b_ * c_ * d_ - 27 * a_ * d_**2 - 27 * e_ * b_**2 - 2 * c_**3
    return I, J
say('Quadratic twists (the note\'s "twist pairs": |a_p| equal at 40 primes).  The Jacobian of y^2 = quartic with')
say('invariants I, J is Y^2 = X^3 - 27 I X - 27 J; its twist by d is (d^2 I, d^3 J).  Exact search: from the twist of')
say('the first curve by d, all curves reached by rational 2-isogenies (Velu: y^2 = x^3 + a x + b, T = (x0, 0),')
say('t = 3x0^2 + a: y^2 = x^3 + (a - 5t) x + (b - 7 x0 t)), up to depth 3, tested for isomorphism with the second')
say('(a\' = u^4 a, b\' = u^6 b with u rational).  d is read off from the signs a_p(b)/a_p(a) and confirmed by the search.')
def sqfree(r):
    r = sp.Rational(r); sgn = 1 if r > 0 else -1; out = 1
    for pp, k in sp.factorint(abs(sp.numer(r)) * sp.denom(r)).items():
        if k % 2: out *= pp
    return sgn * out
def is_iso(c1, c2):
    a1, b1 = c1; a2, b2 = c2
    assert a1 != 0 and b1 != 0, 'j = 0 or 1728 is not handled here'
    if a2 == 0 or b2 == 0: return False
    u2 = (b2 / b1) / (a2 / a1)
    return u2 > 0 and sqfree(u2) == 1 and u2**2 == a2 / a1 and u2**3 == b2 / b1
def two_isog(c):
    a, b = c; X = sp.symbols('X')
    out = []
    for r in sp.Poly(X**3 + a * X + b, X).ground_roots():
        t = 3 * r**2 + a
        out.append((a - 5 * t, b - 7 * r * t))
    return out
def isog_search(c1, c2, depth=3):
    frontier = [(c1, 0)]; seen = [c1]
    while frontier:
        c, dpt = frontier.pop(0)
        if is_iso(c, c2): return dpt
        if dpt < depth:
            for c_ in two_isog(c):
                if not any(is_iso(c_, s_) for s_ in seen): seen.append(c_); frontier.append((c_, dpt + 1))
    return None
for a, b in (('F0*F8', 'L*U'), ('F0*F4', 'F4*L'), ('F0*F4', 'F4*U'), ('F8*F4', 'F4*L'), ('F8*F4', 'F4*U')):
    I1, J1 = IJ(P[a]); I2, J2 = IJ(P[b])
    sg = {}
    for p in primes:
        u, v = ap_c(a, p), ap_c(b, p)
        if u is None or v is None or u == 0: continue
        assert abs(u) == abs(v); sg[p] = u * v // abs(u * v)
    # d with (d/p) = sign at all those p: try squarefree divisors of 2 * (primes of the two discriminants)
    cand = None
    bad = sorted(set(sp.primefactors(sp.discriminant(sp.Poly(list(reversed(P[a])), x)) * sp.discriminant(sp.Poly(list(reversed(P[b])), x)))))
    for r in range(0, len(bad) + 1):
        for sub in itertools.combinations(bad, r):
            for sgn in (1, -1):
                dd = sgn * int(sp.prod(sub)) if sub else sgn
                if all(sp.jacobi_symbol(dd % p, p) == sg[p] for p in sg): cand = dd; break
            if cand is not None: break
        if cand is not None: break
    c1 = (-27 * cand**2 * I1, -27 * cand**3 * J1); c2 = (-27 * I2, -27 * J2)
    dpt = isog_search(c1, c2)
    say('  %-6s %-6s |a_p| equal at the good p < 1000; signs fit d = %s; %s' % (a, b, cand,
        'the twist by d is isomorphic to it' if dpt == 0 else ('isogenous to the twist by d by a chain of %d rational 2-isogenies (exact)' % dpt if dpt else 'NO isogeny found up to degree 8')))
    assert dpt is not None
say()

say('Multiplicities of these classes in the Jacobian of the curve of the five forms: genus3_mult.py.')
if WRITE:
    open(os.path.join(HERE, 'genus3_isog_%d_%d.txt' % (MM, NN)), 'w').write('\n'.join(OUT) + '\n')
