#!/usr/bin/env python3
"""The genus-3 curve X: y^2 = F0 F8 L U of a ratio m : n splits over Q into three elliptic curves (claim (A) of the
separate session's note of 27 November 2026, section 2), with explicit maps, checked symbolically for general m, n.

Written by an AI model (build agent t12); not read by a person.  Experiment, tier D (the written argument in NOTES.md
is tier C).  Nothing here changes a count.

    python3 -I genus3_split.py            prints genus3_split.txt (sympy, exact; then a trace control, pure Python)
    python3 -I genus3_split.py --write    also writes genus3_split.txt

Notation (FOUR_POINTS.md section 1): F0 = m^2 + n^2 t^2, F8 = n^2 + m^2 t^2, G4 = 2S(t^2 + 1), L = GL = G4 - 8mnt,
U = GU = G4 + 8mnt, S = m^2 + n^2; x = t^2, z = x + 1/x; A = m^2 n^2, B = m^4 + n^4, C = S^2, D = 2S^2 - 16 m^2 n^2.

  E1: w^2 = (z^2 - 4)(Az + B)(Cz + D),  E2: w^2 = (z + 2)(Az + B)(Cz + D),  E3: w^2 = (z - 2)(Az + B)(Cz + D).

What is checked (every step exact, sympy over Q(m, n)):
  (1) L, U are (m+n)^2 (t-1)^2 + (m-n)^2 (t+1)^2 and (m-n)^2 (t-1)^2 + (m+n)^2 (t+1)^2 (the note's L, U);
      f = F0 F8 L U is even and palindromic: f(-t) = f(t), t^8 f(1/t) = f(t);  F0 F8 = x (Az + B), L U = 4x (Cz + D).
  (2) sigma (t, y) = (-t, y), tau (t, y) = (1/t, y/t^4), iota (t, y) = (t, -y) are automorphisms of X over Q, of
      order 2, commuting: a group (Z/2)^3.
  (3) Y = y/(2t^2) and z are fixed by G' = <sigma, tau>, and Y^2 = (Az + B)(Cz + D): X/G' is (birational to) a conic.
  (4) phi_i = (z, w_i) with w1 = Y (t^2 - t^-2), w2 = Y (t + 1/t), w3 = Y (t - 1/t) satisfy the equation of E_i, and
      phi_i o g = phi_i exactly for g in K_i, phi_i o g = (z, -w_i) otherwise, with
      K1 = <sigma, iota tau>, K2 = <tau, iota sigma>, K3 = <sigma tau, iota sigma>.
  (5) deg phi_i = 4: z has degree 4 in t, so 8 on X, and degree 2 on E_i.
  (6) In the group ring Z[(Z/2)^3]: sum_i sum_{g in K_i} g = 4 - sum_{g in G'} g once iota is replaced by -1.
  (7) The discriminants: f is squarefree of degree 8 and the three polynomials in z are squarefree of degrees
      4, 3, 3 iff m n (m^2 - n^2)(m^2 + n^2)(m^8 + 14 m^4 n^4 + n^8) != 0 (factorised resultants and discriminants).
The written argument from (2)-(7) to J(X) ~ E1 x E2 x E3 over Q is in NOTES.md.
Then a control (not a proof): #X(F_p) = p + 1 - sum a_p(E_i) for the nine open ratios at every good prime p < 600,
and #X(F_{p^2}) = p^2 + 1 - sum (a_p(E_i)^2 - 2p) at good p < 60, all point counts in pure Python.
"""
import sys, os
import sympy as sp

WRITE = '--write' in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)

m, n, t, y, z = sp.symbols('m n t y z')
S = m**2 + n**2
F0 = m**2 + n**2 * t**2
F8 = n**2 + m**2 * t**2
G4 = 2 * S * (t**2 + 1)
L = G4 - 8 * m * n * t
U = G4 + 8 * m * n * t
f = sp.expand(F0 * F8 * L * U)
A, B, C, D = m**2 * n**2, m**4 + n**4, S**2, 2 * S**2 - 16 * m**2 * n**2
P = {1: (z**2 - 4) * (A * z + B) * (C * z + D), 2: (z + 2) * (A * z + B) * (C * z + D), 3: (z - 2) * (A * z + B) * (C * z + D)}

def zero(e):
    return sp.cancel(sp.together(sp.expand(e))) == 0

say('genus3_split.py: X: y^2 = F0 F8 L U and its three elliptic quotients (sympy %s; exact over Q(m, n))' % sp.__version__)
say()
# (1)
ok1 = zero(L - ((m + n)**2 * (t - 1)**2 + (m - n)**2 * (t + 1)**2)) and zero(U - ((m - n)**2 * (t - 1)**2 + (m + n)**2 * (t + 1)**2))
say('(1) L = (m+n)^2(t-1)^2 + (m-n)^2(t+1)^2 and U = (m-n)^2(t-1)^2 + (m+n)^2(t+1)^2 (the note\'s L, U): %s' % ok1)
ev = zero(f.subs(t, -t) - f); pal = zero(sp.expand(t**8 * f.subs(t, 1 / t)) - f)
say('    f(-t) = f(t): %s;  t^8 f(1/t) = f(t): %s;  deg_t f = %d, leading coefficient %s' % (ev, pal, sp.degree(f, t), sp.factor(sp.LC(sp.Poly(f, t)))))
xx = t**2; zz = xx + 1 / xx
ok_a = zero(F0 * F8 - xx * (A * zz + B)); ok_b = zero(L * U - 4 * xx * (C * zz + D))
say('    F0 F8 = x(Az + B): %s;  L U = 4x(Cz + D): %s   (x = t^2, z = x + 1/x)' % (ok_a, ok_b))
assert ok1 and ev and pal and ok_a and ok_b

# (2) automorphisms as maps (t, y) -> (t', y')
sig = (lambda T, Y_: (-T, Y_))
tau = (lambda T, Y_: (1 / T, Y_ / T**4))
iot = (lambda T, Y_: (T, -Y_))
def comp(g, h): return (lambda T, Y_: g(*h(T, Y_)))
ident = (lambda T, Y_: (T, Y_))
def on_curve(g):
    T2, Y2 = g(t, y)
    # y'^2 - f(t') must vanish when y^2 = f(t): y'^2 is y^2 times a function of t
    e = sp.expand(Y2**2 - f.subs(t, T2))
    return zero(e.subs(y**2, f)) if e.has(y) else zero(e)
def same(g, h):
    a, b = g(t, y), h(t, y)
    return zero(a[0] - b[0]) and zero(a[1] - b[1])
auts = {'sigma': sig, 'tau': tau, 'iota': iot}
ok2 = all(on_curve(g) for g in auts.values()) and all(same(comp(g, g), ident) for g in auts.values())
ok2 = ok2 and same(comp(sig, tau), comp(tau, sig)) and same(comp(sig, iot), comp(iot, sig)) and same(comp(tau, iot), comp(iot, tau))
say('(2) sigma, tau, iota map X to X, have order 2 and commute: %s' % ok2)
assert ok2
elems = {}
for a in range(2):
    for b in range(2):
        for c in range(2):
            g = ident
            if a: g = comp(sig, g)
            if b: g = comp(tau, g)
            if c: g = comp(iot, g)
            elems[(a, b, c)] = g   # exponent vector of (sigma, tau, iota)
say('    the eight elements sigma^a tau^b iota^c are pairwise distinct as maps: %s'
    % all(not same(elems[u], elems[v]) for u in elems for v in elems if u < v))

# (3) the conic
Yf = y / (2 * t**2)
def fn_eq(F, Gm):
    """F(g(t, y)) == G(t, y) for rational functions F of (t, y)"""
    return zero(F - Gm)
def pull(expr, g):
    T2, Y2 = g(t, y)
    return expr.subs({t: T2, y: Y2}, simultaneous=True)
ok3 = all(fn_eq(pull(Yf, g), Yf) and fn_eq(pull(zz, g), zz) for g in (sig, tau))
cone = sp.expand(Yf**2 - (A * zz + B) * (C * zz + D)).subs(y**2, f)
ok3 = ok3 and zero(cone)
say('(3) Y = y/(2t^2) and z are fixed by sigma and tau, and Y^2 = (Az + B)(Cz + D) on X: %s' % ok3)
assert ok3

# (4) the maps
W = {1: Yf * (t**2 - t**-2), 2: Yf * (t + 1 / t), 3: Yf * (t - 1 / t)}
K = {1: [(0, 0, 0), (1, 0, 0), (0, 1, 1), (1, 1, 1)],     # <sigma, iota tau>
     2: [(0, 0, 0), (0, 1, 0), (1, 0, 1), (1, 1, 1)],     # <tau, iota sigma>
     3: [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)]}     # <sigma tau, iota sigma>
for i in (1, 2, 3):
    eq = sp.expand(W[i]**2 - P[i].subs(z, zz)).subs(y**2, f)
    on = zero(eq)
    inv = []
    for e, g in elems.items():
        zi = fn_eq(pull(zz, g), zz)
        if e in K[i]: inv.append(zi and fn_eq(pull(W[i], g), W[i]))
        else: inv.append(zi and fn_eq(pull(W[i], g), -W[i]))
    say('(4) phi_%d = (z, w_%d), w_%d = %s: on E%d: %s; fixed by K%d = %s and (z, w) -> (z, -w) for the other four: %s'
        % (i, i, i, '(y/2t^2)(%s)' % (['t^2 - t^-2', 't + 1/t', 't - 1/t'][i - 1]), i, on, i,
           ['<sigma, iota tau>', '<tau, iota sigma>', '<sigma tau, iota sigma>'][i - 1], all(inv)))
    assert on and all(inv)
# (5) degrees
say('(5) deg phi_i = 4 for each i: z = t^2 + t^-2 is (t^4 + 1)/t^2, of degree %d in t, so of degree 8 on X (deg t = 2);'
    % max(sp.degree(sp.numer(sp.together(zz)), t), sp.degree(sp.denom(sp.together(zz)), t)))
say('    on E_i, z has degree 2; so deg phi_i = 8/2 = 4 = |K_i| and phi_i is the quotient map X -> X/K_i (NOTES.md)')

# (6) the group-ring identity
def mul(u, v): return tuple((a + b) % 2 for a, b in zip(u, v))
ring = {}
for i in (1, 2, 3):
    for g in K[i]: ring[g] = ring.get(g, 0) + 1
# replace iota by -1: g = (a, b, 1) contributes -[ (a, b, 0) ]
red = {}
for g, c in ring.items():
    sgn = -1 if g[2] else 1
    key = (g[0], g[1]); red[key] = red.get(key, 0) + sgn * c
target = {(0, 0): 3, (1, 0): -1, (0, 1): -1, (1, 1): -1}   # 4 - (1 + sigma + tau + sigma tau)
ok6 = all(red.get(k, 0) == v for k, v in target.items()) and all(k in target for k in red)
say('(6) sum_i sum_{g in K_i} g = %s; with iota = -1 it is 3 - sigma - tau - sigma tau = 4 - sum_{G\'} g: %s'
    % (' + '.join('%d*%s' % (c, 'sigma^%d tau^%d iota^%d' % g) for g, c in sorted(ring.items())), ok6))
assert ok6

# (7) discriminants
quads = {'F0': F0, 'F8': F8, 'L': L, 'U': U}
say('(7) discriminants and resultants (factorised):')
FACT = []
for k, q in quads.items():
    d = sp.factor(sp.discriminant(q, t)); say('    disc_t %s = %s' % (k, d)); FACT.append(d)
ks = list(quads)
for a in range(4):
    for b in range(a + 1, 4):
        r = sp.factor(sp.resultant(quads[ks[a]], quads[ks[b]], t)); say('    Res_t(%s, %s) = %s' % (ks[a], ks[b], r)); FACT.append(r)
for i in (1, 2, 3):
    d = sp.factor(sp.discriminant(sp.expand(P[i]), z)); lc = sp.factor(sp.LC(sp.Poly(sp.expand(P[i]), z)))
    say('    E%d: degree %d in z, leading coefficient %s, discriminant %s' % (i, sp.degree(sp.expand(P[i]), z), lc, d)); FACT += [d, lc]
oct_ = m**8 + 14 * m**4 * n**4 + n**8
q1 = sp.expand(m**4 - 2*m**3*n + 2*m**2*n**2 + 2*m*n**3 + n**4); q2 = sp.expand(m**4 + 2*m**3*n + 2*m**2*n**2 - 2*m*n**3 + n**4)
assert sp.expand(q1 * q2 - oct_) == 0
allowed = {m, n, sp.expand(m - n), sp.expand(m + n), sp.expand(S), q1, q2}
for d in FACT:
    c, fl = sp.factor_list(d)
    c = abs(int(c)); assert c & (c - 1) == 0, ('content not a power of 2', c)
    for g, e in fl: assert sp.expand(g) in allowed or sp.expand(-g) in allowed, ('unexpected factor', g)
say('    checked: every content is +-2^k and every irreducible factor is one of m, n, m - n, m + n, m^2 + n^2 and')
say('    the two quartic factors q1 q2 = m^8 + 14m^4n^4 + n^8 (product checked); q1 q2 > 0 for real (m, n) != 0')
say('    (the leading coefficients of F0, F8, L, U are n^2, m^2, 2S, 2S); over Q the condition is m n (m^2 - n^2) != 0')

# trace control ---------------------------------------------------------------------------------------------
OPEN = [(195, 28), (223, 167), (440, 279), (609, 200), (627, 364), (719, 161), (809, 409), (861, 620), (991, 263)]
def primes(N):
    s = bytearray([1]) * N; s[0] = s[1] = 0
    for i in range(2, int(N**0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N) if s[i]]
def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
def coeffs_t(M, N):
    return [int(c) for c in sp.Poly(f.subs({m: M, n: N}), t).all_coeffs()][::-1]   # low to high
def coeffs_z(i, M, N):
    return [int(c) for c in sp.Poly(sp.expand(P[i].subs({m: M, n: N})), z).all_coeffs()][::-1]
def ev_(cs, x, p):
    r = 0
    for c in reversed(cs): r = (r * x + c) % p
    return r
def count(cs, p):
    """#points over F_p of the smooth model of w^2 = poly (degree 3 or 4 or 8, nonzero leading coeff mod p)"""
    deg = len(cs) - 1
    N = sum(1 + chi(ev_(cs, x, p), p) for x in range(p))
    if deg % 2: return N + 1
    return N + 1 + chi(cs[-1], p)
# F_{p^2} = F_p[s]/(s^2 - r), r a non-residue
def count2(cs, p):
    r = next(a for a in range(2, p) if chi(a, p) == -1)
    def mulq(u, v): return ((u[0] * v[0] + r * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)
    def evq(x):
        acc = (0, 0)
        for c in reversed(cs): acc = mulq(acc, x); acc = ((acc[0] + c) % p, acc[1])
        return acc
    def chiq(u):
        nm = (u[0] * u[0] - r * u[1] * u[1]) % p     # the quadratic character of F_{p^2} is chi(norm)
        return 0 if u == (0, 0) else chi(nm, p)
    deg = len(cs) - 1
    N = sum(1 + chiq(evq((a, b))) for a in range(p) for b in range(p))
    if deg % 2: return N + 1
    return N + 1 + 1          # the leading coefficient is in F_p, hence a square in F_{p^2}
say()
say('Control, not a proof: point counts over F_p and F_{p^2} (pure Python; the curves above with m, n substituted).')
tot1 = tot2 = 0
for (M, N) in OPEN:
    cf = coeffs_t(M, N); cz = {i: coeffs_z(i, M, N) for i in (1, 2, 3)}
    bad = 2 * M * N * (M * M - N * N) * (M * M + N * N) * (M**8 + 14 * M**4 * N**4 + N**8)
    good = [p for p in primes(600) if bad % p and p > 3]
    k1 = k2 = 0
    for p in good:
        ap = {i: p + 1 - count(cz[i], p) for i in (1, 2, 3)}
        aX = p + 1 - count(cf, p)
        assert aX == sum(ap.values()), ('trace mismatch', M, N, p)
        k1 += 1
        if p < 60:
            NX2 = count2(cf, p)
            assert NX2 == p * p + 1 - sum(a * a - 2 * p for a in ap.values()), ('F_p^2 mismatch', M, N, p)
            for i in (1, 2, 3):
                assert count2(cz[i], p) == p * p + 1 - (ap[i] ** 2 - 2 * p)
            k2 += 1
    tot1 += k1; tot2 += k2
    say('  (%d,%d): a_p(X) = a_p(E1) + a_p(E2) + a_p(E3) at all %d good primes 5 <= p < 600; #X(F_{p^2}) agrees at %d primes < 60'
        % (M, N, k1, k2))
say('  total: %d prime tests over F_p, %d over F_{p^2}; no mismatch' % (tot1, tot2))
say()
say('Result: the maps phi_1, phi_2, phi_3 exist over Q for all m, n, and (6) with the argument of NOTES.md gives')
say('J(X) ~ E1 x E2 x E3 over Q whenever m n (m^2 - n^2) != 0 (m, n rational).')
if WRITE:
    open(os.path.join(HERE, 'genus3_split.txt'), 'w').write('\n'.join(OUT) + '\n')
