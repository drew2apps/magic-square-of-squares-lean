#!/usr/bin/env python3
"""How often can a rank-1 curve occur in the Jacobian of the curve X5 of the five forms?  (The Demjanenko-Manin
question of the separate session's note, 27 November 2026, section 2: "no rank-1 class occurs twice".)

Written by an AI model (build agent t12); not read by a person.  Experiment, tier D.  Nothing here changes a count.

    python3 -I genus3_mult.py [--write]       genus3_mult.txt   (the nine open ratios)

X5: y_i^2 = F_i(t) for the five forms F0, F8, F4, L, U (FOUR_POINTS.md section 1), genus 49.  Quoted (Kani-Rosen for
(Z/2)^5, and for an extra involution of a hyperelliptic curve): J(X5) ~ prod_I J(C_I), C_I: y^2 = prod_{i in I} F_i,
and each C_I that is even in t (y^2 = h(t^2)) splits as J(y^2 = h(x)) x J(y^2 = x h(x)); each C_I that is palindromic
(t^d h(t + 1/t) for deg 2d) splits by w = t + 1/t (pieces below).  This gives 29 pieces of genus 1 and 10 of genus 2
(genus 49 = 29 + 2 x 10):
  pairs (10):            y^2 = F_i F_j in t;
  {F0,F8,F4}:            y^2 = h(x), y^2 = x h(x), h = F0F8F4 at t^2 = x;
  {F0,F8,L}, {F0,F8,U}:  y^2 = g(w)(w + 2), y^2 = g(w)(w - 2), g = (A(w^2 - 2) + B)(2Sw -+ 8mn)   (each twice);
  {F0,L,U}, {F8,L,U}:    y^2 = F0(x) LU(x), y^2 = x F0(x) LU(x)                                    (each twice);
  {F4,L,U}:              y^2 = F4(x) LU(x), y^2 = x F4(x) LU(x);
  {F0,F8,F4,L}, {..,U}:  y^2 = g4(w) and the genus-2 y^2 = g4(w)(w^2 - 4), g4 = (A(w^2-2)+B) 2Sw (2Sw -+ 8mn) (twice);
  {F0,F8,L,U}:           E1, E2, E3 of the note (genus3_split.py);
  {F0,F4,L,U}, {F8,..}:  y^2 = h(x) and the genus-2 y^2 = x h(x), h = F0 F4 LU at x                 (twice);
  all five:              the genus-2 y^2 = h5(x) and y^2 = x h5(x), h5 = F0 F8 F4 LU at x;
  {F0,F4,L} and its three images under t -> -t, 1/t: one genus-2 curve, four times.
Checked here (a control of the piece list, not of the theorem): at two primes, L_p of each C_I of genus 3 and 4 from
point counts over F_{p^k} equals the product of the L_p of its pieces.
For a curve E (by a model) the multiplicity of E in J(X5) is at most: the number of genus-1 pieces with the same
a_p as E at every good p < 600 (a piece whose a_p differs once is not isogenous to E: exact), plus, for each genus-2
piece, 0 if at some good prime (1 - a_p T + p T^2) does not divide L_p(piece) (exact), else 2.  Lower bounds come
from explicit maps where the script names them.  Demjanenko-Manin with a curve of rank r needs multiplicity > r.
"""
import sys, os, itertools
import sympy as sp

WRITE = '--write' in sys.argv
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = []
def say(s=''):
    OUT.append(s); print(s, flush=True)
OPEN = [(195, 28), (223, 167), (440, 279), (609, 200), (627, 364), (719, 161), (809, 409), (861, 620), (991, 263)]
t, x, w = sp.symbols('t x w')

def chi(a, p):
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)
def ev(cs, v, p):
    r = 0
    for c in reversed(cs): r = (r * v + c) % p
    return r
def good_for(cs, p):
    if cs[-1] % p == 0: return False
    q = sp.Poly(list(reversed([c % p for c in cs])), x, modulus=p)
    return sp.degree(sp.gcd(q, q.diff(x))) == 0
def count1(cs, p):
    N = sum(1 + chi(ev(cs, v, p), p) for v in range(p))
    return N + (1 if (len(cs) - 1) % 2 else 1 + chi(cs[-1], p))
def count2(cs, p):
    r = next(a for a in range(2, p) if chi(a, p) == -1)
    def mul(u, v): return ((u[0] * v[0] + r * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)
    N = 0
    for a in range(p):
        for b in range(p):
            acc = (0, 0)
            for c in reversed(cs): acc = mul(acc, (a, b)); acc = ((acc[0] + c) % p, acc[1])
            N += 1 if acc == (0, 0) else (2 if chi((acc[0] * acc[0] - r * acc[1] * acc[1]) % p, p) == 1 else 0)
    return N + (1 if (len(cs) - 1) % 2 else 2)
def count3(cs, p, F):
    """over F_{p^k} given as (q, exp, log, zech) (Zech logarithms)"""
    q, exp, log, zech = F; Q1 = q - 1
    lc = [None if c % p == 0 else log[c % p] for c in cs]
    v0 = cs[0] % p
    N = 1 if v0 == 0 else (2 if log[v0] % 2 == 0 else 0)
    for lt in range(Q1):
        acc = None
        for c in reversed(lc):
            if acc is not None: acc = (acc + lt) % Q1
            if c is not None:
                if acc is None: acc = c
                else:
                    z = zech[(c - acc) % Q1]; acc = None if z is None else (acc + z) % Q1
        N += 1 if acc is None else (2 if acc % 2 == 0 else 0)
    deg = len(cs) - 1
    if deg % 2: return N + 1
    return N + (2 if log[cs[-1] % p] % 2 == 0 else 0)
def field(p, k):
    q = p**k; T = sp.symbols('T')
    def enc(v):
        r = 0
        for c in reversed(v): r = r * p + c
        return r
    for cand in itertools.product(range(p), repeat=k):
        f = sp.Poly([1] + list(cand), T, modulus=p)
        if not f.is_irreducible: continue
        fc = [int(c) % p for c in f.all_coeffs()]
        for c0 in range(p):
            exp = [0] * (q - 1); v = [1] + [0] * (k - 1); seen = set(); ok = True
            for i in range(q - 1):
                e_ = enc(v)
                if e_ in seen: ok = False; break
                seen.add(e_); exp[i] = e_
                top = v[-1]; ww = [0] + v[:-1]
                ww = [(ww[j] - top * fc[k - j]) % p for j in range(k)]
                v = [(ww[j] + c0 * v[j]) % p for j in range(k)]
            if ok and enc(v) == 1:
                log = [None] * q
                for i, e_ in enumerate(exp): log[e_] = i
                def plus1(e_): return e_ - (e_ % p) + ((e_ % p) + 1) % p
                zech = [None if plus1(exp[i]) == 0 else log[plus1(exp[i])] for i in range(q - 1)]
                return q, exp, log, zech
    raise AssertionError('no field')
def lpoly(cs, p, g, fields=None):
    N = [None, count1(cs, p)]
    if g >= 2: N.append(count2(cs, p))
    for k in range(3, g + 1): N.append(count3(cs, p, fields[k]))
    s = [None] + [p**k + 1 - N[k] for k in range(1, g + 1)]
    c = [1]
    for k in range(1, g + 1):
        tot = -sum(s[i] * c[k - i] for i in range(1, k + 1)); assert tot % k == 0; c.append(tot // k)
    full = c + [0] * g
    for k in range(g + 1, 2 * g + 1): full[k] = p**(k - g) * c[2 * g - k]
    return full
def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b): r[i + j] += u * v
    return r

def pieces(M, N):
    m, n = sp.Integer(M), sp.Integer(N); S = m * m + n * n
    A, B, C, D = m**2 * n**2, m**4 + n**4, S**2, 2 * S**2 - 16 * m**2 * n**2
    Ft = {'F0': m**2 + n**2 * t**2, 'F8': n**2 + m**2 * t**2, 'F4': 2 * S * (t**2 + 1),
          'L': 2 * S * (t**2 + 1) - 8 * m * n * t, 'U': 2 * S * (t**2 + 1) + 8 * m * n * t}
    Fx = {'F0': m**2 + n**2 * x, 'F8': n**2 + m**2 * x, 'F4': 2 * S * (x + 1)}
    LU = 4 * S**2 * (x + 1)**2 - 64 * m**2 * n**2 * x
    gL = (A * (w**2 - 2) + B) * (2 * S * w - 8 * m * n); gU = (A * (w**2 - 2) + B) * (2 * S * w + 8 * m * n)
    g4L = (A * (w**2 - 2) + B) * 2 * S * w * (2 * S * w - 8 * m * n); g4U = (A * (w**2 - 2) + B) * 2 * S * w * (2 * S * w + 8 * m * n)
    P = []   # (name, genus, expr, variable, copies)
    for a, b in itertools.combinations(Ft, 2): P.append((a + '*' + b, 1, Ft[a] * Ft[b], t, 1))
    h = Fx['F0'] * Fx['F8'] * Fx['F4']; P += [('{F0F8F4}/s', 1, h, x, 1), ('{F0F8F4}/is', 1, x * h, x, 1)]
    P += [('{F0F8L}/tau', 1, gL * (w + 2), w, 1), ('{F0F8L}/itau', 1, gL * (w - 2), w, 1),
          ('{F0F8U}/tau', 1, gU * (w + 2), w, 1), ('{F0F8U}/itau', 1, gU * (w - 2), w, 1)]
    P += [('{F0LU}/s', 1, Fx['F0'] * LU, x, 1), ('{F0LU}/is', 1, x * Fx['F0'] * LU, x, 1),
          ('{F8LU}/s', 1, Fx['F8'] * LU, x, 1), ('{F8LU}/is', 1, x * Fx['F8'] * LU, x, 1)]
    P += [('{F4LU}/s', 1, Fx['F4'] * LU, x, 1), ('{F4LU}/is', 1, x * Fx['F4'] * LU, x, 1)]
    P += [('{F0F8F4L}/tau', 1, g4L, w, 1), ('{F0F8F4L}/itau', 2, g4L * (w**2 - 4), w, 1),
          ('{F0F8F4U}/tau', 1, g4U, w, 1), ('{F0F8F4U}/itau', 2, g4U * (w**2 - 4), w, 1)]
    P += [('E1', 1, (w**2 - 4) * (A * w + B) * (C * w + D), w, 1), ('E2', 1, (w + 2) * (A * w + B) * (C * w + D), w, 1),
          ('E3', 1, (w - 2) * (A * w + B) * (C * w + D), w, 1)]
    for f in ('F0', 'F8'):
        h = Fx[f] * Fx['F4'] * LU
        P += [('{%sF4LU}/s' % f, 1, h, x, 1), ('{%sF4LU}/is' % f, 2, x * h, x, 1)]
    h5 = Fx['F0'] * Fx['F8'] * Fx['F4'] * LU
    P += [('{all}/s', 2, h5, x, 1), ('{all}/is', 2, x * h5, x, 1)]
    for a, b in (('F0', 'L'), ('F0', 'U'), ('F8', 'L'), ('F8', 'U')):
        P.append(('{%sF4%s}' % (a, b), 2, Ft[a] * Ft['F4'] * Ft[b], t, 1))
    out = []
    for name, g, e, var, k in P:
        cs = [int(c) for c in sp.Poly(sp.expand(e), var).all_coeffs()][::-1]
        out.append((name, g, cs))
    assert sum(g for _, g, _ in out) == 49
    return out, Ft

def control(M, N, prs):
    """L_p of the genus-3 and genus-4 C_I (direct counts) = product of the L_p of their pieces"""
    PC, Ft = pieces(M, N)
    byname = {nm: (g, cs) for nm, g, cs in PC}
    groups = {('F0', 'F8', 'F4', 'L'): ['{F0F8F4L}/tau', '{F0F8F4L}/itau'], ('F0', 'F8', 'L', 'U'): ['E1', 'E2', 'E3'],
              ('F0', 'F4', 'L', 'U'): ['{F0F4LU}/s', '{F0F4LU}/is'], ('F0', 'F8', 'F4', 'L', 'U'): ['{all}/s', '{all}/is'],
              ('F0', 'F8', 'F4'): ['{F0F8F4}/s', '{F0F8F4}/is'], ('F0', 'F8', 'L'): ['{F0F8L}/tau', '{F0F8L}/itau'],
              ('F0', 'L', 'U'): ['{F0LU}/s', '{F0LU}/is'], ('F4', 'L', 'U'): ['{F4LU}/s', '{F4LU}/is']}
    res = []
    for p in prs:
        fields = {3: field(p, 3), 4: field(p, 4)}
        for I, names in groups.items():
            e = sp.Integer(1)
            for f in I: e *= Ft[f]
            cs = [int(c) for c in sp.Poly(sp.expand(e), t).all_coeffs()][::-1]
            assert good_for(cs, p)
            g = len(I) - 1
            lhs = lpoly(cs, p, g, fields)
            rhs = [1]
            for nm in names:
                gg, cc = byname[nm]; assert good_for(cc, p)
                rhs = pmul(rhs, lpoly(cc, p, gg))
            assert lhs == rhs, ('piece list wrong', I, p)
            res.append((I, p))
    return len(res)

def run():
    say('genus3_mult.py: multiplicities of rank-1 curves in J(X5), X5 the curve of the five forms (genus 49)')
    nc = control(627, 364, [17, 23])
    say('control at (627,364): for 8 curves C_I of genus 2, 3, 4, L_p(C_I) = product of the L_p of the pieces at p = 17, 23'
        ' (%d checks, counts over F_{p^k}, k <= 4)' % nc)
    say()
    prs = [int(p) for p in sp.primerange(5, 600)]
    for (M, N) in OPEN:
        PC, _ = pieces(M, N)
        g1 = [(nm, cs) for nm, g, cs in PC if g == 1]; g2 = [(nm, cs) for nm, g, cs in PC if g == 2]
        APc = {}
        def ap(nm, cs, p):
            if (nm, p) not in APc: APc[(nm, p)] = (p + 1 - count1(cs, p)) if good_for(cs, p) else None
            return APc[(nm, p)]
        # classes of the genus-1 pieces by a_p at good p < 600
        cls = []
        for nm, cs in g1:
            for c in cls:
                rep = c[0]
                if all(ap(nm, cs, p) == ap(rep[0], rep[1], p) for p in prs if ap(nm, cs, p) is not None and ap(rep[0], rep[1], p) is not None):
                    c.append((nm, cs)); break
            else: cls.append([(nm, cs)])
        say('(%d,%d): the 29 genus-1 pieces fall into %d classes by a_p at the good primes < 600:' % (M, N, len(cls)))
        L2 = {}
        for c in cls:
            nm0, cs0 = c[0]
            # genus-2 pieces: exclusion prime
            ex = []
            for nm2, cs2 in g2:
                hit = None
                for p in prs[:40]:
                    a = ap(nm0, cs0, p)
                    if a is None or not good_for(cs2, p): continue
                    if (nm2, p) not in L2: L2[(nm2, p)] = lpoly(cs2, p, 2)
                    l = L2[(nm2, p)]
                    b = l[1] + a
                    if l[2] != 2 * p - a * b: hit = p; break
                ex.append((nm2, hit))
            free = [nm2 for nm2, hp in ex if hp is None]
            bound = len(c) + 2 * len(free)
            say('  %-48s size %d; genus-2 pieces excluded at %s%s => multiplicity <= %d'
                % ('[' + ', '.join(nm for nm, _ in c) + ']', len(c),
                   'p <= %d' % max(hp for _, hp in ex if hp is not None) if any(hp for _, hp in ex) else '-',
                   ('; NOT excluded: %s' % free) if free else '', bound))
        say()

if __name__ == '__main__':
    run()
    if WRITE: open(os.path.join(HERE, 'genus3_mult.txt'), 'w').write('\n'.join(OUT) + '\n')
