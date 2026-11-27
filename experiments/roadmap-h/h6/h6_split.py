"""H6 (part 1): the quotient curves of genus 1 and 2 of the (Z/2)^5-cover y_j^2 = F_j(t), j in {F0, F8, F4, L, U},
for one ratio m : n, and how each Jacobian splits (EXPERIMENT, written by an AI model; not read by a person; nothing
here changes the census).
For each subset S of the five forms with |S| >= 2 the curve C_S: y^2 = prod_{j in S} 4F_j(t) (4F_j is integral; the
factor 4 is a square, so a solution t gives a rational point on every C_S) has genus |S| - 1, and
Jac(cover) ~ prod_S Jac(C_S).  A binary form F of degree 2k is split further by every rational involution of P^1
under which F is exactly invariant: x -> -x (F(x,-z) = F(x,z)) and x <-> z (F(z,x) = F(x,z)), the latter after the
linear change u = x + z, v = x - z.  The quotient by x -> -x of y^2 = P(x^2, z^2) is y^2 = P(w, 1) and y^2 = w P(w, 1);
the Jacobian is isogenous to the product of the two quotients' Jacobians (Kani-Rosen for the group {1, iota}).
Each invariance is an exact polynomial identity, and each quotient is checked by the identity
F(x, z) = P(x^2, z^2) (resp. its transform).  No PARI here: pure Python integers.
Output: the list of pieces (genus 1 quartics or cubics, and genus 2 pieces with no such involution), written to
h6_split_<m>_<n>.json with --write.  usage: h6_split.py m n [--write]   |   h6_split.py --negative"""
import sys, os, json, itertools
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
NAMES = ['F0', 'F8', 'F4', 'L', 'U']

def forms(m, n):
    S = m*m + n*n
    # coefficient lists c[0] + c[1] t + c[2] t^2, times 4 where needed to make them integral (4 is a square)
    return {'F0': [m*m, 0, n*n], 'F8': [n*n, 0, m*m], 'F4': [2*S, 0, 2*S], 'L': [2*S, -8*m*n, 2*S], 'U': [2*S, 8*m*n, 2*S]}

def pmul(a, b):
    r = [0]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r

def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0: a.pop()
    return a

def deg(a): return len(trim(a)) - 1

def genus_of(f, D):
    """f: binary form of degree D (coefficient list of length D+1, c[i] the coefficient of x^i z^(D-i)), D even."""
    return D//2 - 1

def bin_sub(f, D, A):
    """F(a x + b z, c x + d z) for the binary form f of degree D, A = (a, b, c, d); returns coefficients."""
    a, b, c, d = A
    out = [0]*(D + 1)
    for i, ci in enumerate(f):
        if ci == 0: continue
        # x^i z^(D-i) -> (a x + b z)^i (c x + d z)^(D-i)
        p = [1]
        for _ in range(i): p = pmul(p, [b, a])        # list in powers of x: b + a x  (z = 1, homogeneous by degree)
        for _ in range(D - i): p = pmul(p, [d, c])
        for k, v in enumerate(p): out[k] += ci*v
    return out

def even_split(f, D):
    """f(x, -z) = f(x, z)?  Then f = P(x^2, z^2); return the two quotient forms (P as a binary form in (w, z)), each
    padded to even degree: y^2 = P(w) and y^2 = w P(w)."""
    if any(f[i] != 0 for i in range(D + 1) if (D - i) % 2 == 1): return None
    k = D // 2
    P = [f[2*j] for j in range(k + 1)]                 # coefficient of w^j z'^(k-j), w = x^2, z' = z^2
    wP = [0] + P                                        # w * P: degree k + 1
    out = []
    for q, dq in ((P, k), (wP, k + 1)):
        dd = dq + (dq % 2)                              # as a binary form of even degree (a root at infinity if odd)
        out.append((q + [0]*(dd + 1 - len(q)), dd))
    # exact check: F(x, z) = P(x^2, z^2)
    for x, z in ((2, 3), (5, -7), (11, 4)):
        assert sum(c * x**i * z**(D - i) for i, c in enumerate(f)) == sum(c * (x*x)**j * (z*z)**(k - j) for j, c in enumerate(P))
    return out

def swap_split(f, D):
    """f(z, x) = f(x, z)?  Change variables x = u + v, z = u - v (so the swap is v -> -v) and split as above."""
    if any(f[i] != f[D - i] for i in range(D + 1)): return None
    g = bin_sub(f, D, (1, 1, 1, -1))                  # g(u, v) = f(u + v, u - v)... see the order of variables below
    # bin_sub returns coefficients of x^i where x stands for the first variable u and the second for v:
    # f(a u + b v, c u + d v) with (a, b, c, d) = (1, 1, 1, -1): x = u + v, z = u - v
    # we need it as a form in (u, v) with v -> -v the involution; even_split tests the second variable's parity.
    for u, v in ((2, 3), (5, -7)):
        assert sum(c * u**i * v**(D - i) for i, c in enumerate(g)) == sum(c * (u + v)**i * (u - v)**(D - i) for i, c in enumerate(f))
    return even_split(g, D)

def content_sqfree(f):
    """Divide out the largest square dividing the content (a square factor changes nothing over Q)."""
    from math import gcd
    g = 0
    for c in f: g = gcd(g, abs(c))
    s = 1; q = 2
    import sympy
    for p, e in sympy.factorint(g).items(): s *= p**(e // 2)
    return [c // (s*s) for c in f]

def nonsingular(f, D):
    """the binary form has no repeated root (its discriminant as a form is nonzero): test via sympy."""
    import sympy
    x = sympy.Symbol('x')
    poly = sympy.Poly(sum(c * x**i for i, c in enumerate(f)), x)
    dgr = poly.degree()
    if D - dgr >= 2: return False                       # a double root at infinity
    return sympy.discriminant(poly) != 0

def decompose(f, D, path):
    g = genus_of(f, D)
    assert nonsingular(f, D), (path, 'singular')
    if g <= 1: return [dict(path=path, genus=g, form=f, D=D)]
    for name, fun in (('x->-x', even_split), ('x<->z', swap_split)):
        r = fun(f, D)
        if r is not None:
            (q1, d1), (q2, d2) = r
            g1, g2 = genus_of(q1, d1), genus_of(q2, d2)
            assert g1 + g2 == g, (path, name, g1, g2, g)
            return decompose(content_sqfree(q1), d1, path + [name + ':P']) + decompose(content_sqfree(q2), d2, path + [name + ':wP'])
    return [dict(path=path, genus=g, form=f, D=D, simple='no rational involution of these two kinds')]

def run(m, n):
    F = forms(m, n); out = []
    for r in range(2, 6):
        for S in itertools.combinations(NAMES, r):
            f = [1]
            for j in S: f = pmul(f, F[j])
            D = 2*r
            pieces = decompose(content_sqfree(f), D, [])
            assert sum(p['genus'] for p in pieces) == r - 1
            for p in pieces:
                out.append(dict(S='*'.join(S), genus_S=r - 1, path=' '.join(p['path']) or '-', genus=p['genus'],
                                form=p['form'], D=p['D'], note=p.get('simple', '')))
    tot = sum(p['genus'] for p in out); assert tot == 49, tot     # genus of the (Z/2)^5-cover: 2g - 2 = -64 + 10*16
    return out

if __name__ == '__main__':
    argv = sys.argv
    if '--negative' in argv:
        lines = []
        # forged: a form that is not even must not split by x -> -x; a non-palindromic one not by the swap
        f = [1, 1, 0, 0, 1]; lines.append('x->-x on a non-even form: %s' % ('refused' if even_split(f, 4) is None else 'ACCEPTED'))
        f = [1, 2, 0, 3, 5]; lines.append('swap on a non-palindromic form: %s' % ('refused' if swap_split(f, 4) is None else 'ACCEPTED'))
        try: decompose([1, 2, 1], 2, ['forged square']); lines.append('singular form: ACCEPTED')
        except AssertionError: lines.append('singular form (x + 1)^2: refused')
        s = '\n'.join(lines) + '\n'; fn = os.path.join(HERE, 'h6_split_negative.txt')
    else:
        m, n = int(argv[1]), int(argv[2])
        res = run(m, n)
        s = json.dumps(res, indent=0) + '\n'; fn = os.path.join(HERE, 'h6_split_%d_%d.json' % (m, n))
    if '--write' in argv: open(fn, 'w').write(s)
    else:
        sys.stdout.write(s if '--negative' in argv else '%d pieces\n' % len(json.loads(s)))
        if os.path.exists(fn): print('# matches' if open(fn).read() == s else '# DIFFERS from ' + os.path.basename(fn))
