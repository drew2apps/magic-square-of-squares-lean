"""H6 (part 2): ranks of the genus-1 pieces found by h6_split.py, for (627,364) and (991,263), and a test of the
genus-2 pieces (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
Genus 1 piece y^2 = f(x) (f of degree 3 or 4): its Jacobian is y^2 = x^3 - 27 I x - 27 J with I, J the invariants of
the quartic, computed here in Python; PARI's ellfromeqn must give the same minimal model (two constructions).
Pieces are grouped by isogeny class (the set of minimal models of the class, ellisomat).  For each class, PARI's
ellrank (setrand(1), effort 0, 1, 2, 3 until the bounds meet, a time limit per call) on every curve 2-isogenous to it
(ellisomat at p = 2): the rank is the same on isogenous curves, so the class gets the largest lower bound and the
smallest upper bound.  These are PARI's bounds (a finder): an upper bound is a 2-Selmer bound, and a lower bound can
exceed the number of points returned (then it rests on Sha[2] having even rank).  Labels: computed (numerical,
PARI 2.19.1), consistent with PARI.
Genus 2 piece with no involution found: the characteristic polynomial of Frobenius at small good primes, from point
counts over F_p and F_{p^2} in pure Python; if it is irreducible over Q at one good prime, the Jacobian is simple
over Q (exact: a Jacobian isogenous over Q to a product of elliptic curves has a reducible polynomial at every good p).
usage: h6_ranks.py m n [--write]      (reads h6_split_<m>_<n>.json; writes h6_ranks_<m>_<n>.txt)
       h6_ranks.py --negative         (h6_ranks_negative.txt)"""
import sys, os, json
from fractions import Fraction as Fr
import cypari2, sympy
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
HERE = os.path.dirname(os.path.abspath(__file__))
TLIM = int(os.environ.get('H6_TLIM', '900'))

def IJ(f):
    e, d, c, b, a = (f + [0]*5)[:5]                      # f = e + d x + c x^2 + b x^3 + a x^4
    I = 12*a*e - 3*b*d + c*c
    J = 72*a*c*e + 9*b*c*d - 27*a*d*d - 27*e*b*b - 2*c**3
    return I, J

def poly_str(f): return '+'.join('(%d)*x^%d' % (c, i) for i, c in enumerate(f) if c)

def minmodel(ai): return [int(v) for v in pari.ellminimalmodel(pari.ellinit(ai))[:5]]

def jac(f):
    I, J = IJ(f)
    mine = minmodel([0, 0, 0, -27*I, -27*J])
    theirs = minmodel(pari('ellfromeqn(y^2 - (%s))' % poly_str(f)))
    assert mine == theirs, ('Jacobian models differ', f)
    return mine

def item_model(item, j0):
    v = list(item)
    if len(v) == 2: ai = [0, 0, 0, v[0], v[1]]                       # [a4, a6]
    else: ai = [0, 0, 0, -27*v[0], -54*v[1]]; assert pari('ellinit(%s).j' % str(ai)) == v[2]   # [c4, c6, j]
    return tuple(minmodel(ai))

def iso_list(ai, p):
    E = pari.ellinit(ai); out = sorted(set(item_model(it, None) for it in pari.ellisomat(E, p, 1)[0]))
    assert tuple(minmodel(ai)) in out, 'the curve itself is not in its isogeny class'
    return out

def iso_class(ai): return tuple(iso_list(ai, 0))

def two_class(ai): return iso_list(ai, 2)

def rank_one(ai):
    pari('setrand(1)'); pari('E_h6 = ellinit(%s)' % str(list(ai)))
    best = None
    for eff in (0, 1, 2, 3):
        try: r = pari('alarm(%d, ellrank(E_h6, %d))' % (TLIM, eff))
        except cypari2.PariError: return best if best else (0, None, 0, eff, 'time limit')
        lo, hi, npts = int(r[0]), int(r[1]), len(r[3])
        best = (lo, hi, npts, eff, '')
        if lo == hi: break
    return best

def rank_class(ai):
    rows = []
    for c in two_class(ai):
        rows.append((c, rank_one(list(c))))
    los = [r[1][0] for r in rows]; his = [r[1][1] for r in rows if r[1][1] is not None]
    return max(los), (min(his) if his else None), rows

# ---- genus 2: Frobenius polynomial from point counts ----
def chi_p(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1)//2, p) == 1 else -1

def count_Fp(f, p, D):
    N = sum(1 + chi_p(sum(c * pow(x, i, p) for i, c in enumerate(f)), p) for x in range(p))
    dgr = max(i for i, c in enumerate(f) if c % p)
    N += (1 + chi_p(f[dgr], p)) if dgr == D else 1          # points at infinity: degree 6 or 5
    return N

def count_Fp2(f, p, D):
    g = next(a for a in range(2, p) if chi_p(a, p) == -1)   # F_{p^2} = F_p(sqrt g)
    def mul(u, v): return ((u[0]*v[0] + g*u[1]*v[1]) % p, (u[0]*v[1] + u[1]*v[0]) % p)
    def pw(u, e):
        r = (1, 0)
        while e:
            if e & 1: r = mul(r, u)
            u = mul(u, u); e >>= 1
        return r
    q = p*p; N = 0
    for a in range(p):
        for b in range(p):
            x = (a, b); val = (0, 0); xp = (1, 0)
            for c in f:
                val = ((val[0] + c*xp[0]) % p, (val[1] + c*xp[1]) % p); xp = mul(xp, x)
            if val == (0, 0): N += 1
            else: N += 1 + (1 if pw(val, (q - 1)//2) == (1, 0) else -1)
    dgr = max(i for i, c in enumerate(f) if c % p)
    N += 2 if dgr == D else 1                                # every element of F_p is a square in F_{p^2}
    return N

def frob(f, p, D):
    N1, N2 = count_Fp(f, p, D), count_Fp2(f, p, D)
    s1 = p + 1 - N1; s2 = (s1*s1 - (p*p + 1 - N2)) // 2
    assert (s1*s1 - (p*p + 1 - N2)) % 2 == 0
    T = sympy.Symbol('T')
    return sympy.Poly(T**4 - s1*T**3 + s2*T**2 - p*s1*T + p*p, T)

def good(f, p, D):
    x = sympy.Symbol('x'); P = sympy.Poly(sum(c * x**i for i, c in enumerate(f)), x)
    dgr = P.degree()
    return p > 2 and f[dgr] % p != 0 and sympy.discriminant(P) % p != 0

def simple_test(f, D, primes=(3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)):
    for p in primes:
        if not good(f, p, D): continue
        P = frob(f, p, D)
        if len(sympy.factor_list(P.as_expr())[1]) == 1 and sympy.factor_list(P.as_expr())[1][0][1] == 1:
            return 'simple over Q (Frobenius polynomial at p = %d irreducible: %s)' % (p, str(P.as_expr()).replace('**', '^'))
    return 'not shown simple (Frobenius polynomial reducible at every good prime tried)'

def run(m, n):
    pieces = json.load(open(os.path.join(HERE, 'h6_split_%d_%d.json' % (m, n))))
    out = ['H6: quotient curves of the (Z/2)^5-cover for (%d,%d).' % (m, n),
           'Written by an AI model; not read by a person. Experiment, not counted. PARI %s (ellrank as a finder).' % '.'.join(map(str, pari.version())),
           'piece: subset S, splitting path, genus; for genus 1 the isogeny class (c<k>) and its rank bounds [lo,hi] (max lo, min hi over the 2-isogeny class)']
    classes = {}; cls_rows = []
    for p in pieces:
        if p['genus'] == 1:
            f = p['form'][:]; f = f[:max(i for i, c in enumerate(f) if c) + 1]
            ai = jac(f); key = iso_class(ai)
            if key not in classes:
                lo, hi, rows = rank_class(ai)
                classes[key] = ('c%d' % (len(classes) + 1), lo, hi, rows)
                sys.stderr.write('%s %s %s\n' % (classes[key][0], lo, hi))
            c, lo, hi, _ = classes[key]
            out.append('%-14s %-14s genus 1  %s  rank [%s,%s]' % (p['S'], p['path'], c, lo, '?' if hi is None else hi))
        elif p['genus'] == 2:
            f = p['form']; out.append('%-14s %-14s genus 2  %s' % (p['S'], p['path'], simple_test(f, p['D'])))
    out.append('isogeny classes (the 2-isogeny class of the first curve; for each curve [lo, hi, points returned, effort]):')
    for key, (c, lo, hi, rows) in sorted(classes.items(), key=lambda kv: int(kv[1][0][1:])):
        verdict = ('rank 0' if hi == 0 else 'rank 1 (bounds meet)' if (lo == hi == 1) else
                   'rank %d (bounds meet)' % lo if lo == hi else 'bounds [%d,%s]' % (lo, hi))
        out.append('%s  %s  minimal model %s' % (c, verdict, list(key[0])))
        for cv, r in rows: out.append('    %s  %s' % (list(cv), list(r[:4]) + ([r[4]] if r[4] else [])))
    return out

if __name__ == '__main__':
    argv = sys.argv
    if '--negative' in argv:
        lines = []
        # 1 a forged Jacobian: the invariants of a different quartic must not match ellfromeqn's model
        f = [1, 0, 3, 0, 5]; I, J = IJ([1, 0, 3, 0, 7])
        lines.append('forged Jacobian: %s' % ('refused' if minmodel([0, 0, 0, -27*I, -27*J]) != minmodel(pari('ellfromeqn(y^2 - (%s))' % poly_str(f))) else 'ACCEPTED'))
        # 2 positive control: y^2 = (x^3 + x) (x^3 + 2) ... a split genus 2 Jacobian (even sextic) must not be called simple
        f = [1, 0, 1, 0, 3, 0, 2]                                # 1 + x^2 + 3x^4 + 2x^6, even: Jacobian splits
        lines.append('even sextic (split Jacobian): %s' % simple_test(f, 6))
        # 3 a generic sextic should be found simple
        f = [1, 3, 0, 2, 5, 1, 1]
        lines.append('generic sextic: %s' % simple_test(f, 6))
        # 4 rank control: 37a (rank 1) and 11a (rank 0)
        lines.append('37a1 rank bounds %s, 11a1 rank bounds %s' % (rank_class([0, 0, 1, -1, 0])[:2], rank_class([0, -1, 1, -10, -20])[:2]))
        s = '\n'.join(lines) + '\n'; fn = os.path.join(HERE, 'h6_ranks_negative.txt')
    else:
        m, n = int(argv[1]), int(argv[2])
        s = '\n'.join(run(m, n)) + '\n'; fn = os.path.join(HERE, 'h6_ranks_%d_%d.txt' % (m, n))
    if '--write' in argv: open(fn, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(fn): print('# matches' if open(fn).read() == s else '# DIFFERS from ' + os.path.basename(fn))
