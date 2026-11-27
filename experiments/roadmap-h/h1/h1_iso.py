"""H1, part 8: the partner isomorphism, written down and checked (EXPERIMENT, written by an AI model; not read by a
person; nothing here changes the census).
For a ratio (m, n) the four Gaussian forms are, with t = u / v (h1_pairs.forms):
    z_0 = m v + i n u,  z_8 = n v + i m u,  z_L = (m + n)(u - v) + i (m - n)(u + v),  z_U = (m - n)(u - v) + i (m + n)(u + v),
and C is the curve where all four norms are squares.  For the partner (m', n') = ((m + n)/g, (m - n)/g), the script
searches the integer matrices M with entries in {-1, 0, 1} and det M != 0, and the bijections s of {0, 8, L, U}, for
which z'_{s(j)}(M (u, v)) = lam_j z_j(u, v) identically in (u, v), with lam_j in Q(i) (exact Gaussian arithmetic on
the coefficients).  Each lam_j must have square norm.  Then a rational point of C' is M applied to a rational point of
C, and its tuple is d'_{s(j)} = sqf(c_j d_j e), where c_j is the class of lam_j and e is one class for all four forms
(the class of the denominator that brings M (u, v) back to (t', 1)).
CHECKS, per pair and in each direction (C -> C' and C' -> C):
  (a) the identity, for every j;
  (b) every lam_j has square norm;
  (c) TRANSFER: for every known tuple T of C and every class e for which the image lies in the tuple space of C'
      (each d'_k divides the product of its support), the image is NOT an extra tuple of C'.  Every rational point of
      C' has a tuple in known' + extra' (the descent), so if every rational point of C has a known tuple, so does
      every rational point of C'.  This is what turns part 4 to 7's step from one side of a pair to the other into a
      computed statement.
  (d) consistency: the images of known (for some e) include all of known', and of extra include all of extra'.
Output h1_iso.txt.  usage: h1_iso.py [--negative] [--write]   (h1_iso_negative.txt)
The repository is found through MSS_REPO, or else ../../.."""
import sys, os
from itertools import product, permutations
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv; NEG = '--negative' in sys.argv
_argv = sys.argv; sys.argv = [sys.argv[0]]
import h1_rank0_exact as R
sys.argv = _argv
RI = R.RI
OUT = []
def say(s=''): OUT.append(s)

def gm(a, b): return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])
def lin(m, n):
    """z_j = A_j u + B_j v, as Gaussian coefficient pairs (A_j, B_j)"""
    return {'0': ((0, n), (m, 0)), '8': ((0, m), (n, 0)),
            'L': ((m + n, m - n), (-(m + n), m - n)), 'U': ((m - n, m + n), (-(m - n), m + n))}

def check_forms(m, n):
    """lin() against h1_pairs.forms at a few t (the same forms)"""
    from fractions import Fraction as Fr
    for t in (Fr(2), Fr(-3, 5), Fr(7, 3)):
        f = R.HP.forms(m, n, t)
        for j, (A, B) in lin(m, n).items():
            z = (A[0] * t + B[0], A[1] * t + B[1])
            assert z == f[j], 'lin() is not h1_pairs.forms'

def isos(m, n, mp, np_, forced=None):
    F, G = lin(m, n), lin(mp, np_); out = []
    for al, be, ga, de in product((-1, 0, 1), repeat=4):
        if al * de - be * ga == 0: continue
        for perm in ([forced] if forced else permutations('08LU')):
            lam = {}
            for j, k in zip('08LU', perm):
                (Ak, Bk), (Aj, Bj) = G[k], F[j]
                cu = (Ak[0] * al + Bk[0] * ga, Ak[1] * al + Bk[1] * ga); cv = (Ak[0] * be + Bk[0] * de, Ak[1] * be + Bk[1] * de)
                if gm(cu, Bj) != gm(cv, Aj): break                      # (a) fails for this j
                num, den = (cu, Aj) if Aj != (0, 0) else (cv, Bj)
                lam[j] = (num, den)
            else: out.append(((al, be, ga, de), perm, lam))
    return out

def lam_class(num, den):
    """the class c of lam = num / den in Q(i) (lam = c w^2), lam of square norm (b); exact"""
    from fractions import Fraction as Fr
    Nd = den[0] ** 2 + den[1] ** 2
    z = gm(num, (den[0], -den[1])); z = (Fr(z[0], Nd), Fr(z[1], Nd))
    Nz = z[0] ** 2 + z[1] ** 2
    assert R.HP.isq(Nz), 'lam does not have square norm'
    return R.HP.gauss_class(z)

def sqf(x):
    d = 1
    for p in RI.factor_int(x):
        e = 0
        while x % p == 0: x //= p; e += 1
        if e % 2: d *= p
    return d

def space_ok(T, m, n):
    P0 = set([2]) | set(RI.factor_int(m)) | set(RI.factor_int(n)); PL = set([2]) | set(RI.factor_int(abs(m - n))) | set(RI.factor_int(m + n))
    for d, P in zip(T, (P0, P0, PL, PL)):
        if any(q not in P for q in RI.factor_int(d)): return False
    return True

def images(T, s, cls, m2, n2):
    """the images sqf(c_j d_j e) placed at s(j), for every class e that keeps them in C2's tuple space"""
    idx = {'0': 0, '8': 1, 'L': 2, 'U': 3}
    # a prime of e that divides no c_j d_j appears in all four images, so it lies in all four supports of C2
    P0 = set([2]) | set(RI.factor_int(m2)) | set(RI.factor_int(n2)); PL = set([2]) | set(RI.factor_int(abs(m2 - n2))) | set(RI.factor_int(m2 + n2))
    ps = sorted((P0 & PL) | set(q for j in '08LU' for q in RI.factor_int(cls[j] * T[idx[j]])))
    out = set()
    for mask in range(1 << len(ps)):
        e = 1
        for i, p in enumerate(ps):
            if mask >> i & 1: e *= p
        img = [0] * 4
        for j in '08LU': img[idx[s['08LU'.index(j)]]] = sqf(cls[j] * T[idx[j]] * e)
        img = tuple(img)
        if space_ok(img, m2, n2): out.add(img)
    return out

def transfer(m, n, known, extra, m2, n2, known2, extra2, forced=None):
    R.HP.set_base(set(q for x in (2, m, n, abs(m - n), m + n, m2, n2, abs(m2 - n2), m2 + n2) for q in RI.factor_int(x)) | {2})
    sols = isos(m, n, m2, n2, forced)
    if not sols: return None, 'no isomorphism of this shape'
    M, s, lam = sols[0]
    cls = {j: lam_class(*lam[j]) for j in '08LU'}
    imk = set().union(*(images(T, s, cls, m2, n2) for T in known))
    ime = set().union(*(images(T, s, cls, m2, n2) for T in extra))
    bad = sorted(imk & set(extra2))
    ok = not bad
    cons = set(known2) <= imk and set(extra2) <= ime
    txt = 'M = %s, s = (0 8 L U) -> (%s), classes of lam = %s; known -> %d tuples of the space, none extra: %s; consistency (known2 and extra2 reached): %s' % (
        M, ' '.join(s), tuple(cls[j] for j in '08LU'), len(imk), 'yes' if ok else 'NO %s' % bad[:2], 'yes' if cons else 'no')
    return ok, txt

def main():
    say('H1: the partner isomorphism, written down and checked. Written by an AI model; not read by a person. Experiment, not counted.')
    out, part = R.ratios(); seen = set(); good = 0
    for a in sorted(out):
        b = part.get(a)
        if b is None or (b, a) in seen or (a, b) in seen: continue
        seen.add((a, b)); check_forms(*a)
        ka, ea = out[a]; kb, eb = out[b]
        oab, tab = transfer(*a, ka, ea, *b, kb, eb)
        oba, tba = transfer(*b, kb, eb, *a, ka, ea)
        say('%s <-> %s' % (a, b)); say('  %s -> %s: %s' % (a, b, tab)); say('  %s -> %s: %s' % (b, a, tba))
        if oab and oba: good += 1
    say('pairs: %d; the transfer holds in both directions: %d' % (len(seen), good))

def negative():
    say('H1 part 8, negative controls. Written by an AI model; not read by a person.')
    out, part = R.ratios()
    a = (11875080, 6497849); b = part[a]; ka, ea = out[a]; kb, eb = out[b]
    ok, txt = transfer(*a, ka, ea, b[0], b[1] + 2, kb, eb)
    say('forged partner (m\', n\' + 2): %s' % ('refused (%s)' % txt if ok is None else 'ACCEPTED'))
    assert ok is None
    ok, txt = transfer(*a, ka, ea, *b, kb, eb, forced=('L', 'U', '8', '0'))
    say('forged bijection (0 8 L U) -> (L U 8 0): %s' % ('refused (%s)' % txt if ok is None else 'ACCEPTED'))
    assert ok is None
    R.HP.set_base(set(q for x in (2,) + a + b + (a[0] - a[1], a[0] + a[1], b[0] - b[1], b[0] + b[1]) for q in RI.factor_int(x)))
    M, s, lam = isos(*a, *b)[0]; cls = {j: lam_class(*lam[j]) for j in '08LU'}
    planted = sorted(images(ka[0], s, cls, *b))[0]
    ok, txt = transfer(*a, ka, ea, *b, kb, list(eb) + [planted])
    say('an image of a known tuple planted among the partner\'s extra tuples: %s' % ('refused' if ok is False else 'ACCEPTED'))
    assert ok is False

if __name__ == '__main__':
    negative() if NEG else main()
    s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'h1_iso_negative.txt' if NEG else 'h1_iso.txt')
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
