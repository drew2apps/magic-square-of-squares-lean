"""Second descent over K = Q(i) on the phi-coverings of the (627,364) bound-2 curve, edge 1-4 (curve 1 -> curve 4).

Sel^phi(E) has dim 3 = torsion image (dim 1) + 2 excess classes; rank E(K) = number of independent excess classes that
come from K-points.  For an excess class delta the phi-covering is
    C_delta:  N^2 = delta u^4 + a' u^2 e^2 + (b'/delta) e^4,   and  delta N^2 = f1 f2 / delta ... precisely
    (delta N)^2 = delta * f1(u,e) * f2(u,e),  f1 = delta u^2 - A1 e^2,  f2 = delta u^2 - A2 e^2,  A1,2 = a +- 2 beta, beta^2 = b.
A K-point of C_delta gives c in K(S,2) with  f1 = c z1^2,  f2 = delta c z2^2  (the 4-covering D_c).
D_c has good reduction outside S (S: places above 2 N(delta A1 A2 beta)), where a smooth genus-1 curve over F_q has a point
(Hasse), so D_c is everywhere locally soluble iff it is soluble at the places of S.
At each place v of S the set R_v of pairs (class f1, class f2) realised on P^1(K_v) (a zero value of f_i counts as every
class) is computed EXACTLY by recursion on residue discs, with two exact devices:
  * a disc on which v(f_i(x0)) < k (odd v) or v(f_i(x0)) + 5 <= k (v = 1+i) has constant class of f_i;
  * a disc containing a simple root x* of f_i (Hensel: v(f_i(x0)) > 2 v(f_i'(x0)) and k > v(f_i'(x0)) + slack) realises every
    class of f_i together with the constant class of f_j;
  * at split primes p > 1000 the residues where f1 f2 is a unit are not enumerated: the pairs of quadratic characters that
    occur there are read off from Weil's bound (all pairs allowed by the polynomial identities mod p occur once p > 1000).
c survives iff (class_v c, class_v delta c) in R_v for every v in S.
Written by an AI model (Claude); not read by a person.  EXPERIMENT."""
import sys, itertools
from math import isqrt
import cypari2; pari = cypari2.Pari()

def gm(x, y): return (x[0]*y[0] - x[1]*y[1], x[0]*y[1] + x[1]*y[0])
def gp_(x, k):
    r = (1, 0)
    for _ in range(k): r = gm(r, x)
    return r

# ---------- places
def sqrt_m1(p, N):
    g = 2
    while pow(g, (p - 1)//2, p) != p - 1: g += 1
    r = pow(g, (p - 1)//4, p); mod = p
    while mod < p**N:
        mod = min(mod*mod, p**N); r = (r - (r*r + 1) * pow(2*r, -1, mod)) % mod
    return r

SQ8 = set()
for a in range(8):
    for b in range(8):
        if (a + b) % 2:
            z = ((a*a - b*b) % 8, (2*a*b) % 8); SQ8.add(z); SQ8.add(((z[0]+4) % 8, (z[1]+4) % 8))
def ucls2(x):
    return min(((x[0]*s[0] - x[1]*s[1]) % 8, (x[0]*s[1] + x[1]*s[0]) % 8) for s in SQ8)

class Split:
    def __init__(s, p, sign, N=None):
        s.p = p; s.N = N or max(10, 400 // p.bit_length()); s.M = p**s.N
        r = sqrt_m1(p, s.N); s.r = r if sign > 0 else s.M - r; s.name = f"split {p} {'+' if sign>0 else '-'}"
    def emb(s, z): return (z[0] + z[1]*s.r) % s.M
    def val(s, x):
        x %= s.M
        if x == 0: return None
        e = 0
        while x % s.p == 0: x //= s.p; e += 1
        return e
    def cls(s, x):
        x %= s.M; e = s.val(x)
        if e is None or e > s.N - 3: return None
        u = x // s.p**e
        return (e % 2, 1 if pow(u % s.p, (s.p - 1)//2, s.p) == 1 else -1)
    def all_classes(s): return [(e, c) for e in (0, 1) for c in (1, -1)]

class Inert:
    def __init__(s, p, N=8):
        s.p = p; s.N = N; s.M = p**N; s.name = f"inert {p}"
    def emb(s, z): return (z[0] % s.M, z[1] % s.M)
    def val(s, x):
        a, b = x[0] % s.M, x[1] % s.M
        if a == 0 and b == 0: return None
        e = 0
        while a % s.p == 0 and b % s.p == 0: a //= s.p; b //= s.p; e += 1
        return e
    def cls(s, x):
        e = s.val(x)
        if e is None or e > s.N - 3: return None
        a, b = x[0] // s.p**e, x[1] // s.p**e
        return (e % 2, 1 if pow((a*a + b*b) % s.p, (s.p - 1)//2, s.p) == 1 else -1)
    def all_classes(s): return [(e, c) for e in (0, 1) for c in (1, -1)]

class Ram:
    def __init__(s, N=40):
        s.N = N; s.M = 2**N; s.name = "ram 2"   # elements mod 2^N = pi^(2N)
    def emb(s, z): return (z[0] % s.M, z[1] % s.M)
    def val(s, x):
        a, b = x[0] % s.M, x[1] % s.M
        if a == 0 and b == 0: return None
        n = a*a + b*b
        # v_pi(x) = v_2(N x) as long as below precision
        e = (n & -n).bit_length() - 1
        return e
    def cls(s, x):
        a, b = x[0] % s.M, x[1] % s.M
        if a == 0 and b == 0: return None
        e = s.val((a, b))
        if e + 7 > 2*s.N: return None
        for _ in range(e): a, b = (a + b)//2, (b - a)//2      # divide by (1+i): (a+bi)/(1+i) = ((a+b) + (b-a)i)/2
        return (e % 2, ucls2((a % 8, b % 8)))
    def all_classes(s):
        us = sorted(set(ucls2((a, b)) for a in range(8) for b in range(8) if (a + b) % 2))
        assert len(us) == 8
        return [(e, u) for e in (0, 1) for u in us]

def cmul(pl, c1, c2):
    """product of two classes at pl, via representatives"""
    raise NotImplementedError

# ---------- arithmetic in O_v / pi^big, uniform interface
class Ring:
    """O_v elements: split -> int mod M; inert/ram -> Gaussian pair mod M."""
    def __init__(s, pl):
        s.pl = pl; s.gauss = not isinstance(pl, Split)
    def mul(s, x, y):
        if s.gauss: return ((x[0]*y[0] - x[1]*y[1]) % s.pl.M, (x[0]*y[1] + x[1]*y[0]) % s.pl.M)
        return x*y % s.pl.M
    def add(s, x, y):
        if s.gauss: return ((x[0]+y[0]) % s.pl.M, (x[1]+y[1]) % s.pl.M)
        return (x+y) % s.pl.M
    def neg(s, x):
        if s.gauss: return ((-x[0]) % s.pl.M, (-x[1]) % s.pl.M)
        return (-x) % s.pl.M
    def zero(s): return (0, 0) if s.gauss else 0
    def one(s): return (1, 0) if s.gauss else 1
    def is_zero(s, x): return x == s.zero()

def residues(pl):
    """Representatives of O_v / pi, and the uniformiser."""
    if isinstance(pl, Split): return None, pl.p
    if isinstance(pl, Inert): return [(a, b) for a in range(pl.p) for b in range(pl.p)], (pl.p, 0)
    return [(0, 0), (1, 0)], (1, 1)

def slack(pl): return 5 if isinstance(pl, Ram) else 1

def realised(pl, delta, A1, A2):
    """Exact set of pairs (class f1, class f2) on P^1(K_v), with 'ANY' for a zero / all-classes entry.  Returns a set of
    pairs (c1, c2) where c1, c2 are classes (all classes expanded)."""
    R = Ring(pl); E = pl.emb
    d, a1, a2 = E(delta), E(A1), E(A2)
    allc = pl.all_classes()
    out = set()
    sl = slack(pl)
    res, unif = residues(pl)
    piE = E(unif) if not isinstance(pl, Split) else pl.p
    def pi_pow(k):
        if isinstance(pl, Split): return pow(pl.p, k, pl.M)
        return gp_m(piE, k)
    def gp_m(x, k):
        r = R.one()
        for _ in range(k): r = R.mul(r, x)
        return r
    # chart 1: e = 1, u = x;  f_i = d x^2 - a_i ; f_i' = 2 d x
    # chart 2: u = 1, e = x in pi O;  f_i = d - a_i x^2 ; f_i' = -2 a_i x
    def fvals(chart, x):
        x2 = R.mul(x, x)
        if chart == 1:
            return [R.add(R.mul(d, x2), R.neg(a1)), R.add(R.mul(d, x2), R.neg(a2))], [R.mul(R.add(d, d), x)] * 2
        return [R.add(d, R.neg(R.mul(a1, x2))), R.add(d, R.neg(R.mul(a2, x2)))], [R.neg(R.mul(R.add(a1, a1), x)), R.neg(R.mul(R.add(a2, a2), x))]
    def val(z):
        v = pl.val(z); return 10**9 if v is None else v
    count = [0]
    def rec(chart, x0, k):
        count[0] += 1
        if count[0] > 20_000_000: raise RuntimeError("too many discs at " + pl.name)
        fs, ds = fvals(chart, x0)
        vs = [val(f) for f in fs]
        # f(x0 + pi^k y) - f(x0) = f'(x0) pi^k y + lead * pi^(2k) y^2 has valuation >= min(v(f'(x0)) + k, v(lead) + 2k)
        # (k counted in the uniformiser; at 1+i, v counts powers of 1+i).  The class of f is constant on the disc when that
        # bound is >= v(f(x0)) + sl (sl = 1 at odd places, 5 at 1+i).
        leads = [val(d), val(d)] if chart == 1 else [val(a1), val(a2)]
        det = [vs[i] + sl <= min(val(ds[i]) + k, leads[i] + 2*k) for i in range(2)]
        if det[0] and det[1]:
            out.add((pl.cls(fs[0]), pl.cls(fs[1]))); return
        # root device
        for i in range(2):
            if det[i]: continue
            mu = val(ds[i]); vl = leads[i]; v2 = 2 if isinstance(pl, Ram) else 0
            j = 1 - i
            # (a) Hensel: a root x* with v(x* - x0) >= v(f(x0)) - mu > mu;  (b) x* in the disc: v(f(x0)) - mu >= k;
            # (c) f(x) = (x - x*)(f'(x*) + lead (x - x*)) with the second factor of constant class: vl + k >= mu + sl;
            # (d) class f'(x*) = class f'(x0): v2 + vl + v(f(x0)) - mu >= mu + sl.  Then f_i takes every class on the disc.
            if (mu < 10**8 and vs[i] > 2*mu and vs[i] - mu >= k and vl + k >= mu + sl
                    and v2 + vl + vs[i] - mu >= mu + sl and det[j]):
                cj = pl.cls(fs[j])
                for c in allc: out.add((c, cj) if i == 0 else (cj, c))
                return
        # subdivide
        step = pi_pow(k)
        if isinstance(pl, Split):
            if pl.p > 1000 and k == 0:
                generic_split(chart); return
            for t in range(pl.p): rec(chart, (x0 + t*step) % pl.M, k + 1)
        else:
            for t in res: rec(chart, R.add(x0, R.mul(step, (t[0] % pl.M, t[1] % pl.M))), k + 1)
    def generic_split(chart):
        # chart is 1 here (k = 0 means all of Z_p).  Enumerate only residues t mod p where f1(t) f2(t) = 0 mod p; the other
        # residues contribute unit-class pairs given by Weil's bound.
        p = pl.p
        F = [[d % p, 0, (-a1) % p], [d % p, 0, (-a2) % p]]   # coefficients x^2, x, 1
        special = set()
        for co in F:
            A, B, C = co
            if A == 0:
                if B == 0: continue      # constant: no root
                special.add((-C * pow(B, -1, p)) % p)
            else:
                disc = (B*B - 4*A*C) % p
                if disc == 0: special.add((-B * pow(2*A, -1, p)) % p)
                elif pow(disc, (p - 1)//2, p) == 1:
                    sq = int(pari(f"lift(sqrt(Mod({disc},{p})))"))
                    for sg in (sq, -sq): special.add(((-B + sg) * pow(2*A, -1, p)) % p)
        for t in special: rec(chart, t, 1)
        # generic residues: unit values; which character pairs occur?
        def const_times_square(co):
            A, B, C = co
            if A == 0 and B == 0: return (True, C)
            if A != 0 and (B*B - 4*A*C) % p == 0: return (True, A)
            return (False, None)
        s1, s2 = const_times_square(F[0]), const_times_square(F[1])
        # proportional? f1 = lam f2 mod p
        prop = None
        for lam in [None]:
            pass
        A1_, B1_, C1_ = F[0]; A2_, B2_, C2_ = F[1]
        # check rank of [[A1,B1,C1],[A2,B2,C2]] mod p
        minors = [(A1_*B2_ - A2_*B1_) % p, (A1_*C2_ - A2_*C1_) % p, (B1_*C2_ - B2_*C1_) % p]
        if all(mm == 0 for mm in minors):
            # proportional: lam = ratio of a nonzero coefficient
            for x_, y_ in ((A1_, A2_), (B1_, B2_), (C1_, C2_)):
                if y_ % p: prop = x_ * pow(y_, -1, p) % p; break
        chi = lambda z: 1 if pow(z % p, (p - 1)//2, p) == 1 else -1
        pairs = set()
        for e1 in (1, -1):
            for e2 in (1, -1):
                if s1[0] and chi(s1[1]) != e1: continue
                if s2[0] and chi(s2[1]) != e2: continue
                if prop is not None and e1 != e2 * chi(prop): continue
                pairs.add(((0, e1), (0, e2)))
        assert p > 1000
        # Weil: with no identity forcing a character, each allowed pair has (p - O(sqrt p))/4 > 0 residues (p > 1000)
        for pr in pairs: out.add(pr)
    rec(1, R.zero(), 0)
    rec(2, R.zero(), 1)            # e in pi O, u = 1
    return out, count[0]

def gauss_primes(N):
    """Gaussian prime generators above the rational prime factors of N."""
    out = []
    for p in sorted(int(q) for q in pari.factor(N)[0]):
        if p == 2: continue
        if p % 4 == 3: out.append((p, 0))
        else:
            for a in range(1, isqrt(p) + 1):
                b2 = p - a*a; b = isqrt(b2)
                if b*b == b2: g = (a, b); break
            out += [g, (g[0], -g[1])]
    return out

def places_of(primes):
    pls = []
    for p in primes:
        if p == 2: pls.append(Ram())
        elif p % 4 == 3: pls.append(Inert(p))
        else: pls += [Split(p, +1), Split(p, -1)]
    return pls

def nrm(z): return z[0]**2 + z[1]**2

def run(delta, a, beta, label=""):
    A1 = (a[0] + 2*beta[0], a[1] + 2*beta[1]); A2 = (a[0] - 2*beta[0], a[1] - 2*beta[1])
    N = 2 * nrm(delta) * nrm(A1) * nrm(A2) * nrm(beta)
    primes = sorted(set([2]) | set(int(q) for q in pari.factor(N)[0]))
    basis = [(0, 1), (1, 1)] + gauss_primes(N)
    pls = places_of(primes)
    print(f"delta = {delta} {label}: primes {primes}, |basis| = {len(basis)}", flush=True)
    Rs = []
    for pl in pls:
        Rv, cnt = realised(pl, delta, A1, A2)
        Rs.append((pl, Rv))
        print(f"   {pl.name}: {len(Rv)} realised class pairs ({cnt} discs)", flush=True)
    # classes of basis elements and of delta at each place
    def cls_prod(pl, zs):
        z = (1, 0)
        for w in zs: z = gm(z, w)
        return pl.cls(pl.emb(z) if isinstance(pl, Split) else pl.emb(z))
    surv = []
    g = len(basis)
    # precompute per place the class of each basis element as a function; product classes need representatives, so we
    # compute the class of the actual Gaussian integer c (exact, sizes are fine).
    for mask in range(1 << g):
        c = (1, 0)
        for k in range(g):
            if mask >> k & 1: c = gm(c, basis[k])
        dc = gm(delta, c); ok = True
        for pl, Rv in Rs:
            k1 = pl.cls(pl.emb(c)); k2 = pl.cls(pl.emb(dc))
            assert k1 is not None and k2 is not None
            if (k1, k2) not in Rv: ok = False; break
        if ok: surv.append(c)
    print(f"   surviving c: {len(surv)} of {1 << g}" + (f"; e.g. {surv[:4]}" if surv else ""), flush=True)
    return surv

if __name__ == '__main__':
    a = (-1919397480000, 951739973184)
    beta = (833408319744, -547980882500)        # beta^2 = b (curve 1 model of edge 1-4)
    b = (394286379833038533975536, -913383653072318588160000)
    assert gm(beta, beta) == b
    tor = gm((263, 0), (0, 991))                 # the torsion class 263*991 i
    for delta, lab in (((145, 0), "(excess 1)"), ((263, 0), "(excess 2)"), ((145*263, 0), "(excess 1+2)"),
                       ((1, 0), "(control: trivial class, must survive)"), (tor, "(control: torsion class, must survive)")):
        run(delta, a, beta, lab)
