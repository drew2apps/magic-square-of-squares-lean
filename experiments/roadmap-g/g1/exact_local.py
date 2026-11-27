"""The local test of tuple_descent.py at an odd prime p, computed exactly (EXPERIMENT, written by an AI model; not read by a
person).

tuple_descent.Local.allowed() cuts P^1(Z_p) into discs, p children at a time, down to a depth K, and lets every class pass
where the class is not constant at depth K.  That is sound (the set can only be too large) but it costs about p^K steps, so
the repository caps the primes it tests, and a bad prime above the cap is left out (COMPLETE.md, section 7, test 6).

Here the same set is computed with no depth cap and in time independent of p:

  * On a disc x0 + p^k Z_p each form is a linear function f(x) = a + b x of the disc variable.  If v(a) < v(b) the class of
    f is constant on the disc.  Otherwise f = p^v(b) (a' + b' x) and only the residues x of the roots of a' + b' x mod p
    need a smaller disc; on every other residue the class is fixed by Legendre symbols of linear functions of x
    (p = 1 mod 4: the two components in Z_p x Z_p) or of their norms (p = 3 mod 4).
  * Which sign patterns of those symbols occur on the residues that are not roots is decided by listing all residues when
    p is small, and otherwise by Weil's bound |sum_x chi(f(x))| <= (d - 1) sqrt(p) for f with d distinct roots that is not
    a constant times a square: when p is above the bound computed for that system, every pattern occurs.  Functions with
    the same root are grouped first, so no product in the bound is a square.
  * When p = 1 mod 4, a component has a root in Z_p, so the descent toward it never stops.  When all functions that vanish
    at a residue have exactly the same root x* (tested exactly in Z[i]), the disc around x* is described in closed form:
    f = b (x - x*) with x - x* = p^e w, e >= 1 and w a unit, so the classes are those of b p^e w for e of either parity and
    w of either symbol.  Roots that only agree to some depth separate after finitely many steps.  When p = 3 mod 4 a form
    has no root on the line Z_p (it would be a rational zero of the form), so the descent is finite.

The labels are those of tuple_descent.Local (dlabel), and the set is closed under the unit of Z_p as there.  The class is
undefined only where a form vanishes, which no coprime integer point does.  Positive controls: test_exact_local.py.
"""
import sys, os
from math import isqrt, comb
from itertools import product
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'scripts', 'cover_descent'))
import tuple_descent as TD

BRUTE = 5 * 10**6                  # above this, a pattern set must come from the Weil bound (asserted)
STATS = {'brute': 0, 'weil': 0, 'discs': 0, 'tails': 0}

def vp(x, p):
    if x == 0: return None
    v = 0
    while x % p == 0: x //= p; v += 1
    return v
def leg(a, p): return 0 if pow(a % p, (p - 1)//2, p) == 1 else 1          # the convention of tuple_descent: 0 = square
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gconj(a): return (a[0], -a[1])
def gsub(a, b): return (a[0] - b[0], a[1] - b[1])

_QR = {}
def qr_table(p):
    if p not in _QR:
        t = np.zeros(p, dtype=np.int8); x = np.arange(p, dtype=np.int64); t[(x*x) % p] = 1; t[0] = 0; _QR.clear(); _QR[p] = t
    return _QR[p]                                     # 1 for a nonzero square, 0 otherwise

def weil_all(p, degs, excluded):
    """True if Weil's bound shows that every sign pattern of g independent characters occurs on F_p minus `excluded`
    points; degs[i] = number of distinct roots of the i-th function (products of distinct ones have the summed count)."""
    g = len(degs); s = isqrt(p) + 1; total = 0
    for mask in range(1, 1 << g):
        d = sum(degs[i] for i in range(g) if mask >> i & 1)
        total += (d - 1) * s + excluded
    return p - excluded > total

class ExactLocal(TD.Local):
    """tuple_descent.Local with allowed() computed exactly at an odd prime (same labels)."""
    def __init__(self, forms, p, K=None, N=None):
        assert p % 2 == 1
        self.forms, self.p, self.K = forms, p, 1                               # as TD.Local.__init__, without its O(p) search
        self.kind = 1 if p % 4 == 1 else 3                                     # for a square root of -1 (sqrt_m1 below)
        self.N = N or 60
        if self.kind == 1:
            M = p**self.N; r = next(x for x in range(2, min(p, 10**6)) if (x*x + 1) % p == 0) if p < 10**6 else sqrt_m1(p)
            for _ in range(8 + self.N.bit_length()): r = (r - (r*r + 1)*pow(2*r, -1, M)) % M
            assert (r*r + 1) % M == 0
            self.iota, self.M = r, M

    # ---- p = 1 mod 4: 2J functions over Z_p (the components of the forms) -------------------------------------------
    def comp(self, z, s): return (z[0] + s*self.iota*z[1]) % self.M
    def same_root(self, f, g):
        """f = (A, B, s): the function comp_s(A) + comp_s(B) x.  True iff f and g have exactly the same root (exact test in Z[i])."""
        (A1, B1, s1), (A2, B2, s2) = f, g
        if s1 == s2: G = gsub(gmul(A1, B2), gmul(A2, B1))
        else: G = gsub(gmul(A1, gconj(B2)), gmul(gconj(A2), B1))
        return G == (0, 0)

    def allowed(self):
        p = self.p; out = set()
        charts = [([(be, al) for al, be in self.forms], 0)]                       # (t, 1): gamma = beta + alpha t
        charts.append(([(al, (be[0]*p, be[1]*p)) for al, be in self.forms], 1))   # (1, p x): gamma = alpha + beta p x
        for AB, _ in charts:
            (self._disc1 if self.kind == 1 else self._disc3)(AB, out, 0)
        if self.kind == 1: out |= set(tuple((l[0], 1 - l[1]) for l in t) for t in out)
        return out

    def _val(self, a):
        v = vp(a, self.p)
        if v is None or v >= self.N - 2: raise OverflowError('p-adic precision')
        return v

    def _disc1(self, AB, out, depth):
        p = self.p; STATS['discs'] += 1
        assert depth < 400
        J = len(AB); F = []                       # per function: ('c', label) or ('v', vb, bbar, root, (A, B, s))
        for j, (A, B) in enumerate(AB):
            for s in (1, -1):
                a, b = self.comp(A, s), self.comp(B, s); va, vb = self._val(a), self._val(b)
                if va < vb: F.append(('c', (va % 2, leg(a // p**va, p))))
                else:
                    a1, b1 = (a // p**vb) % p, (b // p**vb) % p
                    F.append(('v', vb, b1, (-a1 * pow(b1, -1, p)) % p, (A, B, s)))
        roots = sorted(set(f[3] for f in F if f[0] == 'v'))
        # generic residues: sign pattern of chi(x - r) over the distinct roots r
        pats = self._patterns1(roots)
        for eps in pats:
            e = dict(zip(roots, eps)); lab = []
            for f in F:
                if f[0] == 'c': lab.append(f[1])
                else: lab.append((f[1] % 2, leg(f[2], p) ^ e[f[3]]))
            self._emit1(lab, J, out)
        # residues that are roots
        for x0 in roots:
            G = [f for f in F if f[0] == 'v' and f[3] == x0]
            if all(self.same_root(G[0][4], g[4]) for g in G[1:]):
                STATS['tails'] += 1
                fixed = []
                for f in F:
                    if f[0] == 'c': fixed.append(f[1])
                    elif f[3] == x0: fixed.append(None)
                    else: fixed.append((f[1] % 2, leg((f[2] * (x0 - f[3])) % p, p)))
                for e in (0, 1):
                    for w in (0, 1):
                        lab = [fx if fx is not None else ((f[1] + 1 + e) % 2, leg(f[2], p) ^ w) for fx, f in zip(fixed, F)]
                        self._emit1(lab, J, out)
            else:
                AB2 = [((A[0] + B[0]*x0, A[1] + B[1]*x0), (B[0]*p, B[1]*p)) for A, B in AB]
                self._disc1(AB2, out, depth + 1)

    def _emit1(self, lab, J, out):
        t = []
        for j in range(J):
            if lab[2*j] != lab[2*j + 1]: return
            t.append(lab[2*j])
        out.add(tuple(t))

    def _patterns1(self, roots):
        p, g = self.p, len(roots)
        if g == 0: return [()] if p > 0 else []
        if p - g <= 0: return []
        if weil_all(p, [1]*g, g):
            STATS['weil'] += 1; return list(product((0, 1), repeat=g))
        assert p <= BRUTE, 'p too large to list and below the Weil bound'
        STATS['brute'] += 1
        q = qr_table(p); x = np.arange(p, dtype=np.int64)
        mask = np.ones(p, dtype=bool)
        for r in roots: mask[r] = False
        cols = [(1 - q[(x - r) % p])[mask] for r in roots]           # leg convention: 0 = square
        code = np.zeros(int(mask.sum()), dtype=np.int64)
        for i, c in enumerate(cols): code |= c.astype(np.int64) << i
        return [tuple((int(c) >> i) & 1 for i in range(g)) for c in np.unique(code)]

    # ---- p = 3 mod 4: J Gaussian functions; label (v mod 2, 0), alive iff the norm of the unit part is a square -------
    def _gval(self, z):
        v = min(x for x in (vp(z[0], self.p), vp(z[1], self.p)) if x is not None)
        assert v < 400
        return v
    def _disc3(self, AB, out, depth):
        p = self.p; STATS['discs'] += 1
        assert depth < 400
        F = []
        for A, B in AB:
            va, vb = self._gval(A), self._gval(B)
            if va < vb:
                u = (A[0] // p**va % p, A[1] // p**va % p)
                if leg(u[0]**2 + u[1]**2, p): return                      # dead on the whole disc
                F.append(('c', va % 2))
            else:
                a = (A[0] // p**vb % p, A[1] // p**vb % p); b = (B[0] // p**vb % p, B[1] // p**vb % p)
                nb = (b[0]**2 + b[1]**2) % p; inv = pow(nb, -1, p)
                rho = gmul((-a[0] % p, -a[1] % p), (b[0]*inv % p, -b[1]*inv % p)); rho = (rho[0] % p, rho[1] % p)
                F.append(('v', vb, nb, rho))
        froots = sorted(set(f[3][0] for f in F if f[0] == 'v' and f[3][1] == 0))
        # generic residues: forms with a root in F_p have chi(N) = chi(N(b)) there; the others have chi(N(b) q_rho(x))
        if all(leg(f[2], p) == 0 for f in F if f[0] == 'v' and f[3][1] == 0):
            groups = {}
            ok = True
            for f in F:
                if f[0] == 'v' and f[3][1] != 0:
                    key = (f[3][0], min(f[3][1], p - f[3][1]))           # {rho, conj rho}
                    want = leg(f[2], p)                                  # need chi(q) = chi(N(b))
                    if groups.setdefault(key, want) != want: ok = False
            if ok and self._exists3(list(groups.items()), froots):
                out.add(tuple((f[1], 0) if f[0] == 'c' else (f[1] % 2, 0) for f in F))
        for x0 in froots:
            AB2 = [((A[0] + B[0]*x0, A[1] + B[1]*x0), (B[0]*p, B[1]*p)) for A, B in AB]
            self._disc3(AB2, out, depth + 1)

    def _exists3(self, groups, froots):
        """Is there x in F_p, not in froots, with leg(q_rho(x)) = want for every group (rho, want)?"""
        p, g = self.p, len(groups)
        if g == 0: return p > len(froots)
        if weil_all(p, [2]*g, len(froots)):
            STATS['weil'] += 1; return True
        assert p <= BRUTE, 'p too large to list and below the Weil bound'
        STATS['brute'] += 1
        q = qr_table(p); x = np.arange(p, dtype=np.int64); ok = np.ones(p, dtype=bool)
        for r in froots: ok[r] = False
        for (r0, r1), want in groups:
            val = ((x - r0)**2 + r1*r1) % p                               # (x - rho)(x - conj rho) = (x - r0)^2 + r1^2
            ok &= ((1 - q[val]) == want) & (val != 0)
        return bool(ok.any())

def sqrt_m1(p):
    import cypari2
    P = cypari2.Pari(); return int(P.lift(P.sqrt(P.Mod(-1, p))))

class MixedLocal(TD.Local):
    """The repository's Local at p = 2, the exact one at odd p."""
    def __new__(cls, forms, p, K):
        return TD.Local(forms, p, K) if p == 2 else ExactLocal(forms, p)

def allowed_at(forms, p, K2=12):
    """The allowed set at p: exact at odd p (precision raised if needed), the repository's at p = 2 with depth K2."""
    if p == 2: L = TD.Local(forms, 2, K2); return L, L.allowed()
    N = 60
    while True:
        try: L = ExactLocal(forms, p, N=N); return L, L.allowed()
        except OverflowError: N *= 2

def descend_exact(forms, supports, known, primes):
    """tuple_descent.descend with the local test only at `primes` (exact at odd primes).  Returns (tuples left, primes that cut)."""
    divs = []
    for P in supports:
        ds = [1]
        for q in P: ds += [x*q for x in ds]
        divs.append(ds)
    loc = []
    for p in sorted(primes):
        L, al = allowed_at(forms, p)
        for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('a known point is rejected', p, kt)
        if len(al) < len(L.labels())**len(forms): loc.append((p, L, al))
    return join(divs, loc), [p for p, L, al in loc]

def join(divs, loc):
    """The tuples (d_1, ..., d_J), d_j in divs[j], allowed at every prime of loc (the join of tuple_descent.descend)."""
    J = len(divs); S = [()]
    for j in range(J):
        new = []
        projs = [(L, set(a[:j + 1] for a in al), {d: L.dlabel(d) for d in divs[j]}) for p, L, al in loc]
        cache = [dict() for _ in loc]
        for s in S:
            keys = []
            for idx, (L, proj, lab) in enumerate(projs):
                if s not in cache[idx]: cache[idx][s] = tuple(L.dlabel(d) for d in s)
                keys.append(cache[idx][s])
            for d in divs[j]:
                if all(keys[idx] + (lab[d],) in proj for idx, (L, proj, lab) in enumerate(projs)): new.append(s + (d,))
        S = new
    return sorted(S)
