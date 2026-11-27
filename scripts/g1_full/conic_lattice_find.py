"""FINDER, not certifying (written by an AI model; not read by a person): a point of the conic
x^2 - A y^2 - B z^2 = 0 over a quadratic field L = Q(sqrt D0), A, B integral in L.

A Holzer-type lattice: the coordinates (x, y, z) in Z^6 (integral basis of O_L), cut down at every prime P of A B
so that Q(v) is divisible by a large ideal:
  * odd unramified P of degree 1: a Z_p-minimisation through the embedding y -> r at P (local_split_conditions:
    scaling steps on the radical, then one congruence B(v0, v) = 0 mod p^e for a p-adic isotropic v0);
  * other P with v_P(A) odd, v_P(B) = 0: x = s z mod P (s^2 = B mod P), and symmetrically; v_P(A) = v_P(B) = 1:
    x in P and y = t z mod P; ramified P dividing A and B to order 2: x in P;  P over 2 and the rest: nothing.
The short vectors of the trace form sum_sigma (|x_s|^2 + |A_s| |y_s|^2 + |B_s| |z_s|^2) (random weights over the
embeddings) are enumerated (qfminim), each is tested, and then for pairs v1, v2 of them the binary form
Q(a v1 + b v2) is tested for an L-zero (B(v1, v2)^2 - Q(v1) Q(v2) a square in L: its norm a square in Q, then
nfroots).  PARI does the residue square roots, the lattices and the enumeration; every candidate is tested in exact
arithmetic, and the caller re-checks the point exactly.  Floating point only sets the enumeration radius and the
weights.
"""
import random
from math import isqrt
from fractions import Fraction as Fr
import cypari2
pari = cypari2.Pari()

class Field:
    def __init__(self, D0):
        self.D0 = D0
        self.nf = pari.nfinit(pari('y^2 - (%d)' % D0))
        self.basis = [pari(1), pari.nfbasistoalg(self.nf, pari('[0, 1]~')).lift()]   # 1, w
        w = self.basis[1]
        # w = (u0 + u1 y): coordinates
        self.w0 = Fr(str(pari.polcoef(w, 0, 'y'))); self.w1 = Fr(str(pari.polcoef(w, 1, 'y')))
    def coords(self, el):
        """coordinates (integers) of an integral element (PARI polmod / pol in y) in the basis 1, w."""
        v = pari.nfalgtobasis(self.nf, el)
        return [int(v[0]), int(v[1])]

def _polmod(F, c):
    return pari('Mod(%d + (%d)*(%s), y^2 - (%d))' % (c[0], c[1], F.basis[1], F.D0))

def _vp(n, p):
    if n == 0: return 10 ** 9
    k = 0
    while n % p == 0: n //= p; k += 1
    return k

def _mat_mul(A, B, mod=None):
    C = [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
    return [[x % mod for x in r] for r in C] if mod else C

def _rank_mod_p(M, p):
    M = [[x % p for x in r] for r in M]; rk = 0; rows = len(M); cols = len(M[0])
    for c_ in range(cols):
        piv = next((i for i in range(rk, rows) if M[i][c_]), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c_], -1, p)
        M[rk] = [x * inv % p for x in M[rk]]
        for i in range(rows):
            if i != rk and M[i][c_]:
                f = M[i][c_]; M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk, M

def _kernel_mod_p(M, p):
    rk, R = _rank_mod_p(M, p)
    n = len(M[0]); pivs = []
    for r_ in R[:rk]:
        pivs.append(next(j for j in range(n) if r_[j]))
    free = [j for j in range(n) if j not in pivs]
    ker = []
    for f in free:
        v = [0] * n; v[f] = 1
        for i, pc in enumerate(pivs): v[pc] = (-R[i][f]) % p
        ker.append(v)
    return ker

def _adj3(M):
    m = M
    cof = [[0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            rr = [r for k, r in enumerate(m) if k != i]
            mm = [[x for l, x in enumerate(r) if l != j] for r in rr]
            cof[i][j] = (-1) ** (i + j) * (mm[0][0] * mm[1][1] - mm[0][1] * mm[1][0])
    return [[cof[j][i] for j in range(3)] for i in range(3)]

def _det3(m):
    return (m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]))

def local_split_conditions(a, b, p, vtot):
    """Z_p-minimisation of x^2 - a y^2 - b z^2 (a, b integers: the images of A, B under the embedding at an odd
    unramified degree-1 prime P over p).  Returns integer conditions [(r, j)]: r . (x, y, z) = 0 mod p^j, cutting out
    a sublattice on which the form is divisible by p^(c + e) with index p^(c + e) (Holzer-optimal), or None if the
    form is anisotropic (not locally soluble)."""
    K = 3 * vtot + 8
    pk = p ** K
    G = [[1, 0, 0], [0, (-a) % pk, 0], [0, 0, (-b) % pk]]
    Bm = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]; c = 0
    for _ in range(4 * vtot + 6):
        Gp = _mat_mul(_mat_mul([list(r) for r in zip(*Bm)], G), Bm)
        assert all(x % p ** c == 0 for r in Gp for x in r)
        Gp = [[(x // p ** c) % p ** (K - c) for x in r] for r in Gp]
        e = _vp(_det3(Gp) % p ** (K - c), p)
        if e == 0 or e >= K - c - 2: break
        rk, _ = _rank_mod_p(Gp, p)
        if rk <= 1:
            if rk == 0:
                c += 1; continue
            ker = _kernel_mod_p(Gp, p)            # dim 2: Lambda' = ker + p Z^3
            Nm = [[ker[0][i], ker[1][i], 0] for i in range(3)]
            # complete with p * (a vector outside ker): use p e_k for a k with e_k not in ker mod p
            for k in range(3):
                ek = [1 if i == k else 0 for i in range(3)]
                M3 = [[ker[0][i], ker[1][i], ek[i]] for i in range(3)]
                if _det3(M3) % p: Nm = [[ker[0][i], ker[1][i], p * ek[i]] for i in range(3)]; break
            Bm = _mat_mul(Bm, Nm); c += 1; continue
        # rank 2: radical r, binary part on a complement
        r = _kernel_mod_p(Gp, p)[0]
        comp = [ek for ek in ([1, 0, 0], [0, 1, 0], [0, 0, 1])]
        W = None
        for i in range(3):
            for j in range(i + 1, 3):
                M3 = [[r[k], comp[i][k], comp[j][k]] for k in range(3)]
                if _det3(M3) % p: W = (comp[i], comp[j]); break
            if W: break
        Bf = lambda u, v: sum(u[i] * Gp[i][j] * v[j] for i in range(3) for j in range(3))
        q11, q12, q22 = Bf(W[0], W[0]) % p, Bf(W[0], W[1]) % p, Bf(W[1], W[1]) % p
        disc = (q12 * q12 - q11 * q22) % p
        if disc and pow(disc, (p - 1) // 2, p) != 1:
            # anisotropic binary part: the only F_p-point is the radical; Lambda' = Z_p r + p Z_p^3
            Nm = None
            for i in range(3):
                for j in range(i + 1, 3):
                    if _det3([[r[k], comp[i][k], comp[j][k]] for k in range(3)]) % p:
                        Nm = [[r[k], p * comp[i][k], p * comp[j][k]] for k in range(3)]; break
                if Nm: break
            Bm = _mat_mul(Bm, Nm)
            Gq = _mat_mul(_mat_mul([list(x) for x in zip(*Bm)], G), Bm)
            dv = min(_vp(x % pk, p) for rr in Gq for x in rr)
            if dv < c + 2: return None
            c += 2; continue
        sq = int(pari.lift(pari.sqrt(pari.Mod(disc, p)))) if disc else 0
        if q11:
            s_, t_ = (-q12 + sq) * pow(q11, -1, p) % p, 1
        else:
            s_, t_ = 1, 0
        v0 = [(s_ * W[0][k] + t_ * W[1][k]) % p for k in range(3)]
        # Hensel: v0 <- v0 + tau d with B(v0, d) a unit
        pkc = p ** (K - c)
        d_ = next(dd for dd in ([1, 0, 0], [0, 1, 0], [0, 0, 1]) if Bf(v0, dd) % p)
        for _h in range(12):
            qv = Bf(v0, v0) % pkc
            if qv == 0: break
            bd = Bf(v0, d_) % pkc
            tau = (-qv * pow(2 * bd, -1, pkc)) % pkc
            v0 = [(v0[k] + tau * d_[k]) % pkc for k in range(3)]
        assert Bf(v0, v0) % pkc == 0
        # conditions: phi(v) in Bm Z_p^3 (rows of adj(Bm) mod det) and (G Bm v0) . phi(v) = 0 mod p^(e + c)
        conds = []
        dB = _det3(Bm); kB = _vp(abs(dB), p)
        if kB:
            adj = _adj3(Bm)
            unit = dB // p ** kB
            inv_u = pow(unit % p ** kB, -1, p ** kB) if kB else 1
            for row in adj:
                conds.append(([x * inv_u % p ** kB for x in row], kB))
        GBv = [sum(G[i][j] * sum(Bm[j][k] * v0[k] for k in range(3)) for j in range(3)) % p ** (e + c) for i in range(3)]
        conds.append((GBv, e + c))
        return conds
    # unimodular after the scaling steps: only the lattice Bm
    conds = []
    dB = _det3(Bm); kB = _vp(abs(dB), p)
    if kB:
        adj = _adj3(Bm); unit = dB // p ** kB; inv_u = pow(unit % p ** kB, -1, p ** kB)
        for row in adj: conds.append(([x * inv_u % p ** kB for x in row], kB))
    return conds

def find_point(D0, A, B, tries=40, nvec=3000, seed=1, verbose=False, Rmul=6, npair=400):
    """A, B: PARI t_POLMOD (integral).  Returns (x, y, z) as t_POLMOD with x^2 - A y^2 - B z^2 = 0, or None."""
    F = Field(D0); nf = F.nf
    modp = pari('y^2 - (%d)' % D0)
    Ai, Bi = A.lift(), B.lift()
    nA, nB = abs(int(pari.norm(A))), abs(int(pari.norm(B)))
    ps = sorted(set(int(p) for p in pari.factor(nA)[0]) | set(int(p) for p in pari.factor(nB)[0])) if nA * nB > 1 else []
    # the multiplication by an element of O_L as a 2x2 integer matrix on coordinates
    def mat_of(el):
        cols = []
        for b in F.basis:
            cols.append(F.coords(pari.Mod(el * b, modp)))
        return [[cols[0][0], cols[1][0]], [cols[0][1], cols[1][1]]]
    Lat = pari.matid(6)
    def meet(Lat, phi, J):
        """Lat intersected with {v : phi(v) in J} (phi: 2 x 6 integer rows giving the coordinates of a linear form
        of v in the integral basis; J an ideal)."""
        H = pari.idealhnf(nf, J)
        Hinv = pari.matsolve(H, pari.matid(2))
        R_ = Hinv * pari.matrix(2, 6, phi[0] + phi[1]) * Lat
        d_ = pari.denominator(R_); Ri_ = R_ * d_
        K_ = pari.matkerint(pari.matconcat([Ri_, d_ * pari.matid(2)]))
        Kt_ = pari.matrix(6, int(pari.matsize(K_)[1]), [K_[i, j] for i in range(6) for j in range(int(pari.matsize(K_)[1]))])
        return Lat * pari.mathnf(Kt_)
    I2 = [[1, 0], [0, 1]]; Z2 = [0, 0]
    def lin(cx, cy, cz):
        """rows of v -> coords(cx x + cy y + cz z) for cx, cy, cz in O_L (PARI elements)."""
        Ms = [mat_of(cc) if cc is not None else [[0, 0], [0, 0]] for cc in (cx, cy, cz)]
        return [Ms[0][0] + Ms[1][0] + Ms[2][0], Ms[0][1] + Ms[1][1] + Ms[2][1]]
    def sqrt_high(tgt, P, e):
        """s in O_L with s^2 = tgt mod P^e (P odd), or None if tgt is not a square mod P.  Degree-1 unramified P:
        a rational integer from the p-adic embedding at P; ramified P with e <= 2: s = a + b y modulo p."""
        p_ = int(P[0]); f_ = int(P[3]); ram = int(P[2])
        if p_ == 2: return 'skip'
        if ram == 1 and f_ == 1:
            rts = pari.polrootspadic(pari('y^2 - (%d)' % D0), p_, e + 2)
            for r_ in rts:
                ri = int(pari.lift(pari.truncate(r_)))
                if int(pari.idealval(nf, pari('y') - ri, P)) >= e:
                    val = pari.subst(tgt, 'y', ri)
                    num, den_ = int(pari.numerator(val)), int(pari.denominator(val))
                    a_ = num * pow(den_, -1, p_ ** e) % p_ ** e
                    if a_ % p_ == 0: return 'skip'
                    if pow(a_, (p_ - 1) // 2, p_) != 1: return None
                    sq_ = pari.truncate(pari.sqrt(pari('%d + O(%d^%d)' % (a_, p_, e))))
                    sq_ = int(sq_) % p_ ** e
                    assert (sq_ * sq_ - a_) % p_ ** e == 0
                    return pari(sq_)
            return 'skip'
        if ram == 2 and e <= 2:
            if e == 1: return 'skip'
            # O / P^2 = O / p = F_p[y]/(y^2): tgt = t0 + t1 y, s = a + b y with a^2 = t0, 2 a b = t1
            t0 = pari.polcoef(tgt, 0, 'y'); t1 = pari.polcoef(tgt, 1, 'y')
            t0 = int(pari.numerator(t0)) * pow(int(pari.denominator(t0)), -1, p_) % p_
            t1 = int(pari.numerator(t1)) * pow(int(pari.denominator(t1)), -1, p_) % p_
            if t0 == 0: return 'skip'
            if pow(t0, (p_ - 1) // 2, p_) != 1: return None
            a_ = int(pari.lift(pari.sqrt(pari.Mod(t0, p_))))
            b_ = t1 * pow(2 * a_, -1, p_) % p_
            return pari('%d + %d*y' % (a_, b_))
        return 'skip'
    for p in ps:
        for P in pari.idealprimedec(nf, p):
            eA = int(pari.idealval(nf, Ai, P)) if Ai != 0 else 0
            eB = int(pari.idealval(nf, Bi, P)) if Bi != 0 else 0
            if verbose and (eA or eB): print('    prime over %d (f %d, e %d): v(A) %d, v(B) %d' % (p, int(P[3]), int(P[2]), eA, eB), flush=True)
            if int(P[0]) != 2 and int(P[2]) == 1 and int(P[3]) == 1 and (eA or eB):
                # unramified degree-1 odd prime: Z_p-minimisation through the embedding y -> r at P
                pp = int(P[0]); vt = eA + eB
                rts = pari.polrootspadic(pari('y^2 - (%d)' % D0), pp, 3 * vt + 10)
                rr = None
                for r_ in rts:
                    ri = int(pari.lift(pari.truncate(r_)))
                    if int(pari.idealval(nf, pari('y') - ri, P)) >= 3 * vt + 8: rr = ri; break
                if rr is not None:
                    mod_ = pp ** (3 * vt + 8)
                    def ev(z_):
                        val = pari.subst(z_, 'y', rr)
                        return int(pari.numerator(val)) * pow(int(pari.denominator(val)), -1, mod_) % mod_
                    conds = local_split_conditions(ev(Ai), ev(Bi), pp, vt)
                    if conds is None: return None
                    for (row, j) in conds:
                        Lat = meet(Lat, lin(pari(row[0]), pari(row[1]), pari(row[2])), pari.idealpow(nf, P, j))
                    continue
            if int(P[0]) != 2 and int(P[2]) == 2 and eA >= 2 and eB >= 2:
                # ramified P (P^2 = (p)) dividing A and B to order >= 2: x in P makes Q divisible by p, and
                # Q / p = eps x'^2 - (A/p) y^2 - (B/p) z^2 (x = pi x', pi^2 = eps p); exact for (2, 2)
                Lat = meet(Lat, lin(pari(1), None, None), P)
                if verbose and (eA, eB) != (2, 2): print('    (ramified P over %d: only the (2, 2) part handled)' % p, flush=True)
                continue
            if (eA >= 2 and eB == 0) or (eB >= 2 and eA == 0):
                e_ = max(eA, eB); tgt = Bi if eA else Ai
                s_ = sqrt_high(tgt, P, e_)
                Pe = pari.idealpow(nf, P, e_)
                if s_ is None:
                    if e_ % 2: return None
                    Ph = pari.idealpow(nf, P, e_ // 2)
                    Lat = meet(Lat, lin(pari(1), None, None), Ph)
                    Lat = meet(Lat, lin(None, None, pari(1)) if eA else lin(None, pari(1), None), Ph)
                    continue
                if not isinstance(s_, str):
                    Lat = meet(Lat, lin(pari(1), None, -s_) if eA else lin(pari(1), -s_, None), Pe)
                    continue
            if verbose and (eA or eB) and (eA + eB >= 2): print('    P over %d (f %d, e %d): v(A) %d, v(B) %d' % (p, int(P[3]), int(P[2]), eA, eB), flush=True)
            if eA % 2 == 0 and eB % 2 == 0: continue
            t_shared = None
            if eA % 2 == 1 and eB % 2 == 1:
                rows_el = [(0, 1)]; s = None          # x in P
                if eA == 1 and eB == 1:
                    # then also (A/pi) y^2 + (B/pi) z^2 = 0 mod P: y = t z mod P, t^2 = -(B/pi)/(A/pi) (pi in P - P^2)
                    pi_ = pari.nfbasistoalg(nf, P[1]).lift() if str(pari.type(P[1])) == 't_COL' else pari(P[1])
                    modpr = pari.nfmodprinit(nf, P)
                    if int(pari.idealval(nf, pi_, P)) != 1: pi_ = pi_ + int(P[0])
                    if int(pari.idealval(nf, pi_, P)) == 1:
                        a_ = pari.nfmodpr(nf, pari.lift(pari.Mod(Ai, modp) / pari.Mod(pi_, modp)), modpr)
                        b_ = pari.nfmodpr(nf, pari.lift(pari.Mod(Bi, modp) / pari.Mod(pi_, modp)), modpr)
                        if a_ != 0 and b_ != 0 and pari.issquare(-b_ / a_):
                            t_shared = pari.nfmodprlift(nf, pari.sqrt(-b_ / a_), modpr)
                            t_shared = pari.nfbasistoalg(nf, t_shared).lift() if str(pari.type(t_shared)) == 't_COL' else pari(t_shared)
            else:
                tgt = Bi if eA % 2 == 1 else Ai
                modpr = pari.nfmodprinit(nf, P)
                u = pari.nfmodpr(nf, tgt, modpr)
                if u == 0: continue
                if not pari.issquare(u): return None   # not locally soluble at P: wrong input
                s = pari.nfmodprlift(nf, pari.sqrt(u), modpr)
                s = pari.nfbasistoalg(nf, s).lift() if str(pari.type(s)) == 't_COL' else pari(s)
            H = pari.idealhnf(nf, P)                  # 2x2, columns generate P in the integral basis
            Hinv = pari.matsolve(H, pari.matid(2))
            # phi(v) = coords(x - s z) (or x - s y, or x): 2 x 6 integer matrix
            Mx = [[1, 0], [0, 1]]
            if t_shared is not None:
                # first x in P, then y - t z in P: two conditions, applied one after the other below
                Mt_ = mat_of(-t_shared)
                for phi in ([Mx[0] + [0, 0, 0, 0], Mx[1] + [0, 0, 0, 0]], [[0, 0] + Mx[0] + Mt_[0], [0, 0] + Mx[1] + Mt_[1]]):
                    R = Hinv * pari.matrix(2, 6, phi[0] + phi[1]) * Lat
                    d = pari.denominator(R); Ri = R * d
                    K = pari.matkerint(pari.matconcat([Ri, d * pari.matid(2)]))
                    Kt = pari.matrix(6, int(pari.matsize(K)[1]), [K[i, j] for i in range(6) for j in range(int(pari.matsize(K)[1]))])
                    Lat = Lat * pari.mathnf(Kt)
                continue
            if s is None: phi = [Mx[0] + [0, 0, 0, 0], Mx[1] + [0, 0, 0, 0]]
            else:
                Ms = mat_of(-s)
                if eA % 2 == 1: phi = [Mx[0] + [0, 0] + Ms[0], Mx[1] + [0, 0] + Ms[1]]
                else: phi = [Mx[0] + Ms[0] + [0, 0], Mx[1] + Ms[1] + [0, 0]]
            R = Hinv * pari.matrix(2, 6, phi[0] + phi[1]) * Lat
            d = pari.denominator(R)
            Ri = R * d
            K = pari.matkerint(pari.matconcat([Ri, d * pari.matid(2)]))
            Kt = pari.matrix(6, int(pari.matsize(K)[1]), [K[i, j] for i in range(6) for j in range(int(pari.matsize(K)[1]))])
            Lat = Lat * pari.mathnf(Kt)
    Lat = Lat * pari.qflll(Lat)
    idx = abs(int(pari.matdet(Lat)))
    if verbose: print('  lattice index %d digits (N(m) bound %d digits)' % (len(str(idx)), len(str(nA * nB))), flush=True)
    # exact evaluation of Q(v) from integer coordinates
    def Qv(v):
        x = _polmod(F, v[0:2]); y_ = _polmod(F, v[2:4]); z = _polmod(F, v[4:6])
        return x * x - A * y_ * y_ - B * z * z, (x, y_, z)
    # the trace form, with rational approximations of the embeddings (the enumeration only; Q(v) is exact)
    from math import isqrt
    SC = 10 ** 80
    real = D0 > 0
    w0, w1 = F.w0, F.w1          # w = w0 + w1 sqrt(D0)
    rt = Fr(isqrt(abs(D0) * SC * SC), SC)          # sqrt|D0| to 80 digits
    def ab_of(el):
        c = F.coords(el)          # el = c0 + c1 w = (c0 + c1 w0) + c1 w1 sqrt D0
        return Fr(c[0]) + c[1] * w0, c[1] * w1
    def fsqrt(q):
        q = Fr(q); return Fr(isqrt(int(q * SC * SC)), SC)
    aA, bA = ab_of(Ai); aB, bB = ab_of(Bi)
    if real:
        embs = [1, -1]
        absA = [abs(aA + sg * bA * rt) for sg in embs]; absB = [abs(aB + sg * bB * rt) for sg in embs]
        bas = [[Fr(1), w0 + sg * w1 * rt] for sg in embs]
    else:
        embs = [1]
        nA_ = aA * aA - D0 * bA * bA; nB_ = aB * aB - D0 * bB * bB
        absA = [fsqrt(nA_)]; absB = [fsqrt(nB_)]
    rng = random.Random(seed)
    def rexp(r):    # a rational approximation of exp(r)
        import math
        return Fr(int(math.exp(r) * 10 ** 12), 10 ** 12)
    for trial in range(tries):
        G = [[Fr(0)] * 6 for _ in range(6)]
        for k in range(len(embs)):
            if trial: mu = [Fr(1), rexp(rng.uniform(-3, 3)), rexp(rng.uniform(-3, 3))]
            else: mu = [Fr(1)] * 3
            if real and k == 1 and trial: mu = [1 / m_ for m_ in mu_prev]
            mu_prev = mu
            wts = [mu[0], absA[k] * mu[1], absB[k] * mu[2]]
            for blk in range(3):
                if real:
                    e = bas[k]
                    for i in range(2):
                        for j in range(2): G[2 * blk + i][2 * blk + j] += wts[blk] * e[i] * e[j]
                else:
                    # |x0 + x1 w|^2 = x0^2 + 2 Re(w) x0 x1 + |w|^2 x1^2, w = w0 + w1 sqrt(D0) (D0 < 0)
                    rew, nw = w0, w0 * w0 - D0 * w1 * w1
                    Mh = [[Fr(1), rew], [rew, nw]]
                    for i in range(2):
                        for j in range(2): G[2 * blk + i][2 * blk + j] += wts[blk] * Mh[i][j]
        Lm = [[int(Lat[i, j]) for j in range(6)] for i in range(6)]
        GL = [[sum(Lm[a][i] * G[a][b] * Lm[b][j] for a in range(6) for b in range(6)) for j in range(6)] for i in range(6)]
        m0 = min(GL[i][i] for i in range(6))
        Gi = pari.matrix(6, 6, [int(GL[i][j] * 10 ** 60 / m0) for i in range(6) for j in range(6)])
        Gi = Gi + pari.mattranspose(Gi)
        try:
            U = pari.qflllgram(Gi)
            G2 = pari.mattranspose(U) * Gi * U
            sc = max(1, int(G2[0, 0]) // 10 ** 12)
            G3 = pari.round(G2 / sc)
            if verbose and trial < 2: print('  LLL diag digits', [len(str(G3[i, i])) for i in range(6)], flush=True)
            # the enumeration radius: Rmul times the larger of the first minimum and det^(1/6) (the Holzer scale;
            # for a skew lattice the first minimum is far below it)
            dt = abs(int(pari.matdet(G3)))
            hs = int(round(dt ** (1.0 / 6))) if dt > 0 else int(G3[0, 0])
            bound = max(int(G3[0, 0]), hs) * Rmul
            res = pari.qfminim(G3, bound, nvec, 2)
        except cypari2.PariError as ex:
            if verbose: print('  PARI refused trial %d: %s' % (trial, str(ex)[:80]), flush=True)
            continue
        vs = res[2]
        nr = int(pari.matsize(vs)[0])
        if verbose: print('  trial %d: %d short vectors' % (trial, int(pari.matsize(vs)[1])), flush=True)
        cand = []
        for j in range(int(pari.matsize(vs)[1])):
            u = pari.matrix(nr, 1, [vs[i, j] for i in range(nr)])
            v = Lat * U * u
            vv = [int(v[i, 0]) for i in range(6)]
            q, pt = Qv(vv)
            if q == 0 and any(vv[2:]):
                if verbose: print('  zero found at trial %d, vector %d' % (trial, j), flush=True)
                return pt
            if len(cand) < npair: cand.append((q, pt))
        # binary step: on the O_L-span of two short vectors v1, v2, Q(a v1 + b v2) = a^2 Q1 + 2 a b B12 + b^2 Q2 has
        # an L-zero iff Delta = B12^2 - Q1 Q2 is a square in L (tested first by the norm: a square in Q)
        for i1 in range(len(cand)):
            q1, p1 = cand[i1]
            for i2 in range(i1 + 1, len(cand)):
                q2, p2 = cand[i2]
                B12 = p1[0] * p2[0] - A * p1[1] * p2[1] - B * p1[2] * p2[2]
                Dl = B12 * B12 - q1 * q2
                if Dl == 0: continue
                nD = pari.norm(Dl)
                num, den = int(pari.numerator(nD)), int(pari.denominator(nD))
                if num < 0 or isqrt(num) ** 2 != num or isqrt(den) ** 2 != den: continue
                rts = pari.nfroots(nf, pari('x^2') - pari.lift(Dl))
                if len(rts) == 0: continue
                sq = pari.Mod(pari.lift(pari.nfbasistoalg(nf, rts[0])) if str(pari.type(rts[0])) == 't_COL' else pari.lift(rts[0]), modp)
                al = -B12 + sq; be = q1
                pt = tuple(al * p1[k] + be * p2[k] for k in range(3))
                x, y_, z = pt
                if x * x - A * y_ * y_ - B * z * z == 0 and (y_ != 0 or z != 0):
                    if verbose: print('  zero found on the span of vectors %d, %d (trial %d)' % (i1, i2, trial), flush=True)
                    return pt
    return None
