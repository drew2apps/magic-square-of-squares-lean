#!/usr/bin/env python3
"""PARI port of the repository's rank-1 route for a ratio of square type (EXPERIMENT; the repo's LAST_PAIR.md, Theorem A route:
saturation.py, mw_sieve.py, closure_estimate.py and fibre.py, which need Sage).  Written by an AI model; not reviewed.

What is ported, and what changes:
  * mw_sieve: the same local sets S_p and the same lifting, with elliptic-curve arithmetic over F_p in plain Python integers.
    tau: E2 -> E1 over F_p is (x, y) -> (-x, iota y) followed by the rational isomorphism from the -1 twist of E2 to E1
    (PARI ellisisom).  The two isomorphisms differ by -1, and the residues are taken up to sign, as in the repo.
  * closure_estimate: the same constants and formulas, in exact rationals.  Short-model isomorphism from PARI's ellisisom;
    Kodaira codes and Tamagawa numbers from PARI's elllocalred (the repo takes them from Sage, which calls PARI).
  * fibre: the same polynomial equations over K = Q(i), roots by nfroots.
  * saturation: NOT fully ported.  The repo bounds the index with Sage's Cremona-Siksek lower bound for heights, which PARI does
    not have.  Here PARI's ellsaturation checks that G is saturated at every prime <= SAT; so if G = k G0 + torsion then every
    prime factor of k exceeds SAT, and G0 would have canonical height <= h(G)/SAT'^2.  That last step is a conditional.
usage: route_pari.py m,n d [N p1 p2 p3 p4] [SAT]"""
import sys, json, re, itertools, time
from math import gcd, isqrt
from fractions import Fraction as Fr
sys.path.insert(0, '..')
from sd_pari import pari, nf, I, X, s, family, rational_class, weierstrass, classify_t, fr
sys.path.insert(0, '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/second_descent')
from height_enclosure import ln_bounds, lnplus_hi, alpha as cps_alpha, LN2
pari.allocatemem(4*10**9)

def P(a): return [pari(Fr(x).numerator)/Fr(x).denominator for x in a]
def F(z): return fr(z)
def isprime(q): return q > 1 and all(q % k for k in range(2, isqrt(q) + 1))
def nxt(q):
    q += 4 - (q - 1) % 4 if (q - 1) % 4 else 4
    while not isprime(q): q += 4
    return q

class Setup:
    def __init__(self, m, n, d):
        self.m, self.n, self.d = m, n, d
        fam = family(m, n); (self.Q1, self.r1), (self.Q2, self.r2) = fam['+13'], fam['+24']
        self.ain1, self.A1 = weierstrass(d*self.Q1, self.r1); self.ain2, self.A2 = weierstrass(d*self.Q2, self.r2)
        self.E1 = pari.ellinit(P(self.ain1)); self.E2 = pari.ellinit(P(self.ain2))
        self.E1K = pari.ellinit(P(self.ain1), nf); self.E2K = pari.ellinit(P(self.ain2), nf)
        a = self.ain2; self.E2tw = pari.ellinit(P([0, -a[1], 0, a[3], -a[4]]))
        self.v = pari.ellisisom(self.E2tw, self.E1); assert self.v != 0, 'E2 is not the twist of E1 by -1'
        self.v = [F(c) for c in self.v]
        # G = P + conj(P) at the degenerate solution t = 1 (s = -i), exactly
        val = pari.lift(d*pari.subst(self.Q1, 's', -I)); rts = pari.nfroots(nf, X**2 - val); assert len(rts) > 0
        rr = pari(self.r1.numerator)/self.r1.denominator; Ap = pari(self.A1.numerator)/self.A1.denominator
        w = 1/(-I - rr); Pd = [pari.lift(Ap*w), pari.lift(Ap*rts[0]*w*w)]
        cPd = [pari.subst(c, 'y', -pari('y')) for c in Pd]
        Gk = pari.elladd(self.E1K, Pd, cPd); g = [pari.simplify(pari.lift(c)) for c in Gk]
        assert all(str(c.type()) in ('t_INT', 't_FRAC') for c in g), 'G is not rational'
        self.G = [F(c) for c in g]; assert int(pari.ellorder(self.E1, g)) == 0, 'G has finite order'
        tors = pari.elltors(self.E1); assert int(tors[0]) == 2, 'torsion is not Z/2'
        self.T = [F(c) for c in tors[2][0]]
        self.hG = float(pari.ellheight(self.E1, g))
        rk = pari.ellrank(pari.ellminimalmodel(self.E1), 2); self.rank = [int(rk[0]), int(rk[1])]
        self.rank_class = min(int(pari.ellrank(pari.ellinit(Fc), 2)[1]) for Fc in pari.ellisomat(pari.ellminimalmodel(self.E1), 0, 1)[0])
        self.q1 = [Fr(F(pari.polcoef(d*self.Q1, k, 's'))) for k in range(5)]
        self.q2 = [Fr(F(pari.polcoef(d*self.Q2, k, 's'))) for k in range(5)]

# ---------- arithmetic over F_p ----------
def mk(p, ain):
    md = lambda q: Fr(q).numerator*pow(Fr(q).denominator, -1, p) % p
    return (p, md(ain[1]), md(ain[3]), md(ain[4]))
def add(c, Pt, Qt):
    p, a2, a4, a6 = c
    if Pt is None: return Qt
    if Qt is None: return Pt
    x1, y1 = Pt; x2, y2 = Qt
    if x1 == x2:
        if (y1 + y2) % p == 0: return None
        lam = (3*x1*x1 + 2*a2*x1 + a4)*pow(2*y1, -1, p) % p
    else: lam = (y2 - y1)*pow(x2 - x1, -1, p) % p
    x3 = (lam*lam - a2 - x1 - x2) % p; return (x3, (lam*(x1 - x3) - y1) % p)
def neg(c, Pt): return None if Pt is None else (Pt[0], (-Pt[1]) % c[0])
def red(p, Pt):
    x, y = Fr(Pt[0]), Fr(Pt[1])
    if x.denominator % p == 0: return None
    return (x.numerator*pow(x.denominator, -1, p) % p, y.numerator*pow(y.denominator, -1, p) % p)

def good(S, p):
    if p % 4 != 1: return False
    nd = lambda q: Fr(q).numerator*Fr(q).denominator
    D1 = F(S.E1[11]); D2 = F(S.E2[11]); j = F(S.E1[12])
    if any(Fr(c).denominator % p == 0 for c in S.v) or Fr(S.v[0]).numerator % p == 0: return False
    for q in (D1, D2, nd(S.A1), nd(S.A2), nd(S.r1), nd(S.r2), S.d):
        if Fr(q).numerator % p == 0 or Fr(q).denominator % p == 0: return False
    if j.denominator % p == 0: return False
    jp = j.numerator*pow(j.denominator, -1, p) % p
    return jp not in (0, 1728 % p)

def local(S, p):
    assert good(S, p), 'the prime %d is not allowed' % p
    c1 = mk(p, S.ain1); c2 = mk(p, S.ain2)
    md = lambda q: Fr(q).numerator*pow(Fr(q).denominator, -1, p) % p
    sq = {}
    for y in range(1, (p + 1)//2): sq.setdefault(y*y % p, y)
    sq[0] = 0
    iota = sq[p - 1]
    u, r, ss, t = (md(x) for x in S.v); ui = pow(u, -1, p)
    def tau(Pt):                                           # E2 -> E2tw -> E1 over F_p
        if Pt is None: return None
        x, y = (-Pt[0]) % p, iota*Pt[1] % p
        x2 = (x - r)*ui*ui % p; y2 = (y - ss*(x - r) - t)*ui*ui*ui % p; return (x2, y2)
    Gp = red(p, S.G); Tp = red(p, S.T)
    dlog = {}; Xp = None; a = 0
    while True:
        dlog[Xp] = a; Xp = add(c1, Xp, Gp); a += 1
        if Xp is None: break
    o = a
    q1 = [md(x) for x in S.q1]; q2 = [md(x) for x in S.q2]
    A1, A2, r1, r2 = md(S.A1), md(S.A2), md(S.r1), md(S.r2)
    def pts(A, q, r_, sv):
        if sv is None:
            v = q[4]
            if v not in sq: return []
            yv = sq[v]; return [(0, A*yv % p), (0, (-A*yv) % p)] if yv else [(0, 0)]
        if sv == r_: return [None]
        v = (q[0] + sv*(q[1] + sv*(q[2] + sv*(q[3] + sv*q[4])))) % p
        if v not in sq: return []
        w = pow(sv - r_, -1, p); yv = sq[v]; x = A*w % p; y = A*yv*w*w % p
        return [(x, y), (x, (-y) % p)] if y else [(x, 0)]
    def res(Xp):
        out = set()
        for e in (0, 1):
            Y = Xp if e == 0 else add(c1, Xp, neg(c1, Tp))
            if Y in dlog: out.add(dlog[Y]); out.add((-dlog[Y]) % o)
        return out
    Sp = set()
    for tv in list(range(p)) + [None]:
        s1, s2 = (None, None) if tv is None else ((-iota*tv) % p, iota*tv % p)
        P1, P2, R1, R2 = pts(A1, q1, r1, s1), pts(A1, q1, r1, s2), pts(A2, q2, r2, s1), pts(A2, q2, r2, s2)
        if not (P1 and P2 and R1 and R2): continue
        ra = set(); rb = set()
        for x in P1:
            for y in P2: ra |= res(add(c1, x, y))
        for x in R1:
            for y in R2: rb |= res(tau(add(c2, x, neg(c2, y))))
        for a_ in ra:
            for b_ in rb: Sp.add((a_, b_))
    return o, Sp

def lift(loc, primes, N, say=True):
    (o1, S1), (o2, S2) = loc[primes[0]], loc[primes[1]]
    if gcd(o1, o2) != 1:
        N = min(N, 30000); cand = set((a, b) for (a0, b0) in S1 for a in range(a0, N + 1, o1) for b in range(b0, N + 1, o1))
        if say: print('orders of G modulo %d and %d not coprime: lifting from the first prime only, range %d; pairs: %d' % (primes[0], primes[1], N, len(cand)), flush=True)
        rest = primes[1:]; one = True
    else:
        M = o1*o2; inv = pow(o1, -1, o2); N = min(N, M - 1); cand = set()
        crt = lambda x1, x2: x1 + o1*(((x2 - x1)*inv) % o2)
        B2 = {}; A1s = {}
        for (a2, b2) in S2: B2.setdefault(a2, []).append(b2)
        for (a1, b1) in S1: A1s.setdefault(a1, []).append(b1)
        for a1, b1s in A1s.items():
            for a2, b2s in B2.items():
                a = crt(a1, a2)
                if a > N: continue
                for b1 in b1s:
                    for b2 in b2s:
                        b = crt(b1, b2)
                        if b <= N: cand.add((a, b))
        if say: print('pairs 0 <= a, b <= %d allowed by p = %d and %d: %d' % (N, primes[0], primes[1], len(cand)), flush=True)
        rest = primes[2:]; one = False
    for p in rest:
        o, Sp = loc[p]; cand = set(c for c in cand if (c[0] % o, c[1] % o) in Sp)
        if say: print('  after p = %d: %d' % (p, len(cand)), flush=True)
    return sorted(cand), N, one

# ---------- closure estimate ----------
L2 = Fr(6931471806, 10**10); L3 = Fr(10987, 10**4)
def dec(x): return Fr(int(x.replace('.', '')), 10**len(x.split('.')[1]))
def ulog(x): return abs(int(x)).bit_length()*L2
def hq(q): q = Fr(q); return max(ulog(q.numerator), ulog(q.denominator))
def minmodel(E):
    r = pari('my(v);F=ellminimalmodel(ellinit(%s),&v);[F[1..5],v]' % str([str(c) for c in E[:5]]).replace("'", ''))
    return pari.ellinit(r[0]), [F(c) for c in r[1]]
def chpt(Pt, v):
    u, r, s_, t = v; x, y = Fr(Pt[0]), Fr(Pt[1]); x2 = (x - r)/u**2; return (x2, (y - s_*(x - r) - t)/u**3)
def enclose(Em, xG, k=7):
    a1, a2, a3, a4, a6 = (int(F(c)) for c in Em[:5])
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6; b8 = a1*a1*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3*a3 - a4*a4
    D = -b2*b2*b8 - 8*b4**3 - 27*b6*b6 + 9*b2*b4*b6; c4 = b2*b2 - 24*b4; j = Fr(c4**3, D)
    assert int(Em[11]) == D
    na = Fr(0)
    for pp in pari.factor(abs(D))[0]:
        lr = pari.elllocalred(Em, pp); assert int(lr[0]) >= 1
        assert [int(c) for c in lr[2]] == [1, 0, 0, 0], 'model not minimal at %s' % pp
        na += cps_alpha(int(lr[1]), int(lr[3]))*ln_bounds(int(pp))[1]
    tail = lnplus_hi(Fr(abs(b2), 12)) + (LN2[1] if b2 else 0)
    BUP = na + lnplus_hi(abs(j))/4 + tail + Fr(1946, 1000)
    U = ln_bounds(abs(D))[1]/6 + lnplus_hi(abs(j))/6 + tail + Fr(214, 100)
    x = Fr(xG)
    for _ in range(k):
        num = x**4 - b4*x*x - 2*b6*x - b8; den = 4*x**3 + b2*x*x + 2*b4*x + b6; assert den != 0; x = num/den
    hx = ln_bounds(max(abs(x.numerator), x.denominator))
    return (hx[0] - BUP)/4**k, (hx[1] + U)/4**k, BUP, U
def one_const(S, Q, r):
    ain, A0 = weierstrass(S.d*Q, r); E = pari.ellinit(P(ain)); Em, _ = minmodel(E)
    c4, c6 = F(Em[9]), F(Em[10]); A, B = -27*c4, -54*c6; Es = pari.ellinit(P([0, 0, 0, A, B]))
    v = pari.ellisisom(E, Es); assert v != 0; u, rr, _, _ = (F(c) for c in v)
    al, be = Fr(1)/u**2, -rr/u**2
    HX = hq(al*A0) + hq(r) + hq(be) + 2*L2
    j = F(Es[12]); D = -16*(4*A**3 + 27*B**2); hinf = hq(j)
    c1g = (ulog(D) + hinf)/4 + hq(j)/8 + (ulog(A) + ulog(B))/2 + dec('3.724')
    c2g = (ulog(D) + hinf)/4 + (ulog(A) + ulog(B))/2 + dec('4.015')
    c1q = (ulog(D) + hinf)/4 + hq(j)/8 + ulog(abs(A) + abs(B) + 3)/2 + dec('2.919')
    c2q = (ulog(D) + hinf)/4 + ulog(abs(A) + abs(B) + 3)/2 + dec('3.21')
    Kpt = Fr(3, 2)*HX + ulog(1 + abs(A) + abs(B))/2 + L3/2 + c2g
    return dict(HX=HX, c1g=c1g, c2g=c2g, c1q=c1q, c2q=c2q, Kpt=Kpt)
def closure(S, amax):
    out = {}
    tt = X                                                  # x has the highest variable priority; t is written as x here
    f = [pari.subst(S.d*S.Q1, 's', -I*tt), pari.subst(S.d*S.Q1, 's', I*tt), pari.subst(S.d*S.Q2, 's', -I*tt), pari.subst(S.d*S.Q2, 's', I*tt)]
    def cts(g):
        fa = pari.nffactor(nf, pari.lift(g)); return all(int(e) % 2 == 0 for e in fa[1])
    sqs = []
    for k in (1, 2, 3, 4):
        for cmb in itertools.combinations(range(4), k):
            g = f[cmb[0]]
            for i in cmb[1:]: g = g*f[i]
            if cts(g): sqs.append(''.join(str(i + 1) for i in cmb))
    out['const_times_square_products'] = sqs
    fa = pari.nffactor(nf, pari.lift(f[0]*f[1]*f[2]*f[3])); roots = []
    for g in fa[0]:
        gm = pari.Mod(1, pari('y^2+1'))*g
        k = sum(1 for fi in f if pari.lift(pari.lift(fi % gm)) == 0)
        roots += [k]*int(pari.poldegree(g, 'x'))
    nrad = len(roots); out['distinct_roots'] = nrad; out['roots_in_how_many_f'] = sorted(set(roots))
    out['genus'] = 1 + (8*(-2) + nrad*4)//2
    c = one_const(S, S.Q1, S.r1); c2_ = one_const(S, S.Q2, S.r2)
    mu = 4*c['Kpt'] + 4*c2_['Kpt'] + c['c1g'] + c2_['c1g']
    C1 = dec('72.251'); bounds = {}
    for deg in (96, 192):
        h2C = 2*deg*mu
        for name, c1_, c2__ in (('Q', c['c1q'], c['c2q']), ('gen', c['c1g'], c['c2g'])):
            C2 = C1*(dec('6.019') + 4*c1_); C3 = 4*c2__; bounds[(deg, name)] = C1*h2C*deg + C2*deg*deg + C3
    Em, v = minmodel(S.E1); Gm = chpt(S.G, v)
    lo, hi, BUP, U = enclose(Em, Gm[0]); assert lo <= S.hG <= hi, (float(lo), S.hG, float(hi))
    out['hG_enclosure'] = [float(lo), float(hi)]
    low = Fr(3, 2)*lo*(amax + 1)**2
    for deg in (96, 192):
        top = max(v_ for k, v_ in bounds.items() if k[0] == deg)
        ab = isqrt(int((top/(Fr(3, 2)*lo))))
        out['ab_bound_deg%d' % deg] = ab
    out['mu'] = float(mu); out['hG_check_P'] = None
    return out, Em, Gm, lo, hi

# ---------- fibre ----------
def fibre(S):
    res = {}
    u, r, s_, t = S.v
    def tau_inv(Pt):                                       # E1 -> E2tw (inverse of the change v) -> E2 over K
        x2, y2 = Fr(Pt[0]), Fr(Pt[1]); x = u**2*x2 + r; y = u**3*y2 + s_*u**2*x2 + t
        return [pari(-x.numerator)/x.denominator, -I*(pari(y.numerator)/y.denominator)]
    def solve(ain, A, rr, consts, label):
        a2, a4, a6 = (pari(Fr(z).numerator)/Fr(z).denominator for z in (ain[1], ain[3], ain[4]))
        f = X**3 + a2*X**2 + a4*X + a6; Ap = pari(A.numerator)/A.denominator; rp = pari(rr.numerator)/rr.denominator
        cand = set(); ident = 0
        def add_x(xv):
            xv = pari.lift(pari(xv) + 0*I)
            if xv != 0: cand.add(str(pari.lift(rp + Ap/xv + 0*I)))
        for cst in consts:
            x0, y0 = pari.lift(pari(cst[0]) + 0*I), pari.lift(pari(cst[1]) + 0*I); add_x(x0)
            Dn = 2*rp*X + Ap; hD = (-Ap*X + (a2 + X + x0)*Dn)*(X - x0)**2
            pol = pari.lift((f*Dn - hD) if y0 == 0 else (((f + y0**2)*Dn - hD)**2 - 4*y0**2*f*Dn**2))
            if pol == 0: ident += 1; continue
            for rt in pari.nfroots(nf, pol): add_x(rt)
        ts = {}
        for sv in cand:
            tq = pari.lift(I*pari(sv))
            if pari.poldegree(tq, 'y') >= 1 and pari.polcoef(tq, 1, 'y') != 0: continue
            q = fr(pari.polcoef(tq, 0, 'y')); ts[str(q)] = classify_t(S.m, S.n, q)
        res[label] = dict(consts=len(consts), identically=ident, s_values=len(cand), rational_t=dict(sorted(ts.items())))
    G, T = S.G, S.T
    GT = [F(c) for c in pari.elladd(S.E1, P(G), P(T))]
    c1 = [[pari(Fr(z).numerator)/Fr(z).denominator for z in G], [pari(Fr(z).numerator)/Fr(z).denominator for z in GT]]
    solve(S.ain1, S.A1, S.r1, c1, 'A = +-G or +-G + T')
    solve(S.ain2, S.A2, S.r2, [tau_inv(G), tau_inv(GT)], 'B = +-G or +-G + T')
    return res

if __name__ == '__main__':
    m, n = map(int, sys.argv[1].split(',')); d = int(sys.argv[2])
    N = int(sys.argv[3]) if len(sys.argv) > 3 else 1000000
    primes = [int(x) for x in sys.argv[4:8]] if len(sys.argv) > 7 else [10009, 20089, 30013, 40009]
    SAT = int(sys.argv[8]) if len(sys.argv) > 8 else 2000
    out = dict(pair=[m, n], d=d); t0 = time.time()
    S = Setup(m, n, d); out['rank_E1'] = S.rank; out['rank_upper_isogeny_class'] = S.rank_class; out['hG'] = S.hG
    print('setup: rank bounds of E1 %s, smallest upper bound in its isogeny class %d, h(G) = %.6f' % (S.rank, S.rank_class, S.hG), flush=True)
    Em, v = minmodel(S.E1); Gm = chpt(S.G, v)
    sat = pari.ellsaturation(Em, [P(Gm)], SAT); hsat = float(pari.ellheight(Em, sat[0]))
    out['saturated_up_to'] = SAT; out['saturated'] = abs(hsat - S.hG) < 1e-8*max(1, S.hG)
    print('G saturated at every prime <= %d (PARI ellsaturation): %s (height after saturation %.6f)' % (SAT, out['saturated'], hsat), flush=True)
    loc = {}
    for i in range(len(primes)):
        for _ in range(20):
            if good(S, primes[i]): break
            primes[i] = nxt(primes[i])
    for p in primes:
        t1 = time.time(); loc[p] = local(S, p)
        print('p = %d: order of G %d, pairs allowed %d of %d (%.0f s)' % (p, loc[p][0], len(loc[p][1]), loc[p][0]**2, time.time() - t1), flush=True)
    left, Nc, onep = lift(loc, primes, N)
    extra = 50000
    while (0, 0) in left and len(primes) < 10:
        extra = nxt(extra)
        if not good(S, extra): continue
        loc[extra] = local(S, extra); primes.append(extra); left, Nc, onep = lift(loc, primes, N)
    plant = sorted(set([(2, 1), (Nc, 3), (Nc//2 + 7, Nc - 5), (31337 % Nc, 4242421 % Nc), (loc[primes[0]][0] + 1, Nc - 2)]))
    left2, _, _ = lift({p: (o, Sp | set((a % o, b % o) for (a, b) in plant)) for p, (o, Sp) in loc.items()}, primes, N, say=False)
    out['lift_test'] = all(x in left2 for x in plant) and not (set(left2) - set(plant) - set(left))
    out['primes'] = primes; out['sieve_range'] = Nc; out['one_prime'] = onep; out['left'] = left[:40]; out['n_left'] = len(left)
    print('sieve to %d leaves %d pairs: %s; lifting test passed: %s' % (Nc, len(left), left[:20], out['lift_test']), flush=True)
    ce, Em, Gm, lo, hi = closure(S, Nc); out['closure'] = ce
    print('closure estimate:', json.dumps(ce), flush=True)
    fb = fibre(S); out['fibre'] = fb
    print('fibre:', json.dumps(fb), flush=True)
    b96 = ce['ab_bound_deg96']; inb = [pq for pq in left if max(pq) <= b96]
    out['within_bound'] = inb; out['range_covers_bound'] = b96 <= Nc
    fib_ok = all(v_ in ('no', 'degenerate') for f_ in fb.values() for v_ in f_['rational_t'].values()) and all(f_['identically'] == 0 for f_ in fb.values())
    covered = all(a == 1 or b == 1 for a, b in inb)
    out['closed'] = bool(out['range_covers_bound'] and covered and fib_ok and out['lift_test'] and S.rank_class == 1)
    out['sec'] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)
