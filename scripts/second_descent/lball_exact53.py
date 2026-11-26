#!/usr/bin/env python3
"""L(E,1) > 0 without PARI for the closing curves that lball.py sums (conductor up to 10^14), item N1 of the round-4/5
referee report.  Written by an AI model; not read by a person.

lball.py takes N, w and every a_n of these 53 curves from PARI (its line 22).  lball_exact.py computes N (Tate's
algorithm), the a_n (counting points), w (the sign of f(i/sqrt N)) and an enclosure of L(E,1) in integer interval
arithmetic, but so far only for the uncounted curve of (440,279).  This program applies the same certificate to the
closing curves of second_123.jsonl, smallest conductor first:

 1. The curves: every row of second_123.jsonl with status "closed", the curve named in "by", its "rank0_curve" made
    integral (lball_exact.integral_model), N by lball_exact.local_data (Tate's algorithm, exact); those with N <= 10^14,
    sorted by N.  These are the 53 of lball.txt, in the same order (checked: the pairs are compared with lball.txt).
 2. a_p.  At bad p from Tate's algorithm (+1 split, -1 non-split, 0 additive).  At good p < P0 = 1000 by counting points
    (count_ap: the formula of lball_exact.ap_count in pure Python, checked against it in --aptest).  At good p >= P0 by
    an exact order certificate on E(F_p), found by baby-step giant-step:
      for random points P of E(F_p) (x random, y from a square root that is checked: y^2 + a1 x y + a3 y = x^3 + ...),
      a positive m with m P = O is found by baby-step giant-step over the Hasse interval, m P = O is recomputed, and the
      exact order of P is obtained from m by removing prime factors l of m while (m/l) P = O (m factored by trial
      division).  Let L be the lcm of the orders.  #E(F_p) is a multiple of L (Lagrange) and lies in
      H = [p + 1 - isqrt(4p), p + 1 + isqrt(4p)] (Hasse: (#E - p - 1)^2 <= 4p).  If exactly one multiple of L lies in
      H, it is #E(F_p), and a_p = p + 1 - #E(F_p).  Otherwise another point is taken; after 12 points the prime falls
      back to counting.  Everything is Python integer arithmetic mod p on the given integral model (or, at p | Delta of
      that model with good reduction, on Tate's p-minimal model); the model is asserted to have Delta != 0 mod p.
      What is trusted: Hasse's theorem, Lagrange, and that ec_add is the group law (Silverman III.2.3; it is
      cross-checked against counting points, --aptest).
    a_n from the a_p as in lball_exact.an_from_ap (for the two large curves the same recursion on arrays, an_array).
    No numpy and no PARI: the program runs on Python's integers alone (PARI only for the checks of 5).
 3. w = +1 from F = sum a_n exp(-2 pi n/sqrt N) != 0 (lball_exact, docstring 3(b): w = -1 forces F = 0).  The local
    product of lball_exact is printed beside it (not used; it is missing for (41,1), case at 2 not implemented).
 4. L(E,1): lball_exact.evaluate, with every rounding directed and the tail bounded by |a_n| <= 2n: the enclosure at
    the cut-off y = 1 (what is certified) and at y = 6/5, and the Fricke residual at y = 6/5; the guard (the two
    enclosures meet and the residual contains 0) is a check, not part of the proof.
 5. Checks, not used in 1-4: lball.py's printed enclosure (lball.txt) must meet this one; with PARI available (not
    --nopari), N, w and every a_n with n <= the number of terms used are compared with PARI (ellglobalred, ellrootno,
    ellan of the minimal model).

Certified line: L(E,1) >= lower > 0 with N and the a_n computed here and w = +1 from F != 0.  With modularity, Carayol,
Atkin-Lehner and Kolyvagin (or Kato) this is rank 0, as for (440,279) in lball_exact.py (tier E in LEDGER.md).

usage: lball_exact53.py --list                       the 53 curves, N, w (local), terms needed
       lball_exact53.py [--nopari] [--only i,j,..] [--first K] [--maxN X] [--jsonl out.jsonl]
                                                       certify the selected curves (index in the N-sorted list, 0-based);
                                                       one line per curve on stdout, a JSON line per curve to out.jsonl
       lball_exact53.py --merge a.jsonl b.jsonl ... [--write]   lball_exact53.txt (and _timing.txt, and every record in .jsonl)
       lball_exact53.py --negative [--write]        forged N, w and a_p on three curves through lball_exact.evaluate (its
                                                       guard must refuse each), and
                                                       forged orders #E(F_p) (the order certificate must refuse each)
       lball_exact53.py --aptest [--write]          the order-certificate a_p against counting points
       lball_exact53.py --controls [--write]        67a1 (rank 0) certified; 79a1, 5077a1 (w = -1), 389a1 (rank 2) refused
       lball_exact53.py --apchunk I K J out.bin     the two large curves (I = 51, 52): a_p for every K-th prime from the J-th,
       lball_exact53.py --bigmerge I K prefix       then a_n and L(E,1) from prefix0.bin ... prefix{K-1}.bin (prefix.jsonl),
       lball_exact53.py --bigcheck I prefix 0       and PARI's ellap at every prime used, as a check (added to prefix.jsonl)
CI (proposed): quick:  python3 lball_exact53.py --first 3 --nopari | diff - <(head -n 3 lball_exact53.txt)
                       python3 lball_exact53.py --negative | diff - lball_exact53_negative.txt
                       python3 lball_exact53.py --controls | diff - lball_exact53_controls.txt
               long:   python3 lball_exact53.py --maxN 3e10 --nopari | diff - <(head -n 51 lball_exact53.txt)
                       python3 lball_exact53.py --aptest | diff - lball_exact53_aptest.txt
               large:  for J in 0 1; do python3 lball_exact53.py --apchunk 51 2 $J /tmp/b51_$J.bin & done; wait
                       python3 lball_exact53.py --bigmerge 51 2 /tmp/b51_ --nopari | diff - <(sed -n 52p lball_exact53.txt)
                       (and the same with 52, 3 chunks, line 53; about 20 minutes per chunk)
Outputs: lball_exact53.txt (one line per curve and a summary), lball_exact53.jsonl (every record: both enclosures, F,
the numbers of terms, a_p statistics, the PARI check), lball_exact53_timing.txt, lball_exact53_negative.txt,
lball_exact53_aptest.txt, lball_exact53_controls.txt."""
import sys, os, json, math, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lball_exact as LX

HERE = os.path.dirname(os.path.abspath(__file__))
P0 = 1000                     # below this a_p is counted
TRIES = 12                    # points per prime before falling back to counting


# ---------------------------------------------------------------- the group law on a Weierstrass model mod p
def ec_add(P, Q, A, p):
    """P + Q on y^2 + a1 x y + a3 y = x^3 + a2 x^2 + a4 x + a6 over F_p (A = (a1, ..., a6) reduced mod p); None is O.
    Silverman, Algorithm III.2.3."""
    if P is None: return Q
    if Q is None: return P
    a1, a2, a3, a4, a6 = A
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2 + a1*x2 + a3) % p == 0: return None            # Q = -P
        lam = (3*x1*x1 + 2*a2*x1 + a4 - a1*y1) * pow((2*y1 + a1*x1 + a3) % p, -1, p) % p
    else:
        lam = (y2 - y1) * pow((x2 - x1) % p, -1, p) % p
    nu = (y1 - lam*x1) % p
    x3 = (lam*lam + a1*lam - a2 - x1 - x2) % p
    return (x3, (-(lam + a1)*x3 - nu - a3) % p)

def ec_neg(P, A, p):
    if P is None: return None
    a1, a2, a3, a4, a6 = A
    return (P[0], (-P[1] - a1*P[0] - a3) % p)

def ec_mul(k, P, A, p):
    assert k >= 0
    R = None
    for bit in bin(k)[2:]:
        R = ec_add(R, R, A, p)
        if bit == '1': R = ec_add(R, P, A, p)
    return R

def on_curve(P, A, p):
    a1, a2, a3, a4, a6 = A; x, y = P
    return (y*y + a1*x*y + a3*y - (((x + a2)*x + a4)*x + a6)) % p == 0

def sqrt_mod(g, p):
    """a square root of the nonzero square g mod the odd prime p (Tonelli-Shanks); checked"""
    g %= p
    if p % 4 == 3: r = pow(g, (p + 1)//4, p)
    else:
        q, s = p - 1, 0
        while q % 2 == 0: q //= 2; s += 1
        z = 2
        while pow(z, (p - 1)//2, p) != p - 1: z += 1
        m, c, t, r = s, pow(z, q, p), pow(g, q, p), pow(g, (q + 1)//2, p)
        while t != 1:
            i, t2 = 0, t
            while t2 != 1: t2 = t2*t2 % p; i += 1
            b = pow(c, 1 << (m - i - 1), p)
            m, c, t, r = i, b*b % p, t*b*b % p, r*b % p
    assert r*r % p == g
    return r

def random_point(A, p, rnd):
    a1, a2, a3, a4, a6 = A
    b2 = a1*a1 + 4*a2; b4 = 2*a4 + a1*a3; b6 = a3*a3 + 4*a6; inv2 = (p + 1)//2
    while True:
        x = rnd.randrange(p)
        g = (((4*x + b2)*x + 2*b4)*x + b6) % p          # (2y + a1 x + a3)^2 = g
        if g == 0 or pow(g, (p - 1)//2, p) != 1: continue
        y = (sqrt_mod(g, p) - a1*x - a3) * inv2 % p
        P = (x, y); assert on_curve(P, A, p)
        return P

def primes_upto(M):
    """the primes <= M (sieve of Eratosthenes on a bytearray; no numpy, so the program also runs where numpy is missing)"""
    s = bytearray([1]) * (M + 1); s[0:2] = b'\x00\x00'
    for i in range(2, math.isqrt(M) + 1):
        if s[i]: s[i*i::i] = bytes(len(range(i*i, M + 1, i)))
    return [i for i in range(M + 1) if s[i]] if M < 10**6 else list(_compress(s))

def _compress(s):
    import itertools
    return itertools.compress(range(len(s)), s)

def count_ap(model, p):
    """p + 1 - #E(F_p) by counting points (pure Python; the formula of lball_exact.ap_count: #E(F_p) = p + 1 + sum_x
    chi(4x^3 + b2 x^2 + 2 b4 x + b6), chi the Legendre symbol read off the table of squares)"""
    a1, a2, a3, a4, a6 = model
    if p == 2:
        n = 1 + sum(1 for x in range(2) for y in range(2) if (y*y + a1*x*y + a3*y - x**3 - a2*x*x - a4*x - a6) % 2 == 0)
        return 3 - n
    b2, b4, b6, b8 = LX.b_inv(model); b2 %= p; b4 = 2*b4 % p; b6 %= p
    chi = [-1] * p; chi[0] = 0
    for x in range(1, (p + 1)//2): chi[x*x % p] = 1
    return -sum(chi[(((4*x + b2)*x + b4)*x + b6) % p] for x in range(p))

SMALLP = None
def factor_small(m):
    """the prime factors of 0 < m < 2^32 by trial division"""
    global SMALLP
    if SMALLP is None: SMALLP = primes_upto(1 << 16)
    assert 0 < m < 1 << 32
    out = []
    for l in SMALLP:
        if l*l > m: break
        if m % l == 0:
            out.append(l)
            while m % l == 0: m //= l
    if m > 1: out.append(m)
    return out

def bsgs_multiple(P, A, p, lo, hi):
    """a positive m with m P = O (baby-step giant-step over [lo, hi]); checked by recomputing m P"""
    W = hi - lo; B = math.isqrt(W // 2) + 1
    table = {}; R = P; m = None
    for j in range(1, B + 1):                            # R = jP
        if R is None: m = j; break
        table.setdefault(R[0], (j, R[1])); R = ec_add(R, P, A, p)
    if m is None:
        S = 2*B + 1; step = ec_mul(S, P, A, p); n = lo + B; G = ec_mul(n, P, A, p)
        while n - B <= hi + S:
            if G is None: m = n; break
            if G[0] in table:
                j, yj = table[G[0]]
                m = n - j if G[1] == yj else n + j; break
            G = ec_add(G, step, A, p); n += S
        assert m is not None, ('no multiple found', p)
    assert m > 0 and ec_mul(m, P, A, p) is None
    return m

def exact_order(m, P, A, p):
    """the order of P, from a positive multiple m of it"""
    n = m
    for l in factor_small(m):
        while n % l == 0 and ec_mul(n // l, P, A, p) is None: n //= l
    assert ec_mul(n, P, A, p) is None
    for l in factor_small(n): assert ec_mul(n // l, P, A, p) is not None
    return n

def hasse(p):
    r = math.isqrt(4*p); return p + 1 - r, p + 1 + r       # all n with (n - p - 1)^2 <= 4p

def unique_multiple(L, lo, hi):
    """the multiple of L in [lo, hi] if there is exactly one, else None (and an assertion that there is at least one)"""
    first = -(-lo // L) * L
    assert first <= hi, 'no multiple of the lcm of the orders in the Hasse interval'
    return first if first + L > hi else None

def model_mod(model, p):
    A = tuple(x % p for x in model)
    assert LX.c_inv(model)[2] % p != 0, ('model not good at', p)
    return A

def ap_order(model, p, rnd):
    """(a_p, number of points used) from the order certificate, or (None, TRIES)"""
    A = model_mod(model, p); lo, hi = hasse(p); L = 1
    for k in range(TRIES):
        P = random_point(A, p, rnd)
        o = exact_order(bsgs_multiple(P, A, p, lo, hi), P, A, p)
        L = L * o // math.gcd(L, o)
        n = unique_multiple(L, lo, hi)
        if n is not None: return p + 1 - n, k + 1
    return None, TRIES

def verify_order(model, p, n, rnd, k=TRIES):
    """the certificate for a claimed #E(F_p) = n: n P = O for up to k random points, and n is the only multiple of the
    lcm of their exact orders in the Hasse interval.  True or False (False also if k points do not settle it; used by
    the negative test)"""
    A = model_mod(model, p); lo, hi = hasse(p)
    if not lo <= n <= hi: return False
    L = 1
    for _ in range(k):
        P = random_point(A, p, rnd)
        if ec_mul(n, P, A, p) is not None: return False
        o = exact_order(n, P, A, p); L = L * o // math.gcd(L, o)
        u = unique_multiple(L, lo, hi)
        if u is not None: return u == n
    return False

def ap_one(a, LD, p, st):
    """a_p: Tate at bad p, counting below P0, the order certificate above (counting if it does not settle)"""
    bad = LD['rows']
    if p in bad and bad[p]['red'] != 'good': st['bad'] += 1; return bad[p]['ap']
    model = bad[p]['model'] if p in bad else a
    if p < P0:
        assert LX.c_inv(model)[2] % p != 0, ('model not good at', p)
        x = count_ap(model, p); st['counted'] += 1
    else:
        x, k = ap_order(model, p, random.Random(p)); st['points'] += k
        if x is None: x = count_ap(model, p); st['fallback'] += 1
        else: st['certified'] += 1
    assert x*x <= 4*p, ('Hasse', p, x)
    return x

def new_stats(): return dict(counted=0, certified=0, fallback=0, points=0, bad=0)

def ap_table_fast(a, LD, M):
    st = new_stats(); ap = {p: ap_one(a, LD, p, st) for p in primes_upto(M)}
    return ap, st


# ---------------------------------------------------------------- the curves
def closing_curves():
    rows = []
    for l in open(os.path.join(HERE, 'second_123.jsonl')):
        r = json.loads(l)
        if r['status'] != 'closed': continue
        c = [c for c in r['curves'] if c['curve'] == r['by']][0]
        rows.append((tuple(r['pair']), LX.integral_model(c['rank0_curve'])))
    return rows

def the_53():
    out = []
    for pair, a in closing_curves():
        LD = LX.local_data(a)
        if LD['N'] <= 10**14: out.append((LD['N'], pair, a, LD))
    out.sort(key=lambda t: t[0])
    return out

def lball_txt():
    """pair -> (lower, upper) as printed in lball.txt (9 decimals, from PARI's a_n)"""
    from fractions import Fraction
    d = {}
    for l in open(os.path.join(HERE, 'lball.txt')):
        if 'L(E,1) in' not in l: continue
        m, n = l.split()[:2]; iv = l.split('L(E,1) in [')[1].split(']')[0].split(',')
        d[(int(m), int(n))] = (Fraction(iv[0].strip()), Fraction(iv[1].strip()))
    return d

def dec(x, dig=12, up=False):
    """the integer interval end x / 2^PREC as a decimal rounded outward (down, or up)"""
    return LX.upper_str(x, dig) if up else LX.lower_str(x, dig)


# ---------------------------------------------------------------- one curve
PARI = None
def get_pari():
    """PARI, for the cross-checks only"""
    global PARI
    if PARI is None:
        import cypari2
        PARI = cypari2.Pari(); PARI.allocatemem(2*10**9, silent=True)
    return PARI

def certify(idx, N, pair, a, LD, nopari):
    t = {}; t0 = time.time()
    Mmax = LX.max_terms(N)
    ap, st = ap_table_fast(a, LD, Mmax); t['ap'] = time.time() - t0
    an = LX.an_from_ap(ap, LD, Mmax); t['an'] = time.time() - t0 - t['ap']
    PI = LX.pi_interval(); t1 = time.time()
    ev = LX.evaluate(an, N, 1, PI, Mmax); t['series'] = time.time() - t1
    F = ev['F']; wF = 1 if F[0] > 0 or F[1] < 0 else None
    L1 = ev['L_y1']; L2 = ev['L_y2']
    positive = L1[0] > 0 and wF == 1
    from fractions import Fraction
    lo12 = Fraction(dec(L1[0])); hi12 = Fraction(dec(L1[1], up=True))
    ref = lball_txt().get(pair)
    agree = ref is not None and max(lo12, ref[0]) <= min(hi12, ref[1])
    r = dict(idx=idx, pair=list(pair), N=N, fN=' * '.join('%d^%d' % (p, T['f']) if T['f'] > 1 else str(p) for p, T in LD['rows'].items() if T['f']),
             w_local=LD['w_local'], w_fricke=wF, F=[LX.lower_str(F[0], 3), LX.upper_str(F[1], 3)], terms=list(ev['M']), Mmax=Mmax,
             L_y1=[dec(L1[0]), dec(L1[1], up=True)], L_y2=[dec(L2[0]), dec(L2[1], up=True)], guard=ev['guard'],
             positive=positive, lball=[str(ref[0]), str(ref[1])] if ref else None, agrees_lball=agree, ap_stats=st)
    if not nopari:
        t2 = time.time(); pari = get_pari()
        E = pari.ellminimalmodel(pari.ellinit(list(a)))
        Np = int(pari.ellglobalred(E)[0]); wp = int(pari.ellrootno(E))
        pa = [int(x) for x in pari.ellan(E, Mmax)]
        r['pari'] = dict(N=Np == N, w=wp == wF, an=pa == an[1:Mmax + 1], n=Mmax)
        t['pari'] = time.time() - t2
    t['total'] = time.time() - t0
    r['time'] = {k: round(v, 1) for k, v in t.items()}
    return r

def line(r):
    pq = '%d %d' % tuple(r['pair'])
    return ('%-8s N = %d = %s  w = %s (F in [%s, %s]; local %s)  terms %d  L(E,1) in [%s, %s]  guard %s  lball.py %s  %s' % (
        pq, r['N'], r['fN'], '+1' if r['w_fricke'] == 1 else 'NOT DETERMINED', r['F'][0], r['F'][1],
        {1: '+1', -1: '-1', None: 'not implemented'}[r['w_local']], r['terms'][0], r['L_y1'][0], r['L_y1'][1],
        'passes' if r['guard'] else 'FAILS', 'agrees' if r['agrees_lball'] else 'DISAGREES',
        'certified positive' if r['positive'] and r['guard'] else 'NOT CERTIFIED'))

def run(sel, nopari, jsonl):
    rows = the_53(); fh = open(jsonl, 'a') if jsonl else None
    for i in sel:
        N, pair, a, LD = rows[i]
        r = certify(i, N, pair, a, LD, nopari)
        print(line(r), flush=True)
        if fh: fh.write(json.dumps(r) + '\n'); fh.flush()
        assert r['positive'] and r['guard'] and r['agrees_lball'], ('NOT CERTIFIED or disagreement', pair)
        if 'pari' in r: assert all(r['pari'][k] for k in ('N', 'w', 'an')), ('PARI cross-check disagrees', pair, r['pari'])

def merge(paths, write):
    rs = {}
    for path in paths:
        for l in open(path):
            r = json.loads(l); rs[tuple(r['pair'])] = r
    rows = the_53(); order = [pair for N, pair, a, LD in rows]
    ref = list(lball_txt().keys()); assert ref == order, 'the 53 curves are not those of lball.txt in the same order'
    lines = []; tl = ['ratio  N  terms(max)  a_p: bad/counted/certified/fallback/points  seconds: a_p, a_n, series, PARI, total']
    done = [p for p in order if p in rs]
    for pair in done:
        r = rs[pair]; lines.append(line(r))
        st = r['ap_stats']; tm = r['time']
        tl.append('%d %d  %d  %d(%d)  %d/%d/%d/%d/%d  %s, %s, %s, %s, %s' % (pair[0], pair[1], r['N'], r['terms'][0], r['Mmax'], st['bad'], st['counted'],
                  st['certified'], st['fallback'], st['points'], tm['ap'], tm['an'], tm['series'], tm.get('pari', '-'), tm['total']))
    cert = [p for p in done if rs[p]['positive'] and rs[p]['guard']]
    agree = [p for p in done if rs[p]['agrees_lball']]
    pari_ok = [p for p in done if 'pari' in rs[p] and all(rs[p]['pari'][k] for k in ('N', 'w', 'an'))]
    left = [p for p in order if p not in rs]
    by_an = [p for p in pari_ok if 'primes' not in rs[p]['pari']]
    lines.append('done %d of the 53 closing curves of lball.py (smallest conductor first, up to N = %d): L(E,1) > 0 certified without PARI for %d; '
                 'enclosure meets lball.py\'s for %d; a check, not used: N, w and every a_p used equal to PARI\'s for %d (for %d also every a_n '
                 'used, for the other %d the a_n up to 10^6)' % (len(done), max(rs[p]['N'] for p in done), len(cert), len(agree), len(pari_ok),
                 len(by_an), len(pari_ok) - len(by_an)))
    lines.append('not done here (still on PARI\'s a_n in lball.py): %s' % (', '.join('(%d,%d)' % p for p in left) if left else 'none'))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert len(cert) == len(done) == len(agree)
    if write:
        open(os.path.join(HERE, 'lball_exact53.txt'), 'w').write(out)
        open(os.path.join(HERE, 'lball_exact53_timing.txt'), 'w').write('\n'.join(tl) + '\n')
        open(os.path.join(HERE, 'lball_exact53.jsonl'), 'w').write(''.join(json.dumps(rs[p]) + '\n' for p in done))


# ---------------------------------------------------------------- the two large curves: a_p in K parallel chunks, then one sum
def apchunk(idx, K, J, out):
    """a_p for the primes p <= max_terms(N) of curve idx whose index in the list of primes is J mod K; appended to `out`
    as 32-bit integers (resumes where an earlier run of the same chunk stopped)"""
    from array import array
    N, pair, a, LD = the_53()[idx]; M = LX.max_terms(N)
    mine = primes_upto(M)[J::K]
    done = os.path.getsize(out) // 4 if os.path.exists(out) else 0
    st = new_stats(); t0 = time.time(); fh = open(out, 'ab'); buf = array('i')
    for i in range(done, len(mine)):
        buf.append(ap_one(a, LD, mine[i], st))
        if len(buf) == 20000 or i == len(mine) - 1:
            fh.write(buf.tobytes()); fh.flush(); buf = array('i')
            sys.stderr.write('chunk %d/%d of (%d,%d): %d of %d primes, p = %d, %.0f s\n' % (J, K, pair[0], pair[1], i + 1, len(mine), mine[i], time.time() - t0))
    fh.close()
    json.dump(dict(st, time=round(time.time() - t0, 1), n=len(mine), resumed_at=done), open(out + '.json', 'w'))

def an_array(ap_of, LD, M):
    """a_1..a_M in an array of 64-bit integers, from a_p (ap_of: an array indexed by p, zero elsewhere), by the smallest
    prime factor and the recursion of lball_exact.an_from_ap"""
    from array import array
    bad = {p for p, T in LD['rows'].items() if T['red'] != 'good'}
    spf = array('i', bytes(4*(M + 1)))
    for p in reversed(primes_upto(math.isqrt(M))):
        spf[p*p::p] = array('i', [p]) * len(range(p*p, M + 1, p))
    an = ap_of; an[0] = 0; an[1] = 1
    for n in range(4, M + 1):
        p = spf[n]
        if p == 0: continue                                   # n prime: a_n = a_p already
        m = n // p; q = p
        while m % p == 0: m //= p; q *= p
        if m > 1: an[n] = an[q] * an[m]
        elif p in bad: an[n] = an[p] * an[n // p]
        else: an[n] = an[p] * an[n // p] - p * an[n // (p*p)]
    return an

def bigmerge(idx, K, prefix, nopari):
    from array import array
    N, pair, a, LD = the_53()[idx]; M = LX.max_terms(N); t0 = time.time()
    ps = primes_upto(M); ap_of = array('q', bytes(8*(M + 1))); st = new_stats()
    for J in range(K):
        v = array('i'); v.frombytes(open('%s%d.bin' % (prefix, J), 'rb').read())
        mine = ps[J::K]; assert len(v) == len(mine), ('chunk incomplete', J, len(v), len(mine))
        for p, x in zip(mine, v):
            assert x*x <= 4*p; ap_of[p] = x
        s = json.load(open('%s%d.bin.json' % (prefix, J)))
        for k in st: st[k] += s[k]
    for p, T in LD['rows'].items():                         # bad primes again from Tate's algorithm
        if T['red'] != 'good' and p <= M: assert ap_of[p] == T['ap']
    del ps
    tap = time.time() - t0
    an = an_array(ap_of, LD, M); tan = time.time() - t0 - tap
    PI = LX.pi_interval(); t1 = time.time()
    ev = LX.evaluate(an, N, 1, PI, M); tser = time.time() - t1
    F = ev['F']; wF = 1 if F[0] > 0 or F[1] < 0 else None
    L1 = ev['L_y1']; L2 = ev['L_y2']
    from fractions import Fraction
    lo12 = Fraction(dec(L1[0])); hi12 = Fraction(dec(L1[1], up=True)); ref = lball_txt().get(pair)
    r = dict(idx=idx, pair=list(pair), N=N, fN=' * '.join('%d^%d' % (p, T['f']) if T['f'] > 1 else str(p) for p, T in LD['rows'].items() if T['f']),
             w_local=LD['w_local'], w_fricke=wF, F=[LX.lower_str(F[0], 3), LX.upper_str(F[1], 3)], terms=list(ev['M']), Mmax=M,
             L_y1=[dec(L1[0]), dec(L1[1], up=True)], L_y2=[dec(L2[0]), dec(L2[1], up=True)], guard=ev['guard'],
             positive=L1[0] > 0 and wF == 1, lball=[str(ref[0]), str(ref[1])] if ref else None,
             agrees_lball=ref is not None and max(lo12, ref[0]) <= min(hi12, ref[1]), ap_stats=st,
             time=dict(ap='%d chunks' % K, merge=round(tap, 1), an=round(tan, 1), series=round(tser, 1), total=round(time.time() - t0, 1)))
    print(line(r), flush=True)
    open('%s.jsonl' % prefix, 'w').write(json.dumps(r) + '\n')
    assert r['positive'] and r['guard'] and r['agrees_lball']
    open('%s_an.bin' % prefix, 'wb').write(an.tobytes())

def bigcheck(idx, prefix, nprimes):
    """a check, not used: PARI's ellap at every prime p <= M (or the first nprimes) against the merged a_n"""
    from array import array
    pari = get_pari()
    N, pair, a, LD = the_53()[idx]; M = LX.max_terms(N)
    an = array('q'); an.frombytes(open('%s_an.bin' % prefix, 'rb').read()); assert len(an) == M + 1
    E = pari.ellminimalmodel(pari.ellinit(list(a)))
    assert int(pari.ellglobalred(E)[0]) == N
    ps = primes_upto(M)[:nprimes] if nprimes else primes_upto(M); bad = 0
    for p in ps: bad += int(pari.ellap(E, p)) != an[p]
    m = min(M, 10**6); pa = [int(x) for x in pari.ellan(E, m)]
    same = pa == list(an[1:m + 1])
    wp = int(pari.ellrootno(E))
    print('(%d,%d): PARI ellap at %d primes p <= %d: %d disagreements; PARI ellan for n <= %d: %s; PARI root number %d' % (
        pair[0], pair[1], len(ps), ps[-1], bad, m, 'agree' if same else 'DISAGREE', wp))
    r = json.loads(open('%s.jsonl' % prefix).read())
    r['pari'] = dict(N=True, w=wp == r['w_fricke'], an=bad == 0 and same and not nprimes, n=M, primes=len(ps))
    open('%s.jsonl' % prefix, 'w').write(json.dumps(r) + '\n')
    assert bad == 0 and same


# ---------------------------------------------------------------- negative tests
NEG_TL = 1e-3

def negative(write):
    """forged N, w, a_p on three curves, through LX.evaluate (the function certify() calls): its guard must refuse each;
    forged #E(F_p): the order certificate must refuse.  (Before 26 November 2026 this test used a copy of the overlap test
    of evaluate, without the Fricke residual.)"""
    rows = the_53(); lines = ['lball_exact53.py --negative: forged inputs, each must be refused', '']
    PI = LX.pi_interval(); nref = ntot = 0
    for i in (0, 1, 7):                                        # (21,20), (15,8), (23,7)
        N, pair, a, LD = rows[i]
        bad = sorted(p for p, T in LD['rows'].items() if T['f'])
        lN = bad[-1]; nN = next(q for q in (3, 5, 7, 11, 13) if N % q)
        # every forged input goes through LX.evaluate, the function certify() calls, and is refused by its guard
        # (the two enclosures meet and the Fricke residual contains 0); forged N and w with tails below NEG_TL
        Mmax = max(max(LX.max_terms(n, NEG_TL) for n in (N // lN, N * nN, 2*N)), LX.max_terms(N))
        ap, st = ap_table_fast(a, LD, Mmax); an0 = LX.an_from_ap(ap, LD, Mmax)
        good = [p for p in sorted(ap) if p >= 5 and p not in LD['rows']]
        pmid = next(p for p in good if p > math.isqrt(N) // 4)
        def two(an, Nf, wf, tl):
            ev = LX.evaluate(an, Nf, wf, PI, Mmax, tL=tl)
            return ev['L_y1'], ev['L_y2'], ev['guard'], ev
        def why(ev):
            return ', '.join(x for x, y in (('enclosures disjoint', not ev['overlap']), ('residual excludes 0', not ev['R_has_0'])) if y)
        lines.append('curve (%d,%d), N = %d:' % (pair[0], pair[1], N))
        for tl in (NEG_TL, LX.TARGET):
            Ly1, Ly2, ok, ev = two(an0, N, 1, tl)
            lines.append('  true inputs, tails below %g: L(y=1) %s, L(y=6/5) %s, residual %s: guard %s' % (tl, LX.iv(Ly1, 10), LX.iv(Ly2, 10), LX.iv(ev['R'], 6), 'passes' if ok else 'FAILS'))
            assert ok
        bp = bad[0]; bT = LD['rows'][bp]
        forged = [('N/%d (%d left out)' % (lN, lN), N // lN, 1, None, NEG_TL), ('%dN (%d added)' % (nN, nN), N * nN, 1, None, NEG_TL),
                  ('2N', 2*N, 1, None, NEG_TL), ('w = -1', N, -1, None, NEG_TL),
                  ('a_5 one more' if 5 not in LD['rows'] else 'a_7 one more', N, 1, 'small', LX.TARGET),
                  ('a_%d one more (p near sqrt(N)/4)' % pmid, N, 1, {pmid: ap[pmid] + 1}, LX.TARGET),
                  ('a_%d = %d (bad prime, %s)' % (bp, 1 if bT['ap'] != 1 else -1, bT['red']), N, 1, {bp: 1 if bT['ap'] != 1 else -1}, LX.TARGET)]
        for name, Nf, wf, ov, tl in forged:
            if ov == 'small': q = 5 if 5 not in LD['rows'] else 7; ov = {q: ap[q] + 1}
            an = an0 if ov is None else LX.an_from_ap({**ap, **ov}, LD, Mmax)
            Ly1, Ly2, ok, ev = two(an, Nf, wf, tl)
            lines.append('  forged %-40s (tails below %g) L(y=1) %s, L(y=6/5) %s, residual %s: guard %s' % (name + ':', tl, LX.iv(Ly1, 10), LX.iv(Ly2, 10),
                         LX.iv(ev['R'], 6), 'passes (NOT REFUSED)' if ok else 'refuses (%s)' % why(ev)))
            ntot += 1; nref += not ok
        # forged group orders at the first 10 good primes above 10^5, the true order by counting points: the true one must
        # pass, n + d, d = -4, -2, 2, 4 must be refused
        tested = 0; refused = 0; truepass = 0
        for p in [q for q in primes_upto(101000) if q > 10**5 and q not in LD['rows']][:10]:
            n = p + 1 - count_ap(a, p)
            truepass += verify_order(a, p, n, random.Random(7*p))
            for d in (-4, -2, 2, 4):
                tested += 1; refused += not verify_order(a, p, n + d, random.Random(p + d))
        lines.append('  order certificate at the first 10 good primes above 10^5, #E(F_p) by counting points: the true order accepted at %d of 10; '
                     'forged orders #E(F_p) + d, d = -4, -2, 2, 4: %d of %d refused' % (truepass, refused, tested))
        assert truepass == 10
        ntot += tested; nref += refused
        # two broken group laws (P + Q replaced by -(P + Q); the slope of the chord plus 1): at the same 10 primes the
        # order certificate must fail an assertion or give an a_p different from counting points
        global ec_add
        good_add = ec_add
        def minus_sum(P, Q, A, p):
            R = good_add(P, Q, A, p)
            return None if R is None else (R[0], (-R[1] - A[0]*R[0] - A[2]) % p)
        def slope_plus_one(P, Q, A, p):
            if P is None or Q is None or P[0] == Q[0]: return good_add(P, Q, A, p)
            a1, a2, a3, a4, a6 = A; (x1, y1), (x2, y2) = P, Q
            lam = ((y2 - y1) * pow((x2 - x1) % p, -1, p) + 1) % p; nu = (y1 - lam*x1) % p
            x3 = (lam*lam + a1*lam - a2 - x1 - x2) % p
            return (x3, (-(lam + a1)*x3 - nu - a3) % p)
        ps = [q for q in primes_upto(101000) if q > 10**5 and q not in LD['rows']][:10]
        for gname, law in (('P + Q -> -(P + Q)', minus_sum), ('chord slope + 1', slope_plus_one)):
            ec_add = law; det = 0
            try:
                for p in ps:
                    try: det += ap_order(a, p, random.Random(p))[0] != count_ap(a, p)
                    except AssertionError: det += 1
            finally: ec_add = good_add
            lines.append('  broken group law (%s): detected (an assertion fails or a_p differs from counting) at %d of 10 primes' % (gname, det))
            ntot += 1; nref += det == 10
    lines.append('negative tests: %d forged inputs, %d refused' % (ntot, nref))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert nref == ntot
    if write: open(os.path.join(HERE, 'lball_exact53_negative.txt'), 'w').write(out)


# ---------------------------------------------------------------- controls on curves where the answer is known
# (rank and root number as PARI gives them, ellanalyticrank and ellrootno; lball_exact.evaluate needs N > 57)
CONTROLS = [('67a1', (0, 1, 1, -12, -21), 'certified'),        # rank 0, w = +1
            ('79a1', (1, 1, 1, -2, 0), 'w not determined'),     # rank 1, w = -1: f(i/sqrt N) = 0
            ('389a1', (0, 1, 1, -2, 0), 'not positive'),       # rank 2, w = +1: L(E,1) = 0
            ('5077a1', (0, 0, 1, -7, 6), 'w not determined')]  # rank 3, w = -1

def controls(write):
    """the same certificate on four curves of Cremona's tables: it must certify the rank-0 curve and refuse the others
    (w = -1: F contains 0, so w = +1 is not shown; rank 2 with w = +1: the enclosure of L(E,1) contains 0)"""
    PI = LX.pi_interval(); lines = ['lball_exact53.py --controls: the certificate on curves with known rank', '']; ok = 0
    for label, a, expect in CONTROLS:
        LD = LX.local_data(a); N = LD['N']; M = LX.max_terms(N)
        ap, st = ap_table_fast(a, LD, M); an = LX.an_from_ap(ap, LD, M)
        ev = LX.evaluate(an, N, 1, PI, M); F = ev['F']; L1 = ev['L_y1']
        wF = 1 if F[0] > 0 or F[1] < 0 else None
        got = 'certified' if wF == 1 and L1[0] > 0 and ev['guard'] else ('w not determined' if wF is None else 'not positive')
        ok += got == expect
        lines.append('%-7s N = %d: F in %s, L(E,1) computed with w = +1 in %s, guard %s: %s (expected: %s)' % (
            label, N, LX.iv(F, 6), LX.iv(L1, 9), 'passes' if ev['guard'] else 'fails', got, expect))
    lines.append('controls: %d of %d as expected' % (ok, len(CONTROLS)))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert ok == len(CONTROLS)
    if write: open(os.path.join(HERE, 'lball_exact53_controls.txt'), 'w').write(out)


# ---------------------------------------------------------------- the order certificate against counting points
def aptest(write):
    rows = the_53(); lines = []; total = bad = fb = 0; d2 = 0
    try: import numpy; numpy_ok = True
    except ImportError: numpy_ok = False
    rnd = random.Random(2026)
    for i in (0, 7, 25):
        N, pair, a, LD = rows[i]
        ps = [p for p in primes_upto(30000) if p >= P0]
        ps += sorted(rnd.sample([p for p in primes_upto(2*10**6) if p > 30000], 150))
        n = d = f = 0
        for p in ps:
            model = LD['rows'][p]['model'] if p in LD['rows'] else a
            if p in LD['rows'] and LD['rows'][p]['red'] != 'good': continue
            x, k = ap_order(model, p, random.Random(p))
            if x is None: f += 1; continue
            y = count_ap(model, p)
            if numpy_ok and p < 100000: y2 = LX.ap_count(model, p); d2 += y2 != y
            n += 1; d += x != y
        lines.append('aptest (%d,%d): order certificate against counting points at the %d good primes %d <= p < 30000 and 150 random primes '
                     '30000 < p < 2e6: %d compared, %d disagreements, %d left to counting' % (pair[0], pair[1], len(ps) - 150, P0, n, d, f))
        total += n; bad += d; fb += f
    lines.append('aptest: %d primes compared, %d disagreements; the pure-Python count against lball_exact.ap_count (numpy) at p < 1e5: %s' % (
        total, bad, '%d disagreements' % d2 if numpy_ok else 'not run (no numpy)'))
    out = '\n'.join(lines) + '\n'; sys.stdout.write(out)
    assert bad == 0 and d2 == 0
    if write: open(os.path.join(HERE, 'lball_exact53_aptest.txt'), 'w').write(out)


if __name__ == '__main__':
    args = sys.argv[1:]
    def opt(name, default=None):
        return args[args.index(name) + 1] if name in args else default
    if '--list' in args:
        for i, (N, pair, a, LD) in enumerate(the_53()):
            print(i, '%d %d' % pair, N, LD['w_local'], LX.terms_needed(N, 1, True), LX.max_terms(N))
    elif '--merge' in args:
        merge([x for x in args[args.index('--merge') + 1:] if not x.startswith('--')], '--write' in args)
    elif '--apchunk' in args:
        i = args.index('--apchunk'); apchunk(int(args[i + 1]), int(args[i + 2]), int(args[i + 3]), args[i + 4])
    elif '--bigmerge' in args:
        i = args.index('--bigmerge'); bigmerge(int(args[i + 1]), int(args[i + 2]), args[i + 3], '--nopari' in args)
    elif '--bigcheck' in args:
        i = args.index('--bigcheck'); bigcheck(int(args[i + 1]), args[i + 2], int(args[i + 3]))
    elif '--negative' in args: negative('--write' in args)
    elif '--controls' in args: controls('--write' in args)
    elif '--aptest' in args: aptest('--write' in args)
    else:
        n53 = 53
        if opt('--only'): sel = [int(x) for x in opt('--only').split(',')]
        elif opt('--first'): sel = list(range(int(opt('--first'))))
        elif opt('--maxN'):
            lim = float(opt('--maxN')); sel = [i for i, (N, pair, a, LD) in enumerate(the_53()) if N <= lim]
        else: sel = list(range(n53))
        run(sel, '--nopari' in args, opt('--jsonl'))
