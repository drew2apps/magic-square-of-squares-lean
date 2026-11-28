"""I1: a designed survivor hunt (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
Design.  H3's tables (h3_weil.json) show that at the good primes 17, 41 and 73 the constant label vectors 0000 and 1111
are never allowed (at 17 only the eight odd-weight vectors are, and 8 of its 12 good residues allow nothing), so a
tuple of classes with all d_j in {1, 2} can only survive if 17, 41 and 73 are bad primes.  This hunt makes them bad by
construction: m + i n is a unit times a product of Gaussian primes, always including one above each of 17, 41 and 73
(so all three divide m^2 + n^2), plus K further Gaussian primes drawn from those above the primes p = 1 mod 4 with
p < 400 (random choice, random conjugate).  Then m, n = the absolute values of the real and imaginary parts, m > n.
Of each batch of B such ratios (degenerate ratios and those failing the congruences of modp_classes.txt skipped), the
one with the highest EXCESS (h2_design: log2 of the tuple space plus the H3 score at the good primes <= 1153) is taken
and goes through h2_hunt's pipeline: the repository's class criterion, then the four-class descent with every bad
prime exact and good primes by bisection to 5000 (P = 0: bad primes alone empty the set).  A ratio not emptied by 5000
is then tested at every good prime p = 1 mod 4 with 5000 < p <= 67583 (the corrected complete bound) by the exact
local test; still not emptied: SURVIVOR candidate (to be checked by tuple_bruteforce.py and local points).
Design C (I1_MODE=C): working backwards from H3's tables by CRT, see draw_C.
Design B (environment I1_MODE=B): m + i n a unit times K in 1..3 Gaussian primes above primes = 1 mod 8 below 10^6,
so m^2 + n^2 has few primes, each with 2 a square modulo it (the class criterion's condition there for d = 2).
Output i1_hunt_s<seed>.txt (design B: i1_hunt_B_s<seed>.txt): a count line, then one line per taken ratio: m,n K excess result.
usage: i1_hunt.py seed [TAKE] [--write]   (--write keeps a .part file and resumes after an interruption)
Folders: H2_DIR (default ../../mss-experiments/h2 layout: ../h2), H3_DIR (../h3); the repository: MSS_REPO."""
import sys, os, json, math, random, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
H2 = os.environ.get('H2_DIR', next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-h', 'h2'), os.path.join(HERE, '..', 'h2')) if os.path.isdir(d)), os.path.join(HERE, '..', '..', 'roadmap-h', 'h2'))); sys.path.insert(0, H2)
os.environ.setdefault('H3_DIR', next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-h', 'h3'), os.path.join(HERE, '..', 'h3')) if os.path.isdir(d)), os.path.join(HERE, '..', '..', 'roadmap-h', 'h3')))
argv = sys.argv; sys.argv = [argv[0]]
import h2_hunt as HH
import h2_design as HD
sys.argv = argv
P = HH.P; X = HH.X
B = 1500
MODE = os.environ.get('I1_MODE', 'forced')
if MODE == 'forced':      # design A: 17, 41, 73 forced into m^2 + n^2, K in 1..6 more Gaussian primes above p = 1 mod 4, p < 400
    POOL = [int(q) for q in P.primes([5, 400]) if int(q) % 4 == 1]; FORCED = (17, 41, 73); KR = (1, 6)
elif MODE == 'C':         # design C: CRT backwards from H3's tables (see draw_C)
    POOL = []; FORCED = (17, 41, 73); KR = (1, 1)
else:                     # design B: m^2 + n^2 = a product of K in 1..3 primes = 1 mod 8 (so 2 is a square modulo each), p < 10^6
    POOL = [int(q) for q in P.primes([17, 10**6]) if int(q) % 8 == 1]; FORCED = (); KR = (1, 3)
BOUND = 67583

def gprime(p):
    for a in range(1, int(p ** 0.5) + 1):
        b = math.isqrt(p - a * a)
        if b * b == p - a * a: return (a, b)
GP = {}
def gp(q):
    if q not in GP: GP[q] = gprime(q)
    return GP[q]
def gm(x, y): return (x[0] * y[0] - x[1] * y[1], x[0] * y[1] + x[1] * y[0])

def draw(rnd, K):
    qs = list(FORCED) + rnd.sample([q for q in POOL if q not in FORCED], K)
    z = (1, 0)
    for q in qs:
        a, b = gp(q); z = gm(z, (a, b) if rnd.random() < 0.5 else (a, -b))
    m, n = sorted((abs(z[0]), abs(z[1])), reverse=True)
    return m, n

WEIL = None
CRT_P = (29, 37, 53, 61)
def draw_C(rnd, K):
    """design C: a target pattern t of d_j in {1, 2} (bit j = 1 when d_j = 2); 17, 41, 73 made bad by putting each in a
    random one of m, n, m - n, m + n; at 29, 37, 53, 61 a residue m/n chosen at random among those whose H3 mask allows
    the pattern's label vector there (t at p = 5 mod 8, 0000 at p = 1 mod 8); n random below 10^4, m by CRT, lifted by a
    random multiple of the modulus below 10"""
    global WEIL
    if WEIL is None:
        WEIL = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(os.environ['H3_DIR'], 'h3_weil.json'))).items()}
    t = rnd.randrange(16)
    where = {p: rnd.choice(('m', 'n', 'm-n', 'm+n')) for p in FORCED}
    n = rnd.randrange(1, 10**4)
    for p, w in where.items():
        if w == 'n' and n % p: n *= p
    res = []
    for p, w in where.items():
        if w == 'n': continue
        res.append((p, {'m': 0, 'm-n': n % p, 'm+n': -n % p}[w]))
    for p in CRT_P:
        v = t if p % 8 == 5 else 0
        rs = [r for r, a in WEIL[p].items() if a >> v & 1]
        if n % p == 0: continue
        res.append((p, rnd.choice(rs) * n % p))
    M = 1; m0 = 0
    for p, a in res:
        m0 = m0 + M * ((a - m0) * pow(M, -1, p) % p); M *= p
    m = m0 + M * rnd.randrange(1, 10)
    if m < n: m, n = n, m
    return m, n

def extend(m, n):
    """good p = 1 mod 4, 5000 < p <= BOUND, after the bad primes and good primes <= 5000: the first p that empties"""
    forms, sup, known = HH.K.TD.magic(m, n)
    divs = HH.K.divisors(sup)
    bad = sorted(HH.primes_of(2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8))
    cs = []
    for p in bad + [int(q) for q in P.primes([3, 5000]) if int(q) not in bad]:
        L, al = X.allowed_at(forms, p)
        if len(al) < len(L.labels()) ** len(forms): cs.append((p, L, al))
    for q in P.primes([5001, BOUND]):
        p = int(q)
        if p % 4 != 1 or p in bad: continue
        L, al = X.allowed_at(forms, p)
        if len(al) < len(L.labels()) ** len(forms):
            cs.append((p, L, al))
            if not HH.H4.chunked(divs, cs, True): return p
    return None

def run(seed, take, part):
    rnd = random.Random(seed)
    old = [l.rstrip('\n') for l in open(part)] if part and os.path.exists(part) else []
    pf = open(part, 'a') if part else None
    c = collections.Counter(); lines = []
    for i in range(take):
        K = rnd.randint(*KR); best = None
        for _ in range(B):
            m, n = (draw_C if MODE == 'C' else draw)(rnd, K)
            if m == n or HH.gcd(m, n) != 1: continue
            if HH.G.is_deg(m, n) or HH.inclass(m, n): continue
            sc = HD.score(m, n)
            if sc == float('-inf'): continue
            forms, sup, known = HH.K.TD.magic(m, n)
            ex = sum(math.log2(len(d)) for d in HH.K.divisors(sup)) + sc
            if best is None or ex > best[0]: best = (ex, m, n)
        if best is None: c['taken'] += 1; c['nothing passed the congruences'] += 1; lines.append('- %d none' % K); continue
        ex, m, n = best
        if i < len(old): rec = old[i]
        else:
            sy = HH.timed(HH.CL.system, HH.TLIM_CLASS, m, n)
            if sy == 'TIMEOUT': res = 'timeout class'
            elif sy is None: res = 'no point'
            else:
                ok, dim = HH.CL.solve(*sy)
                if not ok: res = 'no class'
                else:
                    r = HH.timed(HH.need, HH.TLIM_DESC, m, n)
                    if r == 'TIMEOUT': res = 'class dim %d timeout descent' % dim
                    elif r[1] is not None: res = 'class dim %d %.0f %s' % (dim, r[0], r[1])
                    else:
                        e = extend(m, n)
                        res = 'class dim %d %.0f %s' % (dim, r[0], 'SURVIVOR' if e is None else e)
            rec = '%d,%d %d %.1f %s' % (m, n, K, ex, res)
            if pf: pf.write(rec + '\n'); pf.flush()
        res = rec.split(' ', 3)[3]
        c['taken'] += 1
        if res.startswith('class'): c['class'] += 1
        elif res == 'no point': c['no point'] += 1
        elif res == 'no class': c['no class'] += 1
        if 'timeout' in res: c['timeout'] += 1
        lines.append(rec)
    head = [('I1 designed hunt A, seed %d: 17, 41, 73 forced into m^2+n^2, K in 1..6 further Gaussian primes below 400; best of %d by the excess' if MODE == 'forced' else
             'I1 designed hunt B, seed %d: m^2+n^2 a product of K in 1..3 primes = 1 mod 8 below 10^6; best of %d by the excess' if MODE == 'B' else
             'I1 designed hunt C, seed %d: CRT backwards from H3 (17, 41, 73 in m, n or m +- n; residues at 29, 37, 53, 61 allowing a target pattern); best of %d by the excess') % (seed, B),
            'taken %d; no local point %d; no class %d; class kept %d; timed out %d' % (c['taken'], c['no point'], c['no class'], c['class'], c['timeout']),
            'SURVIVORS: %s' % ([l.split()[0] for l in lines if l.endswith('SURVIVOR')] or 'none')]
    return head + lines

if __name__ == '__main__':
    a = [x for x in argv[1:] if not x.startswith('--')]
    seed = int(a[0]); take = int(a[1]) if len(a) > 1 else 200
    f = os.path.join(HERE, 'i1_hunt_%ss%d.txt' % ({'forced': '', 'B': 'B_', 'C': 'C_'}[MODE], seed)); part = f + '.part' if '--write' in argv else None
    t0 = time.time(); out = run(seed, take, part); sys.stderr.write('seed %d %.0f s\n' % (seed, time.time() - t0))
    s = '\n'.join(out) + '\n'
    if '--write' in argv: open(f, 'w').write(s); os.remove(part)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
