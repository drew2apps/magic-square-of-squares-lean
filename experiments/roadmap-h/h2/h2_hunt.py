"""H2, part 1: a targeted survivor hunt, m from 10^4 to 10^12 (EXPERIMENT, written by an AI model; not read by a person;
nothing here changes the census).
Sampler, per dyadic range k (2^k <= m < 2^(k+1), cut to [10^4, 10^12]), seed 4000 + k: random coprime m : n, 0 < n < m,
not of the two degenerate types, that pass the congruences of the repository's modp_classes.txt (no class at those
primes), and are kept only if
  * mn(m^2 - n^2) has at least OM(k) = 6 + k // 4 distinct primes (counted over m, n, m - n, m + n), and
  * the odd part of m^2 + n^2 has at most 2 distinct primes.
(The bias was chosen on a pilot at m ~ 10^6: it keeps a class about three times as often as plain random ratios; G1's
63 kept ratios have many primes in mn(m^2 - n^2) and few in m^2 + n^2.)  Each sampled ratio goes through:
  1. the repository's class criterion (class_selmer_linalg.system / solve; factorisations by PARI, each factor proved
     prime), with a time limit (TLIM_CLASS s); a ratio over the limit is counted as "timed out" and set aside;
  2. for a ratio that keeps a class: the four-class descent with every bad prime exact, then good primes in increasing
     order; P = the smallest cutoff such that bad primes and good primes <= P empty the set (0: bad primes alone),
     by bisection over the good primes to 1500, then 5000; the join is split by d0 (as in H4) to save memory; time
     limit TLIM_DESC s.  Not emptied by 5000: SURVIVOR candidate (to be checked by tuple_bruteforce.py and local points).
Each range stops at KEEP ratios that keep a class (or MAXS sampled ratios).  Later rounds: H2_ROUND, H2_MAXS (see below).
Output h2_hunt_k<k>.txt: a count line, then one line per kept ratio: m,n dim log2|space| P.
usage: h2_hunt.py k [KEEP] [--write]   (with --write a .part file is kept and resumed after an interruption)
G1's, G2's and H4's folders: G1_DIR / G2_DIR / H4_DIR, or else ../../roadmap-g/g1, g2 and ../h4."""
import sys, os, random, time, collections
from math import gcd, log2
HERE = os.path.dirname(os.path.abspath(__file__))
def find(env, name, base=('..', '..', 'roadmap-g')):
    d = os.environ.get(env, os.path.join(HERE, *base, name))
    return d if os.path.isdir(d) else os.path.join(HERE, '..', name)
sys.path.insert(0, find('G1_DIR', 'g1')); sys.path.insert(0, find('G2_DIR', 'g2')); sys.path.insert(0, find('H4_DIR', 'h4', ('..',)))
argv = sys.argv; sys.argv = [argv[0]]
import g1_stress as G
import g2_killsets as K
import h4_split as H4
sys.argv = argv
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
CL = G.CL; P = G.P; X = K.X
TLIM_CLASS, TLIM_DESC, MAXS = 600, 1800, 4000
# later rounds (environment H2_ROUND = r > 0): seed 4000 + k + 1000 r, output h2_hunt_k<k>_r<r>.txt, MAXS from H2_MAXS
ROUND = int(os.environ.get('H2_ROUND', '0')); MAXS = int(os.environ.get('H2_MAXS', MAXS))
CLS = collections.defaultdict(set)
for line in open(os.path.join(X.REPO, 'scripts', 'g1_full', 'modp_classes.txt')):
    p, r = line.split()[:2]; CLS[int(p)].add(int(r))
def inclass(m, n): return any(n % p and (m * pow(n, -1, p)) % p in rs for p, rs in CLS.items())
def primes_of(*xs):
    s = set()
    for x in xs: s |= set(int(q) for q in P.factor(abs(x))[0]) if abs(x) > 1 else set()
    return s

def sample(rnd, lo, hi, om):
    tries = 0
    while True:
        tries += 1
        m = rnd.randrange(lo, hi); n = rnd.randrange(1, m)
        if gcd(m, n) != 1 or G.is_deg(m, n) or inclass(m, n): continue
        if len(primes_of(m, n, m - n, m + n)) < om: continue
        S = m*m + n*n
        while S % 2 == 0: S //= 2
        if len(primes_of(S)) > 2: continue
        return m, n, tries

def need(m, n):
    forms, sup, known = K.TD.magic(m, n)
    divs = K.divisors(sup); space = 1
    for d in divs: space *= len(d)
    bad = sorted(primes_of(2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8))
    for q in bad: assert P.isprime(q)
    cache = {}
    def loc(p):
        if p not in cache:
            L, al = X.allowed_at(forms, p)
            for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, p)
            cache[p] = (p, L, al) if len(al) < len(L.labels())**len(forms) else None
        return cache[p]
    def S_of(ps): return H4.chunked(divs, [c for c in (loc(p) for p in ps) if c is not None], True)
    if not S_of(bad): return log2(space), 0
    for pmax in (1500, 5000):
        good = [int(q) for q in P.primes([3, pmax]) if int(q) not in bad]
        if not S_of(bad + good): break
    else: return log2(space), None
    lo, hi = 0, len(good)
    while lo < hi:
        mid = (lo + hi) // 2
        if S_of(bad + good[:mid]): lo = mid + 1
        else: hi = mid
    return log2(space), good[lo - 1]

def timed(fn, lim, *a):
    alarm(lim)
    try: r = fn(*a)
    except AlarmInterrupt: return 'TIMEOUT'
    finally: cancel_alarm()
    return r

def run(k, keep, part):
    lo, hi = max(2**k, 10**4), min(2**(k + 1), 10**12 + 1)
    om = 6 + k // 4; seed = 4000 + k + 1000 * ROUND; rnd = random.Random(seed)
    c = collections.Counter(); lines = []; done = 0
    if part and os.path.exists(part):
        old = [l.rstrip('\n') for l in open(part)]
        done = sum(1 for l in old if l.startswith('S '))
    pf = open(part, 'a') if part else None
    i = 0
    while c['class'] < keep and c['sampled'] < MAXS:
        m, n, tries = sample(rnd, lo, hi, om); i += 1
        c['tries'] += tries
        rec = None
        if i <= done:
            rec = REPLAY[i - 1]
        else:
            sy = timed(CL.system, TLIM_CLASS, m, n)
            if sy == 'TIMEOUT': rec = 'S %d,%d timeout class' % (m, n)
            elif sy is None: rec = 'S %d,%d no point' % (m, n)
            else:
                ok, dim = CL.solve(*sy)
                if not ok: rec = 'S %d,%d no class' % (m, n)
                else:
                    r = timed(need, TLIM_DESC, m, n)
                    if r == 'TIMEOUT': rec = 'S %d,%d class dim %d timeout descent' % (m, n, dim)
                    else: rec = 'S %d,%d class dim %d %.0f %s' % (m, n, dim, r[0], 'SURVIVOR' if r[1] is None else r[1])
            if pf: pf.write(rec + '\n'); pf.flush()
        c['sampled'] += 1
        kind = rec.split(' ', 2)[2]
        if kind.startswith('class'): c['class'] += 1; lines.append(rec[2:])
        if 'timeout' in kind: c['timeout ' + kind.split('timeout ')[1]] += 1
        elif kind == 'no point': c['no point'] += 1
        elif kind == 'no class': c['no class'] += 1
    head = ['range k = %d: %d <= m < %d (seed %d); bias: at least %d primes in mn(m^2-n^2), at most 2 odd primes in m^2+n^2' % (k, lo, hi, seed, om),
            'sampled %d (draws %d); no local point %d; no class %d; class kept %d; timed out: class %d, descent %d' % (
                c['sampled'], c['tries'], c['no point'], c['no class'], c['class'], c['timeout class'], c['timeout descent']),
            'SURVIVORS: %s' % ([l.split()[0] for l in lines if l.endswith('SURVIVOR')] or 'none')]
    return head + lines

REPLAY = []
if __name__ == '__main__':
    a = [x for x in argv[1:] if not x.startswith('--')]
    k = int(a[0]); keep = int(a[1]) if len(a) > 1 else 80
    f = os.path.join(HERE, 'h2_hunt_k%d%s.txt' % (k, '_r%d' % ROUND if ROUND else '')); part = f + '.part' if '--write' in argv else None
    if part and os.path.exists(part): REPLAY = [l.rstrip('\n') for l in open(part) if l.startswith('S ')]
    t0 = time.time(); out = run(k, keep, part); sys.stderr.write('k=%d %.0f s\n' % (k, time.time() - t0))
    s = '\n'.join(out) + '\n'
    if '--write' in argv: open(f, 'w').write(s); os.remove(part)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
