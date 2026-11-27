"""H2, part 3: a hunt steered by the good-prime tables of H3 (EXPERIMENT, written by an AI model; not read by a person;
nothing here changes the census).
H3 shows that the good primes that can cut the four-class set are exactly the p = 1 mod 4 with 17 <= p <= 1153, and
gives, for each such p and residue r = m/n mod p, the allowed set A_p(r) of the 16 label vectors.  The SCORE of a ratio
is sum over those p that are good for it of log2(|A_p(r)| / 16) (minus infinity if some A_p(r) is empty): the
expected loss, in bits, of a uniformly random tuple at the good primes.  The EXCESS is log2 of the tuple space plus
the score: the expected log2 of the number of tuples left by the good primes, if the tuples' labels were uniform and
independent (a heuristic, used only to choose ratios).  A survivor needs a large excess.
Sampler, per dyadic range k (2^k <= m < 2^(k+1)), seed 6000 + k: random coprime m : n, 0 < n < m, not of the two
degenerate types, passing the congruences of modp_classes.txt; of each batch of B = 200 such ratios the one with the
highest excess is taken (ratios with no local point at some good prime are skipped).  Each taken ratio goes through h2_hunt's pipeline (class criterion, then the four-class
descent with every bad prime exact and good primes by bisection; same time limits).
Output h2_design_k<k>.txt: a count line, then one line per taken ratio: m,n excess result (no point | no class |
class dim d log2|space| P, P = 0 for bad primes alone, SURVIVOR if not emptied by good primes to 5000).
usage: h2_design.py k [TAKE] [--write]   (TAKE ratios per range, default 300; --write keeps a .part file and resumes)
The tables are read from H3's folder (H3_DIR, or else ../h3)."""
import sys, os, json, math, random, time, collections
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
argv = sys.argv; sys.argv = [argv[0]]
import h2_hunt as HH
sys.argv = argv
H3 = os.environ.get('H3_DIR', os.path.join(HERE, '..', 'h3'))
T = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(H3, 'h3_weil.json'))).items()}
PS = sorted(T)
B = 200

def good(m, n, p):
    return all(x % p for x in (m, n, m - n, m + n, m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8))

def score(m, n):
    s = 0.0
    for p in PS:
        if not good(m, n, p): continue
        c = bin(T[p][m * pow(n, -1, p) % p]).count('1')
        if c == 0: return float('-inf')
        s += math.log2(c / 16)
    return s

def run(k, take, part):
    lo, hi = max(2 ** k, 10 ** 4), min(2 ** (k + 1), 10 ** 12 + 1)
    rnd = random.Random(6000 + k)
    old = [l.rstrip('\n') for l in open(part)] if part and os.path.exists(part) else []
    pf = open(part, 'a') if part else None
    c = collections.Counter(); lines = []
    for i in range(take):
        best = None
        for _ in range(B):
            while True:
                m = rnd.randrange(lo, hi); n = rnd.randrange(1, m)
                if HH.gcd(m, n) == 1 and not HH.G.is_deg(m, n) and not HH.inclass(m, n): break
            sc = score(m, n)
            if sc == float('-inf'): continue
            forms, sup, known = HH.K.TD.magic(m, n)
            ex = sum(math.log2(len(d)) for d in HH.K.divisors(sup)) + sc
            if best is None or ex > best[0]: best = (ex, m, n)
        if best is None: continue
        sc, m, n = best
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
                    res = ('class dim %d timeout descent' % dim) if r == 'TIMEOUT' else ('class dim %d %.0f %s' % (dim, r[0], 'SURVIVOR' if r[1] is None else r[1]))
            rec = '%d,%d %.1f %s' % (m, n, sc, res)
            if pf: pf.write(rec + '\n'); pf.flush()
        res = rec.split(' ', 2)[2]
        c['taken'] += 1
        if res.startswith('class'): c['class'] += 1
        elif res == 'no point': c['no point'] += 1
        elif res == 'no class': c['no class'] += 1
        if 'timeout' in res: c['timeout'] += 1
        lines.append(rec)
    kept = [l for l in lines if ' class dim ' in l]
    head = ['range k = %d: %d <= m < %d (seed %d); best of %d by the excess' % (k, lo, hi, 6000 + k, B),
            'taken %d; no local point %d; no class %d; class kept %d; timed out %d' % (c['taken'], c['no point'], c['no class'], c['class'], c['timeout']),
            'SURVIVORS: %s' % ([l.split()[0] for l in kept if l.endswith('SURVIVOR')] or 'none')]
    return head + lines

if __name__ == '__main__':
    a = [x for x in argv[1:] if not x.startswith('--')]
    k = int(a[0]); take = int(a[1]) if len(a) > 1 else 300
    f = os.path.join(HERE, 'h2_design_k%d.txt' % k); part = f + '.part' if '--write' in argv else None
    t0 = time.time(); out = run(k, take, part); sys.stderr.write('k=%d %.0f s\n' % (k, time.time() - t0))
    s = '\n'.join(out) + '\n'
    if '--write' in argv: open(f, 'w').write(s); os.remove(part)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
