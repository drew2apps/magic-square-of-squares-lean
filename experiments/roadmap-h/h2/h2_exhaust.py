"""H2, part 2: every coprime ratio with LO < m <= HI (EXPERIMENT, written by an AI model; not read by a person;
nothing here changes the census).
For each coprime m : n, 0 < n < m, that is not of the two degenerate types (m^2 + n^2 a square or twice a square):
  1. the repository's class criterion (class_selmer.empty, class_selmer_linalg.system / solve), unchanged except that
     its factorisations are PARI's with every factor proved prime (G1's g1_stress sets this);
  2. for a ratio that keeps a class: the four-class descent with every bad prime (exact local test, G1's
     exact_local), then good primes in increasing order: the smallest P such that the bad primes and the good primes
     <= P empty the set (0 when the bad primes alone do; found by bisection over the good primes to 1500, then 5000).
     A ratio not emptied by then is a SURVIVOR candidate.
Output, one line per ratio that keeps a class:  m,n dim log2|space| P (or SURVIVOR)  and a count line per m.  Worker i of N takes the m with (m - LO - 1) % N == i.
usage: h2_exhaust.py LO HI i N [--write]    (h2_ex_<LO>_<i>of<N>.txt; with --write it keeps a .part file and
       resumes from it after an interruption)
G1's and G2's folders: G1_DIR / G2_DIR, or else ../../roadmap-g/g1 and g2 (or ../g1, ../g2)."""
import sys, os, time
from math import gcd, isqrt, log2
HERE = os.path.dirname(os.path.abspath(__file__))
def find(env, name):
    d = os.environ.get(env, os.path.join(HERE, '..', '..', 'roadmap-g', name))
    return d if os.path.isdir(d) else os.path.join(HERE, '..', name)
sys.path.insert(0, find('G1_DIR', 'g1')); sys.path.insert(0, find('G2_DIR', 'g2'))
argv = sys.argv; sys.argv = [argv[0]]
import g1_stress as G                      # sets the proved factorisations of the criterion
import g2_killsets as K
sys.argv = argv
CL, C = G.CL, G.C

def need(m, n):
    """(log2 of the tuple space, the smallest good-prime cutoff P that empties the set, or None)."""
    X = K.X; forms, sup, known = K.TD.magic(m, n)
    divs = K.divisors(sup); space = 1
    for d in divs: space *= len(d)
    bad = set()
    for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): bad |= set(K.factorint_proved(x))
    bad = sorted(bad); cache = {}
    def loc(p):
        if p not in cache:
            L, al = X.allowed_at(forms, p)
            for kt in known: assert tuple(L.dlabel(d) for d in kt) in al, ('known point rejected', m, n, p)
            cache[p] = (p, L, al) if len(al) < len(L.labels())**len(forms) else None
        return cache[p]
    def S_of(ps): return X.join(divs, [c for c in (loc(p) for p in ps) if c is not None])
    if not S_of(bad): return log2(space), 0
    for pmax in (1500, 5000):
        good = [int(q) for q in K.P.primes([3, pmax]) if int(q) not in bad]
        if not S_of(bad + good): break
    else: return log2(space), None
    lo, hi = 0, len(good)
    while lo < hi:
        mid = (lo + hi) // 2
        if S_of(bad + good[:mid]): lo = mid + 1
        else: hi = mid
    return log2(space), good[lo - 1]

def one(m, n):
    sy = CL.system(m, n)
    if sy is None: return None
    ok, dim = CL.solve(*sy)
    if not ok: return False
    ls, P = need(m, n)
    return '%d,%d %d %.0f %s' % (m, n, dim, ls, 'SURVIVOR' if P is None else P)

def main(lo, hi, i, N, part=None):
    out = []; t0 = time.time(); done = set()
    if part and os.path.exists(part):                 # resume: whole values of m already in the part file are kept
        for l in open(part):
            out.append(l.rstrip('\n'))
            if l.startswith('m='): done.add(int(l.split()[0][2:]))
        while out and not out[-1].startswith('m='): out.pop()
        open(part, 'w').write(''.join(x + '\n' for x in out))
    pf = open(part, 'a') if part else None
    for m in range(lo + 1 + i, hi + 1, N):
        if m in done: continue
        start = len(out)
        c = {'coprime': 0, 'degenerate': 0, 'no point': 0, 'no class': 0, 'class': 0}
        for n in range(1, m):
            if gcd(m, n) != 1: continue
            c['coprime'] += 1
            if G.is_deg(m, n): c['degenerate'] += 1; continue
            r = one(m, n)
            if r is None: c['no point'] += 1
            elif r is False: c['no class'] += 1
            else: c['class'] += 1; out.append(r)
        out.append('m=%d %s' % (m, ' '.join('%s=%d' % (k.replace(' ', '_'), v) for k, v in c.items())))
        if pf: pf.write(''.join(x + '\n' for x in out[start:])); pf.flush()
        if m % 50 == i % 50: sys.stderr.write('m=%d %.0f s\n' % (m, time.time() - t0))
    return out

if __name__ == '__main__':
    a = [x for x in argv[1:] if not x.startswith('--')]
    lo, hi, i, N = map(int, a[:4])
    f = os.path.join(HERE, 'h2_ex_%d_%dof%d.txt' % (lo, i, N))
    out = main(lo, hi, i, N, f + '.part' if '--write' in argv else None)
    s = '\n'.join(out) + '\n'
    if '--write' in argv: open(f, 'w').write(s); os.remove(f + '.part')
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
