"""I4: the exchange test on every point that G3's point searches found (EXPERIMENT, written by an AI model; not read by
a person; nothing here changes the census).
The five forms n^2u^2 + m^2v^2, m^2u^2 + n^2v^2, S(u^2 + v^2), S(u^2 + v^2) -+ 4mnuv (S = m^2 + n^2) are symmetric under
(m, n, u, v) <-> (u, v, m, n), so a point (m, n, u, v) of C(m : n) is also a point of C(u : v): a magic square with
first ratio m : n has second ratio u : v, and u : v has to be a ratio of its own.  The point searches of G3's pipeline
report rational values t = v/u: the second descent (census-wide/sd_pari.py, 'rational_t') and Demjanenko-Manin
(census2000/dm_pari.py, 't'), each classified 'no', 'degenerate' or 'SOLUTION' (the same nine-cell test).  For each
reported point with a verdict other than 'no' this script
  1. recomputes the verdict at (m, n, t) in exact rationals;
  2. recomputes it at the exchanged point, (u, v, m, n), i.e. ratio u : v and t' = n/m, and requires the same verdict
     (the nine cells are the same up to the factor u^2/m^2, so a mismatch means a wrong t <-> (u : v) convention or a
     wrong verdict);
  3. sorts the point: 'trivial' when uv(u^2 - v^2) = 0 (t in {0, 1, -1, infinity}: u : v is not a ratio, the point is one
     of the degenerate points of the census); otherwise u : v is reduced to u' > v' > 0 coprime (the census's
     symmetries: signs, and exchanging u and v) and, for a SOLUTION, has to be a valid ratio: not degenerate
     (g1_stress.is_deg), outside every class of the repository's modp_classes.txt, and keeping a class (class_selmer's
     criterion).  A SOLUTION whose u : v fails this is printed as EXCHANGE FAILS.
Output i4_exchange.txt: per chunk file, the number of points by verdict and sort, and every nontrivial point.
usage: i4_exchange.py [DIR] [--negative] [--write]   DIR: a folder of sd_*.jsonl / dm_*.jsonl files (default: this one);
       with DIR the output is printed only.  Folders: MSS_REPO; G1_DIR, G2_DIR, H2_DIR, H4_DIR as in roadmap H."""
import sys, os, json, glob
from fractions import Fraction as Fr
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
H2 = os.environ.get('H2_DIR', os.path.join(REPO, 'experiments', 'roadmap-h', 'h2'))
_hh = None
def HH():
    """h2_hunt (the class criterion, modp_classes.txt, is_deg), loaded only when a SOLUTION is to be tested."""
    global _hh
    if _hh is None:
        os.environ.setdefault('G1_DIR', os.path.join(REPO, 'experiments', 'roadmap-g', 'g1'))
        os.environ.setdefault('G2_DIR', os.path.join(REPO, 'experiments', 'roadmap-g', 'g2'))
        os.environ.setdefault('H4_DIR', os.path.join(REPO, 'experiments', 'roadmap-h', 'h4'))
        sys.path.insert(0, H2); argv = sys.argv; sys.argv = [argv[0]]
        import h2_hunt as h; sys.argv = argv; _hh = h
    return _hh

def issq(q): return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator
def classify(m, n, t):
    """census2000/dm_pari.classify_t and census-wide/sd_pari.classify_t, in exact rationals (the same nine cells)."""
    t = Fr(t); s = Fr(m*m + n*n, 2)
    F0 = m*m + n*n*t*t; F8 = n*n + m*m*t*t; F4 = s*(1 + t*t); L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
    if not all(issq(Fr(v)) for v in (F0, F4, F8, L, U)): return 'no'
    cells = [(m*t + n)**2, (m*t - n)**2, (n*t + m)**2, (n*t - m)**2, F0, F8, F4, L, U]
    return 'SOLUTION' if len(set(cells)) == 9 and all(c != 0 for c in cells) else 'degenerate'

def valid_ratio(a, b):
    """(ok, why) for the ratio a : b, a > b > 0 coprime."""
    h = HH()
    if h.G.is_deg(a, b): return False, 'degenerate ratio'
    if h.inclass(a, b): return False, 'in a class of modp_classes.txt'
    sy = h.CL.system(a, b)
    if sy is None: return False, 'no local point (class criterion)'
    ok, dim = h.CL.solve(*sy)
    return (True, 'keeps a class (dim %d)' % dim) if ok else (False, 'no class')

def exchange(m, n, t, verdict):
    """the exchange test for the point t = v/u of C(m : n) reported with this verdict; returns (sort, detail)."""
    t = Fr(t)
    if classify(m, n, t) != verdict: return 'FAILS', 'verdict %s does not recompute (got %s)' % (verdict, classify(m, n, t))
    if verdict == 'no': return 'no', ''
    u, v = t.denominator, t.numerator                       # t = v/u, (u, v) coprime
    if u * v * (u*u - v*v) == 0: return 'trivial', 'u : v = %d : %d' % (u, v)
    if classify(u, v, Fr(n, m)) != verdict: return 'FAILS', 'exchanged point (%d, %d, %d, %d): %s' % (u, v, m, n, classify(u, v, Fr(n, m)))
    a, b = max(abs(u), abs(v)), min(abs(u), abs(v)); g = gcd(a, b); a, b = a // g, b // g
    if verdict == 'degenerate': return 'nontrivial degenerate', 'u : v = %d : %d (is_deg %s)' % (a, b, HH().G.is_deg(a, b))
    ok, why = valid_ratio(a, b)
    return ('SOLUTION, exchange holds' if ok else 'EXCHANGE FAILS'), 'u : v = %d : %d: %s' % (a, b, why)

def points(path):
    """(m, n, t, verdict) for every point a jsonl file reports, partners included."""
    out = []
    for l in open(path):
        if not l.startswith('{'): continue
        o = json.loads(l)
        for rec in (o, o.get('partner') or {}):
            if 'pair' not in rec: continue
            m, n = rec['pair']
            for t, ver in (rec.get('rational_t') or {}).items(): out.append((m, n, t, ver))
            for t, ver in (rec.get('t') or {}).items(): out.append((m, n, t, ver))
    return out

def report(d):
    out = ['I4 exchange test on the points of G3\'s point searches. Written by an AI model; not read by a person. Experiment, not counted.']
    files = sorted(glob.glob(os.path.join(d, 'sd_*_*.jsonl')) + glob.glob(os.path.join(d, 'dm_*_*.jsonl')))
    tot = {}
    for f in files:
        c = {}; notes = []
        for m, n, t, ver in points(f):
            s, det = exchange(m, n, t, ver)
            c[s] = c.get(s, 0) + 1; tot[s] = tot.get(s, 0) + 1
            if s not in ('no', 'trivial'): notes.append('  (%d,%d) t = %s %s: %s; %s' % (m, n, t, ver, s, det))
        out.append('%s: %s' % (os.path.basename(f), ', '.join('%s %d' % kv for kv in sorted(c.items())) or 'no points'))
        out += notes
    out.append('all files: %s' % (', '.join('%s %d' % kv for kv in sorted(tot.items())) or 'no points'))
    out.append('EXCHANGE FAILS or FAILS: %d' % (tot.get('EXCHANGE FAILS', 0) + tot.get('FAILS', 0)))
    return out

def negative():
    out = ['I4 exchange test, negative and positive controls. Written by an AI model; not read by a person.']
    s, _ = exchange(14744, 3633, '1', 'degenerate')
    out.append('positive control: (14744,3633), t = 1 degenerate (reported by the second descent): %s' % s); assert s == 'trivial'
    s, _ = exchange(14744, 3633, '1', 'SOLUTION')
    out.append('forged: the same point reported as a SOLUTION: %s' % ('refused' if s == 'FAILS' else 'ACCEPTED')); assert s == 'FAILS'
    s, _ = exchange(14744, 3633, '2/3', 'degenerate')
    out.append('forged: t = 2/3 reported as degenerate: %s' % ('refused' if s == 'FAILS' else 'ACCEPTED')); assert s == 'FAILS'
    # a point that is a SOLUTION of the forms but whose u : v is not a valid ratio is what the test exists for; with no
    # such point known, the ratio check is exercised on its own
    for a, b in ((7, 1), (3, 2)):
        ok, why = valid_ratio(a, b)
        out.append('forged: a SOLUTION whose u : v is %d : %d: %s (%s)' % (a, b, 'refused' if not ok else 'ACCEPTED', why)); assert not ok
    ok, why = valid_ratio(308, 185)
    out.append('positive control: u : v = 308 : 185 (I1, design B): %s (%s)' % ('accepted' if ok else 'REFUSED', why)); assert ok
    return out

if __name__ == '__main__':
    a = [x for x in sys.argv[1:] if not x.startswith('--')]
    if '--negative' in sys.argv: out = negative(); f = os.path.join(HERE, 'i4_exchange_negative.txt')
    else: out = report(a[0] if a else HERE); f = None if a else os.path.join(HERE, 'i4_exchange.txt')
    s = '\n'.join(out) + '\n'
    if '--write' in sys.argv and f: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if f and os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
