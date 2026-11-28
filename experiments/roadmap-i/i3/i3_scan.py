"""I3: H1's phantom scan redone with every bad prime (EXPERIMENT, written by an AI model; not read by a person; nothing
here changes the census).
H1's scan (G1's descent) tested the primes of the supports (2mn(m^2 - n^2)) and the primes <= pmax, so a bad prime of
m^2 + n^2 or of the octic m^8 + 14 m^4 n^4 + n^8 above pmax was never tested.  Here the four-class descent tests every
bad prime: the primes of 2mn(m^2 - n^2), of m^2 + n^2 and of the octic (PARI's factorisation, each factor proved
prime), all exactly (G1's exact_local), then the good primes to 150, 1500, 5000, 20000 while more than the known tuples
is left.  For a phantom of H1 that is no longer one, the KILLING prime is given: the first prime, in increasing order
among those H1 did not test, after which no extra tuple is left.
usage: i3_scan.py k [COUNT] [--write]   H1's ratios for range k (G1's deg_ratios, seed 3000 + k, COUNT = 1177)
       i3_scan.py k --new [COUNT] [--write]   new ratios for range k (seed 5000 + k, COUNT = 2942; 17 ranges, 50,014)
       i3_scan.py --pairs [--write]     the 24 ratios of H1's 12 phantom pairs (h1_partners.txt's list)
       i3_scan.py --negative [--write]  controls: the new primes left out (H1's phantoms come back), a forged extra tuple
Folders: G1_DIR (the repository's experiments/roadmap-g/g1)."""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
G1 = os.environ.get('G1_DIR', os.path.join(HERE, '..', '..', 'roadmap-g', 'g1'))
sys.path.insert(0, G1)
argv = sys.argv; sys.argv = [argv[0]]
import g1_stress as G
import exact_local as X
sys.argv = argv
P = G.P
PAIRS = [(204724, 97893), (302617, 106831), (413004, 148603), (561607, 264401), (756492, 255635), (1012127, 500857),
         (853587, 255916), (1109503, 597671), (1992873, 999664), (2992537, 993209), (11875080, 6497849), (18372929, 5377231),
         (32396425, 3694248), (36090673, 28702177), (175775580, 15972419), (191747999, 159803161), (189089400, 156063271),
         (345152671, 33026129), (231641256, 72371617), (304012873, 159269639), (234319635, 115640732), (349960367, 118678903),
         (518804652, 31743595), (550548247, 487061057)]

def proved_primes(x):
    x = abs(x); s = set()
    if x > 1:
        for q in P.factor(x)[0]:
            q = int(q); assert P.isprime(q, 1) if q < 2**64 else P.isprime(q, 2); s.add(q)   # isprime(.,1/2): a certificate
    return s

def bad_primes(m, n):
    return proved_primes(2*m*n*(m*m - n*n)) | proved_primes(m*m + n*n) | proved_primes(m**8 + 14*m**4*n**4 + n**8)

def descent(m, n, every=True):
    forms, sup, known = G.TD.magic(m, n)
    base = set(q for S in sup for q in S)
    bad = bad_primes(m, n) if every else set()
    for pmax in (150, 1500, 5000, 20000):
        primes = sorted(base | bad | set(int(q) for q in P.primes([2, pmax])))
        S, used = X.descend_exact(forms, sup, known, primes)
        if S == known or not S: break
    return S, known, pmax, sorted(bad - base - set(int(q) for q in P.primes([2, pmax])))

def killing(m, n, old_primes, new_primes):
    forms, sup, known = G.TD.magic(m, n)
    ps = sorted(old_primes)
    for q in sorted(new_primes):
        ps.append(q)
        S, used = X.descend_exact(forms, sup, known, ps)
        if set(S) <= set(known): return q
    return None

def one(m, n, every=True):
    sy = G.CL.system(m, n)
    if sy is None: return 'no point', None
    ok, dim = G.CL.solve(*sy)
    if not ok: return 'no class', None
    return 'class dim %d' % dim, descent(m, n, every)

def line(m, n, every=True, h1=False):
    st, d = one(m, n, every)
    if d is None: return '(%d,%d) FAILED: %s (a degenerate ratio must keep its points)' % (m, n, st), False
    S, known, pmax, new = d
    extra = sorted(set(S) - set(known))
    if set(known) - set(S): return '(%d,%d) FAILED: a known tuple was removed' % (m, n), False
    if extra: return '(%d,%d) PHANTOM: %s; every bad prime; good primes to %d; known %s; extra %s' % (m, n, st, pmax, known, extra), True
    if h1:
        forms, sup, known = G.TD.magic(m, n)
        old = set(q for S_ in sup for q in S_) | set(int(q) for q in P.primes([2, 20000]))
        q = killing(m, n, old, new)
        bp = [x for x in new]
        kind = ' '.join('%s' % ('m^2+n^2' if (m*m + n*n) % x == 0 else 'octic') for x in [q]) if q else '-'
        return '(%d,%d) not a phantom: %s; extra tuples removed at %s (%s), a bad prime H1 did not test (untested bad primes: %d)' % (
            m, n, st, q, kind, len(bp)), False
    return None, False

def scan(k, count, base=3000):
    lo, hi = max(2**k, 10**4), min(2**(k + 1), 10**9 + 1)
    rs = G.deg_ratios(lo, hi, count, base + k)
    lines = []; nph = 0
    for m, n in rs:
        s, ph = line(m, n)
        nph += ph
        if s: lines.append(s)
    return ['range k = %d: %d <= m < %d; %d degenerate ratios (seed %d); every bad prime tested' % (k, lo, hi, len(rs), base + k),
            'ratios with tuples beyond those of their degenerate points: %d' % nph] + lines

def pairs():
    out = []; nph = 0
    for m, n in PAIRS:
        s, ph = line(m, n, h1=True); nph += ph; out.append(s)
    return ['H1\'s 12 phantom pairs (24 ratios), every bad prime tested', 'still phantoms: %d of 24' % nph] + out

def negative():
    out = []
    for m, n in ((853587, 255916), (11875080, 6497849)):
        s, ph = line(m, n, every=False)
        out.append('new primes left out: %s (%s)' % (s.split(' ')[0] + ': ' + ('PHANTOM' if ph else 'no phantom'),
                                                   'comes back as in H1, as it must' if ph else 'DOES NOT COME BACK'))
    m, n = 1109503, 597671
    forms, sup, known = G.TD.magic(m, n)
    forged = (1109503, 597671, 6, 1)
    primes = sorted(set(q for S_ in sup for q in S_) | bad_primes(m, n) | set(int(q) for q in P.primes([2, 150])))
    S, used = X.descend_exact(forms, sup, known, primes)
    out.append('forged tuple %s at (%d,%d): %s' % (forged, m, n, 'refused, as it must be' if forged not in S else 'KEPT (WRONG)'))
    return out

if __name__ == '__main__':
    t0 = time.time()
    if '--negative' in argv: lines = negative(); f = 'i3_negative.txt'
    elif '--pairs' in argv: lines = pairs(); f = 'i3_pairs.txt'
    else:
        a = [x for x in argv[1:] if not x.startswith('--')]
        k = int(a[0])
        if '--new' in argv: count = int(a[1]) if len(a) > 1 else 2942; lines = scan(k, count, 5000); f = 'i3_new_k%d.txt' % k
        else: count = int(a[1]) if len(a) > 1 else 1177; lines = scan(k, count); f = 'i3_scan_k%d.txt' % k
    sys.stderr.write('%s %.0f s\n' % (f, time.time() - t0))
    f = os.path.join(HERE, f); s = '\n'.join(lines) + '\n'
    if '--write' in argv: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
