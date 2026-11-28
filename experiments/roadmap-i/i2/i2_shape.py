"""I2: what the bad primes leave (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
For each class-keeping ratio of a list: the four-class tuple space (d_j over the squarefree divisors of the supports),
the set S of tuples allowed at every bad prime (exact local test; the archimedean place of Q(i) is complex and imposes
nothing), and:
  |S|, log2 of the space, COSET = whether S is a coset of a subgroup (tuples as vectors over F_2 of the support primes),
  for S nonempty: the good primes (in increasing order, to 5000) that empty S, the killing prime P = the first that
  empties, and the tuples of S when |S| <= 8;
  KILL = an irredundant killing set (from the bad primes and the good primes <= P that cut S; found by dropping
  primes in decreasing order, P last kept, while the set stays empty) and, for each prime p in it, the smallest number of label coordinates (each class's label
  (v_p mod 2, Legendre bit), or the 2-adic label's entries) whose projection already separates the labels of the
  tuples allowed at the other primes of the set from the allowed labels at p.
Output i2_<name>.txt: one line per ratio:  m,n |S| log2space coset P kill=[p:k ...] [tuples]
usage: i2_shape.py LIST i N [--write]   LIST: a file of 'm,n ...' lines (h2's outputs); worker i of N.
Folders: G1_DIR, G2_DIR (the repository's experiments/roadmap-g), H2_DIR (roadmap H's h2)."""
import sys, os, time, itertools
HERE = os.path.dirname(os.path.abspath(__file__))
H2 = os.environ.get('H2_DIR', next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-h', 'h2'), os.path.join(HERE, '..', 'h2')) if os.path.isdir(d)), os.path.join(HERE, '..', '..', 'roadmap-h', 'h2'))); sys.path.insert(0, H2)
argv = sys.argv; sys.argv = [argv[0]]
import h2_hunt as HH
H4 = HH.H4
sys.argv = argv
P = HH.P; X = HH.X; K = HH.K
P.default('nbthreads', 1)      # PARI's threads cannot get their stacks under run_i2.sh's memory cap and then wait forever

def vec(t, sup):
    v = 0; b = 0
    for d, S in zip(t, sup):
        for q in S:
            if d % q == 0: v |= 1 << b
            b += 1
    return v

def is_coset(vs):
    if not vs: return None
    v0 = vs[0]; basis = []
    for v in vs[1:]:
        w = v ^ v0
        for b in basis: w = min(w, w ^ b)
        if w: basis.append(w)
    return 2 ** len(basis) == len(vs)

def sep(lab_left, allowed):
    """smallest number of label coordinates whose projection separates lab_left from allowed"""
    flat = lambda lab: tuple(x for c in lab for x in (c if isinstance(c, tuple) else (c,)))
    A = [flat(l) for l in lab_left]; B = [flat(l) for l in allowed]
    n = len(A[0])
    for k in range(0, n + 1):
        for cs in itertools.combinations(range(n), k):
            pa = {tuple(a[c] for c in cs) for a in A}
            if not any(tuple(b[c] for c in cs) in pa for b in B): return k
    return None

def shape(m, n):
    forms, sup, known = K.TD.magic(m, n)
    divs = K.divisors(sup); space = 1
    for d in divs: space *= len(d)
    bad = sorted(HH.primes_of(2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8))
    cache = {}
    def loc(p):
        if p not in cache:
            L, al = X.allowed_at(forms, p)
            cache[p] = (p, L, al) if len(al) < len(L.labels()) ** len(forms) else None
        return cache[p]
    def S_of(ps): return H4.chunked(divs, [c for c in (loc(p) for p in ps) if c is not None], False)   # = X.join, split by d0 to save memory
    S = S_of(bad)
    co = is_coset([vec(t, sup) for t in S])
    Pk = 0; chain = list(bad)
    if S:
        left = S
        for q in P.primes([3, 5000]):
            q = int(q)
            if q in bad or loc(q) is None: continue
            _, L, al = loc(q)
            new = [t for t in left if tuple(L.dlabel(d) for d in t) in al]
            if len(new) < len(left): chain.append(q)
            left = new
            if not left: Pk = q; break
        else: Pk = None
    kill = []
    if Pk is not None:
        ks = [p for p in chain if loc(p) is not None]
        for p in sorted(ks, reverse=True):
            if len(ks) > 1 and p != Pk or (p == Pk and Pk == 0):
                trial = [x for x in ks if x != p]
                if not S_of(trial): ks = trial
        for p in ks:
            rest = S_of([x for x in ks if x != p])
            _, L, al = loc(p)
            kill.append((p, sep([tuple(L.dlabel(d) for d in t) for t in rest], al) if rest else 0))
    return len(S), space, co, Pk, kill, (S if 0 < len(S) <= 8 else None)

def fmt(m, n, r):
    nS, space, co, Pk, kill, S = r
    from math import log2
    s = '%d,%d %d %.0f %s %s kill=[%s]' % (m, n, nS, log2(space), {None: '-', True: 'coset', False: 'not-coset'}[co],
                                         'NONE<=5000' if Pk is None else Pk, ' '.join('%d:%s' % kp for kp in sorted(kill)))
    if S: s += ' S=' + ';'.join(','.join(map(str, t)) for t in S)
    return s

def read(path):
    out = []
    for l in open(path):
        t = l.split()
        if not t or ',' not in t[0]: continue
        if 'class' in l or (len(t) >= 4 and t[1].isdigit()) or (len(t) == 2 and not t[1][0].isdigit()):
            try: m, n = map(int, t[0].split(','))
            except ValueError: continue
            if 'no class' in l or 'no point' in l: continue
            out.append((m, n))
    return out

if __name__ == '__main__':
    a = [x for x in argv[1:] if not x.startswith('--')]
    lst, i, N = a[0], int(a[1]), int(a[2])
    rs = [r for j, r in enumerate(read(lst)) if j % N == i]
    name = os.path.splitext(os.path.basename(lst))[0]
    f = os.path.join(HERE, 'i2_%s_%dof%d.txt' % (name, i, N)); part = f + '.part'
    done = [l.rstrip('\n') for l in open(part)] if '--write' in argv and os.path.exists(part) else []
    pf = open(part, 'a') if '--write' in argv else None
    out = []; t0 = time.time()
    for j, (m, n) in enumerate(rs):
        if j < len(done): out.append(done[j]); continue
        try: s = fmt(m, n, shape(m, n))
        except MemoryError: s = '%d,%d MEMORY: over the per-process memory limit of run_i2.sh, not done' % (m, n)
        out.append(s)
        if pf: pf.write(s + '\n'); pf.flush()
    sys.stderr.write('%s %d of %d: %d ratios %.0f s\n' % (name, i, N, len(rs), time.time() - t0))
    s = '\n'.join(out) + '\n'
    if '--write' in argv: open(f, 'w').write(s); os.remove(part)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
