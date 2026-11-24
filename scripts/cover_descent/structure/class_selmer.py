#!/usr/bin/env python3
"""The classes that are the class of a local point at every prime, from the laws alone (COMPLETE.md, section 7).

criterion(m, n) uses no local point: it factors mn(m^2-n^2), m^2+n^2 and the octic, and looks at residue symbols.
  no class at all if:  not (4 | mn or m = ±n mod 8)                      [no 2-adic point; classes modulo 8]
                       5 or 13 divides m^2+n^2 exactly once               [no 5-adic, 13-adic point]
                       m/n mod p is in a row E of class_good_table.txt    [no point at a good prime]
  otherwise the classes are d = 2 * (a product of primes = ±1 mod 8 of mn(m^2-n^2)) with
                       d a square modulo every odd prime of (m^2+n^2)(m^8+14m^4n^4+n^8)     [Theorems G, O]
                       d = ±2 mod 5 if 5 divides mn(m^2-n^2) exactly once                   [the law at 5]
                       d a residue (non-residue) modulo p if m/n mod p is in a row R (N)    [the table at good primes]
  (if 25 or 169 divides m^2+n^2 the existence of a point at 5 or 13 is taken from the test modulo 5^3, 13^2.)
direct(m, n) is the local test of lean_cert.py at every prime below 620, on the same candidates or on all products.

usage: class_selmer.py            the 425 ratios of exceptions1000, criterion against test (all products); a minute
       class_selmer.py --lists    the 21,104 ratios with m <= 5000 of tuple_681, tuple_scan_5000, tuple_trivial_5000
       class_selmer.py --all H    every coprime ratio with m <= H: how many have a class
"""
import sys, os, re, time, random
from itertools import combinations
from math import prod, gcd, isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
import lean_cert as L
try:
    import cypari2; _pari = cypari2.Pari()
    def fac(x): return [int(q) for q in _pari.factor(x)[0]]
    def sqfree(x): return prod(int(q) for q, e in zip(*[list(c) for c in _pari.factor(x)[:2]]) if int(e) % 2)
except ImportError:
    from sympy import factorint
    def fac(x): return sorted(factorint(x))
    def sqfree(x): return prod(q for q, e in factorint(x).items() if e % 2)
TAB = {}
for line in open(os.path.join(HERE, 'class_good_table.txt')):
    if not line.startswith('#'):
        a = line.rstrip('\n').split(' '); TAB[int(a[0])] = tuple(set(int(x) for x in a[i].split(',') if x) for i in (2, 4, 6))
TP = sorted(TAB)
def chi(a, p): return 1 if pow(a % p, (p - 1)//2, p) == 1 else -1
def K(p): return {2: 6, 3: 4, 5: 3}.get(p, 2 if p < 40 else 1)
def empty(m, n):
    """a reason why there is no class at all, or ''"""
    if not ((m*n) % 4 == 0 or (m - n) % 8 == 0 or (m + n) % 8 == 0): return '2'
    N = m*m + n*n
    for q in (5, 13):
        if N % q == 0:
            if N % (q*q): return str(q)
            if not L.survivors(m, n, q**K(q), q): return str(q)
    X = m*n*(m*m - n*n)
    for p in TP:
        if TAB[p][0] and X % p and N % p:
            tau = (m*pow(n, -1, p)) % p
            if tau in TAB[p][0] and (pow(m, 8, p) + 14*pow(m, 4, p)*pow(n, 4, p) + pow(n, 8, p)) % p: return 'good prime %d' % p
    return ''
def criterion(m, n):
    """(the classes, the reason if there is none because no point exists somewhere)"""
    why = empty(m, n)
    if why: return [], why
    X = m*n*(m*m - n*n); NW = (m*m + n*n)*(m**8 + 14*m**4*n**4 + n**8)
    P8 = [q for q in fac(X) if q > 2 and q % 8 in (1, 7)]
    cond = [(q, 1) for q in fac(NW) if q > 2]
    if L.exact5(m, n): cond.append((5, -1))
    for p in TP:
        if (X*NW) % p: 
            tau = (m*pow(n, -1, p)) % p
            if tau in TAB[p][1]: cond.append((p, 1))
            elif tau in TAB[p][2]: cond.append((p, -1))
    out = []
    for k in range(len(P8) + 1):
        for S in combinations(P8, k):
            d = 2*prod(S)
            if all(chi(d, q) == s for q, s in cond): out.append(d)
    return sorted(out), ''
SMALL = [q for q in range(2, 620) if all(q % r for r in range(2, isqrt(q) + 1))]
def direct(m, n, allprods):
    X = m*n*(m*m - n*n); NW = (m*m + n*n)*(m**8 + 14*m**4*n**4 + n**8)
    P = sorted(fac(2*X)); C = L.prods(P)
    if not allprods: C = [d for d in C if d % 2 == 0 and all(q % 8 in (1, 7) for q in P if q > 2 and d % q == 0)]
    C = [d for d in C if all(chi(d, q) == 1 for q in fac(NW) if q >= 620)]
    for p in SMALL:
        if not C: break
        M = p**K(p); surv = L.survivors(m, n, M, p); C = [d for d in C if not L.killed(d, surv, M)]
    return sorted(C)
def is_deg(m, n): N = m*m + n*n; return isqrt(N)**2 == N or isqrt(2*N)**2 == 2*N
def read(fn):
    out = []
    for line in open(os.path.join(HERE, '..', fn)):
        mm = re.match(r'\((\d+),(\d+)\): .*left: (\[.*?\]);', line)
        if mm: out.append((int(mm.group(1)), int(mm.group(2)), eval(mm.group(3))))
    return out
if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    if mode == '':
        src = open(os.path.join(HERE, '..', '..', '..', 'MagicSquareOfSquares/G1Census1000.lean')).read()
        body = src[src.index('def exceptions1000'):]; body = body[:body.index(']')]
        pairs = [tuple(int(x) for x in t) for t in re.findall(r'\((\d+), (\d+)\)', body)]
        same = 0; sizes = {}
        for m, n in pairs:
            a, why = criterion(m, n); b = direct(m, n, True)
            k = ('square or twice a square' if is_deg(m, n) else 'neither', len(a)); sizes[k] = sizes.get(k, 0) + 1
            if a == b: same += 1
            else: print('DIFFERENT (%d,%d): criterion %s %s; test %s' % (m, n, a, why, b))
        print('the 425 ratios of exceptions1000, every squarefree product of the primes of 2mn(m^2-n^2) as a candidate:')
        print('  the criterion and the local test at the primes below 620 leave the same classes for %d ratios of %d' % (same, len(pairs)))
        for k in sorted(sizes): print('  m^2+n^2 %s, %d class%s: %d ratios' % (k[0], k[1], '' if k[1] == 1 else 'es', sizes[k]))
    elif mode == '--lists':
        non = read('tuple_681.txt') + read('tuple_scan_5000.txt'); exc = read('tuple_trivial_5000.txt'); random.seed(11)
        cnt = {}; whys = {}; keep = []; emp = []
        for m, n, left in non:
            a, why = criterion(m, n); cnt[len(a)] = cnt.get(len(a), 0) + 1
            if why: whys[why.split()[0]] = whys.get(why.split()[0], 0) + 1
            (keep if a else emp).append((m, n, a))
        print('m^2+n^2 neither a square nor twice a square: %d ratios; classes left by the criterion: %s' % (len(non), ', '.join('%d: %d ratios' % kv for kv in sorted(cnt.items()))))
        print('  no point at all: %s' % ', '.join('at %s %d' % kv for kv in sorted(whys.items())))
        print('  with a class: ' + ' '.join('(%d,%d)%s' % (m, n, str(a).replace(' ', '')) for m, n, a in keep))
        bad = sum(direct(m, n, False) != a for m, n, a in keep); smp = random.sample(emp, 1500); bad0 = sum(direct(m, n, False) != [] for m, n, a in smp)
        print('  the local test below 620 (candidates with F and E): different for %d of the %d with a class, and for %d of 1500 drawn from the others' % (bad, len(keep), bad0))
        cnt = {}; eq = miss = bad = 0; more = []
        for m, n, left in exc:
            a, why = criterion(m, n); cnt[len(a)] = cnt.get(len(a), 0) + 1
            known = sorted(set(sqfree(prod(t)) for t in left)); miss += (not set(known) <= set(a)); bad += (direct(m, n, False) != a)
            if a == known: eq += 1
            else: more.append('(%d,%d)%s' % (m, n, str(a).replace(' ', '')))
        print('m^2+n^2 a square or twice a square: %d ratios; classes left by the criterion: %s' % (len(exc), ', '.join('%d: %d ratios' % kv for kv in sorted(cnt.items()))))
        print('  exactly the classes of the degenerate solutions: %d; a class of a degenerate solution missing: %d; the local test below 620 different: %d' % (eq, miss, bad))
        print('  with a further class: ' + ' '.join(more))
    elif mode == '--all':
        H = int(sys.argv[2]); step = max(H//4, 1); t0 = time.time()
        tot = {}; 
        for m in range(2, H + 1):
            b = (m - 1)//step
            for n in range(1, m):
                if gcd(m, n) != 1: continue
                r = tot.setdefault(b, {'ratios': 0, 'no point': 0, 'point, no class': 0, 'deg': 0, 'deg, more than its own classes': 0, 'not deg, a class': 0, 'list': []})
                r['ratios'] += 1; a, why = criterion(m, n); dg = is_deg(m, n)
                if why: r['no point'] += 1
                elif not a: r['point, no class'] += 1
                elif dg: r['deg'] += 1
                else: r['not deg, a class'] += 1; r['list'].append((m, n))
        for b in sorted(tot):
            r = tot[b]; print('m in (%d, %d]: coprime ratios %d; no local point somewhere %d; points everywhere but no class %d; a class, m^2+n^2 a square or twice a square %d; a class, neither %d' % (b*step, (b + 1)*step, r['ratios'], r['no point'], r['point, no class'], r['deg'], r['not deg, a class']))
        print('the ratios with a class and m^2+n^2 neither a square nor twice a square, m <= %d: %s' % (min(H, 400), ' '.join('(%d,%d)' % t for b in sorted(tot) for t in tot[b]['list'] if t[0] <= 400)))
