"""I1: the constraints that bite (EXPERIMENT, written by an AI model; not read by a person; nothing here changes the census).
(a) From H3's table h3_weil.json (for each good prime p <= 1153 and residue r = m/n mod p, the 16-bit mask of the label
    vectors (one Legendre bit per class, classes in the order d0, d8, dL, dU) allowed at p): the union of the masks over
    all residues, i.e. the vectors that are allowed at p for SOME ratio, and the number of residues that allow nothing.
    A vector outside the union is never allowed at p, whatever m : n is, so a tuple with that vector at p needs p bad.
(b) The laws at a large bad prime q (q exactly once in one factor, the other factors prime to q): for each factor type
    (m, n, m - n, m + n, m^2 + n^2, and the two quartic factors Q+, Q- of m^8 + 14 m^4 n^4 + n^8) and q in fixed lists
    per class mod 8, six ratios each (fixed seed), the repository's exact local test (exact_local.allowed_at on
    tuple_descent.magic) gives the allowed tuples; printed: the unit-label vectors allowed (all valuations even) and
    the number of allowed tuples with an odd valuation.  Rows group by (type, q mod 8); 'varies' means the six ratios
    did not agree.
usage: i1_laws.py [--negative] [--write]   (i1_laws.txt, i1_laws_negative.txt)
Folders: MSS_REPO (the repository), H3_DIR (roadmap H's h3, default ../h3)."""
import sys, os, json, random, collections
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
H3 = os.environ.get('H3_DIR', next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-h', 'h3'), os.path.join(HERE, '..', 'h3')) if os.path.isdir(d)), os.path.join(HERE, '..', '..', 'roadmap-h', 'h3')))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'roadmap-g', 'g1')); sys.path.insert(0, os.path.join(REPO, 'scripts', 'cover_descent'))
argv = sys.argv; sys.argv = [argv[0]]
import exact_local as X, tuple_descent as TD
sys.argv = argv
OUT = []
def say(s=''): OUT.append(s)
def vec(i): return ''.join(str(i >> j & 1) for j in range(4))       # bit j = class j (d0, d8, dL, dU)

def table():
    return {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(H3, 'h3_weil.json'))).items()}

def part_a(T):
    say('(a) H3 table h3_weil.json: %d good primes <= %d' % (len(T), max(T)))
    for p in sorted(T):
        u = 0
        for a in T[p].values(): u |= a
        empty = sum(1 for a in T[p].values() if a == 0)
        if u != 0xFFFF or empty:
            miss = [vec(i) for i in range(16) if not u >> i & 1]
            say('  p = %d: never allowed %s; residues allowing nothing: %d of %d' % (p, miss or 'none', empty, len(T[p])))
    say('  at every other prime of the table every vector is allowed for some residue, and no residue allows nothing')

def Qp(m, n): return m**4 + 2*m**3*n + 2*m**2*n**2 - 2*m*n**3 + n**4
def Qm(m, n): return m**4 - 2*m**3*n + 2*m**2*n**2 + 2*m*n**3 + n**4
FAC = [('m', lambda m, n: m), ('n', lambda m, n: n), ('m-n', lambda m, n: m - n), ('m+n', lambda m, n: m + n),
       ('S', lambda m, n: m*m + n*n), ('Q+', Qp), ('Q-', Qm)]
QS = [313, 337, 409, 433, 317, 349, 373, 389, 307, 331, 347, 379, 311, 359, 383, 431]   # 4 primes in each class 1, 5, 3, 7 mod 8

def others_prime(m, n, q, ty):
    return all(f(m, n) % q for t, f in FAC if t != ty)

def part_b():
    rnd = random.Random(5)
    say('(b) laws at a bad prime q dividing one factor exactly once (exact local test, 6 ratios per type and q)')
    res = collections.defaultdict(list)
    for q in QS:
        for ty, f in FAC:
            roots = [r for r in range(q) if f(r, 1) % q == 0] if ty != 'n' else [None]
            if not roots: continue
            for k in range(6):
                while True:
                    if ty == 'n': m = rnd.randrange(q + 1, 10**6); n = q * rnd.randrange(1, m // q + 1)
                    else:
                        r0 = rnd.choice(roots); n = rnd.randrange(1, 10**6); m = (r0 * n) % q + q * rnd.randrange(1, 10**5)
                    if m > n > 0 and gcd(m, n) == 1 and f(m, n) % q == 0 and f(m, n) % (q*q) and others_prime(m, n, q, ty): break
                forms, _, _ = TD.magic(m, n)
                L, al = X.allowed_at(forms, q)
                unit = tuple(sorted(''.join(str(b) for (v, b) in a) for a in al if all(v == 0 for v, b in a)))
                odd = sum(1 for a in al if any(v for v, b in a))
                res[(ty, q % 8)].append((unit, odd))
    for (ty, c) in sorted(res, key=lambda k: ([t for t, _ in FAC].index(k[0]), k[1])):
        rs = res[(ty, c)]
        same = len(set(rs)) == 1
        unit, odd = rs[0]
        say('  %-3s q = %d mod 8 (%d ratios): %s' % (ty, c, len(rs), ('unit vectors %s; tuples with an odd valuation %d' % (list(unit), odd)) if same else 'varies'))

def negative():
    say('I1 laws, negative and positive controls. Written by an AI model; not read by a person.')
    T = table()
    u = 0
    for a in T[17].values(): u |= a
    say('forged claim "0000 is allowed at good 17 for some residue": %s' % ('refused' if not u & 1 else 'ACCEPTED')); assert not u & 1
    say('positive control: 0001 is allowed at good 17 for some residue: %s' % ('yes' if u >> 1 & 1 else 'NO')); assert u >> 1 & 1
    # at q = 313 (1 mod 8) dividing Q+ once, the exact local test allows the unit vectors 0000 and 1111 only (part (b))
    rnd = random.Random(1)
    roots = [r for r in range(313) if Qp(r, 1) % 313 == 0]
    while True:
        n = rnd.randrange(1, 10**6); m = (roots[0] * n) % 313 + 313 * rnd.randrange(1, 10**5)
        if m > n and gcd(m, n) == 1 and Qp(m, n) % (313 * 313) and others_prime(m, n, 313, 'Q+'): break
    forms, _, _ = TD.magic(m, n); L, al = X.allowed_at(forms, 313)
    unit = {''.join(str(b) for (v, b) in t) for t in al if all(v == 0 for v, b in t)}
    say('ratio (%d,%d), 313 | Q+ once: forged claim "0110 allowed": %s' % (m, n, 'refused' if '0110' not in unit else 'ACCEPTED')); assert '0110' not in unit
    say('positive control: 1111 allowed: %s' % ('yes' if '1111' in unit else 'NO')); assert '1111' in unit

if __name__ == '__main__':
    if '--negative' in argv: negative(); f = 'i1_laws_negative.txt'
    else:
        say('I1: the constraints that bite. Written by an AI model; not read by a person. Experiment, not counted.')
        part_a(table()); part_b(); f = 'i1_laws.txt'
    s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, f)
    if '--write' in argv: open(f, 'w').write(s)
    else:
        sys.stdout.write(s); print('# matches' if os.path.exists(f) and open(f).read() == s else '# DIFFERS')
