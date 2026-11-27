"""H3: the local image of the four-class cover at good primes p = 1 mod 4, p <= 211 (EXPERIMENT, written by an AI
model; not read by a person; nothing here changes the census).

At a prime p = 1 mod 4 that divides none of m n (m^2 - n^2)(m^2 + n^2)(m^8 + 14 m^4 n^4 + n^8) (a GOOD prime), the
classes d_j of a tuple are units at p, so only their Legendre symbols (d_j / p) matter: a vector in F_2^4, and the set
allowed at p is closed under adding (1, 1, 1, 1) (a unit of Z_p that is not a square).  It depends only on the residue
r = m / n mod p (scaling m, n by a unit multiplies every form by it).  TABLE: for each such p and good residue r, the
allowed set, as a 16-bit mask (bit sum_j e_j 2^j for the vector e of symbols, 0 = square; j = 0, 8, L, U).

Two methods, which must agree on every residue:
  A. G1's exact local test (exact_local.ExactLocal, no depth cap) on the repository's forms tuple_descent.magic(r, 1),
     restricted to the label vectors of units;
  B. an enumeration over P^1(F_p): the forms split as l_j^+ l_j^- over F_p (sqrt(-1) mod p); a point x where no l
     vanishes gives the vector (chi(l_j^+(x)))_j when chi(l_j^+(x)) = chi(l_j^-(x)) for every j; at the root x of one
     l_j^s (the eight roots are distinct at a good prime, asserted), the class d_j is free on that component (Hensel:
     l_j^s = b p^e w with e even and w of either symbol), so the vector has chi(l_j^-s(x)) in place j; then closed
     under (1, 1, 1, 1).
CHECKER: random coprime (m, n) up to 10^12 and a random p of the list that is good for them; G1's exact test on
tuple_descent.magic(m, n) must give the table's set at r = m / n mod p (h3_check.txt).
usage: h3_tables.py [--write]       (h3_tables.json, h3_tables.txt)
       h3_tables.py --check [N]     (h3_check.txt: N random inputs, default 300)
       h3_tables.py --negative      (h3_negative.txt: forged tables that the checker must refuse)
G1's folder is found through G1_DIR, or else ../../roadmap-g/g1; the repository through MSS_REPO, or else ../../.."""
import sys, os, json, random
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
G1 = os.environ.get('G1_DIR', os.path.join(HERE, '..', '..', 'roadmap-g', 'g1'))
if not os.path.isdir(G1): G1 = os.path.join(HERE, '..', 'g1')
sys.path.insert(0, G1)
import exact_local as X
TD = X.TD
WRITE = '--write' in sys.argv
PS = [p for p in range(5, 212) if p % 4 == 1 and all(p % q for q in range(2, int(p ** 0.5) + 1))]

def good(m, n, p):
    return all(x % p for x in (m, n, m - n, m + n, m * m + n * n, m ** 8 + 14 * m ** 4 * n ** 4 + n ** 8))

def mask_of(vecs):
    s = 0
    for v in vecs: s |= 1 << sum(b << j for j, b in enumerate(v))
    return s

def method_A(m, n, p):
    forms, sup, known = TD.magic(m, n)
    L, al = X.allowed_at(forms, p)
    vecs = [tuple(l[1] for l in t) for t in al if all(l[0] == 0 for l in t)]
    return mask_of(vecs)

def method_B(r, p):
    i = X.sqrt_m1(p); assert (i * i + 1) % p == 0
    forms, sup, known = TD.magic(r, 1)
    def lin(f, s):            # l^s(u, v) = (a0 + s i a1) u + (b0 + s i b1) v mod p, for the form ((a0, a1), (b0, b1))
        (a0, a1), (b0, b1) = f
        return ((a0 + s * i * a1) % p, (b0 + s * i * b1) % p)
    ls = [[lin(f, 1), lin(f, -1)] for f in forms]
    pts = [(x, 1) for x in range(p)] + [(1, 0)]
    val = lambda l, P: (l[0] * P[0] + l[1] * P[1]) % p
    roots = []
    for j in range(4):
        for s in range(2):
            a, b = ls[j][s]; assert (a, b) != (0, 0)
            roots.append(((-b) * pow(a, -1, p) % p, 1) if a else (1, 0))
    assert len(set(roots)) == 8, 'two of the eight roots agree: not a good prime'
    chi = lambda x: 0 if pow(x, (p - 1) // 2, p) == 1 else 1
    vecs = set()
    for P in pts:
        zs = [[val(ls[j][s], P) for s in range(2)] for j in range(4)]
        z0 = [(j, s) for j in range(4) for s in range(2) if zs[j][s] == 0]
        if not z0:
            if all(chi(zs[j][0]) == chi(zs[j][1]) for j in range(4)): vecs.add(tuple(chi(zs[j][0]) for j in range(4)))
        else:
            (j0, s0), = z0
            v = []
            ok = True
            for j in range(4):
                if j == j0: v.append(chi(zs[j][1 - s0]))
                elif chi(zs[j][0]) == chi(zs[j][1]): v.append(chi(zs[j][0]))
                else: ok = False
            if ok: vecs.add(tuple(v))
    vecs |= set(tuple(1 - b for b in v) for v in vecs)
    return mask_of(vecs)

def tables():
    T = {}
    for p in PS:
        T[p] = {}
        for r in range(1, p):
            if not good(r, 1, p): continue
            a = method_A(r, 1, p); b = method_B(r, p)
            assert a == b, ('the two methods disagree', p, r, a, b)
            T[p][r] = a
    return T

FULL = (1 << 16) - 1
def summary(T):
    out = ['H3: allowed sets of the four-class cover at good primes p = 1 mod 4, p <= 211 (methods A and B agree on every residue).',
           'Written by an AI model; not read by a person. Experiment, not counted.',
           'p: good residues; residues where every class is allowed; sizes of the allowed sets (of 16 label vectors) that occur']
    for p in PS:
        sz = {}
        for r, a in T[p].items(): c = bin(a).count('1'); sz[c] = sz.get(c, 0) + 1
        out.append('  %3d: %3d good residues, full at %3d; sizes %s' % (p, len(T[p]), sz.get(16, 0), ' '.join('%d:%d' % (k, sz[k]) for k in sorted(sz))))
    out.append('primes with no good residue at all (always bad): %s' % [p for p in PS if not T[p]])
    never = [p for p in PS if T[p] and all(a == FULL for a in T[p].values())]
    out.append('primes with good residues at which no good residue removes anything: %s' % (never if never else 'none'))
    out.append('good residues where the allowed set is empty (no local point at p): %d of %d' % (sum(a == 0 for p in PS for a in T[p].values()), sum(len(T[p]) for p in PS)))
    return out

def check(N, T, seed=7):
    rnd = random.Random(seed); out = []; bad = 0
    for _ in range(N):
        while True:
            m = rnd.randrange(2, 10 ** 12); n = rnd.randrange(1, m); p = rnd.choice(PS)
            if gcd(m, n) == 1 and good(m, n, p): break
        r = m * pow(n, -1, p) % p
        a = method_A(m, n, p)
        ok = a == T[p][r]
        bad += not ok
        out.append('(%d,%d) p = %d r = %d: %s' % (m, n, p, r, 'agrees' if ok else 'DISAGREES'))
    return ['checker: %d random inputs, %d disagree' % (N, bad)] + out

def emit(lines, f):
    s = '\n'.join(lines) + '\n'; f = os.path.join(HERE, f)
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))

if __name__ == '__main__':
    jf = os.path.join(HERE, 'h3_tables.json')
    if '--check' in sys.argv or '--negative' in sys.argv:
        T = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(jf)).items()}
        if '--check' in sys.argv:
            k = sys.argv.index('--check'); N = int(sys.argv[k + 1]) if len(sys.argv) > k + 1 and sys.argv[k + 1].isdigit() else 2000
            emit(check(N, T), 'h3_check.txt')
        else:
            out = []
            for p, r in ((17, 2), (61, 7), (181, 5)):
                if r not in T[p]: r = sorted(T[p])[0]
                F = {q: dict(d) for q, d in T.items()}
                # forge: drop one allowed vector (and its flip), or add one that is not allowed
                a = F[p][r]
                if a == FULL: F[p][r] = a & ~3 & ~(1 << 15) & ~(1 << 14)
                else:
                    miss = next(b for b in range(16) if not a >> b & 1); F[p][r] = a | (1 << miss) | (1 << (15 - miss))
                got = method_A(r, 1, p)
                out.append('forged entry p = %d r = %d: %s' % (p, r, 'refused (differs from the exact test), as it must be' if got != F[p][r] else 'NOT REFUSED'))
            emit(out, 'h3_negative.txt')
    else:
        T = tables()
        if WRITE: json.dump({str(p): {str(r): a for r, a in d.items()} for p, d in T.items()}, open(jf, 'w'), sort_keys=True)
        elif os.path.exists(jf):
            old = json.load(open(jf)); print('# table matches h3_tables.json' if old == {str(p): {str(r): a for r, a in d.items()} for p, d in T.items()} else '# TABLE DIFFERS')
        emit(summary(T), 'h3_tables.txt')
