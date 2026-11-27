"""H3 (part 2): the good primes that can matter at all, for every ratio (EXPERIMENT, written by an AI model; not read by
a person; nothing here changes the census).

WITHDRAWN (review z8, run 7; see NOTES.md): the STATEMENT and PROOF below use the wrong curve.  The label is the
class of the linear factor l_j^+ at points where every F_j is a non-zero square; the curve carrying a label vector at
p = 1 mod 4 has genus 129, good primes above 1153 do remove vectors, and good primes p = 3 mod 4 can empty a set.  The
tables this script computes for p <= 1153 are right.

STATEMENT.  At a good prime p >= 1154, the four-class set loses nothing: every one of the 16 label vectors is allowed.
PROOF (Hasse-Weil, quoted).  For a label vector e, let u be a non-square unit of Z_p and C_e the curve
y_j^2 = u^(e_j) F_j(t), j = 0, 8, L, U.  At a good prime the eight roots of the four forms are distinct mod p and the
leading coefficients are units, so C_e is the smooth (Z/2)^4-cover of P^1 branched at 8 points, of genus
1 + 2^3 (8/2 - 2) = 17, with smooth reduction (p odd, tame, branch points distinct mod p).  Hasse-Weil gives
#C_e(F_p) >= p + 1 - 34 sqrt(p) > 0 once sqrt(p) > 17 + sqrt(288), i.e. p >= 1154.  A point of C_e(F_p) lifts by Hensel to
a Z_p-point, whose t has F_j(t) in u^(e_j) Q_p*^2 (or F_j(t) = 0, where the class is free), so e is allowed.
At a good prime p = 3 mod 4 no tuple is lost either (G2: every divisor has the same label there).
So the four-class set of ANY ratio is cut out by its bad primes, the real place, and the good primes p = 1 mod 4 below 1154.

COMPUTATION.  For each p = 1 mod 4 with 17 <= p < 1154 and each good residue r = m/n mod p, the allowed set, by
h3_tables.method_B (enumeration of P^1(F_p)) vectorised with integer numpy arrays (exact: residues mod p, no floating
point).  Checks: equal to h3_tables.json for p <= 197, and to G1's exact local test (method A) on every good residue
at every larger prime.  Output: per prime, the residues at which some vector is lost, and the largest prime at which
anything is lost.
usage: h3_weil.py [--write]   (h3_weil.txt, h3_weil.json)      h3_weil.py --negative   (h3_weil_negative.txt)"""
import sys, os, json, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv
_argv = sys.argv; sys.argv = [sys.argv[0]]
import h3_tables as H
sys.argv = _argv
LIMIT = 1154
PS = [p for p in range(17, LIMIT) if p % 4 == 1 and all(p % q for q in range(2, int(p ** 0.5) + 1))]
FULL = (1 << 16) - 1

def hw_threshold():
    """the least integer N with N + 1 - 34 sqrt(N) > 0 for all p >= N, in exact integer arithmetic:
    p + 1 > 34 sqrt(p)  <=>  (p + 1)^2 > 1156 p  (both sides positive)"""
    N = 1
    while not ((N + 1) ** 2 > 1156 * N and all((q + 1) ** 2 > 1156 * q for q in range(N, N + 2000))): N += 1
    return N

def method_B_fast(r, p, chi_tab):
    i = H.X.sqrt_m1(p); assert (i * i + 1) % p == 0
    forms, sup, known = H.TD.magic(r, 1)
    ls = []
    for (a0, a1), (b0, b1) in forms:
        ls.append([((a0 + i * a1) % p, (b0 + i * b1) % p), ((a0 - i * a1) % p, (b0 - i * b1) % p)])
    roots = []
    for j in range(4):
        for s in range(2):
            a, b = ls[j][s]; assert (a, b) != (0, 0)
            roots.append(((-b) * pow(a, -1, p) % p, 1) if a else (1, 0))
    assert len(set(roots)) == 8, 'two of the eight roots agree: not a good prime'
    x = np.arange(p, dtype=np.int64)
    Z = [[(ls[j][s][0] * x + ls[j][s][1]) % p for s in range(2)] for j in range(4)]       # points (x, 1)
    Zinf = [[ls[j][s][0] % p for s in range(2)] for j in range(4)]                         # the point (1, 0)
    C = [[chi_tab[Z[j][s]] for s in range(2)] for j in range(4)]                          # 0 square, 1 not, 2 zero
    vecs = set()
    # points where no linear form vanishes
    nz = np.ones(p, dtype=bool)
    for j in range(4):
        for s in range(2): nz &= Z[j][s] != 0
    agree = nz.copy()
    for j in range(4): agree &= C[j][0] == C[j][1]
    code = sum(C[j][0].astype(np.int64) << j for j in range(4))
    for c in np.unique(code[agree]): vecs.add(tuple((int(c) >> j) & 1 for j in range(4)))
    # the eight root points (and the point at infinity, done by the general rule below)
    pts = [int(xx) for xx in np.nonzero(~nz)[0]]
    def chi1(v): return 0 if pow(v, (p - 1) // 2, p) == 1 else 1
    def one(zs):
        z0 = [(j, s) for j in range(4) for s in range(2) if zs[j][s] == 0]
        if not z0:
            if all(chi1(zs[j][0]) == chi1(zs[j][1]) for j in range(4)): return tuple(chi1(zs[j][0]) for j in range(4))
            return None
        (j0, s0), = z0
        v = []
        for j in range(4):
            if j == j0: v.append(chi1(zs[j][1 - s0]))
            elif chi1(zs[j][0]) == chi1(zs[j][1]): v.append(chi1(zs[j][0]))
            else: return None
        return tuple(v)
    for xx in pts:
        t = one([[int(Z[j][s][xx]) for s in range(2)] for j in range(4)])
        if t is not None: vecs.add(t)
    t = one(Zinf)
    if t is not None: vecs.add(t)
    vecs |= set(tuple(1 - b for b in v) for v in vecs)
    return H.mask_of(vecs)

def chi_table(p):
    t = np.ones(p, dtype=np.int64); t[0] = 2
    t[(np.arange(1, p, dtype=np.int64) ** 2) % p] = 0
    return t

def tables():
    T = {}
    for p in PS:
        ct = chi_table(p); T[p] = {}
        for r in range(1, p):
            if H.good(r, 1, p): T[p][r] = method_B_fast(r, p, ct)
    return T

def checks(T, nA=4, seed=11):
    out = []
    old = {int(p): {int(r): a for r, a in d.items()} for p, d in json.load(open(os.path.join(HERE, 'h3_tables.json'))).items()}
    for p in PS:
        if p in old: assert T[p] == old[p], ('differs from h3_tables.json', p)
    out.append('equal to h3_tables.json at every p <= 197: yes')
    rnd = random.Random(seed); nchk = 0
    for p in PS:
        if p in old: continue
        for r in sorted(T[p]):
            a = H.method_A(r, 1, p); assert a == T[p][r], ('method A disagrees', p, r); nchk += 1
    out.append('G1 exact local test (method A) agrees on all %d good residues at the primes 197 < p < %d' % (nchk, LIMIT))
    return out

def summary(T):
    N = hw_threshold()
    out = ['H3 (part 2): good primes p = 1 mod 4 below %d, where some label vector can be lost. Written by an AI model; not read by a person. Experiment, not counted.' % LIMIT,
           'Hasse-Weil threshold (genus 17): every good prime p >= %d keeps all 16 vectors (exact integer check of (p + 1)^2 > 1156 p).' % N]
    assert N <= LIMIT
    last = 0; tot = 0; lost_tot = 0
    for p in PS:
        lost = [r for r, a in T[p].items() if a != FULL]
        tot += len(T[p]); lost_tot += len(lost)
        if lost: last = p
        mean = sum(bin(a).count('1') for a in T[p].values()) / max(1, len(T[p]))
        out.append('  %4d: %4d good residues, %4d lose something; mean allowed %5.2f of 16' % (p, len(T[p]), len(lost), mean))
    out.append('largest prime at which some good residue loses a vector: %d' % last)
    out.append('so the four-class set of every ratio is cut out by its bad primes, the real place and the good primes p = 1 mod 4 with p <= %d (residues losing something: %d of %d)' % (last, lost_tot, tot))
    return out, last

def emit(lines, f):
    s = '\n'.join(lines) + '\n'; f = os.path.join(HERE, f)
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))

if __name__ == '__main__':
    if '--negative' in sys.argv:
        out = []
        # (1) a forged threshold: genus 16 would give a smaller N; the exact check must not accept p = N - 1
        N = hw_threshold(); out.append('threshold %d; p = %d fails (p + 1)^2 > 1156 p: %s' % (N, N - 1, not ((N) ** 2 > 1156 * (N - 1))))
        # (2) a forged allowed set (a lost vector put back) must be refused by method A
        for p, r in ((229, None), (409, None)):
            ct = chi_table(p)
            rs = [r for r in range(1, p) if H.good(r, 1, p)]
            r = next(r for r in rs if method_B_fast(r, p, ct) != FULL)
            a = method_B_fast(r, p, ct); miss = next(b for b in range(16) if not a >> b & 1)
            forged = a | (1 << miss) | (1 << (15 - miss))
            out.append('forged set at p = %d r = %d: %s' % (p, r, 'refused by the exact test' if H.method_A(r, 1, p) != forged else 'NOT REFUSED'))
        # (3) a bad residue must be refused by the root assertion
        try:
            p = 29; r = next(r for r in range(1, p) if not H.good(r, 1, p) and r % p)
            method_B_fast(r, p, chi_table(p)); out.append('bad residue at p = 29: NOT REFUSED')
        except AssertionError:
            out.append('bad residue at p = %d r = %d: refused (two roots agree)' % (p, r))
        emit(out, 'h3_weil_negative.txt')
    else:
        T = tables()
        lines, last = summary(T)
        lines += checks(T)
        jf = os.path.join(HERE, 'h3_weil.json')
        if WRITE: json.dump({str(p): {str(r): a for r, a in d.items()} for p, d in T.items()}, open(jf, 'w'), sort_keys=True)
        emit(lines, 'h3_weil.txt')
