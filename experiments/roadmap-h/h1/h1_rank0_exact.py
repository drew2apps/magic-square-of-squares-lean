"""H1, part 7: the rank-0 quotients of h1_rank0.py, with the rank bound computed exactly (EXPERIMENT, written by an AI
model; not read by a person; nothing here changes the census).
Two of the six quotients of C in h1_rank0.py have a model y^2 = x (x^2 + a x + b) with integer a, b:
  A: y^2 = w f0 f8, x = m^2 n^2 w:   y^2 = x (x + m^4)(x + n^4)          (full rational 2-torsion)
  B: y^2 = w K,     x = S^2 w (S = m^2 + n^2):  y^2 = x (x^2 + (2 S^2 - 16 m^2 n^2) x + S^4)
Rank bound: for A, the complete 2-descent and the Cassels-Tate pairing of the repository's rank_second.py
(full_selmer2, DFull.pairing), rank <= dim Sel^2 - 2 - rank of the pairing, taken over the curves of A's isogeny
class with full rational 2-torsion reached by 2-isogenies (Velu's formula); for B, the exact descent via 2-isogeny
with the exact second descent of rank_indep.py (rank_bound; its conic points from rank_second.find_conic_point_any).  No Magma; PARI only as a finder of conic points (checked).
Torsion: PARI's elltors as a finder; each point checked on the curve and of finite order by the exact group law, and
the number of torsion points equal to the gcd of #E(F_p) over the good odd primes p < 400 (torsion injects into E(F_p)
at a good odd prime), so the list is all of E(Q)_tors.
If the bound is 0: E(Q) is that list, so every point t of C has w = t^2 among its w-values (or t = infinity); the tuples
of C over them are computed exactly (h1_rank0.tuple_at) and an extra tuple that is none of them is the class of no
rational point of C.  The ratios: the phantoms of the scan, (1109503,597671), and their partners (h1_partners.txt); a
ratio and its partner have isomorphic systems (CENSUS9000.md), so one rank-0 curve on either side excludes both.
Output h1_rank0_exact.txt.   usage: h1_rank0_exact.py [--write]   |   h1_rank0_exact.py --negative [--write]
The repository is found through MSS_REPO, or else ../../.."""
import sys, os, re, ast, glob, random
from fractions import Fraction as Fr
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv; NEG = '--negative' in sys.argv
_argv = sys.argv; sys.argv = [sys.argv[0]]
import h1_pairs as HP
sys.argv = _argv
RS, RI = HP.RS, HP.RI
# when PARI's conic point has N = 0 (curve B: b = S^4 is a square), take the repository's own rank_second.find_conic_point_any,
# which moves to another point of the same conic; the point is a hint either way, checked by rank_indep.second_descent
RI.find_conic_point = RS.find_conic_point_any
OUT = []
def say(s=''): OUT.append(s)

def tuple_at(m, n, t):
    if t is None: f = {'0': (Fr(0), Fr(n)), '8': (Fr(0), Fr(m)), 'L': (Fr(m + n), Fr(m - n)), 'U': (Fr(m - n), Fr(m + n))}
    else: f = HP.forms(m, n, Fr(t))
    tp = []
    for j in '08LU':
        z = f[j]; Nz = z[0] ** 2 + z[1] ** 2
        if Nz == 0 or not HP.isq(Nz): return None
        tp.append(HP.gauss_class(z))
    return tuple(tp)

def ratios():
    out = {(1109503, 597671): ([(1109503, 597671, 2, 1), (1195342, 2219006, 2, 1)],
                              [(1109503, 597671, 2, 959), (1109503, 597671, 1918, 1), (1195342, 2219006, 2, 959), (1195342, 2219006, 1918, 1)])}
    for f in sorted(glob.glob(os.path.join(HERE, 'h1_scan_k*.txt'))):
        for l in open(f):
            mm = re.match(r'\((\d+),(\d+)\) PHANTOM: .*known (\[.*\]); extra (\[.*\])', l.strip())
            if mm: out[(int(mm.group(1)), int(mm.group(2)))] = (ast.literal_eval(mm.group(3)), ast.literal_eval(mm.group(4)))
    pf = os.path.join(HERE, 'h1_partners.txt')
    part = {}
    for l in open(pf):
        mm = re.match(r'\((\d+),(\d+)\) -> partner \((\d+),(\d+)\).*: PHANTOM; .*known (\[.*\]); extra (\[.*\])', l.strip())
        if mm:
            a, b = (int(mm.group(1)), int(mm.group(2))), (int(mm.group(3)), int(mm.group(4)))
            part[a] = b
            if b not in out: out[b] = (ast.literal_eval(mm.group(5)), ast.literal_eval(mm.group(6)))
    return out, part

# ---- exact group law and torsion on y^2 = x^3 + a x^2 + b x
def add(a, b, P, Q):
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2:
        if y1 + y2 == 0: return None
        l = (3 * x1 * x1 + 2 * a * x1 + b) / (2 * y1)
    else: l = (y2 - y1) / (x2 - x1)
    x3 = l * l - a - x1 - x2; return (x3, -(y1 + l * (x3 - x1)))
def mul(a, b, k, P):
    R = None
    for _ in range(k): R = add(a, b, R, P)
    return R
def npts(a, b, p):
    c = 1
    for x in range(p):
        v = (x * x * x + a * x * x + b * x) % p
        c += 1 if v == 0 else (2 if pow(v, (p - 1) // 2, p) == 1 else 0)
    return c
def torsion(a, b, drop=0):
    import cypari2
    pari = cypari2.Pari()
    E = pari.ellinit([0, a, 0, b, 0]); T = pari.elltors(E)
    pts = [None]
    for g, d in zip(list(T[2]), [int(x) for x in T[1]]):
        G = (Fr(str(g[0])), Fr(str(g[1])))
        pts = [add(a, b, P, mul(a, b, k, G)) for P in pts for k in range(d)]
    pts = list(dict.fromkeys(pts))
    if drop: pts = pts[:-drop]
    for P in pts:
        if P is not None: assert P[1] ** 2 == P[0] ** 3 + a * P[0] ** 2 + b * P[0], 'a torsion point is not on the curve'
        assert mul(a, b, 12, P) is None or mul(a, b, 8, P) is None or mul(a, b, 10, P) is None, 'a point of infinite order'
    disc = b * b * (a * a - 4 * b); g = 0
    for p in range(3, 400):
        if all(p % q for q in range(2, isqrt(p) + 1)) and disc % p: g = gcd(g, npts(a % p, b % p, p))
    assert len(pts) == g, 'the torsion list (%d points) is not all of E(Q)_tors (gcd of #E(F_p): %d)' % (len(pts), g)
    return pts

def bound_A(m, n):
    return bound_class(m, n)

def bound_full(es):
    bad, loc, sel, kum, enc, places = RS.full_selmer2(es, random.Random(1))
    Ds = [RS.DFull(es, s_) for s_ in sel]
    Mx = [[Ds[i].pairing(b_)[0] for b_ in sel] for i in range(len(sel))]
    for i in range(len(sel)):
        assert Mx[i][i] == 1, 'the pairing is not alternating'
        for j in range(len(sel)): assert Mx[i][j] == Mx[j][i], 'the pairing is not symmetric'
    rho = len(RS.f2_span([sum((0 if v == 1 else 1) << j for j, v in enumerate(r)) for r in Mx]))
    r = len(sel) - 2 - rho; info = 'dim Sel^2 = %d, Cassels-Tate pairing of rank %d' % (len(sel), rho)
    return r, info

def isog(a, b):
    """the 2-isogenous curve of y^2 = x (x^2 + a x + b) by (0, 0): y^2 = x (x^2 - 2 a x + a^2 - 4 b) (Velu; standard)"""
    return -2 * a, a * a - 4 * b

def roots_of(a, b):
    """the rational roots 0, e2, e3 of x (x^2 + a x + b) if x^2 + a x + b splits over Q, else None"""
    D = a * a - 4 * b
    if D < 0 or isqrt(D) ** 2 != D: return None
    r = isqrt(D); assert (-a + r) % 2 == 0
    return tuple(sorted((0, (-a + r) // 2, (-a - r) // 2)))

def shift(es, e):
    """(a, b) of the model with the root e moved to x = 0"""
    u, v = [x - e for x in es if x != e]
    return -(u + v), u * v

def bound_class(m, n, limit=6):
    """the least rank bound over the curves with full rational 2-torsion in A's isogeny class reached by 2-isogenies
    (at most LIMIT curves), each by the complete 2-descent and the pairing; the rank is the same on isogenous curves"""
    M, N = m ** 4, n ** 4
    todo = [roots_of(M + N, M * N)]; seen = []; best = None; info = []
    while todo and len(seen) < limit:
        es = todo.pop(0)
        g = gcd(es[1] - es[0], es[2] - es[0]); u = 1
        for q in range(2, 200):
            while (g // (u * u)) % (q * q) == 0: u *= q
        key = tuple((x - es[0]) // (u * u) for x in es)     # up to x -> x + c and x -> u^2 x
        if key in seen: continue
        seen.append(key)
        try: r, inf = bound_full(es)
        except AssertionError as ex:          # the pairing's local search gave up on this curve: no bound from it
            if 'no local point found' not in str(ex): raise
            info.append('pairing not computed (no local point found at %s)' % ex.args[0][1]); continue
        info.append('rank <= %d (%s)' % (r, inf))
        if best is None or r < best: best = r
        if best == 0: break
        for e in es:
            nx = roots_of(*isog(*shift(es, e)))
            if nx: todo.append(nx)
    if best is None: best = 99
    return best, 'over %d curves of the isogeny class with full 2-torsion: %s' % (len(seen), ', '.join(info))

def bound_B(m, n):
    S = m * m + n * n; a, b = 2 * S * S - 16 * m * m * n * n, S ** 4
    r, rec = RI.rank_bound(a, b, [], [], say=lambda s: None)
    return r, 'descent via 2-isogeny: |Sel(phi)| %d, |Sel(phi^)| %d; second descent leaves 2^%d and 2^%d' % (rec['sel'], rec['selp'], rec['E']['log2'], rec["E'"]['log2'])

def curve(m, n, which):
    S = m * m + n * n
    if which == 'A': return m ** 4 + n ** 4, m ** 4 * n ** 4, Fr(m * m * n * n), 'w f0 f8'
    return 2 * S * S - 16 * m * m * n * n, S ** 4, Fr(S * S), 'w K'

def analyse(m, n, known, extra, drop=0):
    HP.set_base(set(q for x in (2, m, n, m - n, m + n, m * m + n * n) for q in RI.factor_int(x)) | {2})
    for which in 'AB':
        a, b, sc, nm = curve(m, n, which)
        r, info = (bound_A if which == 'A' else bound_B)(m, n)
        if r > 0:
            say('  %s: y^2 = %s: rank <= %d (%s): no conclusion' % (which, nm, r, info)); continue
        pts = torsion(a, b, drop)
        ts = set()
        for P in pts:
            if P is None: ts.add(None); continue
            w = P[0] / sc
            if w >= 0 and HP.isq(w): s_ = HP.rsqrt(w); ts |= {s_, -s_}
        tps = sorted(set(x for x in (tuple_at(m, n, t) for t in ts) if x is not None))
        for k in known: assert k in tps, 'a known tuple is not the tuple of a torsion point: the method is wrong'
        v = [tp for tp in extra if tp not in tps]
        say('  %s: y^2 = %s: rank 0 (%s); %d torsion points (all of them); t over them: %s; tuples on C: %d (the known ones)'
            % (which, nm, info, len(pts), sorted(['inf' if t is None else str(t) for t in ts]), len(tps)))
        return which, len(v) == len(extra)
    return None, False

def main():
    say('H1: rank-0 quotients of C at the phantoms and their partners, rank bounds exact. Written by an AI model; not read by a person. Experiment, not counted.')
    R, part = ratios(); res = {}
    for (m, n) in sorted(R):
        known, extra = R[(m, n)]
        say(''); say('(%d,%d): %d known tuples, %d extra' % (m, n, len(known), len(extra)))
        try: res[(m, n)] = analyse(m, n, known, extra)
        except AssertionError as ex:
            say('  stopped: %s' % str(ex).split(':')[0]); res[(m, n)] = (None, False)
    say(''); say('VERDICT (a ratio and its partner together):')
    seen = set(); nk = 0; npair = 0
    for a in sorted(part):
        b = part[a]; key = tuple(sorted((a, b)))
        if key in seen: continue
        seen.add(key); npair += 1
        kill = [x for x in key if res.get(x, (None, False))[1]]
        if kill: nk += 1
        say('  %s and %s: %s' % (key[0], key[1], ('every extra tuple excluded (rank 0 of y^2 = %s at %s)' % ({'A': 'w f0 f8', 'B': 'w K'}[res[kill[0]][0]], kill[0])) if kill else 'not excluded'))
    say('phantom pairs: %d; excluded by an exact rank-0 quotient: %d' % (npair, nk))

def negative():
    say('H1: rank-0 quotients, forged inputs. Written by an AI model; not read by a person.')
    m, n = 1109503, 597671
    HP.set_base(set(q for x in (2, m, n, m - n, m + n, m * m + n * n) for q in RI.factor_int(x)) | {2})
    a, b, sc, nm = curve(m, n, 'A')
    try: torsion(a, b, drop=1); say('a torsion point withheld: NOT REFUSED')
    except AssertionError: say('a torsion point withheld: refused (the list is shorter than the gcd of #E(F_p))')
    forged = [tuple_at(m, n, Fr(0))]
    try:
        w, ok = analyse(m, n, [(1109503, 597671, 2, 1)], forged)
        say('an "extra" tuple that is the tuple of t = 0: %s' % ('reported as excluded: WRONG' if ok else 'not excluded, as it must be'))
    except AssertionError as ex: say('forged extra: stopped (%s)' % str(ex).split(':')[0])
    try:
        analyse(m, n, [(1109503, 597671, 2, 3)], [])
        say('a known tuple that is not over the torsion: NOT REFUSED')
    except AssertionError: say('a known tuple that is not over the torsion: refused')

if __name__ == '__main__':
    (negative if NEG else main)()
    s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'h1_rank0_exact%s.txt' % ('_negative' if NEG else ''))
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
