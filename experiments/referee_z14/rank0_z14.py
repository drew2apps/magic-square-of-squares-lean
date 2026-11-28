"""z14 task 4: EXPLORATION_3 §4 (rank0.py) recomputed with the repository's descent: Sel^2 by rank_second.full_selmer2
and the Cassels-Tate pairing by t10_pairct.DFullX (rank_second.DFull.pairing), for E_A: y^2 = x(x - A^2)(x - C^2) (= E_08)
and its twist by -1, y^2 = x(x - B^2)(x - C^2) (= E_LU up to 2^2).  bound = dim Sel^2 - 2 - rank(pairing): the rank
bound from Sel^2 and the pairing (E(Q)[2] = (Z/2)^2).  Also explore3's own rank0.py numbers (its code, copied as
e3_rank0.py with paths only changed, run on the same ratios) for comparison, line by line.
Ratios: exactly rank0.py's generator (random.Random(5); 'nd' with H = 3000, N = 150; 'sq' with H = 3000, N = 80), and an
independent draw ('nd', random.Random(6), H = 3000, N = 150).
Written by an AI model (reviewer z14); not read by a person.   usage: rank0_z14.py OUT.jsonl"""
import sys, os, json, time, random
HERE = os.path.dirname(os.path.abspath(__file__))
ARGS = sys.argv[1:]; sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
RS = T.RS
T.PARI.allocatemem(2 * 10 ** 8, 2 * 10 ** 9)          # z14: smaller PARI stack (memory is shared)
from math import gcd, isqrt

def is_sq(x): r = isqrt(x); return r * r == x
def ratios(seed, H, N, kind):
    """rank0.py's generator, verbatim in effect"""
    rng = random.Random(seed); pairs = []
    if kind == 'nd':
        while len(pairs) < N:
            m = rng.randint(2, H); n = rng.randint(1, m - 1)
            if gcd(m, n) != 1 or not ((m * n) % 4 == 0 or (m - n) % 8 == 0 or (m + n) % 8 == 0): continue
            S = m * m + n * n
            if is_sq(S) or (S % 2 == 0 and is_sq(S // 2)): continue
            pairs.append((m, n))
    else:
        for a in range(2, 80):
            for b in range(1, a):
                if gcd(a, b) == 1 and (a - b) % 2 == 1:
                    m, n = max(a * a - b * b, 2 * a * b), min(a * a - b * b, 2 * a * b)
                    if (m * n) % 4 == 0 and m <= H: pairs.append((m, n))
        rng.shuffle(pairs); pairs = pairs[:N]
    return pairs

def rank2(rows):
    B = []
    for v in rows:
        for b in B: v = min(v, v ^ b)
        if v: B.append(v)
    return len(B)

def bound(es):
    bad, loc, sel, kum, enc, places = RS.full_selmer2(es, random.Random(1))
    d = len(sel)
    rows = [0] * d; fails = 0; asym = 0
    D = [T.DFullX(es, a) for a in sel]
    for i in range(d):
        for j in range(i + 1, d):
            h, per = T.safe_pair(D[i], sel[j])
            if h is None: fails += 1; continue
            if h == -1: rows[i] |= 1 << j; rows[j] |= 1 << i
    diag = [T.safe_pair(D[i], sel[i])[0] for i in range(d)]
    return dict(dim=d, rank=rank2(rows), bound=d - rank2(rows) - 2, fails=fails, diag_ok=all(x == 1 for x in diag))

if __name__ == '__main__':
    out = ARGS[0]
    done = set()
    if os.path.exists(out):
        for l in open(out): r = json.loads(l); done.add((r['set'], r['m'], r['n']))
    sets = [('nd5', ratios(5, 3000, 150, 'nd')), ('sq5', ratios(5, 3000, 80, 'sq')), ('nd6', ratios(6, 3000, 150, 'nd'))]
    for name, prs in sets:
        for m, n in prs:
            if (name, m, n) in done: continue
            t0 = time.time()
            A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
            rA = bound((0, A * A, C * C)); rB = bound((0, B * B, C * C))
            rec = dict(set=name, m=m, n=n, EA=rA, Etw=rB, time=round(time.time() - t0, 1))
            with open(out, 'a') as f: f.write(json.dumps(rec, sort_keys=True) + '\n')
            sys.stderr.write('%s %d,%d %s %s %.1f\n' % (name, m, n, rA, rB, time.time() - t0))
