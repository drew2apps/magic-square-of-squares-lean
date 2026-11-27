"""Review x6 (run 6), written by an AI model; not read by a person.  (a) 20 kill sets of roadmap-g/g2 (random, seed 66;
10 of bad primes only, 10 with good primes) recomputed: empty with the repository's tuple_descent.Local at every prime
<= 3000 (exact_local above), and irredundant: after dropping each prime, one tuple left has explicit local points
(mylocal.py) at every remaining prime; (b) 10 ratios that the bad primes alone do not empty: every tuple they leave
gets explicit points at 2 and at every odd bad prime.  Prints only.   usage: python3 -I killsets_witness.py (about 1 minute)"""
import sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mylocal as ML
from mylocal import X, TD, P, witness, divisors
G2 = os.path.join(HERE, '..', 'roadmap-g', 'g2')
CAP = 3000
rnd_unused = None

def parse(f):
    rows = []
    for l in open(f):
        a = l.split()
        if len(a) < 4 or 'SKIPPED' in l: continue
        m, n = map(int, a[0].split(',')); ks = l.split('kill:')[1].split()
        rows.append((m, n, a[3] == 'EMPTY', None if ks[:2] == ['NOT', 'EMPTIED'] else [int(x.split(':')[0]) for x in ks]))
    return rows

def main():
    rnd = random.Random(66); out = []
    rows = parse(os.path.join(G2, 'g2_997_all.txt')) + parse(os.path.join(G2, 'g2_15000_all.txt'))
    multi = [r for r in rows if r[3] and len(r[3]) >= 2]
    pick = rnd.sample([r for r in multi if r[2]], 10) + rnd.sample([r for r in multi if not r[2]], 10)
    for m, n, empty, kill in pick:
        t0 = time.time()
        forms, sup, known = TD.magic(m, n); divs = divisors(sup)
        locs = []; how = []
        for p in kill:
            if p <= CAP:
                L = TD.Local(forms, p, 12 if p == 2 else (6 if p < 50 else 3)); al = L.allowed(); how.append('%d:repo' % p)
            else:
                L, al = X.allowed_at(forms, p); how.append('%d:exact' % p)
            locs.append((p, L, al))
        S = X.join(divs, locs)
        drops = []
        for i, p in enumerate(kill):
            rest = [q for q in kill if q != p]
            ex = [X.allowed_at(forms, q) for q in rest]
            T = X.join(divs, [(q, L, al) for q, (L, al) in zip(rest, ex)])
            ok = None
            for tup in T[:40]:
                ws = [ML.witness(forms, (m, n), q, tup, rnd) for q in rest]
                if all(w is not None for w in ws): ok = tup; break
            drops.append((p, len(T), ok is not None))
        line = '(%d,%d) kill %s (%s): with these primes %s; dropping each prime: %s' % (m, n, kill, 'bad only' if empty else 'with good primes', 'EMPTY' if not S else 'NOT EMPTY %d' % len(S),
              ', '.join('%d -> %d left%s' % (p, k, ', one with explicit points at every other prime' if w else ' (NO WITNESS FOUND)') for p, k, w in drops))
        out.append(line); print(line, '(%.0f s)' % (time.time() - t0), flush=True)
        assert not S
    # (b) witnesses
    fails = [r for r in parse(os.path.join(G2, 'g2_997_all.txt')) + parse(os.path.join(G2, 'g2_15000_all.txt')) if not r[2]]
    for m, n, empty, kill in rnd.sample(fails, 10):
        t0 = time.time()
        forms, sup, known = TD.magic(m, n); divs = divisors(sup)
        bad = sorted(set(int(q) for x in (2*m*n*(m*m-n*n), m*m+n*n, m**8+14*m**4*n**4+n**8) for q in P.factor(x)[0]))
        locs = [(p,) + X.allowed_at(forms, p) for p in bad]
        T = X.join(divs, [l for l in locs if len(l[2]) < len(l[1].labels())**4])
        res = []
        for tup in T:
            miss = [q for q in bad if ML.witness(forms, (m, n), q, tup, rnd) is None]
            res.append(miss)
        nw = sum(not x for x in res)
        line = '(%d,%d) bad primes %s: %d tuples left; with explicit points at 2 and at every odd bad prime: %d%s' % (m, n, bad if len(bad) < 12 else '%d primes up to %d' % (len(bad), bad[-1]), len(T), nw,
              '' if nw == len(T) else '; missing at %s' % sorted(set(q for x in res for q in x)))
        out.append(line); print(line, '(%.0f s)' % (time.time() - t0), flush=True)
    return out

main()
