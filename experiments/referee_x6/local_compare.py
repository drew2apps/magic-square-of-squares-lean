"""Review x6 (run 6), written by an AI model; not read by a person.  exact_local.allowed_at (roadmap-g/g1) against
(1) the repository's tuple_descent.Local at increasing depth, where it can run (the exact set must be inside it), and
(2) the labels of explicit points computed by mylocal.py (they must be inside the exact set; the line says how many of
the exact set they reach).  107 (ratio, prime) cases: every odd bad prime of six ratios (octic primes to 1.9*10^15),
and chosen good and bad primes of four more, up to 2^61 - 1.  Prints only.   usage: python3 -I local_compare.py"""
import sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mylocal as ML
X, TD, P = ML.X, ML.TD, ML.P

def bad_primes(m, n):
    return sorted(set(int(q) for x in (2*m*n*(m*m-n*n), m*m+n*n, m**8+14*m**4*n**4+n**8) for q in P.factor(x)[0]))

def cases():
    c = [((627, 364), 3, 7), ((627, 364), 5, 7), ((627, 364), 7, 6), ((627, 364), 13, 5), ((627, 364), 17, 5),
         ((627, 364), 37, 4), ((627, 364), 263, 3), ((627, 364), 10009, 2), ((627, 364), 10007, 2),
         ((291, 224), 5, 7), ((291, 224), 7, 6), ((291, 224), 67, 4), ((291, 224), 1201, 3), ((291, 224), 29, 5),
         ((1109503, 597671), 5, 7), ((1109503, 597671), 7, 6), ((1109503, 597671), 137, 3), ((1109503, 597671), 593, 3),
         ((1109503, 597671), 7129, 2), ((1109503, 597671), 1871, 2), ((1109503, 597671), 597671, 0),
         ((1109503, 597671), 100049, 2), ((1109503, 597671), 1000033, 0), ((1109503, 597671), 1000000007, 0),
         ((993260025, 766111057), 5, 6), ((993260025, 766111057), 13, 5), ((993260025, 766111057), 3, 7),
         ((993260025, 766111057), 100003, 2), ((993260025, 766111057), 999999937, 0), ((993260025, 766111057), 2**61 - 1, 0)]
    m, n = 993260025, 766111057
    c += [((m, n), q, 0) for q in [q for q in bad_primes(m, n) if q > 10**5][:3]]
    for (m, n) in [(10125, 4459), (626, 1), (7, 4), (627, 364), (1189, 820), (4999, 1802)]:
        for q in bad_primes(m, n):
            if q > 2: c.append(((m, n), q, 6 if q < 10 else (4 if q < 60 else (3 if q < 2000 else (2 if q < 200000 else 0)))))
    return c

lines = []; rnd = random.Random(30)
for (m, n), p, Kmax in cases():
    forms, sup, known = TD.magic(m, n)
    L, E = X.allowed_at(forms, p)
    info = 'repository test not run (p too large)'
    if Kmax:
        res = []
        for K in range(2, Kmax + 1):
            R = TD.Local(forms, p, K).allowed(); res.append((K, len(R), R == E))
            assert E <= R, ('EXACT SET NOT INSIDE THE REPOSITORY SET', m, n, p, K)
            if R == E: break
        info = 'repository: ' + ', '.join('depth %d: %d%s' % (K, k, ' (equal)' if eq else '') for K, k, eq in res)
    S = set(ML.points(forms, p, rnd).keys())
    assert S <= E, ('AN EXPLICIT POINT OUTSIDE THE EXACT SET', m, n, p)
    lines.append('(%d,%d) p=%d (%d mod 4, %s): exact %d tuples; %s; explicit points reach %d of them'
                 % (m, n, p, p % 4, 'bad' if p in bad_primes(m, n) else 'good', len(E), info, len(S)))
    print(lines[-1], flush=True)
print('# %d cases; exact set inside the repository set wherever it ran; every explicit point inside the exact set; exact set equal to the repository set: %d of %d run; reached entirely by explicit points: %d of %d'
      % (len(lines), sum('(equal)' in l for l in lines), sum('repository:' in l for l in lines),
         sum(l.split('exact ')[1].split(' ')[0] == l.split('reach ')[1].split(' ')[0] for l in lines), len(lines)))
