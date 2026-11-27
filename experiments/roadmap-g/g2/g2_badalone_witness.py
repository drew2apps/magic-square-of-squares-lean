"""G2, test 6 made exact: for the ratios that the bad primes alone do not empty (g2_killsets.py), every tuple left is
shown to be locally possible at 2 by an explicit 2-adic point: integers (u, v) = (t, 1) or (1, 2t), 0 <= t < 2^14, at
which the classes key2 of the four forms are the tuple's (or 3 times them: (u, v) is defined up to a 2-adic unit).  At
the odd bad primes the exact test of exact_local.py is used (no prime left out).  So for these ratios the bad primes
alone leave a tuple that no bad prime excludes, and the answer to test 6 of COMPLETE.md section 7 is no.
(EXPERIMENT, written by an AI model; not read by a person.)   usage: g2_badalone_witness.py RATIOS.txt [--write]"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import g2_killsets as G
X, TD = G.X, G.TD

def witnessed_at_2(forms, tuples, T=2**14):
    L = TD.Local(forms, 2, 12); seen = set()
    for c in (0, 1):
        for t in range(T):
            u, v = (t, 1) if c == 0 else (1, 2*t)
            gs = [(al[0]*u + be[0]*v, al[1]*u + be[1]*v) for al, be in forms]
            if any(g == (0, 0) for g in gs): continue
            ks = tuple(TD.key2(g)[1] for g in gs)
            if all(k in TD.D2 for k in ks):
                seen.add(ks); seen.add(tuple(L.times3(k) for k in ks))
    return [tuple(L.dlabel(d) for d in tp) in seen for tp in tuples]

if __name__ == '__main__':
    out = []
    for l in open(sys.argv[1]):
        if 'bad-alone: EMPTY' in l or not l.strip(): continue
        m, n = map(int, l.split()[0].split(','))
        forms, sup, known = TD.magic(m, n)
        bad = set()
        for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): bad |= set(G.factorint_proved(x))
        locs = []
        for p in sorted(bad):
            L, al = X.allowed_at(forms, p)
            if len(al) < len(L.labels())**len(forms): locs.append((p, L, al))
        left = X.join(G.divisors(sup), locs)
        w = witnessed_at_2(forms, left)
        out.append('%d,%d  left by the bad primes: %d  witnessed at 2: %d  %s' % (m, n, len(left), sum(w), 'OK' if all(w) else 'NOT ALL'))
        print(out[-1], flush=True)
    out.append('# %d ratios; all witnessed: %s' % (len(out), all(o.endswith('OK') for o in out)))
    print(out[-1])
    f = os.path.join(HERE, 'g2_badalone_witness_' + os.path.basename(sys.argv[1]))
    if '--write' in sys.argv: open(f, 'w').write('\n'.join(out) + '\n')
    elif os.path.exists(f): assert open(f).read() == '\n'.join(out) + '\n'; print('# matches ' + os.path.basename(f))
