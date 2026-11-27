"""Review x6 (run 6), written by an AI model; not read by a person.  G1's 63 kept ratios (roadmap-g/g2/g1_kept.txt) through
the four-class descent (a) with exact_local at odd primes, (b) with the repository's tuple_descent.Local only, at the
same primes (supports and good primes <= 150) up to CAP = 2*10^5, larger support primes left out.  Prints only.
usage: python3 -I kept63_repo_only.py   (about 10 minutes)"""
import sys, time
import os; HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "roadmap-g", "g1"))
import exact_local as X
TD = X.TD
CAP = 200000
for m, n in [tuple(map(int, l.split(","))) for l in open(os.path.join(HERE, "..", "roadmap-g", "g2", "g1_kept.txt")) if l.strip()]:
    t0 = time.time()
    forms, sup, known = TD.magic(m, n)            # sympy's factorint (the repository's default)
    primes = sorted(set(q for S in sup for q in S) | set(TD.primerange(2, 151)))
    S, used = X.descend_exact(forms, sup, known, primes)
    t1 = time.time()
    # (b) repository's Local only, at primes <= CAP, with the repository's depths
    divs = []
    for Pp in sup:
        ds = [1]
        for q in Pp: ds += [x*q for x in ds]
        divs.append(ds)
    loc = []; skipped = []
    for p in primes:
        if p > CAP: skipped.append(p); continue
        L = TD.Local(forms, p, 12 if p == 2 else (6 if p < 50 else 3)); al = L.allowed()
        if len(al) < len(L.labels())**4: loc.append((p, L, al))
    Sr = X.join(divs, loc)
    print('(%d,%d): tuple space %d; supports %s; exact descent (supports + good primes <= 150): %s (%.0f s); repository Local only, primes <= %d (left out: %s): %s (%.0f s)'
          % (m, n, len(divs[0])*len(divs[1])*len(divs[2])*len(divs[3]), [max(Pp) for Pp in sup], 'nothing left' if not S else S, t1 - t0, CAP, skipped or 'none', 'nothing left' if not Sr else '%d left' % len(Sr), time.time() - t1), flush=True)
