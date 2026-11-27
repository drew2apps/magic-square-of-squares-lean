"""Review z8 of roadmap H (written by an AI model; not read by a person; experiment, not counted).
H1's scan (through G1's g1_stress.one / descent) tests the primes of the supports (2mn, 2(m^2 - n^2)) and the primes up
to 20,000 only.  This lists, for the 24 phantom ratios of H1, the bad primes that run left out (primes of m^2 + n^2 and of
the octic m^8 + 14 m^4 n^4 + n^8 = ((m^2-mn-n^2)^2 + 3m^2n^2)((m^2+mn-n^2)^2 + 3m^2n^2), not in the supports, above
20,000; PARI's factor, each factor proved prime by isprime), and runs G1's exact_local.allowed_at there: the known tuples
must be allowed (asserted); is each extra tuple allowed?  Rests on G1's exact_local.py and the repository's
tuple_descent.py.
usage: badprimes_left.py [--write]   (badprimes_left.txt; about 6 s)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from z8common import emit, G1, load_phantoms, pari, fac_proved
sys.path.insert(0, G1)
import exact_local as X
P = pari(); out = ['z8: the bad primes above 20,000 that H1\'s scan left out, tested with G1\'s exact_local. Written by an AI model; not read by a person.']
cut_pairs = []
for (m, n), (known, extra) in sorted(load_phantoms().items()):
    forms, sup, kn = X.TD.magic(m, n)
    supp = set(q for S in sup for q in S)
    a = (m*m - m*n - n*n)**2 + 3*m*m*n*n; b = (m*m + m*n - n*n)**2 + 3*m*m*n*n
    assert a * b == m**8 + 14*m**4*n**4 + n**8
    S2 = set(fac_proved(P, m*m + n*n)); Oc = set(fac_proved(P, a)) | set(fac_proved(P, b))
    left = sorted(q for q in (S2 | Oc) if q not in supp and q > 20000)
    res = []
    for p in left:
        L, al = X.allowed_at(forms, p)
        for T in known: assert tuple(L.dlabel(d) for d in T) in al, ('a known tuple rejected', m, n, p)
        cut = [T for T in extra if tuple(L.dlabel(d) for d in T) not in al]
        res.append('%d (%s, p mod 4 = %d): %d allowed label tuples, %s' % (p, 'm^2+n^2' if p in S2 else 'octic', p % 4, len(al), ('cuts %d of %d extra tuples' % (len(cut), len(extra))) if cut else 'keeps every extra tuple'))
        if len(cut) == len(extra) and (m, n) not in cut_pairs: cut_pairs.append((m, n))
    out.append('(%d,%d): %d such primes' % (m, n, len(left)))
    out += ['  ' + x for x in res]
out.append('ratios whose every extra tuple is cut at one of these primes: %d %s' % (len(cut_pairs), cut_pairs))
emit(out, 'badprimes_left.txt')
