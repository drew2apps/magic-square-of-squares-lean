"""Reviewer z16 (an AI model; not read by a person): the four-class descent of tuple_descent.magic with the centre
F4 = 2(m^2+n^2)(u^2+v^2) added as a fifth Gaussian linear form, recomputed from the repository's own local code.

gamma4 = (1+i)(m+in)(u+iv) = c*u + (c*i)*v, c = (1+i)(m+in) = (m-n) + i(m+n); N(gamma4) = 2(m^2+n^2)(u^2+v^2) = F4.
Support of the class d4 (Lemma S2 with the determinant -4i(m^2+n^2) of gamma4 and its conjugate): the primes of
2(m^2+n^2).  With --exact-support only 2 and the primes of odd exponent in m^2+n^2 are used (a check: the local tests
at the primes of m^2+n^2 must remove the others, so the two give the same tuples).

The local sets: exact_local.allowed_at (exact at odd p, the repository's TD.Local at 2 with depth 12), at the same bad
primes as G2 (every prime of 2mn(m^2-n^2), m^2+n^2 and m^8+14m^4n^4+n^8).  Factors by PARI, each certified prime by
PARI's isprime (APRCL).

usage: python3 -I z16_five.py ratios FILE [--exact-support] [--four]     one line per ratio "m,n ..." (first word)
       python3 -I z16_five.py one m,n
Output: one JSON line per ratio: m, n, n4 (four-form tuples), n5 (five-form tuples), ext (4-tuples that extend)."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
os.environ['MSS_REPO'] = REPO
sys.path.insert(0, os.path.join(REPO, 'experiments', 'roadmap-g', 'g1'))
import exact_local as X
TD = X.TD
import cypari2
P = cypari2.Pari(); P.allocatemem(6 * 10**8, silent=True)

_FC = {}
def fact(x):
    x = abs(int(x))
    if x <= 1: return {}
    if x not in _FC:
        f = P.factor(x); d = {int(q): int(e) for q, e in zip(f[0], f[1])}
        for q in d: assert P.isprime(q), ('not a certified prime', q)
        _FC[x] = d
    return _FC[x]
TD.factorint = fact                         # TD.magic and TD.sqf factor m, n, m - n, m + n with this

def gamma4(m, n):
    c = TD.gmul((1, 1), (m, n))
    return (c, TD.gmul(c, (0, 1)))

def divisors(S):
    ds = [1]
    for q in S: ds += [x*q for x in ds]
    return ds

def bad_primes(m, n):
    bad = set()
    for x in (2*m*n*(m*m - n*n), m*m + n*n, m**8 + 14*m**4*n**4 + n**8): bad |= set(fact(x))
    return sorted(bad)

def residual(m, n, five=True, exact_support=False, primes=None):
    forms, sup, known = TD.magic(m, n)
    C = m*m + n*n
    if five:
        g = gamma4(m, n); forms = forms + [g]
        # identity check on a few points: the norm of the fifth form is F4
        for u, v in ((1, 0), (0, 1), (3, -7), (11, 5)):
            z = (g[0][0]*u + g[1][0]*v, g[0][1]*u + g[1][1]*v)
            assert z[0]**2 + z[1]**2 == 2*C*(u*u + v*v)
        fC = fact(C)
        sup = sup + [sorted({2} | {q for q, e in fC.items() if (e % 2 or not exact_support)})]
    divs = [divisors(S) for S in sup]
    loc = []
    for p in (primes if primes is not None else bad_primes(m, n)):
        L, al = X.allowed_at(forms, p)
        if len(al) < len(L.labels())**len(forms): loc.append((p, L, al))
    return X.join(divs, loc)

def main():
    if sys.argv[1] == 'one':
        m, n = (int(x) for x in sys.argv[2].split(','))
        T5 = residual(m, n); print(len(T5), T5[:20]); return
    pairs = [tuple(int(x) for x in l.split()[0].split(',')) for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
    exs = '--exact-support' in sys.argv
    for m, n in pairs:
        t0 = time.time()
        T4 = residual(m, n, five=False)
        T5 = residual(m, n, five=True, exact_support=exs)
        ext = sorted(set(t[:4] for t in T5))
        assert set(ext) <= set(T4), 'a 5-tuple whose first four classes are not a 4-tuple'
        print(json.dumps(dict(m=m, n=n, n4=len(T4), n5=len(T5), ext=[list(t) for t in ext], t5=[list(t) for t in T5])), flush=True)
        print('%d,%d %.2f' % (m, n, time.time() - t0), file=sys.stderr, flush=True)

if __name__ == '__main__':
    main()
