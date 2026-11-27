"""G2, looking for the law: candidate residue-symbol features of a ratio, tabulated against whether the bad primes alone
empty its four-class set (from g2_*_all.txt).  A feature (or its negation) that holds for 5 ratios or more and forces one
outcome is listed at the end.  None does.  (EXPERIMENT, written by an AI model; not read
by a person.)
Features, over the primes q = 1 mod 8 dividing m^2+n^2 to an odd power (S8 below):
  has_q8: S8 is non-empty;  q16: some q in S8 is 1 mod 16;
  quartic_<a>_trivial_some / _nontriv_some: for some q in S8, a is (is not) a quartic residue mod q, a in 2, m, n, m-n,
  m+n, mn;  all_X_primes_quad_res / quartic_res: for some q in S8, every prime of mn(m^2-n^2) is a square (a fourth
  power) mod q;  all_pm1mod8_X_quartic_res: the same for the primes of mn(m^2-n^2) that are +-1 mod 8.
usage: g2_features.py g2_<list>_all.txt [--write]     (writes g2_features_<list>.txt)"""
import sys, os, collections
import cypari2
P = cypari2.Pari()
HERE = os.path.dirname(os.path.abspath(__file__))
def fac(x):
    f = P.factor(x); out = {int(q): int(e) for q, e in zip(f[0], f[1])}
    for q in out: assert P.isprime(q)
    return out
def power_res(a, q, k):            # a^((q-1)/k) mod q, 1 iff a is a k-th power residue mod q (a prime to q)
    a %= q
    return None if a == 0 else pow(a, (q - 1)//k, q)
def table(path):
    rows = []
    for l in open(path):
        if not l.strip(): continue
        m, n = map(int, l.split()[0].split(',')); emp = 'bad-alone: EMPTY' in l
        S = fac(m*m + n*n); X = fac(m*n*(m*m - n*n)); F = collections.defaultdict(bool)
        qs = [q for q, e in S.items() if e % 2 and q % 8 == 1]
        F['has_q8'] = bool(qs)
        for q in qs:
            F['q16'] |= q % 16 == 1
            for name, a in (('2', 2), ('m', m), ('n', n), ('m-n', m - n), ('m+n', m + n), ('mn', m*n)):
                r = power_res(a, q, 4)
                if r is None: continue
                F['quartic_%s_trivial_some' % name] |= r == 1
                F['quartic_%s_nontriv_some' % name] |= r != 1
            F['all_X_primes_quartic_res'] |= all(power_res(p, q, 4) == 1 for p in X if p != q)
            F['all_X_primes_quad_res'] |= all(power_res(p, q, 2) == 1 for p in X if p != q)
            F['all_pm1mod8_X_quartic_res'] |= all(power_res(p, q, 4) == 1 for p in X if p % 8 in (1, 7))
        rows.append((emp, F))
    keys = sorted(set(k for r in rows for k in r[1]))
    out = ['# %s: %d ratios; empty = bad primes alone empty the set, left = they do not; F = feature holds' % (os.path.basename(path), len(rows))]
    for k in keys:
        c = collections.Counter((r[0], bool(r[1].get(k))) for r in rows)
        out.append('%-30s empty&F %5d  empty&~F %5d  left&F %4d  left&~F %4d' % (k, c[(True, True)], c[(True, False)], c[(False, True)], c[(False, False)]))
    imp = []
    for k in keys:
        c = collections.Counter((r[0], bool(r[1].get(k))) for r in rows)
        for val in (True, False):           # does F (or not F), where it occurs at least 5 times, force one outcome?
            tot = c[(True, val)] + c[(False, val)]
            if tot >= 5 and 0 in (c[(True, val)], c[(False, val)]):
                imp.append('%s%s => %s (%d ratios)' % ('' if val else 'not ', k, 'empty' if c[(False, val)] == 0 else 'left', tot))
    out.append('# features that, where they hold for 5 ratios or more, force one outcome: %s' % ('; '.join(imp) if imp else 'none'))
    return out
if __name__ == '__main__':
    path = os.path.join(HERE, [a for a in sys.argv[1:] if not a.startswith('--')][0])
    out = '\n'.join(table(path)) + '\n'
    f = os.path.join(HERE, 'g2_features_' + os.path.basename(path).replace('g2_', '').replace('_all', ''))
    if '--write' in sys.argv: open(f, 'w').write(out)
    print(out, end='')
    if '--write' not in sys.argv and os.path.exists(f): assert open(f).read() == out; print('# matches ' + os.path.basename(f))
