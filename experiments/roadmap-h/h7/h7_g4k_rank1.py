"""H7: rank C(Q) = 1 for the quotient C: y^2 = g4 K, by a complete 2-descent on a 2-isogenous curve (EXPERIMENT, written
by an AI model; not read by a person; nothing here changes the census).
F4's g4k_exact.py shows, in exact arithmetic: IF rank C(Q) = 1, every solution t of the ratio is degenerate.  F4's
exact bound (g4k_rank_exact.py, descent via 2-isogeny on C) is 2 to 5, so its closures stay conditional.  Here:
  (1) C is dm_g4k.g4k_curve(m, n), the same curve as F4's; a rational root e of its cubic (PARI's nfroots as a hint,
      checked exactly) moves it to y^2 = X (X^2 + a X + b), X = x - e (an isomorphism over Q);
  (2) E' : y^2 = X (X^2 - 2 a X + a^2 - 4 b) is 2-isogenous to it (Velu's formula, quoted; the script also checks
      that the two curves have the same number of points mod every good prime < 1000, as isogenous curves must), so
      rank E'(Q) = rank C(Q);
  (3) if E' has full rational 2-torsion (a^2 - 4 (a^2 - 4 b) = 16 b a square), the repository's complete 2-descent
      and Cassels-Tate pairing (rank_second.full_selmer2 and DFull.pairing, through h1/h1_rank0_exact.bound_full,
      unchanged) give rank E'(Q) <= dim Sel^2 - 2 - rank of the pairing;
  (4) a point G of C (PARI's ellrank as a hint) is checked on the curve and of infinite order exactly (n G != O for
      n = 1 .. 12, Mazur, quoted), so rank >= 1.
If the bound of (3) is 1, rank C(Q) = 1 and F4's condition holds.  Output h7_g4k_rank1.txt.
usage: h7_g4k_rank1.py [--negative] [--write] [m,n ...]   (h7_g4k_rank1_negative.txt)
Finds the repository through MSS_REPO, or else ../../..; H1's folder through H1_DIR, or else ../h1."""
import sys, os, io, contextlib
from fractions import Fraction as Fr
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
WRITE = '--write' in sys.argv; NEG = '--negative' in sys.argv
ARGS = [a for a in sys.argv[1:] if not a.startswith('--')]
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..', '..'))
H1 = os.environ.get('H1_DIR', os.path.join(HERE, '..', 'h1'))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'census-wide', 'g4k')); sys.path.insert(0, H1)
_argv = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import dm_g4k
    import h1_rank0_exact as R
sys.argv = _argv
import cypari2
pari = cypari2.Pari()
OUT = []
def say(s=''): OUT.append(s)
RATIOS = ['12201,7360', '7735,4968', '6557,324', '6881,6233']

def is_prime(p): return p > 1 and all(p % q for q in range(2, isqrt(p) + 1))
def npts(a, b, p):
    """#E(F_p) for y^2 = x (x^2 + a x + b), by Legendre symbols (exact)"""
    s = 1
    for x in range(p):
        v = (x * (x * x + a * x + b)) % p
        s += 1 if v == 0 else (2 if pow(v, (p - 1) // 2, p) == 1 else 0)
    return s

def shifted(m, n):
    ainv = dm_g4k.g4k_curve(m, n)[0]; _, a2, _, a4, a6 = [int(v) for v in ainv]
    rts = sorted(int(r) for r in pari.nfroots(None, pari('x^3 + (%d)*x^2 + (%d)*x + (%d)' % (a2, a4, a6))))
    e = rts[0]; assert e ** 3 + a2 * e * e + a4 * e + a6 == 0, 'the hint root is not a root'
    return ainv, e, 3 * e + a2, 3 * e * e + 2 * a2 * e + a4

def isogenous_check(a, b, a2, b2):
    disc = b * b * (a * a - 4 * b) * b2 * b2 * (a2 * a2 - 4 * b2)
    for p in range(3, 1000):
        if is_prime(p) and disc % p:
            if npts(a % p, b % p, p) != npts(a2 % p, b2 % p, p): return False
    return True

def infinite_order(a, b, G):
    assert G is not None and G[1] ** 2 == G[0] ** 3 + a * G[0] ** 2 + b * G[0], 'the hint point is not on the curve'
    for k in range(1, 13):
        if R.mul(a, b, k, G) is None: return False
    return True

def run(arg, forge=None):
    m, n = (int(v) for v in arg.split(','))
    ainv, e, a, b = shifted(m, n)
    E = pari.ellinit(ainv); pts = pari.ellrank(E, 2)[3]
    if len(pts) == 0: pts = pari.ellrank(E, 4)[3]
    G = [Fr(int(pari.numerator(pts[0][0])), int(pari.denominator(pts[0][0]))) - e,
         Fr(int(pari.numerator(pts[0][1])), int(pari.denominator(pts[0][1])))] if len(pts) else None
    a2, b2 = R.isog(a, b)
    if forge == 'curve': b2 += 4 * 10 ** 6                  # a forged "isogenous" curve
    if forge == 'point': G = [G[0] + 1, G[1]] if G else [Fr(1), Fr(1)]               # a forged hint point
    say('(%d,%d) g4*K: X = x - (%d); y^2 = X (X^2 + a X + b), a = %d, b = %d' % (m, n, e, a, b))
    if not isogenous_check(a, b, a2, b2):
        say('  E\' is not isogenous to C (point counts differ): refused'); return None
    try: inf = G is not None and infinite_order(a, b, G)
    except AssertionError as ex:
        say('  %s: refused' % ex.args[0]); return None
    es = R.roots_of(a2, b2)
    if es is None:
        say('  E\': y^2 = X (X^2 + %d X + %d) has no full rational 2-torsion; no bound here' % (a2, b2)); return 'none'
    r, info = R.bound_full(es)
    say('  E\' (2-isogenous, same point counts mod every good p < 1000): full rational 2-torsion, roots %s' % (es,))
    say('  complete 2-descent on E\': %s, so rank <= %d' % (info, r))
    say('  a point of infinite order is known (n G != O for n <= 12), so rank >= 1' if inf else '  no point of infinite order known')
    if r == 1 and inf: say('  so rank C(Q) = 1 (computed, exact), the condition of F4\'s g4k_exact.py')
    return r

def main():
    say('H7: rank C(Q) for C: y^2 = g4 K by a complete 2-descent on a 2-isogenous curve. Written by an AI model; not read by a person. Experiment, not counted.')
    res = {a: run(a) for a in (ARGS or RATIOS)}
    say('rank C(Q) = 1, computed (exact), for: %s' % ', '.join('(%s)' % a for a, r in res.items() if r == 1))

def negative():
    say('H7 g4*K rank 1, negative and positive controls. Written by an AI model; not read by a person.')
    say('forged curve: %s' % ('refused' if run('12201,7360', forge='curve') is None else 'ACCEPTED')); assert OUT[-1].endswith('refused')
    say('forged point: %s' % ('refused' if run('12201,7360', forge='point') is None else 'ACCEPTED')); assert OUT[-1].endswith('refused')
    for arg in ('11165,7548', '12204,9853'):                 # PARI finds 2 independent points on g4 K here
        m, n = (int(v) for v in arg.split(','))
        lo = int(pari.ellrank(pari.ellinit(dm_g4k.g4k_curve(m, n)[0]), 2)[0])
        r = run(arg)
        ok = r == 'none' or (r is not None and r >= lo)
        say('positive control %s: PARI finds %d independent points; the exact bound is %s: %s' % (arg, lo, r, 'consistent' if ok else 'NOT consistent'))
        assert ok

if __name__ == '__main__':
    negative() if NEG else main()
    s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'h7_g4k_rank1_negative.txt' if NEG else 'h7_g4k_rank1.txt')
    if WRITE: open(f, 'w').write(s)
    else:
        sys.stdout.write(s)
        if not ARGS and os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
