"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Exact identities behind H6 (experiments/roadmap-h/h6), in sympy, for generic m, n:
  * the five forms of G1_FULL_SYSTEM.md (F0 = m^2 + n^2 t^2, F8 = n^2 + m^2 t^2, F4 = s(1 + t^2) with s = (m^2 + n^2)/2,
    L = F4 - 2mn t, U = F4 + 2mn t, in t = a/b) are, times 1, 1, 4, 4, 4, those of h6_split.forms;
  * sigma: t -> -t fixes F0, F8, F4 and swaps L, U; tau: t -> 1/t (times t^2) swaps F0, F8 and fixes F4, L, U; so
    F0F4L(-t) = F0F4U(t), t^6 F0F4L(1/t) = F8F4L(t), t^6 F0F4L(-1/t) = F8F4U(t): the four genus-2 pieces are one curve;
  * the partner map: with (m', n') = (m + n, m - n) and t' = (t + 1)/(t - 1), (t - 1)^2 F'(t') is L, U, 4F4, 16F0,
    16F8 for F' = F0, F8, F4, L, U (h6_split's forms), all constants squares: the (Z/2)^5-covers of a ratio and of
    its partner are isomorphic over Q; (627,364) -> (991,263).
Rests on: sympy's polynomial arithmetic.  usage: python3 -I h6_ident.py [--write]   (h6_ident.txt; a few seconds)"""
import os, sys
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
O = Out('h6_ident.txt')
sys.path.insert(0, H6)
import h6_split
O('Review z9: identities for H6 (generic m, n). Written by an AI model; not read by a person.')
t, m, n = sp.symbols('t m n')
def mine(m, n):
    s = (m * m + n * n) / sp.Integer(2); F4 = s * (1 + t ** 2)
    return {'F0': m * m + n * n * t ** 2, 'F8': n * n + m * m * t ** 2, 'F4': F4, 'L': F4 - 2 * m * n * t, 'U': F4 + 2 * m * n * t}
def theirs(m, n): return {k: sum(c * t ** i for i, c in enumerate(v)) for k, v in h6_split.forms(m, n).items()}
F, G = mine(m, n), theirs(m, n)
O('h6_split.forms = (1, 1, 4, 4, 4) x the forms of G1_FULL_SYSTEM.md: %s' % all(sp.expand(G[k] - c * F[k]) == 0 for k, c in zip(('F0', 'F8', 'F4', 'L', 'U'), (1, 1, 4, 4, 4))))
F = G
sig = all(sp.expand(F[a].subs(t, -t) - F[b]) == 0 for a, b in (('F0', 'F0'), ('F8', 'F8'), ('F4', 'F4'), ('L', 'U'), ('U', 'L')))
tau = all(sp.expand(sp.expand(t ** 2 * F[a].subs(t, 1 / t)) - F[b]) == 0 for a, b in (('F0', 'F8'), ('F8', 'F0'), ('F4', 'F4'), ('L', 'L'), ('U', 'U')))
O('sigma: t -> -t fixes F0, F8, F4 and swaps L, U: %s' % sig)
O('tau: t -> 1/t (times t^2) swaps F0, F8 and fixes F4, L, U: %s' % tau)
P = lambda *ks: sp.expand(sp.prod([F[k] for k in ks]))
A = P('F0', 'F4', 'L')
O('F0F4L(-t) = F0F4U(t): %s' % (sp.expand(A.subs(t, -t) - P('F0', 'F4', 'U')) == 0))
O('t^6 F0F4L(1/t) = F8F4L(t): %s' % (sp.expand(sp.expand(t ** 6 * A.subs(t, 1 / t)) - P('F8', 'F4', 'L')) == 0))
O('t^6 F0F4L(-1/t) = F8F4U(t): %s' % (sp.expand(sp.expand(t ** 6 * A.subs(t, -1 / t)) - P('F8', 'F4', 'U')) == 0))
Gp = theirs(m + n, m - n)
O('partner (m + n, m - n), t\' = (t + 1)/(t - 1): (t - 1)^2 F\'(t\') / F_k(t), constant in t:')
for a in ('F0', 'F8', 'F4', 'L', 'U'):
    e = sp.cancel(sp.together((t - 1) ** 2 * Gp[a].subs(t, (t + 1) / (t - 1))))
    hits = [(b, sp.factor(sp.cancel(e / F[b]))) for b in F if sp.cancel(e / F[b]).free_symbols <= {m, n}]
    O('  %s -> %s' % (a, ', '.join('%s times %s (a square: %s)' % (b, q, sp.sqrt(q).is_rational) for b, q in hits)))
O('627 + 364 = %d, 627 - 364 = %d, gcd(991, 263) = %d; 627^2 + 364^2 = %s, 991^2 + 263^2 = %s' % (
    627 + 364, 627 - 364, sp.gcd(991, 263), sp.factorint(627 ** 2 + 364 ** 2), sp.factorint(991 ** 2 + 263 ** 2)))
O.finish()
