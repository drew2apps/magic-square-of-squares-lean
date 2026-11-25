#!/usr/bin/env python3
"""The exact route of exact_last_pair.py on other ratios of square type (a test of the program, not new closures).

    python3 exact_last_pair_route.py m,n d [model=0,a2,0,a4,a6] [p1 p2 ...]      prints; with --write also appends to exact_last_pair_route.txt

Written by an AI model (Claude); not read by a person.  For a ratio m : n with m^2 + n^2 a square and the class d of
its degenerate solution (route_all.txt), every step of exact_last_pair.py is run on the curves of that ratio: models,
G and the degenerate points, torsion, the exact rank bound, the height enclosure, saturation to 199 and the exact lower
bound for heights, Lemma 2, the bound of Lemmas 4 and 5, the sieve at the given primes (the orders of G modulo the
first two must be coprime), the fibre.  The model used for E1 is the one given (for example Sage's minimal model, as printed in
saturation_*.txt; setup() checks that it is isomorphic over Q to the model of weierstrass(), and S6b that it is
minimal), or the model of weierstrass() if it has integer coefficients.  These ratios are closed in the census by other means; agreement with
route_all.txt (Sage and PARI) tests the exact programs.
"""
import os, sys, json
WRITE = '--write' in sys.argv   # outputs and hint files are written only with --write, so CI can compare stdout with them
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import exact_last_pair as X
import rank_indep as RI

def route(m, n, d, primes, em=None):
    out = []
    def say(s): print(s, flush=True); out.append(s)
    S, D, Q1, r1, Q2, r2 = X.family(m, n)
    E1w, A1 = X.weierstrass([d * c for c in Q1], r1)
    co = (E1w.a2, E1w.a4, E1w.a6)
    if em is not None: X.EM1 = tuple(em)                 # a stated model (checked to be isomorphic over Q by setup)
    else:
        assert all(Fr(c).denominator == 1 for c in co), 'the model of weierstrass() is not integral; give a model'
        X.EM1 = (0, int(co[0]), 0, int(co[1]), int(co[2]))
    # G from the degenerate solution, then the full setup with it
    s = X.GQ(0, -1); dQ1 = [d * c for c in Q1]
    Y = X.gsqrt(X.peval([X.gq(c) for c in dQ1], s)); assert Y is not None
    P = X.psi(A1, X.gq(r1), s, Y); Pc = (X.gq(P[0]).conj(), X.gq(P[1]).conj())
    Gk = E1w.add(P, Pc); Gw0 = (X.gq(Gk[0]).re, X.gq(Gk[1]).re)
    ip = X.iso_params(E1w, X.Ell(X.EM1[1], X.EM1[3], X.EM1[4])); assert ip is not None, 'the model is not E1'
    G = X.to_model(Gw0, *ip)                              # G on the stated model
    St = X.setup(m, n, d, G)
    Em = St['Em']
    T, tq = X.torsion(Em)
    r, _ = X.rank_part(Em, G, say=lambda s: None)
    lo, hi, _, _ = X.enclose(Em, G, 7)
    sat = X.saturation(Em, G, T, 199)
    LX = X.lam_exact(Em, hi / 210 ** 2)
    X.lemma2(St)
    L = X.lemma45(St, lo)
    Gw = X.from_model(G, St['u1'], St['rr1']); Tw = X.from_model(T, St['u1'], St['rr1'])
    loc = {p: X.sieve_local(St, p, Gw, Tw) for p in primes}
    left, Ncov, n12 = X.lift(loc, list(primes), max(30000, 2 * L['ab']))
    ts, _ = X.fibre(St, [Gw, St['E1w'].add(Gw, Tw)])
    sols = [t for t in ts if X.is_solution(St, t)]
    say('%d,%d d=%d: rank <= %d; hhat(G) in [%.4f, %.4f]; saturated to 199; heights >= hhat(G)/210^2 (m = %d, search to %d); '
        'a, b <= %d; sieve (%s) to %d leaves %s (within the bound: %s); fibre: rational t %s, solutions %s'
        % (m, n, d, r, float(lo), float(hi), LX['m'], LX['X0'], L['ab'], ' '.join(str(p) for p in primes), Ncov, left,
           [c for c in left if max(c) <= L['ab']], [str(t) for t in ts], [str(t) for t in sols]))
    within = [c for c in left if max(c) <= L['ab']]
    ok = r == 1 and within == [(1, 1)] and Ncov >= L['ab'] and sorted(sols) == [Fr(-1), Fr(1)]
    say('  the analogue of Theorem A holds for %d,%d by the exact route: %s' % (m, n, ok))
    return ok, out

if __name__ == '__main__':
    m, n = (int(v) for v in sys.argv[1].split(',')); d = int(sys.argv[2])
    em = None
    args = [a for a in sys.argv[3:] if a != '--write']
    if args and args[0].startswith('model='):
        em = [int(v) for v in args[0][6:].split(',')]; args = args[1:]
    primes = [int(v) for v in args] or [10009, 20089, 30013, 40009]
    hf = os.path.join(here, 'exact_last_pair_route_hints.json')
    if os.path.exists(hf):
        for k, v in json.load(open(hf)).items(): RI.HINTS.setdefault(k, {}).update(v)
    ok, out = route(m, n, d, primes, em)
    if WRITE: json.dump(RI.HINTS, open(hf, 'w'), indent=0, sort_keys=True)
    print('\n'.join(out))
    if WRITE:   # appends one ratio's lines; the committed file is three runs: 187,84 / 440,279 / 420,341
        with open(os.path.join(here, 'exact_last_pair_route.txt'), 'a') as f: f.write('\n'.join(out) + '\n')
    sys.exit(0 if ok else 1)
