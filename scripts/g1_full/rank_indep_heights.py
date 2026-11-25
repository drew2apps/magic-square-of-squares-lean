#!/usr/bin/env python3
"""The generator of each curve of sigma_evidence.txt without Sage's lower bound for the canonical height.

    python3 rank_indep_heights.py              output: rank_indep_heights.txt
    python3 rank_indep_heights.py --negative   a curve of rank 2 must be refused; output rank_indep_heights_negative.txt
    (add --write to rewrite the output files and the hints; without it the script only prints, so CI can diff its output)

Written by an AI model (Claude); not read by a person.

sigma_exact.py proves saturation of G at every prime up to an index bound K that comes from Sage's lower bound lam
(Cremona-Siksek) for the heights of points of infinite order: the one floating-point input of the 18 closures of
G1_FULL_SYSTEM.md, Result 4.  Here lam is replaced, curve by curve, by the exact argument of exact_last_pair.py (S6b):
on the minimal model (the a-invariants sent to Magma, magma_sigma*.json, checked to be the same curve and minimal at
every bad prime), every point P has 24 P in E0 at every prime (multiplicative types I_n with n | 24, additive types
with a component group of order <= 4); for such a point h(x) - hhat <= log+|j|/4 + log+|b2/12| + log 2 + 1.946
(Silverman's archimedean bound as quoted by Cremona-Prickett-Siksek, every non-archimedean term 0 on E0: the
statement height_enclosure.py uses); a search of the small x = a/e^2 finds no point.  So hhat(P) >= hhat(G)/210^2 for
every P of infinite order; with saturation at every prime <= 199 (proved here again, exactly), and IF rank E(Q) = 1,
G generates E(Q) modulo torsion (G = k G0 + t with k > 1 would need k >= 211 and hhat(G0) < hhat(G)/210^2).
hhat(G) is enclosed by Silverman's bound and exact doubling.
The rank: the script calls the exact descent of rank_indep.py and states "G generates" without condition only where
that bound is 1 (four curves); for the other eight the exact bound is 3, and the statement is printed as conditional
on rank 1, with the source of rank 1 named from sigma_evidence.txt (PARI's ellrank, Magma's
TwoPowerIsogenyDescentRankBound, eclib).  generates() refuses unless the exact bound is 1; --negative checks this on a
curve of rank 2 (found by a reviewer) on which every other step passes.
The height steps use no Sage, PARI or Magma; the statements quoted are Cremona-Prickett-Siksek (J. Number Theory 116
(2006)) Theorems 2 and 15, Silverman (1990), and the Kodaira-Neron bounds (tier E of LEDGER.md).
"""
import os, sys
WRITE = '--write' in sys.argv   # outputs and hint files are written only with --write, so CI can compare stdout with them
from fractions import Fraction as Fr
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..', 'second_descent'))
import sigma_exact as SE
import rank_indep as RI
import exact_last_pair as X

HINTS_FILE = os.path.join(here, 'rank_indep_hints.json')

def height_steps(Em, G, T, lmax=199):
    """Torsion Z/2, the enclosure of hhat(G), saturation at every prime <= lmax and the exact lower bound for heights.
    Together with rank 1 these prove that G generates E(Q) modulo torsion; alone they do not."""
    g = 0
    for q in range(3, 600):
        if RI.is_prime_td(q) and Em.D % q:
            g = X.gcd(g, X.count_mod(Em, q))
    assert g == 2, ('torsion', g)
    lo, hi, B, U = X.enclose(Em, G, 7)
    assert lo > 0
    X.saturation(Em, G, T, lmax)
    nxt = lmax + 1
    while not RI.is_prime_td(nxt): nxt += 1
    need = hi / (nxt - 1) ** 2
    L = X.lam_exact(Em, need)
    return dict(lo=lo, hi=hi, m=L['m'], Binf=L['Binf'], X0=L['X0'], need=need, nxt=nxt)

def exact_rank(Em, G):
    e, a, b = RI.shift_to_two_torsion(Em.a2, Em.a4, Em.a6)
    GE = (G[0] - e, G[1])
    ap, bp = -2 * a, a * a - 4 * b
    r, _ = RI.rank_bound(a, b, [GE], RI.find_points([0, ap, 0, bp, 0]), say=lambda s: None)
    return r

def generates(Em, G, T, lmax=199):
    """G generates E(Q) modulo torsion: refused unless the exact rank bound is 1."""
    r = exact_rank(Em, G)
    assert r == 1, ('exact rank bound %d: that G generates would need rank 1, which is not proved here' % r)
    return height_steps(Em, G, T, lmax)

def curve(m, n, cv, G):
    M = SE.model(m, n, 2 if cv == 'E+' else -2)
    Es = X.Ell(M['a2'], M['a4'], M['a6'])
    ai = RI.magma_ainvs()[((m, n), cv)]
    assert ai[0] == 0 and ai[2] == 0
    Em = X.Ell(ai[1], ai[3], ai[4])
    ip = X.iso_params(Es, Em); assert ip is not None, 'the model of Magma is not the same curve over Q'
    Gm = X.to_model(G, *ip); assert Em.on(Gm)
    e, _, _ = RI.shift_to_two_torsion(Em.a2, Em.a4, Em.a6); T = (Fr(e), Fr(0))
    return Em, Gm, T

def rank_sources():
    out = {}
    for line in open(os.path.join(here, 'sigma_evidence.txt')):
        if line.startswith('#') or not line.strip(): continue
        w = line.split(); f = dict(x.split('=', 1) for x in w[2:])
        m, n = (int(v) for v in w[0].split(','))
        out[((m, n), w[1])] = f['rank']
    return out

def negative(say):
    """A reviewer's curve of rank >= 2 (two points with height-matrix determinant 238.06, by PARI): every height step
    passes for its point G, so the conclusion "G generates" must come from the rank and be refused."""
    Em = X.Ell(0, -41159003824935, -713237279244058)
    G = (Fr(-33478094, 9), Fr(272194294168, 27)); assert Em.on(G)
    e, _, _ = RI.shift_to_two_torsion(Em.a2, Em.a4, Em.a6); T = (Fr(e), Fr(0))
    hs = height_steps(Em, G, T)
    say('rank-2 curve y^2 = x^3 - 41159003824935 x - 713237279244058, G = (-33478094/9, 272194294168/27): the height steps '
        'all pass (hhat(G) in [%.4f, %.4f], saturated to 199, search to %d), as the reviewer found' % (float(hs['lo']), float(hs['hi']), hs['X0']))
    try:
        generates(Em, G, T)
        say('  NOT REFUSED: "G generates" on a curve of rank 2'); return False
    except AssertionError as ex:
        say('  refused: %s' % str(ex)[:120]); return True

if __name__ == '__main__':
    import json
    if os.path.exists(HINTS_FILE):
        for k, v in json.load(open(HINTS_FILE)).items(): RI.HINTS.setdefault(k, {}).update(v)
    out = []
    def say(s): print(s, flush=True); out.append(s)
    if '--negative' in sys.argv:
        ok = negative(say)
        say('negative tests: %d of 1 refused' % int(ok))
        if WRITE: json.dump(RI.HINTS, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
        if WRITE: open(os.path.join(here, 'rank_indep_heights_negative.txt'), 'w').write('\n'.join(out) + '\n')
        sys.exit(0 if ok else 1)
    say('rank_indep_heights.py: the generators of sigma_evidence.txt without Sage\'s lower bound (exact), and the rank')
    src = rank_sources()
    n_unc = n_cond = 0
    for (mn, cv), G in RI.evidence().items():
        m, n = mn
        Em, Gm, T = curve(m, n, cv, G)
        r = exact_rank(Em, Gm)
        hs = height_steps(Em, Gm, T)
        base = ('%d,%d %s: hhat(G) in [%.6f, %.6f]; torsion Z/2; saturated at every prime <= 199; minimal model, m = %d, '
                'h(x) - hhat <= %.3f on E0, no point with x = a/e^2, |a|, e^2 <= %d; so hhat(P) >= %.3e for every P of '
                'infinite order' % (m, n, cv, float(hs['lo']) - 5e-7, float(hs['hi']) + 5e-7, hs['m'], float(hs['Binf']), hs['X0'], float(hs['need'])))
        if r == 1:
            generates(Em, Gm, T); n_unc += 1
            say(base + '; exact rank bound 1 (rank_indep.py): G GENERATES E(Q) modulo torsion, no Sage, PARI or Magma')
        else:
            n_cond += 1
            say(base + '; exact rank bound %d only: IF rank E(Q) = 1, G generates; rank 1 rests on %s' % (r, src[(mn, cv)].replace('_', ' ')))
    say('rank_indep_heights: %d curves; G generates unconditionally (exact rank 1) for %d; for %d the conclusion is '
        'conditional on rank 1 (from PARI and Magma, or from Magma alone; eclib gives no bound on these); Sage\'s lam is used for none' % (n_unc + n_cond, n_unc, n_cond))
    if WRITE: json.dump(RI.HINTS, open(HINTS_FILE, 'w'), indent=0, sort_keys=True)
    if WRITE: open(os.path.join(here, 'rank_indep_heights.txt'), 'w').write('\n'.join(out) + '\n')
