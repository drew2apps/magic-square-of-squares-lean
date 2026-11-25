#!/usr/bin/env python3
"""The inputs of sigma_exact.py that need Sage, PARI or eclib: ranks, the generator, the index bound.

    python3 sigma_certify.py m,n:E+ m,n:E- ...      one line of sigma_evidence.txt per argument

For each curve (the model of dm_sigma.py):
  rank  PARI's ellrank on the isogeny class (tight_bounds of dm.py): upper bound 1.  Also eclib's mwrank
        (two-descent; on a 2-isogenous curve when it fails on the curve itself), recorded whatever it gives; its
        second descent fails on some of these curves (an NTL error).  The second independent upper bound used for
        counting is Magma's TwoPowerIsogenyDescentRankBound (magma_sigma2.py, magma_sigma2_raw.txt).
  G     PARI's point, saturated by ellsaturation up to the index bound; also checked against eclib's
        generator: the two have the same canonical height and differ by a torsion point or a sign.
  K     the index bound: if G = k G1 + torsion then hhat(G1) = hhat(G)/k^2 >= lam, where lam is the lower bound
        of Cremona and Siksek for the canonical height of a point of infinite order (Sage's height_function().min,
        on the minimal model), so k <= sqrt(hhat(G)/lam).  K = floor(sqrt(1.01 hhat(G)/lam)) + 1.
        sigma_exact.py then proves p-saturation for every prime p <= K itself.
  amax  the bound dm_sigma.py listed to (dm_sigma_24.jsonl); sigma_exact.py recomputes its own bound and
        checks it is not larger.
What is trusted here: the rank programs (PARI, and Magma or eclib, which share no code with it), and Sage's lam."""
import sys, os, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sage.all__sagemath_schemes import EllipticCurve, QQ
from dm import tight_bounds
import dm_sigma


MAGMA = {}
def read_magma(path):
    import re
    for mm, nn, c, r in re.findall(r'<line>BOUNDS (\d+) (\d+) (E[+-]) (\d+)</line>', open(path).read()):
        MAGMA[((int(mm), int(nn)), c)] = int(r)


def eclib(Em, iso):
    """eclib's mwrank in a child process with a time limit (its second descent fails on some of these curves)."""
    import subprocess
    code = ('from sage.libs.eclib.interface import mwrank_EllipticCurve as M\n'
            'mw = M(%s, verbose=False)\nb = mw.rank_bound()\nprint(b)\n'
            'print(" ".join(str(c) for c in mw.gens()[0]) if b == 1 else "")' % str([int(c) for c in Em.ainvs()]))
    try:
        r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=300)
        out = r.stdout.split('\n')
        if r.returncode or not out[0].strip(): return 'failed (%s)' % (r.stderr.strip().split('\n')[-1][:30] if r.stderr else 'no output'), None
        b = int(out[0])
        if b != 1: return 'upper bound %d' % b, None
        c = [QQ(v) for v in out[1].split()]
        Z = c[2] if len(c) == 3 else 1
        return 'upper bound 1', iso(Em([c[0]/Z, c[1]/Z]))
    except subprocess.TimeoutExpired:
        return 'no answer within 300 s', None


def certify(m, n, cv, rec):
    s = 2 if cv == 'E+' else -2
    Z, q = dm_sigma.quartic(m, n, s); z0 = -s
    d = [q.derivative(k)(z0) for k in range(5)]; A = d[1]
    E = EllipticCurve(QQ, [0, d[2]//2, 0, A*d[3]//6, A*A*d[4]//24])
    Em = E.minimal_model(); iso = Em.isomorphism_to(E); Ep = Em.pari_curve()
    lo, hi, _ = tight_bounds(E)
    magma = MAGMA.get(((m, n), cv))
    assert lo == 1 and (hi == 1 or magma == 1), (lo, hi, magma)
    assert magma == 1, 'no Magma bound'
    rk = 'PARI[%d,%d];Magma_TwoPowerIsogenyDescentRankBound:%d' % (lo, hi, magma)
    pts = []
    for effort in range(6):
        pts = Ep.ellrank(effort)[3]
        if len(pts) >= 1: break
    P0 = Em([QQ(pts[0][0]), QQ(pts[0][1])]); h0 = float(P0.height())
    lam = float(Em.height_function().min(0.01, 20)); assert lam > 0
    k0 = int(math.floor(math.sqrt(h0/lam))) + 1
    sat = Ep.ellsaturation([pts[0]], k0)
    G = iso(Em([QQ(sat[0][0]), QQ(sat[0][1])]))
    hG = float(G.height())
    K = int(math.floor(math.sqrt(1.01*hG/lam))) + 1
    ecl, mwgen = eclib(Em, iso)
    if mwgen is not None:
        assert abs(float(mwgen.height()) - hG) < 1e-9*hG, 'eclib generator has another height'
        assert any((G - e*mwgen).has_finite_order() for e in (1, -1)), 'eclib generator differs'
        gen = 'PARI, saturated; equal to eclib\'s generator up to sign and torsion'
    else:
        gen = 'PARI, saturated'
    x, y = G.xy()
    return ('%d,%d %s G=%s:%s K=%d amax=%d rank=%s;eclib:%s lam=%.6g gen=%s'
            % (m, n, cv, x, y, K, rec['amax'], rk, ecl.replace(' ', '_'), lam, gen.replace(' ', '_')))


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    rec = {}
    for f in ('dm_sigma_306.jsonl', 'dm_sigma_24.jsonl'):
        for l in open(os.path.join(here, f)):
            r = json.loads(l); rec[(tuple(r['pair']), r['curve'])] = r
    read_magma(os.path.join(here, 'magma_sigma2_raw.txt'))
    read_magma(os.path.join(here, 'magma_sigma3_raw.txt'))
    for arg in sys.argv[1:]:
        pr, cv = arg.split(':'); m, n = (int(v) for v in pr.split(','))
        r = rec[((m, n), cv)]; assert r['status'] == 'closed'
        print(certify(m, n, cv, r), flush=True)
