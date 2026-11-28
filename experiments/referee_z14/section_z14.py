"""z14 task 3: E_mix has a section of infinite order: the class of (iota t0) - (t0) on C_0L, t0 = m/n (t10_pairct.Pair's
psi(iota t0)), is a section of the family over the r-line; one specialization of infinite order proves the section has
infinite order (torsion sections specialize to torsion points).  PARI's ellorder (0 = infinite order) and
elltors are used as a finder and checked by multiplying out (the point times the torsion exponent is not O).
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I section_z14.py | diff - section_z14.txt"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
PARI = T.PARI
print('# z14 section check. Written by an AI model (reviewer z14); not read by a person.')
for m, n in [(2, 1), (3, 2), (4, 1), (7, 4), (1234, 567)]:
    Pr = T.Pair(m, n, '0', 'L')
    P = Pr.iota_base()
    E = PARI.ellinit([0, Pr.A[0], 0, Pr.A[1], Pr.A[2]])
    Pp = [PARI(str(P[0])), PARI(str(P[1]))]
    assert PARI.ellisoncurve(E, Pp)
    o = int(PARI.ellorder(E, Pp)); tors = PARI.elltors(E)
    N = int(tors[0]); Q = PARI.ellmul(E, Pp, N)
    print((m, n), 'psi(iota t0) =', (str(P[0])[:30], '...'), 'ellorder', o, '| #E(Q)_tors', N, tors[1], '| N P = O:', str(Q) == '[0]')
