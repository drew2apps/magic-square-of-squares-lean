#!/usr/bin/env python3
"""z11: the whole Cassels-Tate pairing matrix on Sel^2 of three pair curves of (5623,1017), my code (z7's pairing).
Written by an AI model (review agent z11); not read by a person.  Backs t10's "the pairing trivial on Sel^2" of E_0L
and the bound rank E_LU <= 1.  PARI as in z11_pair.py (no ellrank here).
    python3 -I pmat.py      (prints pmat.txt)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import Out
from z11_pair import PairCurve, Pairing, CT
say = Out('pmat.txt')
say('pmat.py (review z11; written by an AI model): pairing matrices <b, D_c> on a basis of Sel^2')
for m, n, jk in [(5623, 1017, '08'), (5623, 1017, '0L'), (5623, 1017, 'LU')]:
    C = PairCurve(m, n, jk); PR = Pairing(C.es)
    M = [[PR.pair(b, c)[0] for b in PR.sel] for c in PR.sel]
    rk = CT.f2rank(M)
    sym = all(M[i][j] == M[j][i] for i in range(len(M)) for j in range(len(M))) and all(M[i][i] == 1 for i in range(len(M)))
    say('(%d,%d) E_%s: dim Sel^2 %d, pairing matrix of rank %d (symmetric and alternating: %s), so rank E(Q) <= %d'
        % (m, n, jk, PR.dim, rk, sym, PR.dim - 2 - rk))
    for r in M: say('    ' + ' '.join('%+d' % v for v in r))
say.finish()
