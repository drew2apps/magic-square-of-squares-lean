"""Reviewer z13 (revision 3 follow-up): the five curves of the note's h2_check.py, recomputed with z13_indep.run
(rank_second's full_selmer2 and DFull pairing against the reviewer's own formula (2)), because h2_check.py imports a
module 'pipeline' from /tmp paths that were not supplied.  Counts: pairs (a, b), a = (1, d, d) in Sel^2, b over all of
Sel^2.  Written by an AI model (reviewer z13), not read by a person."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = sys.argv[:1]
sys.path.insert(0, HERE)
import z13_indep as Z
for m, n, cv in ((725, 497, 'LU'), (725, 497, '08'), (535, 23, '08'), (465, 416, 'LU'), (465, 416, '08')):
    r = Z.run(m, n, cv, 'h2_check')
    print((m, n), cv, '(A,C)=', (r['A'], r['C']), 'H2' if r['H2'] else 'H2 fails', 'dim Sel2', r['dim'],
          'agree', r['cmp']['agree'], 'disagree', r['cmp']['disagree'], 'undefined', r['cmp']['undefined'],
          "S' = {b_r = 1}", r['Sprime_eq'], 'Theta agree/disagree', r['theta']['agree'], r['theta']['disagree'])
