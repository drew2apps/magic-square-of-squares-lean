#!/usr/bin/env python3
"""Review x4 (run 6): the torsion step of the (6557,324) closure without PARI: census_lvalue_exact.torsion_step (torsion_exact.py's
certificate + global.py's own classify_t/analyse, with rank 0 assumed from the L-value) on w*f0*f8 of (6557,324).  No PARI.  Written by an AI model; not read by a person.   python3 tors6557.py  (prints only)"""
import sys, os
R = os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
sys.path.insert(0, R + '/scripts/g1_full'); sys.path.insert(0, R + '/scripts/second_descent')
import census_lvalue_exact as CE
G = CE.global_py()
a = CE.model(CE.poly(6557, 324, 'w*f0*f8'))
print('model', a)
rec = CE.torsion_step((6557, 324), 'w*f0*f8', a, G)
import torsion_exact as TE
print(CE.torsion_line(rec, len(TE.good_odd_primes(TE.Curve(a[1], a[3], a[4]), 1000))))
print('closes (rank 0 assumed):', rec.get('closes'), 'alert:', rec.get('alert'))
assert rec['cert']['ok'] and rec['closes'] and not rec['alert'] and list(a) == [0, 1848516467034577, 0, 20370457151276745574584576, 0]
