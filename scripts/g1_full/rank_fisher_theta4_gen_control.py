#!/usr/bin/env python3
"""External controls for rank_fisher_theta4_gen.py: rank_fisher_theta4_control.py (four rank-0 curves whose Theta_4
is predicted by the BSD order of Sha: A, B non-zero; C, D zero) re-run with the GENERAL module's pieces in place of
rank_fisher_theta4's: setup_gen (no (72,65)-specific assertion), choose_Xi_lin (level-4 twists by F_2 linear
algebra), theta4_value (with the extra f0-degeneracy primes in S).  The stored conic points of
rank_fisher_theta4_control_hints.json are used (re-checked exactly by rank_fisher5.analyse).

    python3 rank_fisher_theta4_gen_control.py [--write]      (rank_fisher_theta4_gen_control.txt)

Written by an AI model (Claude); not read by a person.  Experiment, not counted.  The predictions are conditional
on the BSD formula for the order of Sha (rank 0 itself is proved), as in rank_fisher_theta4_control.py.
"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
WRITE = '--write' in sys.argv
if WRITE: sys.argv.remove('--write')
import rank_fisher_theta4 as T
import rank_fisher_theta4_gen as TG
import rank_fisher_theta4_control as CT

def choose_Xi_wrap(M, D, C, say=None, places=None, images=None):
    sols, K, Rs = TG.choose_Xi_lin(M, D, C, places, say=lambda s: None)
    return sols, K, Rs

def theta4_wrap(M, D, C, Xi, ys, say=None, label='', B0=None, kmul=None):
    out, S, _ = TG.theta4_value(M, D, C, Xi, ys, B0, label=label, say=say)
    return out, S

T.setup = TG.setup_gen
T.choose_Xi = choose_Xi_wrap
T.theta4_for = theta4_wrap
CT.OUT = TG.OUT          # one transcript
CT.say = TG.say
if __name__ == '__main__':
    CT.WRITE = False
    TG.say('rank_fisher_theta4_gen_control.py: the four external controls of rank_fisher_theta4_control.py through the '
           'general module (setup_gen, choose_Xi_lin, theta4_value); its transcript follows')
    CT.main()
    if WRITE: open(os.path.join(here, 'rank_fisher_theta4_gen_control.txt'), 'w').write('\n'.join(TG.OUT) + '\n')
