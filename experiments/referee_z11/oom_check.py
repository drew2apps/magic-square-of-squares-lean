#!/usr/bin/env python3
"""z11: t10_997_a.txt was cut from the log of a run killed for lack of memory after 79 ratios.  Rerun t10's committed
script (../conductor-t10/t10_pairct.py, without --write: it only prints) on ratios 35..40 and 74..78 of
g2_997_all.txt and compare each ratio's block line by line with t10_997_a.txt.  Written by an AI model (review agent
z11); not read by a person.  PARI: t10's code as it is (its stderr and its stack message are dropped).
    python3 -I oom_check.py      (prints oom_check.txt)"""
import sys, os, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _z11common import T10, G2, REPO, Out
say = Out('oom_check.txt')
say('oom_check.py (review z11; written by an AI model): t10_pairct.py rerun against t10_997_a.txt')
def blocks(lines):
    out = {}; cur = None
    for l in lines:
        if l.startswith('ratio ('): cur = l.split()[1].rstrip(':'); out[cur] = [l]; continue
        if l.startswith('SUMMARY') or l == '': cur = None
        if cur: out[cur].append(l)
    return out
ref = blocks(open(os.path.join(T10, 't10_997_a.txt')).read().split('\n'))
allok = True
for lo, hi in ((35, 41), (74, 79)):
    env = dict(os.environ, MSS_REPO=REPO)
    r = subprocess.run([sys.executable, '-I', os.path.join(T10, 't10_pairct.py'), '--list', os.path.join(G2, 'g2_997_all.txt'),
                        '--from', str(lo), '--to', str(hi), '--out', 'z11rerun'], cwd=T10, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-500:]
    got = blocks([l for l in r.stdout.split('\n') if not l.startswith('PARI stack')])
    for k, v in got.items():
        same = ref.get(k) == v; allok &= same
        say('  ratio %s (index %d..%d): %d lines, identical to t10_997_a.txt: %s' % (k, lo, hi - 1, len(v), same))
say('ALL RERUN BLOCKS IDENTICAL' if allok else 'A BLOCK DIFFERS')
say.finish()
