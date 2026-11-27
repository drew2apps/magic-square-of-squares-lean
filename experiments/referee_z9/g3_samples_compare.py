"""Review z9 of roadmap H (run 7), written by an AI model; not read by a person.  Experiment, not counted.
Compares the reruns of g3_samples.sh (in DIR) with G3's committed records (experiments/roadmap-g/g3): identical apart
from "sec", or the fields that differ.  Rests on: the reruns (PARI 2.17.4 here; the records were made with 2.19.1).
usage: python3 -I g3_samples_compare.py DIR [--write]   (g3_samples_compare.txt)"""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import *
D = [a for a in sys.argv[1:] if not a.startswith('--')][0]
O = Out('g3_samples_compare.txt')
O('Review z9: G3 block 13000 sample reruns against the committed records. Written by an AI model; not read by a person.')
def strip(o):
    o = dict(o); o.pop('sec', None)
    if o.get('partner'): o['partner'] = strip(o['partner'])
    return o
def dif(a, b): return 'identical apart from sec' if strip(a) == strip(b) else 'DIFFERENT in %s' % ', '.join('%s: %s here, %s committed' % (k, a.get(k), b.get(k)) for k in sorted(set(a) | set(b)) if k != 'sec' and a.get(k) != b.get(k))
ref = {}
for f in glob.glob(os.path.join(G3, 'sd_*_*.jsonl')):
    for l in open(f):
        if l.startswith('{'): o = json.loads(l); ref[tuple(o['pair'])] = o
for l in open(os.path.join(D, 'sample_sd.jsonl')):
    if l.startswith('{'): o = json.loads(l); O('second descent %s %s: %s' % (o['pair'], o['status'], dif(o, ref[tuple(o['pair'])])))
r0 = {}
for f in glob.glob(os.path.join(G3, 'r0_*_out*.json')): r0.update(json.load(open(f)))
for k, v in sorted(json.load(open(os.path.join(D, 'r0_out.json'))).items()):
    O('rank 0 %s closes by %s: %s' % (k, [c['curve'] for c in v if c['closes']] or 'none', 'identical' if v == r0.get(k) else 'DIFFERENT'))
dm = {}
for f in glob.glob(os.path.join(G3, 'dm_*_*.jsonl')):
    for l in open(f):
        if l.startswith('{'): o = json.loads(l); dm.setdefault((tuple(o['pair']), o['curve'], tuple(o.get('of', []))), o)
for l in open(os.path.join(D, 'dm_out.jsonl')):
    if l.startswith('{'):
        o = json.loads(l); r = dm.get((tuple(o['pair']), o['curve'], tuple(o.get('of', []))))
        O('Demjanenko-Manin %s %s (of %s) %s: %s' % (o['pair'], o['curve'], o.get('of'), o['status'], dif(o, r) if r else 'NO COMMITTED RECORD'))
O.finish()
