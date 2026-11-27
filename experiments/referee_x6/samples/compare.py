"""Review x6 (run 6), written by an AI model; not read by a person.  The reruns in this folder (samples.sh) against the
records committed in roadmap-g/g3 and g4: identical apart from "sec", or the differing fields.  Prints only.
usage: python3 -I compare.py [DIR]   (DIR: where the reruns are; default this folder)"""
import json, glob, os, sys
H = os.path.dirname(os.path.abspath(__file__)); D = sys.argv[1] if len(sys.argv) > 1 else H
G3 = os.path.join(H, '..', '..', 'roadmap-g', 'g3'); G4 = os.path.join(H, '..', '..', 'roadmap-g', 'g4')
def strip(o):
    o = dict(o); o.pop('sec', None)
    if 'partner' in o: o['partner'] = strip(o['partner'])
    return o
ref = {}
for f in glob.glob(os.path.join(G3, 'sd_*_*.jsonl')):
    for l in open(f):
        if l.startswith('{'): o = json.loads(l); ref[tuple(o['pair'])] = o
for l in open(os.path.join(D, 'sample_sd.jsonl')):
    if l.startswith('{'):
        o = json.loads(l); r = ref[tuple(o['pair'])]
        print('second descent', o['pair'], 'identical apart from sec' if strip(o) == strip(r) else 'DIFFERENT')
r0 = {}
for f in glob.glob(os.path.join(G3, 'r0_*_out*.json')): r0.update(json.load(open(f)))
for k, v in json.load(open(os.path.join(D, 'r0_out.json'))).items():
    print('rank 0', k, 'identical' if v == r0.get(k) else 'DIFFERENT')
dm = {}
for f in glob.glob(os.path.join(G3, 'dm_*_*.jsonl')):
    for l in open(f):
        if l.startswith('{'): o = json.loads(l); dm.setdefault((tuple(o['pair']), o['curve'], tuple(o.get('of', []))), o)
for l in open(os.path.join(D, 'dm_out.jsonl')):
    if l.startswith('{'):
        o = json.loads(l); r = dm.get((tuple(o['pair']), o['curve'], tuple(o.get('of', []))))
        print('Demjanenko-Manin', o['pair'], o['curve'], o['status'], 'identical apart from sec' if r and strip(o) == strip(r) else 'DIFFERENT')
a = [json.loads(l) for l in open(os.path.join(D, 'log_12416.txt')) if l.startswith('{')][-1]
b = [json.loads(l) for l in open(os.path.join(G3, 'route', 'log_12416.txt')) if l.startswith('{')][-1]
print('route (12416,5313): closed %s; %s' % (a.get('closed'), 'identical apart from sec' if strip(a) == strip(b) else 'DIFFERENT in %s' % [k for k in a if a.get(k) != b.get(k)]))
g4 = {}
for f in sorted(glob.glob(os.path.join(G4, 'g4_*_3.json'))):
    for r in json.load(open(f)): g4[tuple(r['ratio'])] = r
for l in open(os.path.join(D, 'g4_rerun.jsonl')):
    if not l.startswith('{'): continue
    r = json.loads(l); q = g4[tuple(r['ratio'])]
    dif = [c for c in r['curves'] if r['curves'][c] != q['curves'][c]]
    print('G4 search', r['ratio'], 'hits %d (committed %d), solutions %s (committed %s);' % (len(r['hits']), len(q['hits']), r['solutions'], q['solutions']),
          'every curve identical' if not dif else 'curves that differ (PARI 2.17.4 against 2.19.1): %s' % dif)
