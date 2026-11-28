"""z15: move the timing fields (cpu, wall) of the raw t10 run records t10_*.jsonl into t10_times.txt, so that the
committed records are deterministic and diffable.  Run once after the runs.  Written by an AI model (reviewer z15),
not read by a person."""
import os, json
HERE = os.path.dirname(os.path.abspath(__file__))
out = ['# t10_pairct CPU and wall-clock seconds per ratio (machine shared, load ~7 on 2 cores; not deterministic)']
for f in sorted(os.listdir(HERE)):
    if not (f.startswith('t10_') and f.endswith('.jsonl')): continue
    recs = [json.loads(l) for l in open(os.path.join(HERE, f))]
    for r in recs:
        if 'cpu' in r: out.append('%s %d,%d cpu %.1f wall %s' % (f, r['m'], r['n'], r.pop('cpu'), r.pop('wall', '?')))
    with open(os.path.join(HERE, f), 'w') as fo:
        for r in recs: fo.write(json.dumps(r, sort_keys=True) + '\n')
open(os.path.join(HERE, 't10_times.txt'), 'w').write('\n'.join(out) + '\n')
print(len(out) - 1, 'timing lines moved')
