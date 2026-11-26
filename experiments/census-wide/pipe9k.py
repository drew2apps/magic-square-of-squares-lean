"""Stages after the second descent for 7000 < m <= 9000 (same programs as the 7000 run).  EXPERIMENT."""
import json, sys, subprocess, math
G = '/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full/global.py'
def partner(m, n):
    a, b = m + n, abs(m - n)
    if a % 2 == 0 and b % 2 == 0: a, b = a // 2, b // 2
    g = math.gcd(a, b); a, b = a // g, b // g
    return [max(a, b), min(a, b)]
stage = sys.argv[1]
if stage == 'rank0':
    op = []
    for f in ['s9_0.jsonl','s9_1.jsonl','s9_2.jsonl','s9r_0.jsonl','s9r_1.jsonl','s9r_2.jsonl']:
        for l in open(f):
            if l.startswith('{'):
                o = json.loads(l)
                if o['status'] == 'open': op.append(o['pair'])
    pm = {str(tuple(p)): partner(*p) for p in op}
    json.dump({'pairs': op, 'partner': pm}, open('open9k_map.json', 'w'))
    U = sorted(set(tuple(p) for p in op) | set(tuple(v) for v in pm.values()))
    print(len(op), 'open,', len(U), 'with partners')
    procs = []
    for k in range(3):
        json.dump({'open': [list(u) for u in U[k::3]]}, open(f'g9k_in{k}.json', 'w'))
        procs.append(subprocess.Popen(f'python3 {G} g9k_in{k}.json g9k_out{k}.json > g9k{k}.txt 2>&1', shell=True))
    for p in procs: p.wait()
    closed = set()
    for k in range(3):
        for key, lst in json.load(open(f'g9k_out{k}.json')).items():
            if any(c.get('closes') for c in lst): closed.add(tuple(int(t) for t in key.split(',')))
    left = [p for p in op if tuple(p) not in closed and tuple(pm[str(tuple(p))]) not in closed]
    json.dump(left, open('after_rank0_9k.json', 'w')); print('after rank 0:', len(left))
if stage == 'dm':
    left = json.load(open('after_rank0_9k.json')); pm = json.load(open('open9k_map.json'))['partner']
    U = sorted(set(tuple(p) for p in left) | set(tuple(pm[str(tuple(p))]) for p in left))
    procs = []
    for k in range(3):
        json.dump({'open': [list(u) for u in U[k::3]]}, open(f'dm9k_in{k}.json', 'w'))
        procs.append(subprocess.Popen(f'python3 ../census2000/dm_pari.py dm9k_in{k}.json dm9k{k}.jsonl > dm9k{k}.log 2>&1', shell=True))
    for p in procs: p.wait()
    closed = set(); alerts = []
    for k in range(3):
        for l in open(f'dm9k{k}.jsonl'):
            if l.startswith('{'):
                o = json.loads(l)
                if o.get('status') == 'closed': closed.add(tuple(o['pair']))
                if 'alert' in str(o.get('status')): alerts.append(o)
    left2 = [p for p in left if tuple(p) not in closed and tuple(pm[str(tuple(p))]) not in closed]
    json.dump(left2, open('after_dm_9k.json', 'w')); print('after DM:', len(left2), left2, 'alerts', len(alerts))
