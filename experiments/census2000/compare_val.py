"""Compare dm_pari.py closures on the 157 pairs with the repo's records (dm_log.jsonl: E04, the repo's dm.py;
dm2 lines of dm_closures.txt: D1; dm_sigma_306.jsonl and dm_sigma_24.jsonl: E+, E-)."""
import json, collections
R='/tmp/claude-0/-home-user/9ec7a8d3-42dc-53c8-8518-1a663ac809f1/scratchpad/repo/scripts/g1_full/'
mine = collections.defaultdict(set); mine_rank = {}
for l in open('val157.jsonl'):
    r = json.loads(l); key = tuple(r['of'])
    if r.get('status') == 'closed': mine[key].add({'E':'E','D1':'D1','Sig+':'E+','Sig-':'E-'}[r['curve'].replace('-partner','')])
    if r.get('status') == 'ALERT': print('ALERT', r)
repo = collections.defaultdict(set)
for l in open(R+'dm_log.jsonl'):
    r = json.loads(l)
    if r.get('status') == 'closed': repo[tuple(int(v) for v in str(r['of']).split(','))].add('E') if isinstance(r['of'], str) else repo[tuple(r['of'])].add('E')
for l in open(R+'dm_closures.txt'):
    if l.startswith('#') or not l.strip(): continue
    a = l.split()
    if len(a) < 3: continue
    if a[2] == 'D1': repo[tuple(map(int, a[0].split(',')))].add('D1')
    if a[2] == 'E04': repo[tuple(map(int, a[0].split(',')))].add('E')
for f in ('dm_sigma_306.jsonl', 'dm_sigma_24.jsonl'):
    for l in open(R+f):
        r = json.loads(l)
        if r.get('status') == 'closed': repo[tuple(r['pair'])].add(r['curve'])
pairs = [tuple(p) for p in json.load(open('val157.json'))['open']]
agree = dis = 0
for p in pairs:
    for c in ('E', 'D1', 'E+', 'E-'):
        a, b = c in mine[p], c in repo[p]
        if c in ('E+', 'E-') and not any(c2 in repo[p] for c2 in ('E+','E-')) and b is False:
            pass
        if a == b: agree += 1
        else: dis += 1; print('differ', p, c, 'mine' if a else '', 'repo' if b else '')
print('agree', agree, 'differ', dis)
print('pairs closed by some curve: mine', sum(bool(mine[p]) for p in pairs), 'repo', sum(bool(repo[p] & {'E','D1','E+','E-'}) for p in pairs))
