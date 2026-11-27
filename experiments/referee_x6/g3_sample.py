"""Review x6 (run 6), written by an AI model; not read by a person.  10 rows of CENSUS_TABLE_G3.tsv against the raw stage records (my own parser), plus every OPEN/special row.
Prints only.  usage: python3 -I g3_sample.py"""
import sys, json, os, glob, random, math
D = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'roadmap-g', 'g3')
rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(D, 'CENSUS_TABLE_G3.tsv'))][1:]
rows = [r for r in rows if int(r[0]) <= 13000]   # run 7 (conductor): the table now runs to 14000; this check keeps the range x6 reviewed
def partner(m, n):
    a, b = m + n, abs(m - n)
    if a % 2 == 0 and b % 2 == 0: a, b = a//2, b//2
    g = math.gcd(a, b); return (max(a, b)//g, min(a, b)//g)
sd = {}
for f in glob.glob(os.path.join(D, 'sd_*_*.jsonl')):
    for l in open(f):
        if l.startswith('{'):
            o = json.loads(l); sd.setdefault(tuple(o['pair']), []).append(o)
            if 'partner' in o: sd.setdefault(tuple(o['partner']['pair']), []).append(o['partner'])
r0 = {}
for f in glob.glob(os.path.join(D, 'r0_*_out*.json')):
    for k, v in json.load(open(f)).items(): r0[tuple(map(int, k.split(',')))] = v
dm = {}
for f in glob.glob(os.path.join(D, 'dm_*_*.jsonl')):
    for l in open(f):
        if l.startswith('{'):
            o = json.loads(l); dm.setdefault(tuple(o['pair']), []).append(o)
rnd = random.Random(2026)
pick = rnd.sample(rows, 10) + [r for r in rows if r[3] in ('open', 'special')]
for r in pick:
    m, n = int(r[0]), int(r[1]); q = partner(m, n); s = m*m + n*n
    typ = 'square' if math.isqrt(s)**2 == s else ('twice' if s % 2 == 0 and math.isqrt(s//2)**2 == s//2 else 'OTHER')
    ev = []
    for who in ((m, n), q):
        for o in sd.get(who, []): ev.append('sd%s: %s by %s bounds %s t %s' % (who, o['status'], o.get('by'), o.get('rank_bounds'), o.get('rational_t')))
        for c in r0.get(who, []):
            if c.get('closes'): ev.append('r0%s: closes on %s' % (who, c.get('curve')))
        for o in dm.get(who, []):
            if o.get('status') == 'closed': ev.append('dm%s: closed on %s rank %s ncand %s' % (who, o.get('curve'), o.get('rank'), o.get('ncand')))
    print('ROW %s,%s %s [%s] step=%s | %s | %s' % (m, n, typ, 'type ok' if typ == r[2] else 'TYPE MISMATCH', r[3], r[4][:90], r[5][:60]))
    for e in ev[:6]: print('    ' + e[:200])
