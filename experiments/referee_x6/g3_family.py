"""Review x6 (run 6), written by an AI model; not read by a person.  The G3 family for 11000 < m <= 13000 by brute force
over all pairs (numpy; not via Pythagorean parameters), against roadmap-g/g3/fam_11000_30000.json.  Prints only.
usage: python3 -I g3_family.py   (about 1 minute)"""
import sys, os, json, numpy as np, collections, math
HERE = os.path.dirname(os.path.abspath(__file__)); REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
cls = collections.defaultdict(set)
for line in open(os.path.join(REPO, 'scripts', 'g1_full', 'modp_classes.txt')):
    p, r = line.split()[:2]; cls[int(p)].add(int(r))
fam = json.load(open(os.path.join(HERE, '..', 'roadmap-g', 'g3', 'fam_11000_30000.json')))['open']
mine = []
for m in range(11001, 13001):
    n = np.arange(1, m, dtype=np.int64); s = m*m + n*n
    r = np.floor(np.sqrt(s.astype(np.float64))).astype(np.int64)
    sq = np.zeros(len(n), bool)
    for d in (-1, 0, 1): sq |= (r + d)**2 == s
    h = s // 2; r2 = np.floor(np.sqrt(h.astype(np.float64))).astype(np.int64); tw = np.zeros(len(n), bool)
    for d in (-1, 0, 1): tw |= ((r2 + d)**2 == h) & (s % 2 == 0)
    for nn in n[sq | tw]:
        nn = int(nn)
        if math.gcd(m, nn) != 1: continue
        if any(nn % p and (m * pow(nn, -1, p)) % p in rs for p, rs in cls.items()): continue
        mine.append((m, nn))
ref = sorted(tuple(p) for p in fam if 11000 < p[0] <= 13000)
print('brute force: %d; fam_11000_30000.json in 11000 < m <= 13000: %d; identical: %s' % (len(mine), len(ref), sorted(mine) == ref))
