#!/usr/bin/env python3
"""The Magma transcripts against the output of the second descent.

magma_second_raw.txt holds, verbatim, the requests sent to the public Magma calculator and its answers.  For every ratio
that second_123.jsonl closes, this checks that some request contains the closing curve with exactly the coefficients of the
JSON record (as MinimalModel(EllipticCurve([...]))), and that Magma gives the upper bound 0 for its rank: RankBounds = "0 0",
or TwoPowerIsogenyDescentRankBound = 0 (that function needs a single point of order 2; RankBounds can be weaker, and for some
curves it did not finish within the 60 seconds of the public calculator).
Magma is not run here; this only reads the transcripts.  usage: magma_check.py [raw file [jsonl file]]"""
import sys, re, json

raw = open(sys.argv[1] if len(sys.argv) > 1 else 'magma_second_raw.txt').read()
recs = [json.loads(l) for l in open(sys.argv[2] if len(sys.argv) > 2 else 'second_123.jsonl')]
found = {}; versions = set()
for block in raw.split('=== ')[1:]:
    inp, resp = block.split('--- response\n', 1); inp = inp.split('--- input\n', 1)[1]
    versions.update(re.findall(r'<version>(.*?)</version>', resp))
    curves = re.findall(r'E := MinimalModel\(EllipticCurve\(\[(.*?)\]\)\);\nprint "pair (\d+),(\d+) curve (\S+) twist (-?\d+)";', inp)
    lines = [l.replace('&gt;', '>') for l in re.findall(r'<line>(.*?)</line>', resp, re.S)]
    for coeffs, m, n, name, tw in curves:
        tag = 'pair %s,%s curve %s twist %s' % (m, n, name, tw)
        if tag not in lines: continue
        i = lines.index(tag); seg = lines[i + 1:]
        nxt = [k for k, l in enumerate(seg) if l.startswith('pair ')]
        seg = seg[:nxt[0]] if nxt else seg
        rb = None; tp = None
        for k, l in enumerate(seg):
            if l == 'RankBounds:' and k + 1 < len(seg): rb = seg[k + 1].strip()
            if l.startswith('TwoPowerIsogenyDescentRankBound:'):                     # several lines when it was called with each point of order 2: the least bound
                v = l.split(':', 1)[1].split()[0]; tp = v if tp is None or (v.isdigit() and tp.isdigit() and int(v) < int(tp)) else tp
        old = found.get((int(m), int(n)))                 # a curve may have been sent twice: keep every answer
        if old: rb = rb if rb not in (None, '') else old[3]; tp = tp if tp is not None else old[4]
        found[(int(m), int(n))] = (coeffs, name, int(tw), rb, tp)
ok = 0; tp0 = 0; nrb = 0; both = 0; bad = []; weaker = []
closed = [r for r in recs if r['status'] == 'closed']
for r in closed:
    p = tuple(r['pair']); c = [c for c in r['curves'] if c['curve'] == r['by']][0]
    f = found.get(p)
    if f is None or f[0] != ', '.join(c['rank0_curve']) or f[1] != c['curve'] or f[2] != c['rank0_twist']: bad.append(p); continue
    rb0 = f[3] == '0 0'; t0 = f[4] == '0'
    if not (rb0 or t0) or (f[4] is not None and not t0): bad.append(p); continue       # one of the two must give 0, and the isogeny descent must not contradict
    ok += 1; nrb += rb0; tp0 += t0; both += rb0 and t0
    if f[3] not in ('0 0', None, ''): weaker.append((p, f[3]))
print('closing curves in second_123.jsonl: %d; in the Magma transcripts with the same coefficients and an upper bound 0: %d' % (len(closed), ok))
print('RankBounds = 0 0 for %d; TwoPowerIsogenyDescentRankBound = 0 for %d (it needs a single point of order 2); both for %d' % (nrb, tp0, both))
print('RankBounds weaker than 0 (the isogeny descent gives 0 there):', weaker)
print('Magma version in the transcripts:', ', '.join(sorted(versions)), '; missing or different:', bad)
if '--list' in sys.argv: sys.exit(0)
assert not bad and ok == len(closed)
