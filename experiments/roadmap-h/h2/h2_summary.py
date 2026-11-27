"""H2: summary of the hunt (h2_hunt_k*.txt), the excess-guided design (h2_design_k*.txt) and the exhaustive run (h2_ex_15000_*of3.txt) (EXPERIMENT, written by an
AI model; not read by a person).  Reads the committed outputs only (and .part files with --partial).
  * hunt: per range, ratios sampled and kept, survivors; the largest good-prime cutoff P needed, against log2 m and
    against log2 of the tuple space;
  * exhaustive: per block of 1000 in m, the counts and the cutoffs.
usage: h2_summary.py [--partial] [--write]   (h2_summary.txt)"""
import sys, os, glob, re, collections
from math import log2
HERE = os.path.dirname(os.path.abspath(__file__))
PART = '--partial' in sys.argv
def files(pat):
    fs = sorted(glob.glob(os.path.join(HERE, pat)))
    if PART: fs += sorted(glob.glob(os.path.join(HERE, pat + '.part')))
    return fs
out = ['H2 summary. Written by an AI model; not read by a person. Experiment, not counted.%s' % (' (PARTIAL: includes unfinished .part files)' if PART else '')]
# hunt
rows = []; per = []
for f in files('h2_hunt_k*.txt'):
    k = int(re.search(r'k(\d+)', os.path.basename(f)).group(1)); samp = kept = 0; surv = []
    rd = re.search(r'_r(\d+)', os.path.basename(f)); rd = int(rd.group(1)) if rd else 0
    for l in open(f):
        l = l.rstrip('\n')
        if l.startswith('S '): l = l[2:]
        if l.startswith('sampled'): samp = int(l.split()[1])
        mm = re.match(r'(\d+),(\d+) class dim (\d+) (\d+) (\S+)$', l)
        if mm:
            m, n, d, ls, P = mm.groups(); rows.append((k, int(m), int(ls), P, rd)); kept += 1
            if P == 'SURVIVOR': surv.append((int(m), int(n)))
        elif l.startswith(('S ',)) or re.match(r'\d+,\d+ ', l): samp = max(samp, 0)
    if f.endswith('.part'): samp = sum(1 for l in open(f) if l.startswith('S '))
    per.append((k, rd, samp, kept, surv, f.endswith('.part')))
out.append('HUNT (biased sampler, m from 10^4 to 10^12): ranges %d, ratios sampled %d, keeping a class %d' % (len(per), sum(p[2] for p in per), sum(p[3] for p in per)))
out.append('SURVIVORS: %s' % ([s for p in per for s in p[4]] or 'none'))
out.append('k round  sampled  kept  [(part) = unfinished]  largest P  P values (count)')
for k, rd, samp, kept, surv, part in sorted(per):
    Ps = [r[3] for r in rows if r[0] == k and r[4] == rd]
    num = [int(P) for P in Ps if P.isdigit()]
    c = collections.Counter(Ps)
    out.append('%2d %d %7d %5d %s  %s  %s' % (k, rd, samp, kept, '(part)' if part else '', max(num) if num else '-', ' '.join('%s:%d' % (P, c[P]) for P in sorted(c, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else 0)))))
out.append('largest P needed by log2 of the tuple space (hunt):')
by = collections.defaultdict(list)
for k, m, ls, P, _rd in rows:
    if P.isdigit(): by[ls].append(int(P))
for ls in sorted(by): out.append('  %2d: ratios %d, need a good prime %d, largest P %d' % (ls, len(by[ls]), sum(1 for P in by[ls] if P), max(by[ls])))
# design (excess-guided)
dsg = []; dsurv = []
for f in files('h2_design_k*.txt'):
    k = int(re.search(r'k(\d+)', os.path.basename(f)).group(1)); tk = kp = 0; Ps = []
    for l in open(f):
        mm = re.match(r'(\d+),(\d+) \S+ class dim \d+ \d+ (\S+)$', l.strip())
        if mm:
            kp += 1
            if mm.group(3) == 'SURVIVOR': dsurv.append((int(mm.group(1)), int(mm.group(2))))
            elif mm.group(3).isdigit(): Ps.append(int(mm.group(3)))
        if re.match(r'\d+,\d+ ', l): tk += 1
    dsg.append((k, tk, kp, max(Ps) if Ps else '-', f.endswith('.part')))
out.append('DESIGN (best of 200 random ratios by the excess, per range): ranges %d, ratios taken %d, keeping a class %d' % (len(dsg), sum(d[1] for d in dsg), sum(d[2] for d in dsg)))
out.append('SURVIVORS: %s' % (dsurv or 'none'))
for k, tk, kp, P, part in sorted(dsg): out.append('  k=%2d taken %3d kept %3d largest P %s%s' % (k, tk, kp, P, ' (part)' if part else ''))
# exhaustive
blk = collections.defaultdict(collections.Counter); bP = collections.defaultdict(list); exsurv = []
for f in files('h2_ex_15000_*of3.txt'):
    for l in open(f):
        if l.startswith('m='):
            ws = l.split(); m = int(ws[0][2:]); b = (m - 1) // 1000 * 1000
            for w in ws[1:]: kk, v = w.split('='); blk[b][kk] += int(v)
            blk[b]['m values'] += 1
        else:
            mm = re.match(r'(\d+),(\d+) (\d+) (\d+) (\S+)$', l.strip())
            if mm:
                m = int(mm.group(1)); b = (m - 1) // 1000 * 1000; P = mm.group(5)
                if P == 'SURVIVOR': exsurv.append((m, int(mm.group(2))))
                else: bP[b].append(int(P))
out.append('EXHAUSTIVE 15000 < m <= 30000 (every coprime ratio not of the degenerate types): blocks %d' % len(blk))
out.append('SURVIVORS: %s' % (exsurv or 'none'))
out.append('block       m done  coprime  degenerate  no point  no class  class  need a good prime  largest P')
for b in sorted(blk):
    c = blk[b]; Ps = bP[b]
    out.append('%5d-%5d %4d %9d %6d %10d %8d %6d %6d %6s' % (b, b + 1000, c['m values'], c['coprime'], c['degenerate'], c['no_point'], c['no_class'], c['class'], sum(1 for P in Ps if P), max(Ps) if Ps else '-'))
s = '\n'.join(out) + '\n'; f = os.path.join(HERE, 'h2_summary.txt')
if '--write' in sys.argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f) and not PART: print('# matches' if open(f).read() == s else '# DIFFERS from h2_summary.txt')
