"""H2, control: h2_exhaust.py against G2 on 14900 < m <= 15000 (EXPERIMENT, written by an AI model; not read by a
person; nothing here changes the census).
G2 (roadmap-g/g2) lists every class-keeping ratio with m <= 15000 (class_selmer_15000_list.txt) and, for each, whether
the bad primes alone empty its four-class set (g2_15000_all.txt).  h2_exhaust.py, run on 14900 < m <= 15000
(h2_ex_14900_0of1.txt), must keep exactly the same ratios, and its P must be 0 exactly where G2 says
"bad-alone: EMPTY"; where a good prime is needed, P must not exceed the largest good prime of G2's irredundant kill set
(g2_15000_need_good.txt).  With --negative, a forged copy of the H2 output (one ratio dropped; a P changed from 0; a P raised above G2's kill set) is refused.
usage: h2_control.py [--negative] [--write]   (h2_control.txt, h2_control_negative.txt)
G2's folder: G2_DIR, or else ../../roadmap-g/g2 (or ../g2)."""
import sys, os, re
HERE = os.path.dirname(os.path.abspath(__file__))
G2 = os.environ.get('G2_DIR') or next((d for d in (os.path.join(HERE, '..', '..', 'roadmap-g', 'g2'), os.path.join(HERE, '..', 'g2')) if os.path.isdir(d)), None)
LO, HI = 14900, 15000
NEG = '--negative' in sys.argv

def g2_sets():
    keep = set(); bad = {}
    lst = os.path.join(G2, 'class_selmer_15000_list.txt')
    if not os.path.exists(lst):  # in the repository the list is kept once, with the class criterion's outputs
        lst = os.path.join(HERE, '..', '..', '..', 'scripts', 'cover_descent', 'structure', 'class_selmer_15000_list.txt')
    for l in open(lst):
        w = l.split()
        if w:
            m, n = map(int, w[0].split(','))
            if LO < m <= HI: keep.add((m, n))
    for l in open(os.path.join(G2, 'g2_15000_all.txt')):
        mm = re.match(r'(\d+),(\d+)\s+(\d+)\s+bad-alone: (\S+)', l)
        if mm:
            m, n = int(mm.group(1)), int(mm.group(2))
            if LO < m <= HI: bad[(m, n)] = (int(mm.group(3)), mm.group(4) == 'EMPTY')
    kill = {}
    for l in open(os.path.join(G2, 'g2_15000_need_good.txt')):
        mm = re.match(r'(\d+),(\d+)\s.*kill: (.*)$', l.strip())
        if mm and LO < int(mm.group(1)) <= HI:
            kill[(int(mm.group(1)), int(mm.group(2)))] = max(int(w[:-2]) for w in mm.group(3).split() if w.endswith(':g'))
    return keep, bad, kill

def h2_lines():
    out = {}
    for l in open(os.path.join(HERE, 'h2_ex_14900_0of1.txt')):
        mm = re.match(r'(\d+),(\d+) (\d+) (\d+) (\S+)$', l.strip())
        if mm: out[(int(mm.group(1)), int(mm.group(2)))] = (int(mm.group(4)), mm.group(5))
    return out

def compare(h2, keep, bad, kill):
    res = []
    ok = set(h2) == keep
    res.append('class-keeping ratios: H2 %d, G2 %d, same set: %s' % (len(h2), len(keep), 'yes' if ok else 'NO'))
    if not ok:
        res.append('  only in H2: %s' % sorted(set(h2) - keep)[:5]); res.append('  only in G2: %s' % sorted(keep - set(h2))[:5])
    bad_ok = True; nz = 0
    for r, (ls, P) in sorted(h2.items()):
        if r not in bad: bad_ok = False; res.append('  %s not in g2_15000_all.txt' % (r,)); continue
        g_ls, g_empty = bad[r]
        if g_ls != ls or (P == '0') != g_empty: bad_ok = False; res.append('  %s: H2 log2 %d P %s; G2 log2 %d bad-alone %s' % (r, ls, P, g_ls, 'EMPTY' if g_empty else 'not empty'))
        if P != '0': nz += 1
    res.append('log2 of the tuple space and "P = 0 exactly when the bad primes alone empty the set": %s' % ('agree on every ratio' if bad_ok else 'DISAGREE'))
    # G2's irredundant kill set lies in the bad primes and the good primes <= q (its largest good prime), so P <= q
    kq = all(P.isdigit() and P != '0' and int(P) <= kill[r] for r, (ls, P) in h2.items() if r in kill) and all(r in kill for r, (ls, P) in h2.items() if P != '0')
    eq = sum(1 for r, (ls, P) in h2.items() if r in kill and P.isdigit() and int(P) == kill[r])
    res.append('P <= the largest good prime of G2\'s irredundant kill set, on every ratio that needs one: %s (equal on %d of %d)' % ('yes' if kq else 'NO', eq, len(kill)))
    bad_ok = bad_ok and kq
    res.append('ratios that need a good prime: %d; survivors: %d' % (nz, sum(1 for _, P in h2.values() if P == 'SURVIVOR')))
    return ok and bad_ok, res

out = ['H2 control: h2_exhaust.py against G2 on %d < m <= %d. Written by an AI model; not read by a person. Experiment, not counted.' % (LO, HI)]
keep, bad, kill = g2_sets(); h2 = h2_lines()
if NEG:
    forged = dict(h2); r0 = sorted(forged)[0]; del forged[r0]
    ok1, _ = compare(forged, keep, bad, kill)
    out.append('forged: ratio %s dropped: %s' % (r0, 'refused' if not ok1 else 'ACCEPTED'))
    forged = dict(h2); r1 = next(r for r in sorted(forged) if forged[r][1] == '0'); forged[r1] = (forged[r1][0], '29')
    ok2, _ = compare(forged, keep, bad, kill)
    out.append('forged: P of %s set to 29: %s' % (r1, 'refused' if not ok2 else 'ACCEPTED'))
    forged = dict(h2); r2 = next(r for r in sorted(forged) if r in kill); forged[r2] = (forged[r2][0], str(kill[r2] + 4))
    ok3, _ = compare(forged, keep, bad, kill)
    out.append('forged: P of %s set above G2\'s kill set: %s' % (r2, 'refused' if not ok3 else 'ACCEPTED'))
    assert not ok1 and not ok2 and not ok3
else:
    ok, res = compare(h2, keep, bad, kill); out += res
    out.append('RESULT: %s' % ('h2_exhaust agrees with G2' if ok else 'DISAGREES'))
s = '\n'.join(out) + '\n'; f = os.path.join(HERE, 'h2_control_negative.txt' if NEG else 'h2_control.txt')
if '--write' in sys.argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
