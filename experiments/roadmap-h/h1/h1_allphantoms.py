"""H1, part 4: every phantom of the scan (h1_scan_k*.txt) on its pair curves (EXPERIMENT, written by an AI model; not
read by a person; nothing here changes the census).
For each phantom ratio and each pair (j, k) of the four forms F_0, F_8, F_L, F_U, the class of d_j d_k is compared with
those of the known tuples; a pair curve E_jk can only exclude an extra tuple whose class differs from every known one
(a known tuple is the class of a rational point).  Each such E_jk goes through h1_pairs.analyse: complete 2-descent,
Cassels-Tate pairing (exact), points from PARI as a finder (pooled over pair curves with the same integral model,
each checked on the curve), and a tuple is EXCLUDED when either the rank bounds meet and its class is not that of a
rational point, or its class is not that of any element of the kernel of the pairing (see h1_pairs.py).  A pair curve
that does not have full rational 2-torsion, or where the root of the descent map is not identified, is reported as
stopped and decides nothing.
Output h1_allphantoms.txt: per phantom, per useful pair curve, the analysis; then a verdict per extra tuple.
usage: h1_allphantoms.py [--write]   (PARI points kept as hints in h1_allphantoms_hints.json)
The repository is found through MSS_REPO, or else ../../.."""
import sys, os, re, ast, glob, json
from fractions import Fraction as Fr
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv
_argv = sys.argv; sys.argv = [sys.argv[0]]
import h1_pairs as HP
sys.argv = _argv
HF = os.path.join(HERE, 'h1_allphantoms_hints.json')
if os.path.exists(HF):
    H = json.load(open(HF)); HP.HINTS.setdefault('points', {}).update(H.get('points', {}))
    HP.RS.HINTS['conic'].update(H.get('rs_conic', {}))
PAIRS = [('0', '8'), ('0', 'L'), ('0', 'U'), ('8', 'L'), ('8', 'U'), ('L', 'U')]

def phantoms():
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, 'h1_scan_k*.txt'))):
        for l in open(f):
            mm = re.match(r'\((\d+),(\d+)\) PHANTOM: class dim (\d+).*known (\[.*\]); extra (\[.*\])', l.strip())
            if mm: out.append((int(mm.group(1)), int(mm.group(2)), ast.literal_eval(mm.group(4)), ast.literal_eval(mm.group(5))))
    return sorted(out)

def main():
    say = HP.say
    say('H1: every phantom of the scan on its pair curves. Written by an AI model; not read by a person. Experiment, not counted.')
    verdicts = []
    for m, n, known, extra in phantoms():
        S = m * m + n * n
        say(''); say('(%d,%d), %s type: known %s; extra %s' % (m, n, 'square' if isqrt(S) ** 2 == S else 'twice-square', known, extra))
        HP.set_base(set(q for x in (2, m, n, m - n, m + n, S) for q in HP.RI.factor_int(x)) | {2})
        useful = []
        for j, k in PAIRS:
            kc = set(HP.sqcl(Fr(t[HP.IDX[j]] * t[HP.IDX[k]])) for t in known)
            if any(HP.sqcl(Fr(t[HP.IDX[j]] * t[HP.IDX[k]])) not in kc for t in extra): useful.append((j, k))
        say('  pair curves on which some extra tuple has a class no known tuple has: %s' % ' '.join('E_' + j + k for j, k in useful))
        # pool PARI's points over the useful pair curves that share an integral model
        prs = {}
        for jk in useful:
            try: prs[jk] = HP.Pair(m, n, *jk)
            except AssertionError as ex: say('  E_%s%s: stopped (%s)' % (jk[0], jk[1], str(ex).split(':')[0]))
        groups = {}
        for jk, p in prs.items(): groups.setdefault((p.es, p.t0), []).append(jk)
        for key, jks in groups.items():
            if len(jks) < 2: continue
            pool = []
            for jk in jks:
                for P in HP.find_points(prs[jk]):
                    if P not in pool: pool.append(P)
            for jk in jks: HP.HINTS['points']['%d,%d %s%s' % (m, n, jk[0], jk[1])] = [[str(x), str(y)] for x, y in pool]
            say('  %s have the same integral model and base point; %d points pooled' % (' '.join('E_' + a + b for a, b in jks), len(pool)))
        res = {}
        for jk in prs:
            try:
                r, _, _ = HP.analyse(m, n, jk[0], jk[1], known + extra)
            except AssertionError as ex:
                say('  E_%s%s: stopped (%s)' % (jk[0], jk[1], str(ex).split(':')[0])); continue
            for tp, v in r.items(): res.setdefault(tp, []).append((jk[0] + jk[1], v))
        for tp in known:
            assert not any(v.startswith('EXCLUDED') for _, v in res.get(tp, [])), 'a tuple of a rational point is excluded: the method is wrong'
        for tp in extra:
            ex = [jk for jk, v in res.get(tp, []) if v.startswith('EXCLUDED')]
            verdicts.append(((m, n), tp, ex))
    say(''); say('VERDICT per extra tuple:')
    for mn, tp, ex in verdicts:
        say('  %s %s: %s' % (mn, tp, ('killed by E_%s' % ' and E_'.join(ex)) if ex else 'not excluded by its pair curves'))
    kr = sorted(set(mn for mn, _, _ in verdicts)); kd = [mn for mn in kr if all(ex for m2, _, ex in verdicts if m2 == mn)]
    say('phantom ratios: %d; all extra tuples killed by pair curves: %d %s' % (len(kr), len(kd), kd))

if __name__ == '__main__':
    main()
    s = '\n'.join(HP.OUT) + '\n'; f = os.path.join(HERE, 'h1_allphantoms.txt')
    if WRITE:
        open(f, 'w').write(s)
        json.dump({'points': HP.HINTS.get('points', {}), 'rs_conic': HP.RS.HINTS['conic']}, open(HF, 'w'), indent=0, sort_keys=True)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
