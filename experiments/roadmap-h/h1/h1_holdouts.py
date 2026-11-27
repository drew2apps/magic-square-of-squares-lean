"""H1, part 8: the pair curves on the partners of the phantoms that parts 4 and 7 leave open (EXPERIMENT, written by an
AI model; not read by a person; nothing here changes the census).
Parts 4 (h1_allphantoms.txt, the pair curves) and 7 (h1_rank0_exact.txt, exact rank-0 quotients) leave some partner
pairs with an extra tuple not excluded on either side.  Part 4 ran only on the scan's own ratios; this runs the same
analysis (h1_allphantoms: h1_pairs.analyse on every useful pair curve, PARI's points pooled and each checked) on the
partner of each such ratio, with the known and extra tuples that h1_partners.py lists for it.
The open pairs are read from h1_rank0_exact.txt ('not excluded') and dropped when h1_allphantoms.txt kills every extra
tuple of one side.
Output h1_holdouts.txt.   usage: h1_holdouts.py [--write]   (PARI points kept as hints in h1_holdouts_hints.json)
The repository is found through MSS_REPO, or else ../../.."""
import sys, os, re, ast, json
from fractions import Fraction as Fr
from math import isqrt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv
_argv = sys.argv; sys.argv = [sys.argv[0]]
import h1_allphantoms as HA
sys.argv = _argv
HP = HA.HP
HF = os.path.join(HERE, 'h1_holdouts_hints.json')
if os.path.exists(HF):
    H = json.load(open(HF)); HP.HINTS.setdefault('points', {}).update(H.get('points', {}))
    HP.RS.HINTS['conic'].update(H.get('rs_conic', {}))

def open_pairs():
    killed = set()
    for l in open(os.path.join(HERE, 'h1_allphantoms.txt')):
        mm = re.match(r'phantom ratios: \d+; all extra tuples killed by pair curves: \d+ (\[.*\])', l.strip())
        if mm: killed = set(ast.literal_eval(mm.group(1)))
    out = []
    for l in open(os.path.join(HERE, 'h1_rank0_exact.txt')):
        mm = re.match(r'\((\d+), (\d+)\) and \((\d+), (\d+)\): not excluded', l.strip())
        if mm:
            a = (int(mm.group(1)), int(mm.group(2))); b = (int(mm.group(3)), int(mm.group(4)))
            if a not in killed and b not in killed: out.append((a, b))
    return out

def partners():
    d = {}
    for l in open(os.path.join(HERE, 'h1_partners.txt')):
        mm = re.match(r'\((\d+),(\d+)\) -> partner \((\d+),(\d+)\): PHANTOM;.*known (\[.*\]); extra (\[.*\])', l.strip())
        if mm: d[(int(mm.group(3)), int(mm.group(4)))] = (ast.literal_eval(mm.group(5)), ast.literal_eval(mm.group(6)))
    return d

def main():
    say = HP.say
    say('H1: pair curves on the partners of the open phantom pairs. Written by an AI model; not read by a person. Experiment, not counted.')
    pr = partners(); verdicts = []
    pairs = open_pairs(); say('open pairs after parts 4 and 7: %s' % pairs)
    for a, b in pairs:
        for m, n in (a, b):
            if (m, n) not in pr: continue          # a scan ratio: done in part 4
            known, extra = pr[(m, n)]
            S = m * m + n * n
            say(''); say('(%d,%d), %s type: known %s; extra %s' % (m, n, 'square' if isqrt(S) ** 2 == S else 'twice-square', known, extra))
            HP.set_base(set(q for x in (2, m, n, m - n, m + n, S) for q in HP.RI.factor_int(x)) | {2})
            useful = []
            for j, k in HA.PAIRS:
                kc = set(HP.sqcl(Fr(t[HP.IDX[j]] * t[HP.IDX[k]])) for t in known)
                if any(HP.sqcl(Fr(t[HP.IDX[j]] * t[HP.IDX[k]])) not in kc for t in extra): useful.append((j, k))
            say('  pair curves on which some extra tuple has a class no known tuple has: %s' % ' '.join('E_' + j + k for j, k in useful))
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
                say('  %s have the same integral model and base point; %d points pooled' % (' '.join('E_' + x + y for x, y in jks), len(pool)))
            res = {}
            for jk in prs:
                try: r, _, _ = HP.analyse(m, n, jk[0], jk[1], known + extra)
                except AssertionError as ex:
                    say('  E_%s%s: stopped (%s)' % (jk[0], jk[1], str(ex).split(':')[0])); continue
                for tp, v in r.items(): res.setdefault(tp, []).append((jk[0] + jk[1], v))
            for tp in known:
                assert not any(v.startswith('EXCLUDED') for _, v in res.get(tp, [])), 'a tuple of a rational point is excluded: the method is wrong'
            for tp in extra:
                verdicts.append(((m, n), tp, [jk for jk, v in res.get(tp, []) if v.startswith('EXCLUDED')]))
    say(''); say('VERDICT per extra tuple:')
    for mn, tp, ex in verdicts:
        say('  %s %s: %s' % (mn, tp, ('killed by E_%s' % ' and E_'.join(ex)) if ex else 'not excluded by its pair curves'))
    kr = sorted(set(mn for mn, _, _ in verdicts)); kd = [mn for mn in kr if all(ex for m2, _, ex in verdicts if m2 == mn)]
    say('partner ratios run: %d; all extra tuples killed by pair curves: %d %s' % (len(kr), len(kd), kd))

if __name__ == '__main__':
    main()
    s = '\n'.join(HP.OUT) + '\n'; f = os.path.join(HERE, 'h1_holdouts.txt')
    if WRITE:
        open(f, 'w').write(s)
        json.dump({'points': HP.HINTS.get('points', {}), 'rs_conic': HP.RS.HINTS['conic']}, open(HF, 'w'), indent=0, sort_keys=True)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
