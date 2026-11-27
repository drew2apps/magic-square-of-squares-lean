"""H1: one phantom ratio of the random scan (h1_scan_k28.txt), on the pair curves (EXPERIMENT, written by an AI
model; not read by a person; nothing here changes the census).

(349960367,118678903) keeps, besides the tuples of its degenerate points t = 0 and t = infinity,
    (1, mn, 2, 1) and (2mn, 2, 2, 1),
the known ones times (m, m, 1, 1): d0 and d8 change by the class m, dL and dU do not.  So D_08 is unchanged and
D_0L, D_0U, D_8L, D_8U all change by m.  Each of those four pair curves is run through h1_pairs.analyse (complete
2-descent, Cassels-Tate pairing, the root of the descent map), exactly as for (1109503,597671); see h1_pairs.py for
what is computed and what it rests on.  A tuple is killed when one pair curve has its rank bounds equal and the
tuple's class is not the class of a rational point.
usage: h1_phantom.py [--write]   (h1_phantom.txt; PARI points kept as hints in h1_phantom_hints.json)
The repository is found through MSS_REPO, or else ../../.."""
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv
sys.argv = [sys.argv[0]]
import h1_pairs as HP
HF = os.path.join(HERE, 'h1_phantom_hints.json')
if os.path.exists(HF):
    H = json.load(open(HF)); HP.HINTS.setdefault('points', {}).update(H.get('points', {}))
    HP.RS.HINTS['conic'].update(H.get('rs_conic', {}))
M, N = 349960367, 118678903
KNOWN = [(237357806, 699920734, 2, 1), (349960367, 118678903, 2, 1)]
EXTRA = [(1, 41532912449037401, 2, 1), (83065824898074802, 2, 2, 1)]
PAIRS = [('0', 'L'), ('0', 'U'), ('8', 'L'), ('8', 'U')]

def main():
    say = HP.say
    say('H1: the extra tuples of the scan phantom (%d,%d). Written by an AI model; not read by a person. Experiment, not counted.' % (M, N))
    # the four pair curves have the same integral model here (they are isomorphic through t -> -t and t -> 1/t), so
    # the points PARI finds on each are points of all four: they are pooled (each checked on the curve in analyse)
    prs = {jk: HP.Pair(M, N, *jk) for jk in PAIRS}
    es = set(p.es for p in prs.values())
    if len(es) == 1:
        pool = []
        for jk, p in prs.items():
            for P in HP.find_points(p):
                if P not in pool: pool.append(P)
        for jk in PAIRS: HP.HINTS['points']['%d,%d %s%s' % (M, N, jk[0], jk[1])] = [[str(x), str(y)] for x, y in pool]
        say('the four pair curves have the same integral model; %d points pooled from PARI on the four' % len(pool))
    allres = {}
    for j, k in PAIRS:
        try:
            res, _, _ = HP.analyse(M, N, j, k, KNOWN + EXTRA)
        except AssertionError as ex:
            say('  E_%s%s: stopped (%s)' % (j, k, str(ex).split(':')[0])); continue
        for tp, v in res.items(): allres.setdefault(tp, []).append((j + k, v))
    say('')
    say('VERDICT per tuple:')
    for tp in KNOWN + EXTRA:
        ex = [jk for jk, v in allres.get(tp, []) if v.startswith('EXCLUDED')]
        kind = 'known (degenerate point)' if tp in KNOWN else 'extra'
        say('  %s %s: %s' % (tp, kind, ('killed by E_%s' % ' and E_'.join(ex)) if ex else 'not excluded by the four pair curves'))
        if tp in KNOWN: assert not ex, 'a tuple of a rational point is excluded: the method is wrong'

if __name__ == '__main__':
    main()
    s = '\n'.join(HP.OUT) + '\n'; f = os.path.join(HERE, 'h1_phantom.txt')
    if WRITE:
        open(f, 'w').write(s)
        json.dump({'points': {k: v for k, v in HP.HINTS.get('points', {}).items() if k.startswith('%d,%d' % (M, N))},
                   'rs_conic': HP.RS.HINTS['conic']}, open(HF, 'w'), indent=0, sort_keys=True)
    else:
        sys.stdout.write(s)
        if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from ' + os.path.basename(f))
