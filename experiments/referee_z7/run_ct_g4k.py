#!/usr/bin/env python3
"""Review z7 of roadmap H, H7 claim 1: rank of F4's quotient C: y^2 = g4 K on the 2-isogenous curve E', independently.
Written by an AI model (review agent z7, 27 November 2026); not read by a person.  Experiment, not counted.

For each ratio: C is dm_g4k.g4k_curve(m, n) (experiments/census-wide/g4k, imported, unchanged); its rational root
e = -c3 = -2 S^3 is checked; the shift X = x - e gives y^2 = X(X^2 + aX + b) with a = -32 m^2 n^2 S, b = 64 m^2 n^2 S^4
(checked against g4k_curve's a-invariants); E' = Velu: Y^2 = X(X^2 - 2aX + a^2 - 4b) = X(X - 16Smn(m-n)^2)(X + 16Smn(m+n)^2)
(velu_check.py proves the formula as a polynomial identity).  Then my own Sel^2 (my2desc.py) and Cassels-Tate pairing
(ct_mine.py) on E', and the kernel of the pairing against the images of T1, T2 and of the points of E' in
z7_points.json (found by PARI's ellrank with --find; each checked on E' exactly here).  Each line is compared with
H7's committed outputs (h7_g4k_rank1.txt, h7_g4k_rank1_negative.txt, h7_b13000_rank1.txt): same roots of E', same
dim Sel^2, same pairing rank.
Rests on: as ct_mine.py (quoted Cassels-Tate theory), Velu's formula and isogeny invariance of the rank (quoted).
PARI: inside ct_mine.py (see there); with --find, ellrank to find points (a finder only).

    python3 run_ct_g4k.py            prints run_ct_g4k.txt
    python3 run_ct_g4k.py --find     rewrite the E' part of z7_points.json (PARI ellrank), then exit
"""
import sys, os, io, re, json, contextlib
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
REPO = os.environ.get('MSS_REPO', os.path.join(HERE, '..', '..'))
H7 = os.environ.get('H7_DIR', os.path.join(HERE, '..', 'roadmap-h', 'h7'))
sys.path.insert(0, os.path.join(REPO, 'experiments', 'census-wide', 'g4k'))
with contextlib.redirect_stdout(io.StringIO()):
    import dm_g4k
import ct_mine as C
PTS = os.path.join(HERE, 'z7_points.json')
RATIOS = ['12201,7360', '7735,4968', '6881,6233', '13348,3795', '6557,324', '11165,7548', '12204,9853', '13899,1180']

def curves(m, n):
    S = m * m + n * n
    _, a2, _, a4, a6 = [int(v) for v in dm_g4k.g4k_curve(m, n)[0]]
    e = -2 * S ** 3
    assert e ** 3 + a2 * e * e + a4 * e + a6 == 0, 'e is not a root of C'
    a, b = 3 * e + a2, 3 * e * e + 2 * a2 * e + a4
    assert a == -32 * m * m * n * n * S and b == 64 * m * m * n * n * S ** 4
    es = tuple(sorted((0, 16 * S * m * n * (m - n) ** 2, -16 * S * m * n * (m + n) ** 2)))
    s = 8 * m * n * S * S; assert s * s == b
    assert sorted((0, a + 2 * s, a - 2 * s)) == list(es)     # the roots of X^2 - 2aX + a^2 - 4b are a +- 2 sqrt(b)
    return a, b, es

def h7_records():
    rec = {}
    for f in ('h7_g4k_rank1.txt', 'h7_g4k_rank1_negative.txt', 'h7_b13000_rank1.txt'):
        cur = None
        for l in open(os.path.join(H7, f)):
            mm = re.match(r'\((\d+),(\d+)\) g4\*K', l)
            if mm: cur = (int(mm.group(1)), int(mm.group(2))); continue
            mm = re.search(r'roots \(([-\d, ]+)\)', l)
            if mm and cur: rec.setdefault(cur, {})['roots'] = tuple(int(v) for v in mm.group(1).split(','))
            mm = re.search(r'dim Sel\^2 = (\d+), Cassels-Tate pairing of rank (\d+)', l)
            if mm and cur: rec.setdefault(cur, {})['sel'] = (int(mm.group(1)), int(mm.group(2)))
    return rec

def find():
    import cypari2
    pari = cypari2.Pari(); pari.setrand(1)
    data = json.load(open(PTS)) if os.path.exists(PTS) else {}
    out = {}
    for arg in RATIOS:
        m, n = (int(v) for v in arg.split(','))
        a, b, es = curves(m, n)
        E = pari.ellinit([0, -2 * a, 0, a * a - 4 * b, 0])
        pts = []
        for seed in range(1, 6):                     # keep the largest list of independent points PARI finds
            pari.setrand(seed)
            for eff in (1, 2, 3):
                r = pari.ellrank(E, eff)
                if len(r[3]) > len(pts): pts = [[str(P[0]), str(P[1])] for P in r[3]]
                if len(pts) >= int(r[1]): break
            if pts and len(pts) >= int(r[1]): break
        out[arg] = pts
    data['g4k_Eprime'] = out
    json.dump(data, open(PTS, 'w'), indent=1, sort_keys=True)
    print('wrote', PTS)

def main():
    print('# run_ct_g4k.py (review z7; written by an AI model): F4\'s g4*K curve C, its 2-isogenous curve E\' with full 2-torsion,')
    print('# my complete 2-descent and Cassels-Tate pairing on E\', compared with H7\'s committed outputs')
    rec = h7_records(); pts = json.load(open(PTS))['g4k_Eprime']
    for arg in RATIOS:
        m, n = (int(v) for v in arg.split(','))
        a, b, es = curves(m, n)
        xs = []
        for x, y in pts.get(arg, []):
            x, y = Fr(x), Fr(y)
            assert y * y == x * (x * x - 2 * a * x + a * a - 4 * b), 'a stored point is not on E\''
            xs.append(x)
        print('(%d,%d): C = dm_g4k.g4k_curve, root e = -2S^3, a = %d, b = %d; E\' roots %s' % (m, n, a, b, es))
        res = C.analyse(es, seed=3, xs=xs, label='   mine:')
        r = rec.get((m, n))
        if r is None: print('   H7: no record'); continue
        same = r['roots'] == es and r['sel'] == (res['dim'], res['rho'])
        print('   H7 (%s): roots and dim Sel^2 = %d, pairing rank %d: %s' % ('committed', r['sel'][0], r['sel'][1], 'agree' if same else 'DIFFER'))
        assert same

if __name__ == '__main__':
    find() if '--find' in sys.argv else main()
