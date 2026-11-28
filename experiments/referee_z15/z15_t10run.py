"""z15 (referee of EXPLORATION_4): run the repository's t10 pairing test (t10_pairct.py of experiments/conductor-t10,
copied here byte for byte as t10_pairct_copy.py, sha256 96c9d7c5...) on given ratios and write one JSON line per ratio.
Written by an AI model (reviewer z15), not read by a person.  Tier D experiment.

usage: python3 z15_t10run.py OUTFILE.jsonl m,n [m,n ...]
       python3 z15_t10run.py OUTFILE.jsonl @RATIOFILE        (one 'm,n' per line)
No hints are read (the copy lives in this directory, which has no t10_hints*.json), so PARI finds every point anew.
Each line: m, n, residual tuples with their least killing good prime (<= 5000), the t10 verdict per pair curve
('killed' | 'outside Sel2' | 'image' | 'kernel' | 'undecided' | 'no base'), per curve dim Sel^2, rank of the point
images and PARI's ellrank bounds (parsed from t10's own output line), and the CPU time.  Lines already in OUTFILE are
skipped (resumable)."""
import os, sys, json, time, re
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('MSS_REPO', os.path.abspath(os.path.join(HERE, '..', '..')))
outf = sys.argv[1]
args = sys.argv[2:]
sys.argv = [sys.argv[0]]          # t10_pairct reads flags from sys.argv at import: give it none
import io, contextlib
sys.path.insert(0, HERE)
import t10_pairct_copy as T

ratios = []
for a in args:
    if a.startswith('@'):
        ratios += [tuple(map(int, l.split()[0].split(','))) for l in open(a[1:]) if l.strip() and l[0] != '#']
    else:
        ratios.append(tuple(map(int, a.split(','))))
done = set()
if os.path.exists(outf):
    for l in open(outf):
        r = json.loads(l); done.add((r['m'], r['n']))
pat = re.compile(r'  E_(\w\w): .*dim Sel2 = (\d+), point images (\d+); PARI ellrank \[(-?\d+),(-?\d+)\]')
with open(outf, 'a') as fo:
    for m, n in ratios:
        if (m, n) in done: continue
        T.OUT.clear()
        t0 = time.process_time(); w0 = time.time()
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                left, kill, res, known = T.analyse_ratio(m, n)
            curves = {}
            for l in T.OUT:
                mm = pat.match(l)
                if mm: curves[mm.group(1)] = dict(sel=int(mm.group(2)), img=int(mm.group(3)), pari=[int(mm.group(4)), int(mm.group(5))])
            rec = dict(m=m, n=n, known=[list(t) for t in known],
                       tuples=[[list(t), kill.get(t), {jk: res[t][jk] for jk in res[t]}] for t in sorted(left)],
                       curves=curves, cpu=round(time.process_time() - t0, 1), wall=round(time.time() - w0, 1))
        except Exception as ex:
            rec = dict(m=m, n=n, error=repr(ex)[:300], cpu=round(time.process_time() - t0, 1))
        fo.write(json.dumps(rec, sort_keys=True) + '\n'); fo.flush()
        sys.stderr.write('%d,%d done %.1f s cpu\n' % (m, n, rec['cpu']))
