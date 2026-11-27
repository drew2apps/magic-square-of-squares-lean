"""H2: the ratios that a time limit set aside, rerun with no time limit (EXPERIMENT, written by an AI model; not read
by a person; nothing here changes the census).
Reads every h2_hunt_k*.txt, h2_design_k*.txt and h2_ex_*.txt line of a ratio that keeps a class and whose descent timed out, and runs h2_hunt.need
(the four-class descent: every bad prime exact, then good primes by bisection to 1500, then 5000) with no limit.
Output h2_retry.txt: one line per ratio, m,n log2|space| P (P = 0 for bad primes alone; SURVIVOR if not emptied).
usage: h2_retry.py [--write]"""
import sys, os, re, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
WRITE = '--write' in sys.argv; argv = sys.argv; sys.argv = [argv[0], '13']
import h2_hunt as HH
sys.argv = argv
todo = []
for f in sorted(glob.glob(os.path.join(HERE, 'h2_hunt_k*.txt')) + glob.glob(os.path.join(HERE, 'h2_design_k*.txt')) + glob.glob(os.path.join(HERE, 'h2_ex_*.txt'))):
    for l in open(f):
        mm = re.match(r'(?:S )?(\d+),(\d+) .*class dim .*timeout descent', l.strip())
        if mm and (int(mm.group(1)), int(mm.group(2))) not in todo: todo.append((int(mm.group(1)), int(mm.group(2))))
out = ['H2: the ratios set aside by a time limit, rerun with no limit. Written by an AI model; not read by a person. Experiment, not counted.',
       'ratios set aside: %d' % len(todo)]
for m, n in todo:
    ls, P = HH.need(m, n)
    out.append('%d,%d %d %s' % (m, n, round(ls), 'SURVIVOR' if P is None else P))
s = '\n'.join(out) + '\n'; f = os.path.join(HERE, 'h2_retry.txt')
if WRITE: open(f, 'w').write(s)
else:
    sys.stdout.write(s)
    if os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS from h2_retry.txt')
