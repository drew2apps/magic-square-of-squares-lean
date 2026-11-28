"""I1 negative and positive controls for i1_hunt.py (EXPERIMENT, written by an AI model; not read by a person).
  (1) forged survivor: a line of i1_hunt_B_s1.txt that reads 'class dim 0 .. 0' (bad primes alone empty the set) is
      claimed as a SURVIVOR; recomputing h2_hunt.need refuses the claim (P = 0 again);
  (2) forged local table: for the degenerate ratio (7,1), a table at p = 5 with its known tuple removed is refused by
      need's own check (a known point is never rejected);
  (3) positive control: H1's phantom ratio (204724,97893), whose set is not empty because of the phantom tuples (I3),
      is not emptied by bad primes and good primes to 5000 (need returns None), so the pipeline does not empty
      every set.
usage: i1_negative.py [--write]   (i1_negative.txt)"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
argv = sys.argv; sys.argv = [argv[0]]
import i1_hunt as I
sys.argv = argv
HH = I.HH; OUT = []
def say(s=''): OUT.append(s)
say('I1 negative and positive controls. Written by an AI model; not read by a person. Experiment, not counted.')
line = [l for l in open(os.path.join(HERE, 'i1_hunt_B_s1.txt')) if ' class dim 0 ' in l and l.rstrip().endswith(' 0')][0].split()
m, n = (int(v) for v in line[0].split(','))
r = HH.need(m, n)
say('forged survivor (%d,%d): recomputed, bad primes alone empty the set (P = %s): %s' % (m, n, r[1], 'claim refused' if r[1] == 0 else 'CLAIM NOT REFUSED'))
assert r[1] == 0
# (2) a degenerate ratio with known tuples: a forged local table (its known tuple removed at p = 5) goes through need
m2, n2 = 7, 1
forms, sup, known = HH.K.TD.magic(m2, n2)
real = HH.X.allowed_at
def forged(fs, p):
    L, al = real(fs, p)
    if p == 5: al = set(al) - {tuple(L.dlabel(d) for d in known[0])}
    return L, al
HH.X.allowed_at = forged
try: HH.need(m2, n2); verdict = 'ACCEPTED'
except AssertionError as ex: verdict = 'refused (%s)' % ex.args[0][0]
HH.X.allowed_at = real
say('degenerate ratio (%d,%d) with %d known tuples, forged table at p = 5 without its known tuple: %s' % (m2, n2, len(known), verdict))
assert verdict.startswith('refused')
r = HH.need(204724, 97893)
say('positive control (204724,97893), an H1 phantom: need returns P = %s (%s)' % (r[1], 'not emptied by good primes to 5000, as the phantom tuples require' if r[1] is None else 'EMPTIED'))
assert r[1] is None
s = '\n'.join(OUT) + '\n'; f = os.path.join(HERE, 'i1_negative.txt')
if '--write' in argv: open(f, 'w').write(s)
else:
    sys.stdout.write(s); print('# matches' if os.path.exists(f) and open(f).read() == s else '# DIFFERS')
