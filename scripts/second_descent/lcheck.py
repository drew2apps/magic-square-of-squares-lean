"""Second check of the rank-0 curves used by second3.py: L(E,1) is not 0 (then the rank is 0, by Kolyvagin and modularity).
Independent of the 2-descent in ranks2.py.  usage: lcheck.py LIMIT second_123.jsonl   (curves of conductor up to LIMIT; time grows like sqrt(N))"""
import json, sys, time
import cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
rows = [json.loads(l) for l in open(sys.argv[2])]
out = []
for r in rows:
    if r['status'] != 'closed': continue
    c = [c for c in r['curves'] if c['curve'] == r['by']][0]; E = pari.ellinit([pari(a) for a in c['rank0_curve']]); N = pari.ellglobalred(E)[0]
    out.append((int(N), r['pair'], c['rank0_curve']))
out.sort()
print('conductors: min %.3g median %.3g max %.3g' % (out[0][0], out[len(out)//2][0], out[-1][0]), flush=True)
lim = float(sys.argv[1]); done = 0
for N, pair, ai in out:
    if N > lim: break
    t0 = time.time(); E = pari.ellinit([pari(a) for a in ai]); L = float(pari.lfun(E, 1)); om = float(E.omega()[0])
    done += 1; print(pair[0], pair[1], 'N = %.3g' % N, 'L(E,1) = %.6f' % L, 'L/Omega = %.6f' % (L/om), flush=True)
print('checked', done, 'of', len(out))
