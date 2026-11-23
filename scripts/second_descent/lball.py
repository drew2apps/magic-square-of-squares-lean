#!/usr/bin/env python3
"""L(E,1) > 0 in ball arithmetic, for curves over Q of root number +1 (a second program beside lproof.py).

For E/Q of conductor N and root number +1 (modularity):  L(E,1) = 2 sum_{n >= 1} (a_n / n) q^n,  q = exp(-2 pi / sqrt(N)).
|a_n| <= d(n) sqrt(n) <= 2 n, so the tail after M terms is at most 4 q^(M+1) / (1 - q); M is taken with this below 10^-9.
Every real number is a ball of Arb (python-flint), with 128 bits: pi, the square root, the exponential, the powers q^n by
repeated multiplication, and the sum.  Arb's balls contain the exact result of every operation, so the printed interval for
L(E,1) is an enclosure, and "positive" means that its lower end is above 0.  No assumption on rounding is made; this is what
lproof.py, which assumes a bound for PARI's rounding, could not say.  The a_n are exact integers (PARI's ellan).
L(E,1) > 0 gives rank 0 (Kolyvagin, with modularity).  What is trusted: PARI's conductor, root number and a_n, Arb, and the
two theorems.

usage: lball.py LIMIT second_123.jsonl          the closing curves of conductor up to LIMIT
       lball.py curve a1,a2,a3,a4,a6 [label]    one curve"""
import sys, json, math
import cypari2
from flint import arb, ctx
ctx.prec = 128
pari = cypari2.Pari(); pari.allocatemem(6*10**9)

def prove(ai):
    E = pari.ellminimalmodel(pari.ellinit(ai)); N = int(pari.ellglobalred(E)[0]); w = int(pari.ellrootno(E))
    if w != 1: return dict(conductor=N, root_number=w)
    x = 2*arb.pi()/arb(N).sqrt(); q = (-x).exp(); xf = 2*math.pi/math.sqrt(N)
    M = int(math.log(4e9/(1 - math.exp(-xf)))/xf) + 2
    tail = 4*q**(M + 1)/(1 - q); assert tail < arb(10)**-9
    an = pari.ellan(E, M).python(); s = arb(0); p = arb(1)
    for n in range(1, M + 1):
        p *= q; a = an[n - 1]
        if a: s += p*a/n
    lo = 2*s - tail; hi = 2*s + tail                      # balls; L(E,1) lies between the lower end of lo and the upper end of hi
    return dict(conductor=N, root_number=w, terms=M, lower=float(lo.lower()), upper=float(hi.upper()), positive=bool(lo.lower() > 0))

def show(label, o):
    if o.get('root_number') != 1: return '%s N = %.3g  root number %s: not applicable' % (label, o['conductor'], o.get('root_number'))
    return '%s N = %.3g  terms %d  L(E,1) in [%.9f, %.9f]  %s' % (label, o['conductor'], o['terms'], math.floor(o['lower']*1e9)/1e9, math.ceil(o['upper']*1e9)/1e9, 'positive' if o['positive'] else 'NOT SHOWN POSITIVE')

if __name__ == '__main__':
    if sys.argv[1] == 'curve':
        print(show(sys.argv[3] if len(sys.argv) > 3 else '', prove([pari(x) for x in sys.argv[2].split(',')]))); sys.exit(0)
    lim = float(sys.argv[1]); rows = []
    for l in open(sys.argv[2]):
        r = json.loads(l)
        if r['status'] != 'closed': continue
        c = [c for c in r['curves'] if c['curve'] == r['by']][0]
        N = int(pari.ellglobalred(pari.ellinit([pari(a) for a in c['rank0_curve']]))[0]); rows.append((N, r['pair'], c['rank0_curve']))
    rows.sort(); done = 0; allpos = True
    for N, pair, ai in rows:
        if N > lim: break
        o = prove([pari(a) for a in ai]); done += 1; allpos = allpos and o.get('positive', False)
        print(show('%d %d' % tuple(pair), o), flush=True)
    print('closing curves of conductor up to %g: %d of %d; L(E,1) > 0 in ball arithmetic for all of them: %s' % (lim, done, len(rows), allpos))
