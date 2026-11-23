#!/usr/bin/env python3
"""L(E,1) > 0 with an explicit error bound, for curves over Q of root number +1.

For E/Q of conductor N and root number +1 (modularity):  L(E,1) = 2 sum_{n >= 1} (a_n / n) q^n,  q = exp(-2 pi / sqrt(N)).
|a_n| <= d(n) sqrt(n) <= 2 n, so the tail after M terms is at most 4 q^(M+1) / (1 - q).  M is taken with q^M < 10^-25.
The a_n are exact integers (ellan).  The printed error ("+-") is this tail bound; it is at most 5.6e-19 in lproof.txt, and
it is the error that matters.

Rounding.  The sum is computed by PARI with 57 significant digits (a mantissa of 192 bits), with q^n by repeated
multiplication.  Assume every operation used here (Pi, sqrt, exp, *, /, +) returns its result with relative error at most
e = 10^-55; one unit in the last place is 2^-191 < 4e-58, so this allows for more than 200 units per operation.  It is the
documented behaviour of PARI's real arithmetic, not interval arithmetic, and it is the one thing assumed.  Then:
  q has relative error at most 6e (the exponent 2 pi / sqrt(N) is below 2 and carries at most 4e, exp adds e, with room);
  the computed q^n has relative error at most 8 n e (n factors and n multiplications; n e is below 10^-46 here);
  the term a_n q^n / n, at most 2 q^n in absolute value, has error at most 2 q^n (8 n + 2) e;
  each of the M additions adds at most e times a partial sum, which is below 2 q / (1 - q).
So the computed sum differs from the exact partial sum by at most
  rnd = e (16 q / (1-q)^2 + 4 q / (1-q) + 2 M q / (1-q)),
and 2 * rnd is added to the error of L.  For the largest conductor used (7.7e13: 1/(1-q) about 1.4e6, M about 8e7) this is
below 10^-40.  The program computes 2*rnd for each curve and calls L positive only if  2*sum - tail - 2*rnd > 0, tested by
PARI at the same precision; the smallest L in lproof.txt is 0.0266, so nothing is near the limit.
L(E,1) > 0 gives rank 0 (Kolyvagin, with modularity).  What is trusted: PARI's conductor, root number and a_n, the assumption
on its real arithmetic, and the two theorems.

usage: lproof.py LIMIT second_123.jsonl          the closing curves of conductor up to LIMIT
       lproof.py curve a1,a2,a3,a4,a6 [label]    one curve"""
import sys, json, time
import cypari2
pari = cypari2.Pari(); pari.allocatemem(6*10**9)
pari('default(realprecision, 57)')
pari('''lsum(E, N) = { my(q, M, an, p, s, tail); q = exp(-2*Pi/sqrt(N)); M = ceil(25*log(10)/(2*Pi/sqrt(N)));
  an = ellan(E, M); p = 1.; s = 0.; for(n = 1, M, p *= q; s += an[n]*p/n); tail = 4*q^(M+1)/(1-q);
  my(e = 1e-55, rnd = 2*e*(16*q/(1-q)^2 + 4*q/(1-q) + 2*M*q/(1-q))); [2*s, tail, M, rnd, 2*s - tail - rnd > 0] }''')

def prove(ai):
    E = pari.ellinit(ai); E = pari.ellminimalmodel(E); N = pari.ellglobalred(E)[0]; w = int(pari.ellrootno(E))
    if w != 1: return dict(conductor=str(N), root_number=w)
    v = pari("lsum")(E, N); L = float(v[0]); tail = float(v[1]); M = int(v[2])
    return dict(conductor='%.3g' % float(N), root_number=w, terms=M, L='%.12f' % L, tail_bound='%.1e' % tail, positive=bool(int(v[4]) == 1),      # v[3] is the rounding bound 2*rnd of the docstring
                L_over_Omega='%.9f' % (L/float(E.omega()[0])))

if __name__ == '__main__':
    if sys.argv[1] == 'curve':
        r = prove([pari(x) for x in sys.argv[2].split(',')]); print(sys.argv[3] if len(sys.argv) > 3 else '', json.dumps(r)); sys.exit(0)
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
        print(pair[0], pair[1], 'N = %s  terms %d  L(E,1) = %s +- %s  L/Omega = %s  %s' % (o['conductor'], o['terms'], o['L'], o['tail_bound'], o['L_over_Omega'], 'positive' if o['positive'] else 'NOT SHOWN POSITIVE'), flush=True)
    print('closing curves of conductor up to %g: %d of %d; L(E,1) > 0 with the error bound for all of them: %s' % (lim, done, len(rows), allpos))
