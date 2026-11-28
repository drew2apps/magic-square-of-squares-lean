#!/usr/bin/env python3
"""Exact check of the sieve hits: X = e_j + b U^2/W^2 on E: Y^2 = X (X - e2)(X - e3); canonical heights with PARI over nf."""
import sys, json, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
pari("nf=nfinit(w^2+1); e2=273990441250; e3=Mod(-102929464375+118967496648*w,w^2+1); ee=[0,e2,e3];")
pari("E=ellinit([0,-(e2+e3),0,e2*e3,0],nf);")
tasks = {l.split()[0]: l.split()[1:] for l in open('tasks.txt')}
seen = {}; 
for fn in sys.argv[1:]:
    for l in open(fn):
        if not l.startswith('HIT'): continue
        _, tag, a, b, c, d = l.split(); j = int(tag[1]); bre, bim = tasks[tag][:2]
        X = pari(f"lift(ee[{j}] + Mod({bre}+{bim}*w,w^2+1)*Mod(({a})+({b})*w,w^2+1)^2/Mod(({c})+({d})*w,w^2+1)^2)")
        key = str(X)
        if key in seen: continue
        seen[key] = (tag, a, b, c, d)
pts = []
for key, src in seen.items():
    pari(f"X=Mod({key},w^2+1); f=X*(X-e2)*(X-e3);")
    rts = pari("nfroots(nf, y^2 - f)")
    if len(rts) == 0: print('NOT A POINT (sieve false positive)', src); continue
    pari("P=[X, nfroots(nf, y^2-f)[1]];")
    assert pari("ellisoncurve(E,P)") == 1
    h = float(pari("ellheight(E,P)"))
    print(src, 'X =', key, ' Y =', pari("lift(P[2])"), ' on curve: yes  canonical height %.6f' % h, 'TORSION' if h < 1e-9 else '')
    if h > 1e-9: pts.append((h, key, str(pari("lift(P[2])")), src))
pts.sort()
json.dump(pts, open('points.json', 'w'), indent=1)
print(len(pts), 'non-torsion points')
