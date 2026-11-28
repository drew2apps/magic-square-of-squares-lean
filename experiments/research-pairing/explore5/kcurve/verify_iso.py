#!/usr/bin/env python3
"""Exact check of hits on E' = E/<T3>: y'^2 = x'(x'^2 + A1 x' + d2^2), x' = d U^2/W^2; push to E by the dual isogeny
X = e3 + (y'/(2x'))^2; heights and 2-descent classes on E with PARI over nf."""
import sys, json, cypari2
pari = cypari2.Pari(); pari.allocatemem(2*10**9)
pari("nf=nfinit(w^2+1); e2=273990441250; e3=Mod(-102929464375+118967496648*w,w^2+1); A1=2*(e2-2*e3);")
pari("E=ellinit([0,-(e2+e3),0,e2*e3,0],nf); Ep=ellinit([0,A1,0,e2^2,0],nf);")
tasks = {l.split()[0]: l.split()[1:] for l in open('tasks_iso.txt')}
seen = {}
for fn in sys.argv[1:]:
    for l in open(fn):
        if not l.startswith('HIT'): continue
        _, tag, a, b, c, d = l.split(); bre, bim = tasks[tag][:2]
        xp = str(pari(f"lift(Mod({bre}+{bim}*w,w^2+1)*Mod(({a})+({b})*w,w^2+1)^2/Mod(({c})+({d})*w,w^2+1)^2)"))
        seen.setdefault(xp, (tag, a, b, c, d))
pts = []
for xp, src in seen.items():
    pari(f"xp=Mod({xp},w^2+1); f=xp*(xp^2+A1*xp+e2^2); rt=nfroots(nf, y^2-f);")
    if pari("#rt") == 0: print('NOT A POINT', src); continue
    pari("Pp=[xp,rt[1]]; X=e3+(Pp[2]/(2*xp))^2; g=X*(X-e2)*(X-e3); P=[X, nfroots(nf,y^2-g)[1]];")
    assert pari("ellisoncurve(Ep,Pp)") == 1 and pari("ellisoncurve(E,P)") == 1
    hp = float(pari("ellheight(Ep,Pp)")); h = float(pari("ellheight(E,P)"))
    print(src, "x' =", xp, " hhat(E') = %.6f  hhat(E) = %.6f" % (hp, h), 'TORSION' if h < 1e-9 else '')
    if h > 1e-9: pts.append(dict(src=src, xprime=xp, yprime=str(pari("lift(Pp[2])")), X=str(pari("lift(P[1])")), Y=str(pari("lift(P[2])")), h_Eprime=hp, h_E=h))
pts.sort(key=lambda r: r['h_E'])
json.dump(pts, open('points_iso.json', 'w'), indent=1)
print(len(pts), 'non-torsion points')
