import sys, os, json, time
sys.argv = [sys.argv[0]] + sys.argv[1:]
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'explore_models.py')).read().split('for lab, src, m, n, d in CURVES:')[0])
out = {}
sel = sys.argv[1:]
for lab, src, m, n, d in CURVES:
    if sel and lab not in sel: continue
    t0 = time.time()
    if src == 'E1':
        Q, r = SD.family(m, n)['+13']; ain, A = SD.weierstrass(d * Q, r)
    else:
        ain = [Fr(int(v)) for v in [c[0] for c in DM.all_curves(m, n) if c[4] == src][0]]
    E0 = pari.ellinit([pari(x.numerator)/x.denominator for x in ain]); Em = pari.ellminimalmodel(E0)
    mi = [int(v) for v in Em[:5]]
    a1_, a2_, a3_, a4_, a6_ = mi
    e, a, b = RI.shift_to_two_torsion(a2_, a4_, a6_)
    ap, bp = -2*a, a*a - 4*b
    pts = RI.find_points([0, a, 0, b, 0])
    ptsp = RI.find_points([0, ap, 0, bp, 0])
    rb, rec = RI.rank_bound(a, b, pts, ptsp, say=print)
    out[lab] = dict(a=str(a), b=str(b), G=[[str(x) for x in P] for P in pts], Gp=[[str(x) for x in P] for P in ptsp], bound=rb,
                    E_open=rec['E']['open'], Ep_open=rec["E'"]['open'], E_known=rec['E']['known'], Ep_known=rec["E'"]['known'],
                    E_log2=rec['E']['log2'], Ep_log2=rec["E'"]['log2'], S=len(rec['E']['selmer']), Sp=len(rec["E'"]['selmer']),
                    bad=rec['bad'], sec=round(time.time()-t0,1))
    print(lab, json.dumps(out[lab])[:3000], flush=True)
    json.dump(out, open('explore_rank_%s.json' % ('all' if not sel else sel[0].replace('(','_').replace(')','').replace(',','_')), 'w'), indent=1)
json.dump(RI.HINTS, open('explore_ri_hints.json','w'))
