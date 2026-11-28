import json,sys,collections,glob
from sympy import factorint
def sqf(x):
    r=1
    for q,e in factorint(abs(x)).items():
        if e%2: r*=q
    return r
for f in sorted(glob.glob('/tmp/ghost/samp_*.jsonl')):
    tried=None; c=collections.Counter(); surv=[]
    for l in open(f):
        r=json.loads(l)
        if 'done' in r: tried=r['tried']; continue
        if 'error' in r: c['errors']+=1; continue
        if 'skipped' in r: c['skipped (over 400 tuples)']+=1; continue
        m,n=r['m'],r['n']; c['ratios with residual tuples']+=1; tried_last=r['tried']
        for t,k in r['tuples']:
            c['tuples']+=1
            a=sqf(m*n*t[0]*t[1]); b=sqf((m*m-n*n)*t[2]*t[3])
            gh=(a,b)==(sqf(2*m*n),1) or (a,b)==(1,sqf(2*(m*m-n*n)))
            c['ghost-shape tuples']+=gh
            if 'outside' in k.values(): c['outside Sel2']+=1
            if not any(v==1 for v in k.values()):
                c['survivors']+=1; c['survivors not of ghost shape']+=(not gh); surv.append((m,n,t,gh,factorint(m*m+n*n)))
    print(f,'ratios tried:',tried if tried else '>= %d (still running)'%tried_last)
    for k in ('ratios with residual tuples','tuples','ghost-shape tuples','survivors','survivors not of ghost shape','outside Sel2','errors','skipped (over 400 tuples)'): print('   ',k,c[k])
    print('    ratios with a survivor:',len(set(s[:2] for s in surv)))
    for s in surv: print('      ',s)
