import json, collections, sys, re
sys.path.insert(0,'/tmp/t10')
recs=[json.loads(l) for f in ('out_997.jsonl','out_15000_0.jsonl','out_15000_1.jsonl') for l in open(f)]
CV=['08','LU','0L','0U','8L','8U']
# 1. validation against the pilot's kill lines (201 ratios)
import os; cwd=os.getcwd(); os.chdir('/tmp/t10'); from parse import load; R=load(); os.chdir(cwd)
pil=collections.defaultdict(set)
for c in R:
    for l in c['lines']: pil[(c['ratio'],l['tup'])].add(c['cv'])
pil_ratios={c['ratio'] for c in R}
agree=collections.Counter(); ntup=0
for r in recs:
    key=(r['m'],r['n'])
    if key not in pil_ratios: continue
    for t,d in r['tuples']:
        ntup+=1
        for cv in CV:
            mine=(d[cv]==1); theirs=cv in pil[(key,tuple(t))]
            agree[(mine,theirs)]+=1
print('pilot ratios in my run:',sum((r['m'],r['n']) in pil_ratios for r in recs),'tuples',ntup,'; (mine killed, pilot killed) counts over (tuple,curve):',dict(agree))
# 2. the whole population
st=collections.Counter(); surv=[]; outside=0; T=0; per=collections.Counter(); expl=0
for r in recs:
    alive=[]
    for t,d in r['tuples']:
        T+=1
        if any(v=='outside' for v in d.values()): outside+=1
        k=[cv for cv in CV if d[cv]==1]
        for cv in k: per[cv]+=1
        if any(cv in ('08','LU') for cv in k): expl+=1
        if not k: alive.append(t)
    st['ratios']+=1
    if alive: st['ratios with a tuple that no pair curve kills']+=1; surv.append((r['m'],r['n'],alive,r['dims']))
print('ratios',st['ratios'],'residual tuples',T,'classes outside Sel2:',outside)
print('killed per curve:',dict(per)); print('tuples killed on E_08 or E_LU:',expl,'of',T)
print('ratios with a surviving tuple:',st['ratios with a tuple that no pair curve kills'])
for s in surv: print('  ',s[0],s[1],'tuples',s[2],'dims',s[3])
json.dump(surv,open('survivors.json','w'))
