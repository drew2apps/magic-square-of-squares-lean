"""3(b) on all 937 ratios with residual 4-tuples (m <= 15000) and on the sampler's ratios: how many keep a residual
5-tuple once the centre form gamma_4 is added (bad primes only), and what happens to the pairing survivors."""
import sys,json,collections
sys.argv=['x']; sys.path.insert(0,'/tmp/rn'); sys.path.insert(0,'/tmp/ghost')
from f4_form import residual5, residual, sqf
files=sys.argv[1:] if False else ['/tmp/rn/out_997.jsonl','/tmp/rn/out_15000_0.jsonl','/tmp/rn/out_15000_1.jsonl']
import glob
files+=sorted(glob.glob('/tmp/ghost/samp_*.jsonl'))
tot=collections.Counter(); out=open('/tmp/ghost/f4_all.jsonl','w')
for f in files:
    tag='list' if '/tmp/rn/' in f else 'sample'
    for l in open(f):
        r=json.loads(l)
        if 'tuples' not in r: continue
        m,n=r['m'],r['n']
        T5,_=residual5(m,n); ext=set(t[:4] for t in T5)
        surv=[tuple(t) for t,k in r['tuples'] if not any(v==1 for v in k.values())]
        tot[(tag,'ratios')]+=1; tot[(tag,'4-tuples')]+=len(r['tuples']); tot[(tag,'4-tuples that extend to a 5-tuple')]+=len(ext)
        tot[(tag,'ratios with a residual 5-tuple')]+=bool(T5)
        tot[(tag,'pairing survivors')]+=len(surv); tot[(tag,'pairing survivors that extend')]+=sum(t in ext for t in surv)
        tot[(tag,'ratios with a pairing survivor')]+=bool(surv); tot[(tag,'ratios with a pairing survivor that extends')]+=any(t in ext for t in surv)
        tot[(tag,'tuples alive after pairing AND F4')]+=sum(t in ext for t in surv)
        out.write(json.dumps(dict(m=m,n=n,n4=len(r['tuples']),n5=len(T5),ext=sorted(ext),surv=surv))+'\n'); out.flush()
for k,v in sorted(tot.items()): print(k,v)
