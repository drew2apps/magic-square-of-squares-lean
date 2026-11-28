import sys,re; sys.argv=['x']; sys.path.insert(0,'/tmp/rn'); sys.path.insert(0,'/tmp/ghost')
from pipeline import *
from f4_form import residual5
src=open('/mnt/user-data/outputs/fable_27nov_scripts/genus3_ranks.py').read()
pairs=[(int(a),int(b)) for a,b in re.findall(r'\((\d+)\s*,\s*(\d+)\)',src)]
pairs=[p for p in dict.fromkeys(pairs) if p[0]>p[1]>0][:12]
for m,n in pairs:
    T4,_=residual(m,n); T5,_=residual5(m,n)
    print((m,n),'four-form residual',len(T4),'| five-form residual',len(T5),'| 4-tuples that extend',len(set(t[:4] for t in T5)))
