import sys,time,random,collections,json
sys.path.insert(0,'/tmp/mix'); sys.path.insert(0,'/tmp/t10')
from parse import load
from cassels import *
R=load(); rng=random.Random(20261127); T0=time.time()
out={}
def run(lines,seed=1):
    cnt=collections.Counter(); bad=[]
    for c,l in lines:
        a=triple(l['cls']); b=triple(l['b'])
        try: v,info=cassels(c['e'],a,b,seed=seed)
        except Exception as ex: v,info=None,{'fail':repr(ex)}
        cnt[v]+=1
        if v!=-1: bad.append((c['ratio'],c['cv'],l['tup'],a,b,v,info))
    return cnt,bad
# (a)
for cv in ('08','LU'):
    L=[(c,l) for c in R if c['cv']==cv for l in c['lines']]
    S=rng.sample(L,300); cnt,bad=run(S)
    print('(a)',cv,'sample 300 of',len(L),dict(cnt),'bad',bad[:3],'t=%.0f'%(time.time()-T0),flush=True); out['a_'+cv]={str(k):v for k,v in cnt.items()}
# (b) all mixed lines
M=[(c,l) for c in R if c['cv'] in ('0L','0U','8L','8U') for l in c['lines']]
cnt,bad=run(M)
print('(b) mixed lines',len(M),dict(cnt),'bad',bad[:5],'t=%.0f'%(time.time()-T0),flush=True); out['b_mixed']={str(k):v for k,v in cnt.items()}; out['b_bad']=[str(x) for x in bad[:20]]
percv=collections.Counter((c['cv']) for c,l in M); print('  per curve',dict(percv))
# (c) controls on 200 mixed lines (+ 100 from 08/LU)
S=rng.sample(M,200)+rng.sample([(c,l) for c in R if c['cv'] in('08','LU') for l in c['lines']],100)
cc=collections.Counter()
for c,l in S:
    a=triple(l['cls']); b=triple(l['b']); e=c['e']
    cc[('<a,a>',cassels(e,a,a)[0])]+=1
    cc[('<b,b>',cassels(e,b,b)[0])]+=1
    for i in range(3): cc[('<a,T%d>'%(i+1),cassels(e,a,tors_image(e,i))[0])]+=1
    for i in range(3): cc[('<T%d,b>'%(i+1),cassels(e,tors_image(e,i),b)[0])]+=1
    cc[('<b,a> (symmetry, expect -1)',cassels(e,b,a)[0])]+=1
    ab=tuple(sqf(x*y) for x,y in zip(a,b))
    cc[('<a,ab> (expect -1)',cassels(e,a,ab)[0])]+=1
print('(c) controls on 300 lines (200 mixed + 100 of 08/LU):')
for k in sorted(cc,key=str): print('   ',k,cc[k])
out['c']={str(k):v for k,v in cc.items()}
# (d) reseed on 200 mixed lines, 3 other seeds, compare local values prime by prime? (local values may differ; total must not)
S=rng.sample(M,200); dd=collections.Counter(); locdiff=0
for c,l in S:
    a=triple(l['cls']); b=triple(l['b']); e=c['e']
    vals=[]; locs=[]
    for seed in (1,777,424242,99991):
        v,info=cassels(e,a,b,seed=seed,detail=True); vals.append(v); locs.append(info)
    dd[tuple(vals)]+=1
    if any(locs[0]!=x for x in locs[1:]): locdiff+=1
print('(d) 200 mixed lines x seeds (1,777,424242,99991): value tuples',dict(dd),'; lines where some local factor changed with the seed:',locdiff)
out['d']={str(k):v for k,v in dd.items()}; out['d_localdiff']=locdiff
print('partial factorizations',STATS,'total time %.0f s'%(time.time()-T0))
out['stats']=STATS
# first-component statistics of b on mixed kill lines
fc=collections.Counter()
for c,l in M:
    m,n=c['ratio']; A=m*m-n*n;B=2*m*n;C=m*m+n*n; W2=2*(C*C+A*B)*(C*C-A*B)
    b1=l['b'][0]   # component at e1 = e_r
    if b1==1: fc['b_r = 1']+=1
    elif b1==2: fc['b_r = 2']+=1
    elif b1>0 and W2%b1==0: fc['b_r | 2R+R-, not 1,2']+=1
    elif b1<0 and W2%b1==0: fc['b_r<0, |b_r| divides 2R+R-']+=1
    else: fc['other']+=1
    if b1>0 and W2%b1==0: fc['total: positive divisor of 2R+R- (incl 1,2)']+=1
print('b first component (at e_r=e1) over',len(M),'mixed kill lines:',dict(fc)); out['b_first']=dict(fc)
json.dump(out,open('/tmp/mix/task1_results.json','w'),indent=1)
