import re,glob,collections
from math import gcd,isqrt
from sympy import factorint as _fi
def sqf(n):
    if n==0: return 0
    s=1 if n>0 else -1; r=1
    for p,e in _fi(abs(int(n))).items():
        if e%2: r*=p
    return s*r
def load():
    R=[]
    for f in sorted(glob.glob('t10_*_[ab].txt')):
        cur=None; cv=None
        for l in open(f):
            mm=re.match(r'ratio \((\d+),(\d+)\)',l)
            if mm: cur=(int(mm[1]),int(mm[2])); continue
            mm=re.match(r'\s+E_(\w\w): base t0 = (\S+); e = \(([-\d, ]+)\); e_r = e(\d), e\* = e(\d).*dim Sel2 = (\d+), point images (\d+)',l)
            if mm:
                cv=dict(ratio=cur,cv=mm[1],t0=mm[2],e=tuple(int(x) for x in mm[3].split(',')),er=int(mm[4]),es=int(mm[5]),dim=int(mm[6]),pim=int(mm[7]),lines=[]); R.append(cv); continue
            mm=re.match(r'\s+tuple \(([-\d, ]+)\): class \(([-\d]+), ([-\d]+)\): (.*)',l)
            if mm and cv is not None:
                tup=tuple(int(x) for x in mm[1].split(',')); cls=(int(mm[2]),int(mm[3])); rest=mm[4]
                k=re.match(r'killed b = \(([-\d]+), ([-\d]+)\), local -1 at \[([\d, ]*)\]',rest)
                cv['lines'].append(dict(tup=tup,cls=cls,b=(int(k[1]),int(k[2])) if k else None,loc=[int(x) for x in k[3].split(',')] if k and k[3] else [],rest=rest if not k else ''))
    return R
if __name__=='__main__':
    R=load()
    print('curves',len(R),'tuple lines',sum(len(c['lines']) for c in R),'kill lines',sum(1 for c in R for l in c['lines'] if l['b']))
    print('non-kill kinds',collections.Counter(l['rest'][:40] for c in R for l in c['lines'] if not l['b']).most_common(6))
    okE=collections.Counter(); okC=collections.Counter(); pos=collections.Counter()
    idx={'0':0,'8':1,'L':2,'U':3}
    for c in R:
        m,n=c['ratio']; S=m*m+n*n; e=c['e']; cvn=c['cv']
        R1=(  (m*m-2*m*n-n*n)**2+S*S)//2; R2=((m*m+2*m*n-n*n)**2+S*S)//2
        er=e[c['er']-1]; es=e[c['es']-1]; em=[x for i,x in enumerate(e) if i+1 not in (c['er'],c['es'])][0]
        d=(em-er,es-er)   # differences from e_r
        pos[(cvn,c['er'],c['es'])]+=1
        if cvn=='08': okE[(cvn, d==((m*m-n*n)**2,S*S))]+=1
        elif cvn=='LU': okE[(cvn, d, )] if False else okE.update([(cvn, d==(16*m*m*n*n,4*S*S), d==(4*S*S,16*m*m*n*n), d==(4*(m*m-n*n)**2,4*S*S))])
        else: okE[(cvn, d in ((2*R1,2*R2),(2*R2,2*R1)), d==(2*R1,2*R2))]+=1
        for l in c['lines']:
            D=sqf(l['tup'][idx[cvn[0]]]*l['tup'][idx[cvn[1]]])
            k1=sqf(l['cls'][0]); k2=sqf(l['cls'][1]*D)
            cand={'08':(1,sqf(m*n)),'LU':(1,sqf(m*m-n*n)),'0L':(2,sqf(m*n*(m+n))),'0U':(2,sqf(m*n*(m-n))),'8L':(2,sqf(m*n*(m-n))),'8U':(2,sqf(m*n*(m+n)))}[cvn]
            okC[(cvn,k1==cand[0],abs(k2)==abs(cand[1]) , k2>0)]+=1
    print('e-formulas',dict(okE)); print('positions',dict(pos)); print('class formulas (curve, c_r ok, |c*| ok, sign+)',dict(okC))
