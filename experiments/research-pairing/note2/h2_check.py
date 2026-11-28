# Reviewer points (2) and (4): formula (2) against the general Cassels engine, with and without hypothesis H2.
import sys; sys.argv=['x']; sys.path.insert(0,'/tmp/rn'); sys.path.insert(0,'/tmp/z12'); sys.path.insert(0,'/tmp/mix')
from pipeline import Curve, cassels, triple as tr3, sqf
import closed_formula as CF
from math import gcd
import itertools, collections
def run(m,n,cv):
    A,C=CF.triple(m,n,cv)                      # the note's (A,C) for this curve
    e=(0,A*A,C*C); E=Curve(e); B=E.basis; d=len(B)
    h2 = (A%2==0) or ((C-A)%32==0) or ((C+A)%32==0)
    res=collections.Counter()
    # all Selmer elements
    def mul(x,y): return (sqf(x[0]*y[0]), sqf(x[1]*y[1]))
    els=[(1,1)]
    for b in B: els+= [mul(x,b) for x in els]
    for a in els:
        if a[0]!=1: continue
        dl=a[1]
        for b in els:
            ref=cassels(e,tr3(a),tr3(b),seed=5)[0]
            b3=tr3(b)
            try: val=CF.closed(A,C,dl,b3[0],b3[1])
            except Exception as ex: res['error']+=1; continue
            res[val==ref]+=1
    print((m,n),cv,'(A,C)=',(A,C),'H2' if h2 else 'H2 fails','dim Sel2',d,'agree/disagree',dict(res))
run(725,497,'LU'); run(725,497,'08'); run(535,23,'08'); run(465,416,'LU'); run(465,416,'08')
