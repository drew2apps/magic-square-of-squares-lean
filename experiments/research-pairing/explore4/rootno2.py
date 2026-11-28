# root numbers on all 937 ratios whose tuples survive the bad primes, split by "has a tuple no pair curve kills"
import json, collections, cypari2
P=cypari2.Pari(); P.allocatemem(2*10**9)
def w(e):
    E=P.ellinit([0,-(e[0]+e[1]+e[2]),0,e[0]*e[1]+e[0]*e[2]+e[1]*e[2],-e[0]*e[1]*e[2]]); return int(P.ellrootno(E))
rows=[json.loads(l) for f in ('out_997.jsonl','out_15000_0.jsonl','out_15000_1.jsonl') for l in open(f)]
c=collections.Counter(); cm=collections.Counter()
for r in rows:
    m,n=r['m'],r['n']; A,B,C=m*m-n*n,2*m*n,m*m+n*n
    c[(w((0,A*A,C*C)),w((0,B*B,C*C)))]+=1; cm[w((0,2*(C*C-A*B),2*(C*C+A*B)))]+=1
print('937 ratios with residual tuples: (w_08,w_LU)',dict(c),' w_mix',dict(cm))
