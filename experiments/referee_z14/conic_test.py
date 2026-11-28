"""z14: explore3's conic_pairing (the tangent form of the conic of w^2 = d t^4 + 8C^2 t^2 + 16A^2B^2/d; function text copied
verbatim from task23.py, untrusted) against the repository's DFull pairing (bilinear extension of the basis matrix stored
in mix_*.jsonl), for pairs (d, beta) of V on the 134 test curves (up to 6 pairs per curve, fixed seed), soluble or not.
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I conic_test.py | grep -v '^PARI' | diff - conic_test.txt"""
import sys, os, json, glob, random, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from e3_cassels import *
pari.allocatemem(10 ** 8, 4 * 10 ** 8)
import mpmath
from fractions import Fraction as Fr
from math import gcd
def conic_pairing(m,n,d,beta,seed=1,detail=False):
    """prod_v (ell(t^2,z^2,w), beta)_v, ell = tangent form of w^2 = d T^2 + 8C^2 T Z + (16A^2B^2/d) Z^2 at a qfsolve point,
    (t,z=1,w) a local point of w^2 = d t^4 + 8 C^2 t^2 + 16A^2B^2/d found by an independent search."""
    rg=random.Random(seed)
    A=m*m-n*n;B=2*m*n;C=m*m+n*n; W=(C*C+A*B)*(C*C-A*B)
    c0=(16*A*A*B*B)//d; assert c0*d==16*A*A*B*B
    G=pari.matrix(3,3,[d,4*C*C,0,4*C*C,c0,0,0,0,-1])
    s=pari.qfsolve(G)
    if s.type()!='t_COL': return None,{'fail':'conic insoluble'}
    den=1
    for x in s: den=den*int(x.denominator())//gcd(den,int(x.denominator()))
    pt=[int(x*den) for x in s]; g=gcd(gcd(pt[0],pt[1]),pt[2]); T0_,Z0_,w0_=[x//g for x in pt]
    assert d*T0_**2+8*C*C*T0_*Z0_+c0*Z0_**2==w0_**2
    L=[d*T0_+4*C*C*Z0_, 4*C*C*T0_+c0*Z0_, -w0_]; g=gcd(gcd(L[0],L[1]),L[2]); L=[x//g for x in L]
    P={2}|primes_of(A)|primes_of(B)|primes_of(C*C+A*B)|primes_of(C*C-A*B)|primes_of(d)|primes_of(beta)
    for x in L: P|=primes_of(x)
    q=lambda t: d*t**4+8*C*C*t*t+c0
    tot=1; loc={}
    # real place
    if d>0: t=Fr(rg.randrange(1,10**6),rg.randrange(1,10**3))
    else:
        t0=mpmath.sqrt(mpmath.mpf(4*C*C)/abs(d)); t=Fr(int(t0*10**6),10**6)
    qv=q(t)
    if qv<=0: return None,{'fail':'real'}
    w=mpmath.sqrt(mpmath.mpf(qv.numerator)/qv.denominator)*rg.choice([1,-1])
    val=L[0]*mpmath.mpf(t.numerator)**2/t.denominator**2+L[1]+L[2]*w
    if abs(val)<mpmath.mpf(10)**(-60)*(abs(L[1])+abs(L[2]*w)): return None,{'fail':'real cancel'}
    tv=-1 if (val<0 and beta<0) else 1
    loc[-1]=tv; tot*=tv
    for p in sorted(P):
        N=prec_for(p); vmax=vp(16*A*A*B*B,p); done=False
        for attempt in range(60000):
            if p==2: u=rg.randrange(1,2**12)
            else: u=rg.randrange(1,min(p*p,10**12))
            j=rg.randint(-(vmax+4),vmax+4)
            t=Fr(u)*Fr(p)**j; qv=q(t)
            if qv==0 or not issq_int(qv.numerator*qv.denominator,p): continue
            w=psqrt(qv,p,N)*rg.choice([1,-1])
            val=L[0]*pari('(%d)/(%d)'%((t*t).numerator,(t*t).denominator))+L[1]+L[2]*w
            if val==0 or int(val.padicprec(p))-int(pari.valuation(val,p))<4: continue
            tv=int(pari.hilbert(val,beta,p)); done=True; break
        if not done: return None,{'fail':('no local point on quartic',p)}
        loc[p]=tv; tot*=tv
    return tot,(loc if detail else {})

def sqf_(x): return sqf(x) if abs(x) > 1 else x
R = {}
for f in sorted(glob.glob(os.path.join(HERE, 'mix_*.jsonl'))):
    for l in open(f): r = json.loads(l); R[(r['m'], r['n'])] = r
cnt = collections.Counter(); bad = []
for key in sorted(R):
    r = R[key]; m, n = key; sel = [tuple(x) for x in r['sel']]; M = r['M']; k = len(sel)
    def elt(mask):
        c1 = c2 = 1
        for i in range(k):
            if mask >> i & 1: c1 *= sel[i][0]; c2 *= sel[i][1]
        return (sqf_(c1), sqf_(c2))
    V = [mask for mask in range(1, 2 ** k) if elt(mask)[0] == 1]
    prs = [(x, y) for x in V for y in V if x != y]
    rng = random.Random(m * 7919 + n)
    for x, y in (rng.sample(prs, min(6, len(prs))) if prs else []):
        d, beta = elt(x)[1], elt(y)[1]
        s = 0
        for i in range(k):
            for j in range(k):
                if x >> i & 1 and y >> j & 1 and M[i][j] == -1: s ^= 1
        want = -1 if s else 1
        v, info = conic_pairing(m, n, d, beta, seed=5)
        key2 = 'agree' if v == want else ('fail' if v is None else 'DISAGREE')
        cnt[key2] += 1; cnt['pairing %d' % want] += 1
        if key2 != 'agree' and len(bad) < 5: bad.append((m, n, d, beta, v, want, str(info)[:80]))
print('# z14 conic_test. Written by an AI model (reviewer z14); not read by a person.')
print('explore3 conic_pairing vs DFull on %d curves:' % len(R), dict(sorted(cnt.items())))
print('first disagreements/failures:', bad)
