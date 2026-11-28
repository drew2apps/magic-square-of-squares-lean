"""I8, a "neither" ratio: an exact upper bound for the rank of E+ (s = 2) or E- (s = -2) with the repository's rank_indep.rank_bound
(descent via 2-isogeny with a second descent; no PARI there). PARI only finds the point G, which is checked on the curve.
EXPERIMENT, written by an AI model; not read by a person.  usage: ri.py m n s     Folder: MSS_REPO (default: this repository)."""
import sys,os,time
from fractions import Fraction as Fr
REPO = os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
G1=os.path.abspath(os.path.join(REPO,'scripts','g1_full')); sys.path.insert(0,G1); os.chdir(G1)
argv=sys.argv; sys.argv=[argv[0]]
import rank_indep as RI, sigma_exact as SX
import cypari2; P=cypari2.Pari(); P.default('nbthreads',1)
m,n,s=int(argv[1]),int(argv[2]),int(argv[3])
M=SX.model(m,n,s); E=SX.Curve(M['a2'],M['a4'],M['a6'])
pe=P.ellinit([0,M['a2'],0,M['a4'],M['a6']]); r=P.ellrank(pe,1)
G=[Fr(str(c)) for c in r[3][0]] if len(r[3]) else None
print('PARI point',G is not None); assert G and E.on(tuple(G))
e,a,b=RI.shift_to_two_torsion(M['a2'],M['a4'],M['a6']); GE=(G[0]-e,G[1]); ap,bp=-2*a,a*a-4*b
t0=time.time(); rb,rec=RI.rank_bound(a,b,[GE],RI.find_points([0,ap,0,bp,0]),say=print)
print('rank <=',rb)
