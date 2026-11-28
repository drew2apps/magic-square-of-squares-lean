"""I8, a "neither" ratio: E+ and E- of G1_FULL_SYSTEM.md (sigma_exact.model) with PARI's ellrootno, ellrank (effort 1) and elltors.
EXPERIMENT, written by an AI model; not read by a person. Finder output only: the ranks are "consistent with PARI", not proved.
usage: epm.py m n     Folder: MSS_REPO (default: the repository this file sits in)."""
import sys,os,time
REPO = os.environ.get('MSS_REPO', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
sys.path.insert(0,os.path.join(REPO,'scripts','g1_full'))
import sigma_exact as SX
import cypari2; P=cypari2.Pari(); P.allocatemem(2*10**9); P.default('nbthreads',1)
m,n=map(int,sys.argv[1:3])
for s,name in ((2,'E+'),(-2,'E-')):
    M=SX.model(m,n,s); E=P.ellinit([0,M['a2'],0,M['a4'],M['a6']])
    t0=time.time()
    r=P.ellrank(E,int(sys.argv[3]) if len(sys.argv)>3 else 1)
    print(name,'a2,a4,a6 =',M['a2'],M['a4'],M['a6']); print('  rootno',P.ellrootno(E),'ellrank [lo,hi]=',r[0],r[1],'s=',r[2],'pts',len(r[3]),'tors',P.elltors(E)[1],); sys.stdout.flush()
