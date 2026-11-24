import sys, os, importlib.util
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare as C
from mytuple import *
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tuple_descent.py')).read()
l1="        if self.kind == 1: out |= set(tuple((l[0], 1 - l[1]) for l in t) for t in out)\n"
l2="        if self.kind == 2: out |= set(tuple(self.times3(l) for l in t) for t in out)\n"
h ="            out |= set(product(*[[r[2]] if (r[1] or r[2] is not None) else every for r in res]))"
assert l1 in src and l2 in src and h in src
muts={'no unit closure (split)':src.replace(l1,''),'no unit closure (p=2)':src.replace(l2,''),
      'hint used but wrong component':src.replace("(out[1] if out[0] is None else out[0])","(out[0] if out[0] is None else (out[0][0], 1-out[0][1]))")}
for name,code in muts.items():
    open(os.path.join(__import__('tempfile').gettempdir(), '_mut_tuple_descent.py'),'w').write(code)
    spec=importlib.util.spec_from_file_location('m',os.path.join(__import__('tempfile').gettempdir(), '_mut_tuple_descent.py')); M=importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
    C.TD=M; bad=0; tot=0; first=None
    for sysname,(forms,sup) in (('magic 627,364',magic(627,364)),('magic 15,8',magic(15,8)),('rd 26/225',ratdist(26,225)),('rd 1/8',ratdist(1,8))):
        for p in (2,5,13,17,29,37,41):
            U,L,N,how=local_best(forms,p); R=C.repo_allowed(forms,p); tot+=1
            if not set(L)<=R:
                bad+=1
                if first is None: t=sorted(set(L)-R)[0]; first=(sysname,p,t,L[t])
    print("mutation '%s': certified tuple rejected in %d of %d (system,prime) cases; first: %s"%(name,bad,tot,first))
