"""Does any test notice a wrong hint in rule R4 (LAST_PAIR.md, section 14)?

mut.py answers no on the four systems it tries: there the hint is never the deciding rule.  It is the deciding rule on
the system of the perfect cuboid (section 14, "What the exhaustion found"), so the mutation is repeated here on the
three cuboid bases of compare.py, at every split prime up to 61 and every split prime of the support.  Two wrong
hints: the class with the residue symbol reversed, and the hint dropped (every class allowed, the program as it was
before the exhaustion found it too weak).  Reported for each: the number of (base, prime) pairs where the allowed set
changes, where a tuple with an explicit witness is rejected, and where a tuple that the exhaustion of mytuple.py removes
is allowed (the comparison `R inside U` of compare.py).
"""
import sys, os, importlib.util, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare as C
from mytuple import *
TD0 = C.TD
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'tuple_descent.py')).read()
a = "(out[1] if out[0] is None else out[0])"
assert src.count(a) == 1
muts = {'hint with the residue symbol reversed': src.replace(a, "(lambda h: h if h is None else (h[0], 1 - h[1]))(out[1] if out[0] is None else out[0])"),
        'hint dropped': src.replace(a, "None")}
cases = []
for (m, n) in ((15, 8), (3, 2), (2, 1)):
    f3, s3, kn = TD0.cuboid(m, n)
    forms = [tuple(map(tuple, f)) for f in f3] + [tuple(map(tuple, f3[2]))]                 # as in compare.py: the third form twice
    ps = sorted(set(q for q in range(5, 62) if q % 4 == 1 and all(q % r for r in range(2, q))) | set(q for s in s3 for q in s if q % 4 == 1))
    cases += [((m, n), forms, p) for p in ps]
base = {}
for (mn, forms, p) in cases:
    C.TD = TD0; base[(mn, p)] = C.repo_allowed(forms, p)
for name, code in muts.items():
    path = os.path.join(tempfile.gettempdir(), '_mut_hint_tuple_descent.py'); open(path, 'w').write(code)
    spec = importlib.util.spec_from_file_location('mh', path); M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
    C.TD = M; changed = []; rejected = []; larger = 0; notinU = 0
    for (mn, forms, p) in cases:
        R = C.repo_allowed(forms, p)
        if R != base[(mn, p)]: changed.append((mn, p, len(base[(mn, p)]), len(R)))
        if R > base[(mn, p)]: larger += 1
        U, L, N, how = local_best(forms, p)
        if not set(L) <= R: rejected.append((mn, p))
        if not R <= set(U): notinU += 1
    print("mutation '%s': %d (base, prime) pairs; allowed set changes at %d, strictly larger at %d, a tuple with a witness is rejected at %d, a tuple that the exhaustion removes is allowed at %d" % (name, len(cases), len(changed), larger, len(rejected), notinU))
    print("   changes (base, prime, size before, size after):", changed[:12])
    print("   rejected witnesses at:", rejected[:12])
