"""Review x6 (run 6), written by an AI model; not read by a person.  Does the check of local_compare.py have teeth?  Five
mutations of exact_local.py (made in memory from its source; the file is not changed), each of the kind that would make
a local set too small: the tail takes one parity of e only; no closure under a non-square unit; listing keeps one sign
pattern only; the dead test at p = 3 mod 4 reversed; the symbol of b' forgotten.  For each, the labels of explicit
points (mylocal.py) are compared with the mutated sets on 24 (ratio, prime) cases: a mutation is caught if some point
falls outside.  The unmutated code must have none.  Prints only.   usage: python3 -I local_mutation.py"""
import sys, os, random, types
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mylocal as ML
X, TD = ML.X, ML.TD
SRC = open(X.__file__).read()
MUT = [('tail: one parity of e', 'for e in (0, 1):', 'for e in (0,):'),
       ('no unit closure', 'if self.kind == 1: out |= set(tuple((l[0], 1 - l[1]) for l in t) for t in out)\n        return out', 'return out'),
       ('listing keeps one pattern', 'for c in np.unique(code)]', 'for c in np.unique(code)[:1]]'),
       ('dead test at p = 3 mod 4 reversed', 'if leg(u[0]**2 + u[1]**2, p): return', 'if not leg(u[0]**2 + u[1]**2, p): return'),
       ("symbol of b' forgotten", 'lab.append((f[1] % 2, leg(f[2], p) ^ e[f[3]]))', 'lab.append((f[1] % 2, e[f[3]]))')]
CASES = [((627, 364), p) for p in (3, 5, 7, 11, 13, 17, 29, 37, 263)] + [((1109503, 597671), p) for p in (5, 7, 137, 593, 1871, 7129)] + \
        [((10125, 4459), p) for p in (3, 5, 7, 13)] + [((993260025, 766111057), p) for p in (3, 5, 13, 100003, 1670213)]
def module(src):
    M = types.ModuleType('exact_local_mut'); M.__file__ = X.__file__; exec(compile(src, X.__file__, 'exec'), M.__dict__); return M
def caught(M):
    rnd = random.Random(5); bad = 0
    for (m, n), p in CASES:
        forms, sup, known = TD.magic(m, n)
        try: L, E = M.allowed_at(forms, p)
        except Exception as ex: bad += 1; continue                 # a mutation that crashes is caught too
        S = set(ML.points(forms, p, rnd).keys())
        bad += not S <= E
    return bad
print('unmutated: %d of %d cases with an explicit point outside the set' % (caught(module(SRC)), len(CASES)))
for name, a, b in MUT:
    assert SRC.count(a) == 1, name
    k = caught(module(SRC.replace(a, b)))
    print('mutation "%s": caught on %d of %d cases%s' % (name, k, len(CASES), '' if k else '  (NOT CAUGHT)'))
