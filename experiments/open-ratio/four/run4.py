"""Driver: cyclic 4-isogeny descents over K = Q(i) for every entry of a gen4.gp data file.
usage: python3 run4.py data.gp [known_primes_comma_list] [seed]
Written by an AI model (Claude); not read by a person.  EXPERIMENT."""
import sys, random, math
sys.path.insert(0, '/mnt/project-files/mss-experiments/open-ratio/four')
import fourdesc as fd
from fractions import Fraction as F
pari = fd.pari
nf = pari('nfinit(y^2+1)')

def load(path):
    pari('read("%s")' % path)
    return pari('D4')

def gl(z): return fd.gq(z)

def bad_primes(cs, known):
    E = pari.ellinit([0, cs[0], 0, cs[1], cs[2]], nf)
    n = abs(int(pari.idealnorm(nf, E.disc() if hasattr(E, 'disc') else pari('(%s).disc' % E))))
    out = {2}
    for p in sorted(known):
        while n % p == 0: n //= p; out.add(p)
    if n > 1:
        f = pari.factor(n, 10**7)
        for q in f[0]:
            q = int(q)
            if not pari.isprime(q): raise RuntimeError("unfactored part %d" % q)
            out.add(q)
    return out

def tamagawa(cs, pl_list):
    E = pari.ellinit([0, cs[0], 0, cs[1], cs[2]], nf)
    res = {}
    for pl in pl_list:
        prs = pari.idealprimedec(nf, pl.p)
        u = pl.unif
        g = pari('%s + %s*y' % (u[0], u[1]))
        for pr in prs:
            if int(pari.idealval(nf, g, pr)) > 0 or len(prs) == 1:
                res[pl.name] = int(pari.elllocalred(E, pr)[3]); break
    return res

def run_entry(e, known, seed, verbose=True):
    j, tgt = int(e[0]), int(e[1])
    cA = [gl(z) for z in e[2]]; QA = [gl(z) for z in e[3]]
    cB = [gl(z) for z in e[4]]; TB = [gl(z) for z in e[5]]
    torsA = [[gl(z) for z in P] for P in e[6]]; torsB = [[gl(z) for z in P] for P in e[7]]
    S = sorted(bad_primes(e[2], known) | bad_primes(e[4], known))
    places = []; globs = [(F(0), F(1)), (F(1), F(1))]
    for p in S:
        if p == 2: places.append(fd.Ram()); continue
        if p % 4 == 3: places.append(fd.Inert(p)); globs.append((F(p), F(0)))
        else:
            a, b = fd.gprime(p)
            places += [fd.Split(p, (a, b)), fd.Split(p, (a, -b))]
            globs += [(F(a), F(b)), (F(a), F(-b))]
    mB = fd.Kummer4(cB, TB, torsB)     # delta_psi on B(K_v), uses T'
    mA = fd.Kummer4(cA, QA, torsA)     # delta_psihat on A(K_v), uses Q
    tA = tamagawa(e[2], places); tB = tamagawa(e[4], places)
    rng = random.Random(seed)
    genB, genA = [], []
    for pl in places:
        target = 256 if pl.p == 2 else 16
        (gB, gA), (sB, sA) = fd.sample_images(pl, [mB, mA], target, rng)
        ok = sB * sA == target
        expect = None
        if pl.p != 2:
            expect = F(4 * tB[pl.name], tA[pl.name])
            ok = ok and sB == expect
        if verbose: print(f"   {pl.name}: |im psi| = {sB}, |im psihat| = {sA}, product target {target}"
                          + (f", expected |im psi| = 4 c(B)/c(A) = {expect}" if expect is not None else "") + ("" if ok else "   <-- NOT EXACT"), flush=True)
        if not ok: return None
        genB.append(gB); genA.append(gA)
    dl = [len(pl.dims) for pl in places]
    sizeB, invB = fd.selmer(globs, places, genB, dl)
    sizeA, invA = fd.selmer(globs, places, genA, dl)
    bound = math.log(sizeB * sizeA / 16, 4)
    return dict(j=j, tgt=tgt, S=S, selpsi=sizeB, selpsi_inv=invB, selpsihat=sizeA, selpsihat_inv=invA, bound=bound)

if __name__ == '__main__':
    path = sys.argv[1]
    known = [int(x) for x in sys.argv[2].split(',')] if len(sys.argv) > 2 and sys.argv[2] else []
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    D = load(path)
    for k in range(len(D)):
        e = D[k]
        print(f"entry {k+1}: psi: E{int(e[0])} -> E{int(e[1])}", flush=True)
        r = run_entry(e, known, seed)
        if r is None: print("   local image not found exactly; no bound"); continue
        print(f"   S = {r['S']}  |Sel^psi| = {r['selpsi']} (invariants {r['selpsi_inv']}),  |Sel^psihat| = {r['selpsihat']} (invariants {r['selpsihat_inv']})")
        print(f"   rank bound  r <= log_4(|Sel^psi||Sel^psihat|/16) = {r['bound']}", flush=True)
