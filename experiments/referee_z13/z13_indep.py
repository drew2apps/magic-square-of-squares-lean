"""Reviewer z13: independent test of PAIRING_FORMULA_2.md on ratios that are NOT in the pilot.
Written by an AI model (reviewer z13), not read by a person.

For each (m, n, cv): the curve y^2 = x(x - A^2)(x - C^2) of the triple of E_08 / E_LU, shifted as t10_pairct.PairData
does (middle root 0); Sel^2 by rank_second.full_selmer2 (seed 1); the Cassels-Tate pairing by t10_pairct.DFullX
(rank_second.DFull) for a = (1, delta, delta) in Sel^2 against b over a basis of Sel^2; then, by bilinearity of
the pairing, over all of Sel^2.  Compared with:
  formula (2) (z13_formula.formula2) on every (a, b), b over all of Sel^2;
  S' (z13_formula.Sprime) with {b in Sel^2 : b_r = 1};
  Theta (z13_formula.Theta) with the pairing on S' x S'.
a runs over all of {a in Sel^2 : a_r = 1} when it has at most 2^AMAX elements, else over a basis of it plus
EXTRA further elements (and rank_second's values on those extras are checked against bilinearity).
usage: z13_indep.py OUT.jsonl JOBFILE [i0 i1]     (JOBFILE: lines 'm n cv category')
"""
import sys, os, json, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
ARGS = sys.argv[1:]; sys.argv = sys.argv[:1]
sys.path.insert(0, os.path.join(HERE, '..', 'conductor-t10'))
import t10_pairct as T
RS, RI = T.RS, T.RI
sys.path.insert(0, HERE)
import z13_formula as F

AMAX = 7
EXTRA = 10

def span(basis):
    els = [(1, 1)]
    for b in basis:
        els = els + [(RI.sqfree_part(x[0] * b[0]), RI.sqfree_part(x[1] * b[1])) for x in els]
    return els          # element i = product of basis[k] over the bits k of i

def run(m, n, cv, cat):
    t0 = time.time()
    A, B, C = F.triple(m, n, cv)
    es = (0, A * A, C * C)
    esh = tuple(e - es[1] for e in es)            # as PairData: the middle root moved to 0
    bad, loc, sel, kum, enc, places = RS.full_selmer2(esh, random.Random(1))
    allsel = span(sel)
    assert len(set(allsel)) == len(allsel)
    idx = {c: i for i, c in enumerate(allsel)}
    # the part with b_r = 1 and a basis of it
    one = [c for c in allsel if c[0] == 1]
    S1 = sorted(c[1] for c in one)
    Sp = sorted(F.Sprime(A, B, C))
    rec = dict(m=m, n=n, cv=cv, cat=cat, A=A, B=B, C=C, H2=F.H2(A, C), dim=len(sel), sel=sel,
               dimS1=len(one).bit_length() - 1, Sprime_eq=(S1 == Sp), S1=S1, Sprime=Sp)
    # basis of 'one' (as indices into allsel, i.e. bit vectors over sel)
    def red(v, B_):
        for b in B_: v = min(v, v ^ b)
        return v
    obasis = []
    for c in one:
        v = idx[c]
        if red(v, obasis): obasis.append(red(v, obasis)); obasis.sort(reverse=True)
    # which a get rank_second evaluations
    if len(one) <= 2 ** AMAX: adirect = [idx[c] for c in one]
    else:
        rng = random.Random(1000 * m + n)
        rest = [idx[c] for c in one if idx[c] not in obasis and idx[c] != 0]
        adirect = list(obasis) + rng.sample(rest, min(EXTRA, len(rest)))
    rs = {}           # a index -> list of values on the basis of Sel^2
    fails = []
    for ai in adirect:
        a = allsel[ai]
        if ai == 0: rs[ai] = [1] * len(sel); continue     # a = (1,1,1): the pairing is 1 (and DFull needs a covering)
        D = T.DFullX(esh, a)
        vals = []
        for b in sel:
            h, per = T.safe_pair(D, b)
            vals.append(h)
        if None in vals: fails.append((a, vals)); continue
        rs[ai] = vals
    # bilinear extension: value(a_i, b_j) for every a in 'one' and every b in Sel^2
    # first in a: express a as a combination of obasis (bit vectors over sel) -- reduce with tracking
    piv = []          # (vector, combination of obasis elements as a bitmask over obasis positions)
    ob = sorted(set(obasis), reverse=True)
    def combo(v):
        c = 0
        for i, b in enumerate(ob):
            if v ^ b < v: v ^= b; c |= 1 << i
        assert v == 0
        return c
    def rs_val(ai, bj):
        """rank_second's pairing <allsel[ai], allsel[bj]> by bilinearity from the basis values (a over ob, b over sel)."""
        c = combo(ai) if ai not in rs else None
        tot = 1
        for k in range(len(sel)):
            if not (bj >> k & 1): continue
            if c is None: tot *= rs[ai][k]
            else:
                for i, b in enumerate(ob):
                    if c >> i & 1: tot *= rs[b][k]
        return tot
    ok_basis = all(b in rs for b in ob)
    rec['n_a'] = len(one); rec['n_a_direct'] = len(rs); rec['rs_fail'] = [str(x) for x in fails]
    # consistency of rank_second's direct values with bilinearity in a (for the a evaluated directly, not in ob)
    lin_bad = 0; lin_n = 0
    if ok_basis:
        for ai in rs:
            if ai in ob or ai == 0: continue
            c = combo(ai)
            for k in range(len(sel)):
                e = 1
                for i, b in enumerate(ob):
                    if c >> i & 1: e *= rs[b][k]
                lin_n += 1; lin_bad += (e != rs[ai][k])
    rec['lin_n'] = lin_n; rec['lin_bad'] = lin_bad
    # compare formula (2) with rank_second on every (a, b), b over all of Sel^2
    cnt = dict(agree=0, disagree=0, undefined=0, agree_direct=0, disagree_direct=0)
    ex = []
    if ok_basis:
        for c in one:
            ai = idx[c]; dl = c[1]
            for bj, b in enumerate(allsel):
                r = rs_val(ai, bj)
                f = F.formula2(A, B, C, dl, b[0], b[1])
                direct = ai in rs and bin(bj).count('1') == 1
                if f is None: cnt['undefined'] += 1; key = None
                elif f == r: cnt['agree'] += 1; key = 'agree'
                else: cnt['disagree'] += 1; key = 'disagree'
                if direct and key: cnt[key + '_direct'] += 1
                if key == 'disagree' and len(ex) < 5: ex.append(dict(a=c, b=b, rs=r, f=f))
    rec['cmp'] = cnt; rec['ex'] = ex
    # Theta on S' x S' (where S' = S1)
    th = dict(agree=0, disagree=0)
    thex = []
    if ok_basis and rec['Sprime_eq']:
        for c in one:
            for d in one:
                r = rs_val(idx[c], idx[d])
                t = F.Theta(A, B, C, c[1], d[1])
                if t == r: th['agree'] += 1
                else:
                    th['disagree'] += 1
                    if len(thex) < 5: thex.append(dict(delta=c[1], beta=d[1], rs=r, theta=t))
    rec['theta'] = th; rec['theta_ex'] = thex
    rec['time'] = round(time.time() - t0, 2)
    return rec

if __name__ == '__main__':
    out, jobf = ARGS[0], ARGS[1]
    jobs = [l.split() for l in open(jobf) if l.strip() and not l.startswith('#')]
    i0, i1 = (int(ARGS[2]), int(ARGS[3])) if len(ARGS) > 3 else (0, len(jobs))
    done = set()
    if os.path.exists(out):
        for l in open(out):
            d = json.loads(l); done.add((d['m'], d['n'], d['cv']))
    for m, n, cv, cat in jobs[i0:i1]:
        m, n = int(m), int(n)
        if (m, n, cv) in done: continue
        try:
            r = run(m, n, cv, cat)
        except Exception as e:
            import traceback; traceback.print_exc()
            r = dict(m=m, n=n, cv=cv, cat=cat, error=repr(e))
        with open(out, 'a') as f: f.write(json.dumps(r) + '\n')
        print(m, n, cv, cat, 'dim', r.get('dim'), 'S1', r.get('dimS1'), 'eq', r.get('Sprime_eq'), r.get('cmp'), r.get('theta'), 'lin', r.get('lin_bad'), 't', r.get('time'), flush=True)
