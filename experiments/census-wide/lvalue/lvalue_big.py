#!/usr/bin/env python3
"""run5/f5 copy of scripts/g1_full/census_lvalue_big.py at 24d6ef6, changed only so that the ratio and the curve come from
the environment (MSS_PAIR=m,n and MSS_CURVE=name; default (700,529), w*f0*f8, the original's).  The original's text follows.

(700,529): L(E,1) > 0 without PARI on its closing curve w*f0*f8 (N = 86,826,771,232,530, about 9.2e7 terms), in
chunks, the way lball_exact53.py --apchunk / --bigmerge did its two large curves.  Written by an AI model; not read by
a person.

The curve is census_lvalue_exact.py's (pmul, poly, model are copied here so that this file needs only lball_exact.py
and lball_exact53.py; census_lvalue_exact.py --merge asserts that the record's model, curve and N are its own).
Certificate, as in census_lvalue_exact.py: N by Tate (lball_exact.local_data); a_p by lball_exact53.ap_one (Tate at bad
p, counting below 1000, the exact order certificate above, counting as fallback); a_n by lball_exact53.an_array; w = +1
from F = sum a_n exp(-2 pi n/sqrt N) != 0; L(E,1) by lball_exact.evaluate (directed rounding, tail |a_n| <= 2n); the
guard (enclosures at y = 1 and 6/5 meet, Fricke residual contains 0) asserted.

usage: census_lvalue_big.py --list                    N, w (local), terms
       census_lvalue_big.py --apchunk K J out.bin     a_p for every K-th prime from the J-th (resumes an interrupted run)
       census_lvalue_big.py --bigmerge K prefix       a_n and L(E,1) from prefix0.bin ... prefix{K-1}.bin -> prefix.jsonl
                                                       (a record for census_lvalue_exact.py --merge) and prefix_an.bin
       census_lvalue_big.py --bigcheck prefix         a check, not used: PARI's ellap at every prime, ellan(10^6),
                                                       ellglobalred and ellrootno against the merged a_n (added to prefix.jsonl)"""
import sys, os, json, math, time
HERE = os.path.join(os.environ.get('MSS_REPO', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))), 'scripts', 'g1_full')
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'second_descent'))
import lball_exact as LX
import lball_exact53 as L53

PAIR = tuple(int(x) for x in os.environ.get('MSS_PAIR', '700,529').split(',')); NAME = os.environ.get('MSS_CURVE', 'w*f0*f8')

def pmul(a, b):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i + j] += x * y
    return r

def poly(m, n, name):
    W = [0, 1]; f0 = [m*m, n*n]; f8 = [n*n, m*m]; S = m*m + n*n; g4 = [2*S, 2*S]
    K = [S*S, 2*S*S - 16*m*m*n*n, S*S]
    fac = {'w': W, 'f0': f0, 'f8': f8, 'g4': g4, 'K': K}
    P = [1]
    for x in name.split('*'): P = pmul(P, fac[x])
    return P

def model(P):
    if len(P) == 5: assert P[0] == 0; c = [P[4], P[3], P[2], P[1]]
    else: c = P
    c0, c1, c2, c3 = c
    return (0, c2, 0, c1 * c3, c0 * c3 * c3)

def the_curve():
    a = model(poly(PAIR[0], PAIR[1], NAME)); LD = LX.local_data(a)
    return LD['N'], a, LD

def fN(LD):
    return ' * '.join('%d^%d' % (p, T['f']) if T['f'] > 1 else str(p) for p, T in LD['rows'].items() if T['f'])

def apchunk(K, J, out):
    from array import array
    N, a, LD = the_curve(); M = LX.max_terms(N)
    mine = L53.primes_upto(M)[J::K]
    done = os.path.getsize(out) // 4 if os.path.exists(out) else 0
    st = L53.new_stats(); t0 = time.time(); fh = open(out, 'ab'); buf = array('i')
    for i in range(done, len(mine)):
        buf.append(L53.ap_one(a, LD, mine[i], st))
        if len(buf) == 20000 or i == len(mine) - 1:
            fh.write(buf.tobytes()); fh.flush(); buf = array('i')
            sys.stderr.write('chunk %d/%d of (%d,%d): %d of %d primes, p = %d, %.0f s\n' % (J, K, PAIR[0], PAIR[1], i + 1, len(mine), mine[i], time.time() - t0))
    fh.close()
    json.dump(dict(st, time=round(time.time() - t0, 1), n=len(mine), resumed_at=done), open(out + '.json', 'w'))

def bigmerge(K, prefix):
    from array import array
    N, a, LD = the_curve(); M = LX.max_terms(N); t0 = time.time()
    ps = L53.primes_upto(M); ap_of = array('q', bytes(8*(M + 1))); st = L53.new_stats(); tch = []
    for J in range(K):
        v = array('i'); v.frombytes(open('%s%d.bin' % (prefix, J), 'rb').read())
        mine = ps[J::K]; assert len(v) == len(mine), ('chunk incomplete', J, len(v), len(mine))
        for p, x in zip(mine, v):
            assert x*x <= 4*p; ap_of[p] = x
        s = json.load(open('%s%d.bin.json' % (prefix, J))); tch.append(s['time'])
        for k in st: st[k] += s[k]
    for p, T in LD['rows'].items():                         # bad primes again from Tate's algorithm
        if T['red'] != 'good' and p <= M: assert ap_of[p] == T['ap']
    del ps
    tap = time.time() - t0
    an = L53.an_array(ap_of, LD, M); tan = time.time() - t0 - tap
    PI = LX.pi_interval(); t1 = time.time()
    ev = LX.evaluate(an, N, 1, PI, M); tser = time.time() - t1
    F = ev['F']; wF = 1 if F[0] > 0 or F[1] < 0 else None
    L1 = ev['L_y1']; L2 = ev['L_y2']
    r = dict(idx=None, pair=list(PAIR), curve=NAME, model=list(a), N=N, fN=fN(LD), w_local=LD['w_local'], w_fricke=wF,
             F=[LX.lower_str(F[0], 3), LX.upper_str(F[1], 3)], terms=list(ev['M']), Mmax=M,
             L_y1=[L53.dec(L1[0]), L53.dec(L1[1], up=True)], L_y2=[L53.dec(L2[0]), L53.dec(L2[1], up=True)],
             guard=ev['guard'], positive=(L1[0] > 0 and wF == 1), ap_stats=st,
             time=dict(ap='%d chunks (%s s)' % (K, ', '.join('%.0f' % x for x in tch)), merge=round(tap, 1), an=round(tan, 1),
                       series=round(tser, 1), total=round(time.time() - t0, 1)))
    print(json.dumps(r), flush=True)
    open('%s.jsonl' % prefix, 'w').write(json.dumps(r) + '\n')
    assert r['positive'] and r['guard'], 'NOT CERTIFIED'
    open('%s_an.bin' % prefix, 'wb').write(an.tobytes())

def bigcheck(prefix):
    from array import array
    pari = L53.get_pari()
    N, a, LD = the_curve(); M = LX.max_terms(N)
    an = array('q'); an.frombytes(open('%s_an.bin' % prefix, 'rb').read()); assert len(an) == M + 1
    E = pari.ellminimalmodel(pari.ellinit(list(a)))
    Np = int(pari.ellglobalred(E)[0])
    ps = L53.primes_upto(M); bad = 0; t0 = time.time()
    for p in ps: bad += int(pari.ellap(E, p)) != an[p]
    m = min(M, 10**6); pa = [int(x) for x in pari.ellan(E, m)]
    same = pa == list(an[1:m + 1]); wp = int(pari.ellrootno(E))
    print('(%d,%d): PARI conductor %s; ellap at %d primes p <= %d: %d disagreements; ellan for n <= %d: %s; root number %d; %.0f s' % (
        PAIR[0], PAIR[1], 'agrees' if Np == N else 'DISAGREES', len(ps), ps[-1], bad, m, 'agree' if same else 'DISAGREE', wp, time.time() - t0))
    r = json.loads(open('%s.jsonl' % prefix).read())
    r['pari'] = dict(N=Np == N, w=wp == r['w_fricke'], an=bad == 0 and same, n=M, primes=len(ps))
    open('%s.jsonl' % prefix, 'w').write(json.dumps(r) + '\n')
    assert Np == N and bad == 0 and same and wp == r['w_fricke']

if __name__ == '__main__':
    args = sys.argv[1:]
    if '--list' in args:
        N, a, LD = the_curve(); M = LX.max_terms(N)
        print(PAIR, NAME, 'model', a, 'N', N, '=', fN(LD), 'w_local', LD['w_local'], 'terms (y=1)', LX.terms_needed(N, 1, True), 'Mmax', M)
    elif '--apchunk' in args:
        i = args.index('--apchunk'); apchunk(int(args[i + 1]), int(args[i + 2]), args[i + 3])
    elif '--bigmerge' in args:
        i = args.index('--bigmerge'); bigmerge(int(args[i + 1]), args[i + 2])
    elif '--bigcheck' in args:
        i = args.index('--bigcheck'); bigcheck(args[i + 1])
    else: print(__doc__)
