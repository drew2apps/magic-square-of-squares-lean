#!/usr/bin/env python3
"""F5: L(E,1) > 0 for a closing curve too large for census_lvalue_big.py's in-memory arrays, by streaming a_n.
Written by an AI model (Claude); not read by a person.  Experiment, not counted.

The arithmetic of the L-value is the repository's, unchanged in substance: lball_exact.py at 24d6ef6 (PREC = 192 bit
fixed point, directed rounding; q-intervals from pi_interval, rate, exp_neg; the tail bound of series; the guard of
evaluate: the enclosures at y = 1 and y = 6/5 meet and the Fricke residual contains 0).  series() is restated here so
that it can run over a_n in segments; the per-term operations are the same lines.  N, the bad a_p and w_local come from
lball_exact.local_data (Tate's algorithm, exact).

What differs from census_lvalue_big.py:
  * a_n is made in segments of 2^22 by a numpy sieve (exact 64-bit integers; |a_n| <= d(n) sqrt(n) < 2^40 here):
    smallest-prime-power factors from the primes <= sqrt(M), a_{p^e} by the recursion a_{p^(e+1)} = a_p a_{p^e} -
    p a_{p^(e-1)} (good p) or a_p^e (bad p), and the cofactor above sqrt(M), which is prime, looked up in the a_p table;
  * a_p at good p is PARI's ellap on PARI's minimal model (so the result is CONDITIONAL on PARI's ellap), with the
    Hasse bound asserted.  --apcheck compares a sample of these with the repository's own order certificate
    (lball_exact53.ap_one) and every a_n up to 2*10^5 with PARI's ellan.

    MSS_PAIR=m,n MSS_CURVE=name python3 lvalue_stream.py --primes DIR            primes <= Mmax (numpy sieve)
                                               --ap K J DIR            a_p for primes[J::K]
                                               --series A|B|C DIR      the series of one group
                                               --merge DIR [--write]   the record (DIR/record.json with --write)
                                               --apcheck DIR           cross-checks (not part of the certificate)
                                               --negative DIR          the merge with w forged to -1 must fail the guard
Groups: A = (y = 1, a_n/n, M1) and (y = 1, a_n, MF); B = (6/5, a_n/n, Mb) and (5/6, a_n/n, Mb);
C = (6/5, a_n, MR) and (5/6, a_n, MR): the six series of lball_exact.evaluate."""
import sys, os, json, math, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lvalue_big as LB
LX, L53 = LB.LX, LB.L53
PREC, ONE = LX.PREC, LX.ONE
SEG = int(os.environ.get('MSS_SEG', 1 << 22))       # segment length; a test sets it small

TEST = os.environ.get('MSS_STREAM_TEST') == '1'     # the repository's own curve (lball_exact.CURVE): a positive control

def the_curve():
    if TEST: a = tuple(LX.CURVE); LD = LX.local_data(a); return LD['N'], a, LD
    return LB.the_curve()

def setup():
    N, a, LD = the_curve(); M = LX.max_terms(N)
    yn, yd = LX.Y2
    Ms = dict(M1=LX.terms_needed(N, 1, True), Mb=LX.terms_needed(N, yd / yn, True),
              MF=LX.terms_needed(N, 1, False, 0.1), MR=LX.terms_needed(N, yd / yn, False, 0.1))
    return N, a, LD, M, Ms

def sieve(M):
    s = np.ones(M // 2 + 1, dtype=bool); s[0] = False             # s[i] <-> 2i + 1
    for p in range(3, math.isqrt(M) + 1, 2):
        if s[p // 2]: s[p * p // 2::p] = False
    odd = (2 * np.nonzero(s)[0] + 1).astype(np.int64)
    odd = odd[odd <= M]
    return np.concatenate([np.array([2], dtype=np.int64), odd]).astype(np.int32)

def cmd_primes(d):
    N, a, LD, M, Ms = setup(); t = time.time()
    P = sieve(M); np.save(os.path.join(d, 'primes.npy'), P)
    json.dump(dict(M=M, count=int(len(P)), last=int(P[-1]), sec=round(time.time() - t, 1), **Ms), open(os.path.join(d, 'primes.json'), 'w'))

def cmd_ap(K, J, d):
    import cypari2
    pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9)
    N, a, LD, M, Ms = setup()
    Em = pari.ellinit(pari.ellminimalmodel(pari.ellinit(list(a))))
    assert int(pari.ellglobalred(Em)[0]) == N, 'PARI conductor differs from Tate'
    P = np.load(os.path.join(d, 'primes.npy'), mmap_mode='r')[J::K]
    bad = {p: T['ap'] for p, T in LD['rows'].items() if T['red'] != 'good'}
    out = np.zeros(len(P), dtype=np.int32); t = time.time()
    for i, p in enumerate(P.tolist()):
        x = bad[p] if p in bad else int(pari.ellap(Em, p))
        assert x * x <= 4 * p, ('Hasse', p, x)
        out[i] = x
        if i % 2000000 == 0: sys.stderr.write('ap %d/%d: %d of %d, p = %d, %.0f s\n' % (J, K, i, len(P), p, time.time() - t))
    np.save(os.path.join(d, 'ap_%d_%d.npy' % (K, J)), out)
    json.dump(dict(n=len(P), sec=round(time.time() - t, 1)), open(os.path.join(d, 'ap_%d_%d.json' % (K, J)), 'w'))

def load_ap(d):
    P = np.load(os.path.join(d, 'primes.npy'), mmap_mode='r')
    fs = sorted(f for f in os.listdir(d) if f.startswith('ap_') and f.endswith('.npy'))
    K = int(fs[0].split('_')[1]); assert len(fs) == K, 'missing a_p chunks'
    AP = np.zeros(len(P), dtype=np.int32)
    for J in range(K): AP[J::K] = np.load(os.path.join(d, 'ap_%d_%d.npy' % (K, J)))
    return np.asarray(P), AP

def an_segments(P, AP, LD, Mtop):
    """yield (start, a_start..a_end) as Python lists, for start = 1, 1 + SEG, ...; exact"""
    bad = {p for p, T in LD['rows'].items() if T['red'] != 'good'}
    r = math.isqrt(Mtop)
    small = P[:np.searchsorted(P, r, side='right')].astype(np.int64)
    tab = []
    for i, p in enumerate(small.tolist()):
        ap = int(AP[i]); e_max = 1; q = p
        while q * p <= Mtop: q *= p; e_max += 1
        t = [1, ap]
        for e in range(2, e_max + 1): t.append(ap * t[-1] if p in bad else ap * t[-1] - p * t[-2])
        tab.append(np.array(t, dtype=np.int64))
    L = 1
    while L <= Mtop:
        R = min(L + SEG, Mtop + 1)
        n = np.arange(L, R, dtype=np.int64); rem = n.copy(); an = np.ones(R - L, dtype=np.int64)
        for i, p in enumerate(small.tolist()):
            s0 = (-L) % p
            if s0 >= R - L: continue
            idx = np.arange(s0, R - L, p)
            e = np.zeros(len(idx), dtype=np.int64); sub = rem[idx]
            while True:
                mk = sub % p == 0
                if not mk.any(): break
                sub = np.where(mk, sub // p, sub); e += mk
            rem[idx] = sub; an[idx] *= tab[i][e]
        big = rem > 1
        if big.any():
            q = rem[big]; j = np.searchsorted(P, q)
            assert np.all(P[np.minimum(j, len(P) - 1)] == q), 'cofactor above sqrt(M) not in the prime table'
            an[big] *= AP[j].astype(np.int64)
        assert np.abs(an).max() < 1 << 50
        yield L, an.tolist()
        L = R

def series_group(specs, P, AP, LD):
    """specs: list of (q, divide, M); the same per-term lines as lball_exact.series, run side by side"""
    Mtop = max(M for _, _, M in specs)
    st = [[q[0], q[1], ONE, ONE, 0, 0, divide, M] for q, divide, M in specs]   # qlo, qhi, plo, phi, lo, hi, divide, M
    t = time.time()
    for L, seg in an_segments(P, AP, LD, Mtop):
        for s in st:
            qlo, qhi, plo, phi, lo, hi, divide, M = s
            if L > M: continue
            k = min(len(seg), M - L + 1); n = L - 1
            for c in seg[:k]:
                n += 1
                plo = (plo * qlo) >> PREC; phi = -((-(phi * qhi)) >> PREC)
                if c == 0: continue
                d = n if divide else 1
                if c > 0: lo += (c * plo) // d; hi += -((-(c * phi)) // d)
                else: lo += (c * phi) // d; hi += -((-(c * plo)) // d)
            s[2], s[3], s[4], s[5] = plo, phi, lo, hi
        if (L // SEG) % 16 == 0: sys.stderr.write('series: n = %d of %d, %.0f s\n' % (L, Mtop, time.time() - t))
    res = []
    for qlo, qhi, plo, phi, lo, hi, divide, M in st:
        qM1 = LX.mul_hi(phi, qhi); om = ONE - qhi; assert om > 0
        tail = LX.cdiv(2 * qM1 * ONE, om) if divide else LX.cdiv(2 * qM1 * ((M + 1) * ONE - M * qlo) * ONE, om * om)
        res.append((lo - tail, hi + tail))
    return res, round(time.time() - t, 1)

def qs(N):
    PI = LX.pi_interval(); yn, yd = LX.Y2
    return LX.q_interval(N, 1, 1, PI), LX.q_interval(N, yn, yd, PI), LX.q_interval(N, yd, yn, PI)

def cmd_series(g, d):
    N, a, LD, M, Ms = setup(); P, AP = load_ap(d); q1, qa, qb = qs(N)
    specs = dict(A=[(q1, True, Ms['M1']), (q1, False, Ms['MF'])], B=[(qa, True, Ms['Mb']), (qb, True, Ms['Mb'])],
                 C=[(qa, False, Ms['MR']), (qb, False, Ms['MR'])])[g]
    res, sec = series_group(specs, P, AP, LD)
    json.dump(dict(group=g, M=[M_ for _, _, M_ in specs], sums=[[str(x), str(y)] for x, y in res], sec=sec),
              open(os.path.join(d, 'series_%s.json' % g), 'w'))

def combine(d, w=1):
    N, a, LD, M, Ms = setup()
    S = {}
    for g in 'ABC':
        r = json.load(open(os.path.join(d, 'series_%s.json' % g)))
        S[g] = [(int(x), int(y)) for x, y in r['sums']]
    S1, F = S['A']; Sa, Sb = S['B']; Fa, Fb = S['C']
    yn, yd = LX.Y2; k, kd = yn * yn, yd * yd
    out = {}
    out['L_y1'] = (S1[0] + w * S1[0], S1[1] + w * S1[1]) if w == 1 else (S1[0] - S1[1], S1[1] - S1[0])
    out['L_y2'] = (Sa[0] + (Sb[0] if w == 1 else -Sb[1]), Sa[1] + (Sb[1] if w == 1 else -Sb[0]))
    out['F'] = F
    if w == 1: R = (Fb[0] - LX.cdiv(k * Fa[1], kd), Fb[1] - (k * Fa[0]) // kd)
    else: R = (Fb[0] + (k * Fa[0]) // kd, Fb[1] + LX.cdiv(k * Fa[1], kd))
    out['R'] = R
    out['overlap'] = max(out['L_y1'][0], out['L_y2'][0]) <= min(out['L_y1'][1], out['L_y2'][1])
    out['R_has_0'] = R[0] <= 0 <= R[1]
    out['guard'] = out['overlap'] and out['R_has_0']
    return N, a, LD, M, Ms, out

def cmd_merge(d, write):
    N, a, LD, M, Ms, ev = combine(d)
    F = ev['F']; wF = 1 if F[0] > 0 or F[1] < 0 else None
    L1, L2 = ev['L_y1'], ev['L_y2']
    r = dict(pair=None if TEST else list(LB.PAIR), curve='lball_exact.CURVE' if TEST else LB.NAME, model=list(a), N=N, fN=LB.fN(LD), w_local=LD['w_local'], w_fricke=wF,
             F=[LX.lower_str(F[0], 3), LX.upper_str(F[1], 3)], terms=[Ms['M1'], Ms['Mb'], Ms['MF'], Ms['MR']], Mmax=M,
             L_y1=[L53.dec(L1[0]), L53.dec(L1[1], up=True)], L_y2=[L53.dec(L2[0]), L53.dec(L2[1], up=True)],
             guard=ev['guard'], positive=(L1[0] > 0 and wF == 1), ap='PARI ellap on the minimal model (conditional)')
    txt = json.dumps(r) + '\n'; print(txt, end='')
    if write: open(os.path.join(d, 'record.json'), 'w').write(txt)
    else:
        f = os.path.join(d, 'record.json')
        if os.path.exists(f): assert open(f).read() == txt, 'differs from record.json'; print('# matches record.json')
    assert r['positive'] and r['guard'], 'NOT CERTIFIED'

def cmd_negative(d, write):
    lines = []
    N, a, LD, M, Ms, ev = combine(d, w=-1)
    lines.append('the same six sums with w forged to -1: guard %s (overlap %s, Fricke residual contains 0: %s) -> %s'
                 % (ev['guard'], ev['overlap'], ev['R_has_0'], 'refused, as it must be' if not ev['guard'] else 'NOT REFUSED'))
    ok = not ev['guard']
    S = json.load(open(os.path.join(d, 'series_A.json')))
    lo, hi = int(S['sums'][0][0]), int(S['sums'][0][1])
    lines.append('the y = 1 sum S1 (L = 2 S1 when w = +1) lies in [%s, %s]: its width is small against its value' % (LX.lower_str(lo, 12), LX.upper_str(hi, 12)))
    lines.append('negative checks: %s' % ('as required' if ok else 'FAILED'))
    txt = '\n'.join(lines) + '\n'; print(txt, end='')
    if write: open(os.path.join(d, 'negative.txt'), 'w').write(txt)

def cmd_negative_small(write):
    """the repository's forged-N controls (lball_exact.negative), run through this file's streamed series on the
    repository's own curve: each forged conductor must fail the guard; the true one must pass"""
    assert TEST, 'run with MSS_STREAM_TEST=1'
    N, a, LD, M, Ms = setup(); PI = LX.pi_interval(); yn, yd = LX.Y2
    runs = []
    for name, Nf, tl in (('true N', N, LX.TARGET), ('N/2', N // 2, LX.NEG_TL), ('2N', 2 * N, LX.NEG_TL)):
        M1 = LX.terms_needed(Nf, 1, True, tl); Mb = LX.terms_needed(Nf, yd / yn, True, tl)
        MR = LX.terms_needed(Nf, yd / yn, False, 0.1)
        q1 = LX.q_interval(Nf, 1, 1, PI); qa = LX.q_interval(Nf, yn, yd, PI); qb = LX.q_interval(Nf, yd, yn, PI)
        runs.append((name, [(q1, True, M1), (qa, True, Mb), (qb, True, Mb), (qa, False, MR), (qb, False, MR)]))
    P = sieve(max(Mx for _, sp in runs for _, _, Mx in sp))
    ap = LX.ap_table(a, LD, int(P[-1])); AP = np.array([ap[int(p)] for p in P], dtype=np.int32)
    lines = []; ok = True
    for name, sp in runs:
        (S1, Sa, Sb, Fa, Fb), _ = series_group(sp, P, AP, LD)
        L1 = (2 * S1[0], 2 * S1[1]); L2 = (Sa[0] + Sb[0], Sa[1] + Sb[1]); k, kd = yn * yn, yd * yd
        R = (Fb[0] - LX.cdiv(k * Fa[1], kd), Fb[1] - (k * Fa[0]) // kd)
        guard = max(L1[0], L2[0]) <= min(L1[1], L2[1]) and R[0] <= 0 <= R[1]
        want = name == 'true N'; ok &= (guard == want)
        lines.append('%-7s: guard %s (%s)' % (name, guard, 'as required' if guard == want else 'WRONG'))
    lines.append('forged-N controls through the streamed series: %s' % ('as required' if ok else 'FAILED'))
    txt = '\n'.join(lines) + '\n'; print(txt, end='')
    if write: open(os.path.join(HERE, 'test_lball', 'negative_small.txt'), 'w').write(txt)
    assert ok

def cmd_apcheck(d, write):
    import random, cypari2
    pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9)
    N, a, LD, M, Ms = setup(); P, AP = load_ap(d)
    lines = []
    rng = random.Random(1); idx = sorted(rng.sample(range(len(P)), 20000)); st = L53.new_stats(); bad = 0; t = time.time()
    for i in idx:
        p = int(P[i]); x = L53.ap_one(a, LD, p, st)
        bad += (x != int(AP[i]))
    lines.append('a_p at 20000 random primes <= %d: the repository\'s own order certificate (lball_exact53.ap_one) and PARI agree at %d of 20000 (%s)'
                 % (M, 20000 - bad, ', '.join('%s %d' % kv for kv in sorted(st.items()))))
    sys.stderr.write('[apcheck %.0f s]\n' % (time.time() - t))
    T = 200000
    E = pari.ellinit(list(a)); v = [int(x) for x in pari.ellan(E, T)]
    seg = next(an_segments(P, AP, LD, M))[1][:T]
    lines.append('a_n for n <= %d: streamed a_n and PARI ellan agree: %s' % (T, seg == v))
    txt = '\n'.join(lines) + '\n'; print(txt, end='')
    if write: open(os.path.join(d, 'apcheck.txt'), 'w').write(txt)

if __name__ == '__main__':
    A = sys.argv[1:]; W = '--write' in A
    if A[0] == '--primes': cmd_primes(A[1])
    elif A[0] == '--ap': cmd_ap(int(A[1]), int(A[2]), A[3])
    elif A[0] == '--series': cmd_series(A[1], A[2])
    elif A[0] == '--merge': cmd_merge(A[1], W)
    elif A[0] == '--negative': cmd_negative(A[1], W)
    elif A[0] == '--apcheck': cmd_apcheck(A[1], W)
    elif A[0] == '--negative-small': cmd_negative_small(W)
