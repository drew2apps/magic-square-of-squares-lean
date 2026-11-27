#!/usr/bin/env python3
"""F5 positive control for lvalue_stream.py: on the repository's own curve (lball_exact.CURVE), the streamed a_n and the
streamed series give exactly lball_exact.an_from_ap and lball_exact.evaluate's integer enclosures, with the default
segment length and with a short one (4093) that puts many segment boundaries inside the sums.
Written by an AI model (Claude); not read by a person.   python3 stream_selftest.py [--write]  (selftest.txt)"""
import os, sys, importlib
os.environ['MSS_STREAM_TEST'] = '1'
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
lines = []
for seg in (1 << 22, 4093):
    os.environ['MSS_SEG'] = str(seg)
    import lvalue_stream as S; S = importlib.reload(S); LX = S.LX
    N, a, LD, M, Ms = S.setup()
    P = S.sieve(M); ap = LX.ap_table(a, LD, M); AP = np.array([ap[int(p)] for p in P], dtype=np.int32)
    an = LX.an_from_ap(ap, LD, M); full = [0]
    for L, sg in S.an_segments(P, AP, LD, M): full += sg
    q1, qa, qb = S.qs(N)
    specs = [(q1, True, Ms['M1']), (q1, False, Ms['MF']), (qa, True, Ms['Mb']), (qb, True, Ms['Mb']), (qa, False, Ms['MR']), (qb, False, Ms['MR'])]
    got, _ = S.series_group(specs, P, AP, LD)
    ref = [LX.series(an, q, dv, Mx) for q, dv, Mx in specs]
    lines.append('segment %d: a_n identical to lball_exact.an_from_ap: %s; the six sums identical to lball_exact.series: %s'
                 % (seg, full == an, got == ref))
    ok = full == an and got == ref
    assert ok
txt = 'stream_selftest on lball_exact.CURVE, N = %d, M = %d\n' % (N, M) + '\n'.join(lines) + '\n'
sys.stdout.write(txt)
f = os.path.join(HERE, 'selftest.txt')
if '--write' in sys.argv: open(f, 'w').write(txt)
elif os.path.exists(f): assert open(f).read() == txt; print('# matches selftest.txt')
