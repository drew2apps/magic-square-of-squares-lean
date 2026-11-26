#!/usr/bin/env python3
"""Which census pairs can gen_census_closure.py reach?  For each of the 425 exceptions with m <= 1000
(MagicSquareOfSquares/G1Census1000.lean, `exceptions1000`): PARI's rank bounds for the curve
w·f0·f8 (Y^2 = X(X + m^4)(X + n^4), and its isogeny class, as in scripts/g1_full/global.py), and,
where the upper bound is 0, whether the run-3 local test (gen_rank_zero_certs.cert) leaves at most
4 classes, within a time cap.  Guidance only (PARI decides nothing that Lean relies on).

Written by an AI model; not read by a person.
usage: python3 -I census_closure_scan.py CAP_SECONDS PATH/G1Census1000.lean [MMAX] > census_closure_scan.txt
"""
import os
import re
import signal
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gen_census_closure as GC  # noqa: E402
import cypari2  # noqa: E402

pari = cypari2.Pari()
pari.allocatemem(10 ** 9)


def exceptions1000(path):
    src = open(path).read()
    line = src[src.index("def exceptions1000"):].split("\n", 1)[0]
    pairs = [(int(a), int(b)) for a, b in re.findall(r"\((\d+), (\d+)\)", line)]
    assert len(pairs) == 425
    return pairs


def rank_hi(m, n):
    E = pari.ellinit([0, m ** 4 + n ** 4, 0, m ** 4 * n ** 4, 0])
    rk = pari.ellrank(E)
    lo, hi = int(rk[0]), int(rk[1])
    if lo == 0 and hi > 0:
        for item in pari.ellisomat(pari.ellminimalmodel(E), 0, 1)[0]:
            hi = min(hi, int(pari.ellrank(pari.ellinit(item))[1]))
            if hi == 0:
                break
    return lo, hi


class Cap(Exception):
    pass


def _alarm(sig, frm):
    raise Cap()


def main():
    cap = int(sys.argv[1])
    mmax = int(sys.argv[3]) if len(sys.argv) > 3 else 1000
    signal.signal(signal.SIGALRM, _alarm)
    pairs = [pq for pq in exceptions1000(sys.argv[2]) if pq[0] <= mmax]
    n0 = nacc = 0
    for (m, n) in pairs:
        lo, hi = rank_hi(m, n)
        if hi != 0:
            print(f"{m} {n} rank [{lo},{hi}]", flush=True)
            continue
        n0 += 1
        t = time.time()
        signal.alarm(cap)
        try:
            c = GC.G3.cert(m, n)
            verdict = "ACCEPT" if c is not None else "refused"
        except Cap:
            verdict = "timeout"
        finally:
            signal.alarm(0)
        dt = time.time() - t
        extra = ""
        if verdict == "ACCEPT":
            nacc += 1
            _, info = GC.section(m, n)
            extra = f" points {info['N']} p {info['p']}"
        print(f"{m} {n} rank [0,0] cert {verdict} {dt:.1f}s{extra}", flush=True)
    print(f"m <= {mmax}: rank 0 on w*f0*f8: {n0} of {len(pairs)}; "
          f"accepted by the run-3 local test within {cap}s: {nacc}")


if __name__ == "__main__":
    main()
