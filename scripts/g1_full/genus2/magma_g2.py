#!/usr/bin/env python3
"""Experiment: the genus-2 quotient H: y^2 = F0 G4 GL (up to squares, y^2 = (n^2 t^2 + m^2)(t^2 + 1)(S(t^2 + 1) - 4mnt))
of the G1 curve, for the open pairs. Asks the public Magma calculator for RankBound(Jacobian(H)); when the bound
is at most 1 also for the rational points of small height. One pair per request, a long pause between requests.
A pair and its partner give the same H, so only one of them is asked: 23 requests for the 36 open pairs
(13 couples of partners that are both open, and 10 pairs whose partner has m > 1000).
Output: one line per request on stdout (g2.txt), requests and verbatim responses in g2_raw.txt.
g2_grh.txt is the output of the same script with SetClassGroupBounds("GRH") put first in the Magma input."""
import re, sys, time, urllib.parse, urllib.request
URL = 'http://magma.maths.usyd.edu.au/xml/calculator.xml'
OPEN = "72,65 91,60 195,28 209,120 223,167 329,89 420,341 425,168 440,279 572,555 608,105 609,200 627,364 680,111 713,503 719,161 741,580 760,39 761,79 791,569 799,721 809,409 861,620 920,129 945,248 991,263 137,7 151,31 593,257 724,665 740,539 820,799 828,173 840,83 889,320 923,757".split()
def partner(m, n): return (m + n, m - n) if (m + n) % 2 else ((m + n) // 2, (m - n) // 2)
done = set(); pause = float(sys.argv[1]) if len(sys.argv) > 1 else 120
for pr in OPEN:
    m, n = (int(x) for x in pr.split(','))
    if (m, n) in done or partner(m, n) in done: continue
    done.add((m, n))
    code = ('P<t>:=PolynomialRing(Rationals()); m:=%d; n:=%d; S:=m^2+n^2; '
            'f:=(n^2*t^2+m^2)*(t^2+1)*(S*(t^2+1)-4*m*n*t); H:=HyperellipticCurve(f); J:=Jacobian(H); '
            'r:=RankBound(J); print "G2", m, n, r; '
            'if r le 1 then print "PTS", RationalPoints(H:Bound:=5000); end if;') % (m, n)
    for attempt in range(6):
        try: out = urllib.request.urlopen(URL, urllib.parse.urlencode({'input': code}).encode(), timeout=150).read().decode()
        except Exception as ex: out = 'ERROR %s' % ex
        if 'offline' not in out: break
        time.sleep(1800)
    open('g2_raw.txt','a').write('=== %s\n--- input\n%s\n--- response\n%s\n' % (pr, code, out))
    lines = re.findall(r'<line>([^<]*)</line>', out)
    print(pr, '|', ' '.join(lines)[:600] if lines else out[:200].replace('\n', ' '), flush=True)
    time.sleep(pause)
