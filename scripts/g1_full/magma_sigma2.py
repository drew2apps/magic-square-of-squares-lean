#!/usr/bin/env python3
"""usage: python3 magma_sigma2.py magma_sigma2.json [pause]   (appends to magma_sigma2_raw.txt)
Third software for the upper rank bounds of the curves of F0 F8 GL (G1_FULL_SYSTEM.md, Result 4): Magma's
TwoPowerIsogenyDescentRankBound (descent by 2-isogenies, T. Fisher) through the public Magma calculator, as
magma_check.py does for the curves of dm_evidence.txt. Every curve here has a rational point of order 2.
Two curves per request, a pause between requests. Sends the coefficients to a server of the University of Sydney."""
import json, re, sys, time, urllib.parse, urllib.request
URL = 'http://magma.maths.usyd.edu.au/xml/calculator.xml'
curves = json.load(open(sys.argv[1])); pause = float(sys.argv[2]) if len(sys.argv) > 2 else 45
raw = open(sys.argv[1].replace('.json', '_raw.txt'), 'a')
for i in range(0, len(curves), 2):
    part = curves[i:i + 2]
    code = ''.join('E:=MinimalModel(EllipticCurve(%s)); G,mp:=TorsionSubgroup(E); T:=mp((Order(G.1) div 2)*G.1); '
                   'hi:=TwoPowerIsogenyDescentRankBound(E,T); print "BOUNDS", %d, %d, "%s", hi; '
                   % (str(c['ainvs']), c['pair'][0], c['pair'][1], c['curve']) for c in part)
    for attempt in range(6):
        try: out = urllib.request.urlopen(URL, urllib.parse.urlencode({'input': code}).encode(), timeout=150).read().decode()
        except Exception as ex: out = 'ERROR %s' % ex
        if 'offline' not in out: break
        time.sleep(1800)
    raw.write('--- input\n%s\n--- response\n%s\n' % (code, out)); raw.flush()
    lines = re.findall(r'<line>([^<]*)</line>', out)
    print(' | '.join(l for l in lines if l.strip()) if lines else out[:300].replace('\n', ' '), flush=True)
    time.sleep(pause)
