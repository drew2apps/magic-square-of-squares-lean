#!/usr/bin/env python3
"""usage: python3 magma_sigma.py magma_sigma.json [pause]  (writes magma_sigma_raw.txt)
Experiment (review of 7b4ff7a): rank bounds from the public Magma calculator for the six sigma-quotient curves
whose PARI bounds are [1, 3]. One curve per request, a long pause between requests (as magma_g2.py does)."""
import json, re, sys, time, urllib.parse, urllib.request
URL = 'http://magma.maths.usyd.edu.au/xml/calculator.xml'
curves = json.load(open(sys.argv[1]))
pause = float(sys.argv[2]) if len(sys.argv) > 2 else 120
raw = open(sys.argv[1].replace('.json', '_raw.txt'), 'a')
for c in curves:
    code = ('E:=EllipticCurve(%s); r1,r2:=RankBounds(E); print "RB", r1, r2; '
            'if r1 lt r2 then a,b:=MordellWeilShaInformation(E : Silent:=true); print "MWSI", a, b; end if;') % str(c['ainvs'])
    for attempt in range(6):
        try: out = urllib.request.urlopen(URL, urllib.parse.urlencode({'input': code}).encode(), timeout=150).read().decode()
        except Exception as ex: out = 'ERROR %s' % ex
        if 'offline' not in out: break
        time.sleep(1800)
    raw.write('=== %s %s\n--- input\n%s\n--- response\n%s\n' % (c['pair'], c['curve'], code, out)); raw.flush()
    lines = re.findall(r'<line>([^<]*)</line>', out)
    print(c['pair'], c['curve'], '|', ' '.join(lines)[:600] if lines else out[:300].replace('\n', ' '), flush=True)
    time.sleep(pause)
