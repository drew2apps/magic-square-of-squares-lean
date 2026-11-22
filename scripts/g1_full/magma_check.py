#!/usr/bin/env python3
"""Third software for the upper rank bounds: Magma's TwoPowerIsogenyDescentRankBound (descent by
2-isogenies, T. Fisher), through the public Magma calculator (http://magma.maths.usyd.edu.au/calc/,
60 seconds per request). Every curve here has a rational point of order 2.

Without arguments: every line of dm_evidence.txt (curve E04, D1 or E0L of the pair the line was run on);
prints 'of run curve rank <= r' and a summary. With arguments 'CURVE m,n ...': those curves only.
Options: --batch=K curves per request (default 3), --pause=S seconds between requests (default 45; the
calculator refuses clients that connect too often), --skip=FILE to reuse answers from an earlier output.
The E0L curve is built in Magma from the quartic y^2 = F0(t) GL(t) and its point t = m/n.
This sends the curve coefficients to a server of the University of Sydney. No dependencies."""
import os, re, sys, time, urllib.parse, urllib.request

URL = 'http://magma.maths.usyd.edu.au/xml/calculator.xml'

def code(kind, m, n, tag):
    pre = 'm:=%d; n:=%d; S:=m^2+n^2; ' % (m, n)
    if kind == 'E04':
        c = 'E:=EllipticCurve([0,2*S^2,0,4*S^2*m^2*n^2,0]); '
    elif kind == 'D1':
        c = ('D2:=(m^2-n^2)^2; c3:=4*S^2*m^2*n^2; c2:=-16*S*m^3*n^3; c1:=4*S^2*D2; c0:=-16*S*m*n*D2; '
             'E:=EllipticCurve([0,c1,0,c0*c2,c0^2*c3]); ')
    else:
        c = ('P<t>:=PolynomialRing(Rationals()); C:=HyperellipticCurve((n^2*t^2+m^2)*(2*S*t^2-8*m*n*t+2*S)); '
             'E:=EllipticCurve(C, C![m/n, 2*m*(m^2-n^2)/n]); ')
    return pre + c + ('E:=MinimalModel(E); G,mp:=TorsionSubgroup(E); T:=mp((Order(G.1) div 2)*G.1); '
                      'hi:=TwoPowerIsogenyDescentRankBound(E,T); print "BOUNDS", %d, hi;\n' % tag)

def ask(batch):
    data = urllib.parse.urlencode({'input': ''.join(code(k, m, n, i) for i, (k, m, n) in enumerate(batch))}).encode()
    try: out = urllib.request.urlopen(URL, data, timeout=150).read().decode()
    except Exception as ex: out = 'ERROR %s' % ex
    got = {int(a): int(b) for a, b in re.findall(r'BOUNDS (\d+) (\d+)', out)}
    ver = re.search(r'<version>([^<]*)</version>', out)
    return [got.get(i) for i in range(len(batch))], (ver.group(1) if ver else None), ('offline' in out)

if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    opts = dict(a[2:].split('=', 1) for a in sys.argv[1:] if a.startswith('--'))
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    K = int(opts.get('batch', 3)); pause = float(opts.get('pause', 45))
    if args: jobs = [('-', a, args[0]) for a in args[1:]]
    else: jobs = [tuple(l.split()[:3]) for l in open(os.path.join(here, 'dm_evidence.txt')) if l.strip() and not l.startswith('#')]
    seen = {}
    if 'skip' in opts and os.path.exists(opts['skip']):
        for l in open(opts['skip']):
            mt = re.match(r'\S+ (\S+) (\S+) rank <= (\d+)', l)
            if mt: seen[(mt.group(2), mt.group(1))] = int(mt.group(3))
    todo = []
    for of, run, kind in jobs:
        if (kind, run) not in seen and (kind, run) not in todo: todo.append((kind, run))
    version = None
    for i in range(0, len(todo), K):
        part = todo[i:i + K]
        for attempt in range(4):
            res, ver, offline = ask([(k,) + tuple(int(x) for x in r.split(',')) for k, r in part])
            version = ver or version
            if not offline: break
            time.sleep(600)                      # refused: wait ten minutes
        for key, b in zip(part, res): seen[key] = b
        time.sleep(pause)
    ok = 0; bad = []
    for of, run, kind in jobs:
        b = seen.get((kind, run))
        print(of, run, kind, ('rank <= %d' % b) if b is not None else 'no answer', flush=True)
        if b == 1: ok += 1
        else: bad.append((of, run, kind, b))
    print('Magma %s: %d lines, rank <= 1 in %d; other: %s' % (version, len(jobs), ok, bad))
