#!/usr/bin/env python3
"""The route of LAST_PAIR.md on every ratio of square type whose curve E1 has rank exactly 1 (EXPERIMENT).

For each ratio of two_map_scope.txt with rank bounds [1, 1]: saturation.py (the degenerate point generates), mw_sieve.py with
four primes (the pairs (a, b) left, and the range covered), closure_estimate.py with that range (the bound for a, b from the
height bound, and whether the sieve covers it), fibre.py (the rational t with a = 1 or b = 1).  One line per ratio.
All of these ratios except (627,364) are closed by other means, so this is a test of the route, not a new result.
usage: route_all.py [jobs [scope.txt] [m,n ...]]   (scope.txt: an output of two_map_scope.py; default two_map_scope.txt)"""
import re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
PR = [10009, 20089, 30013, 40009]
def run(args):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    return [l for l in r.stdout.splitlines() if 'PARI stack' not in l], r.stderr
def isprime(q): return q > 1 and all(q % k for k in range(2, int(q**0.5) + 1))
def nxt(q):
    q += 4 - (q - 1) % 4 if (q - 1) % 4 else 4
    while not isprime(q): q += 4
    return q
def one(item):
    m, n, d = item; mn = '%d,%d' % (m, n); out = '(%d,%d), d = %d: ' % (m, n, d)
    sat, err = run(['saturation.py', mn, str(d)])
    gen = bool(sat) and sat[-1].endswith('True')
    pr = list(PR)
    for attempt in range(8):                              # a prime that is not admissible for this ratio is replaced by the next one
        sv, err = run(['mw_sieve.py', mn, str(d), '1000000'] + [str(q) for q in pr])
        bad = re.findall(r'the prime (\d+) is not allowed', err)
        if not bad: break
        pr[pr.index(int(bad[0]))] = nxt(int(bad[0]))
    if not sv or not sv[-1].startswith('left:'): return out + 'generator: %s; sieve failed' % gen
    parse = lambda sv: [tuple(int(x) for x in pq) for pq in re.findall(r'\((\d+), (\d+)\)', sv[-1])]
    left = parse(sv); extra = 50000
    while (0, 0) in left and len(pr) < 10:                # t = 0 and infinity can be points modulo all the primes used without being
        extra = nxt(extra)                                # points over K; then (0, 0) is left, and one more prime is added
        sv2, err = run(['mw_sieve.py', mn, str(d), '1000000'] + [str(q) for q in pr + [extra]])
        if sv2 and sv2[-1].startswith('left:'): pr.append(extra); sv = sv2; left = parse(sv)
    rng = [int(x) for l in sv for x in re.findall(r'^pairs 0 <= a, b <= (\d+)', l)][0]
    one_prime = any('first prime only' in l for l in sv)
    ce, err = run(['closure_estimate.py', mn, str(d), str(rng)])
    b96 = [int(x) for l in ce for x in re.findall(r'deg C <= 96: every solution has a, b <= (\d+)', l)]
    if not b96: return out + 'generator: %s; closure_estimate failed' % gen
    inb = [pq for pq in left if max(pq) <= b96[0]]
    fb, err = run(['fibre.py', mn, str(d)])
    ts = sorted(set(re.findall(r"'(-?[0-9/]+)': '", ' '.join(fb))))
    return out + 'generator: %s; primes %s; height bound: a, b <= %d; sieve to %d%s leaves %d pair(s), within the bound: %s; range covers the bound: %s; fibre over a = 1 or b = 1: rational t = %s' % (
        gen, ' '.join(str(q) for q in pr), b96[0], rng, ' (from one prime)' if one_prime else '', len(left), inb, b96[0] <= rng, ', '.join(ts))
scope = [a for a in sys.argv[2:] if a.endswith('.txt')] or ['two_map_scope.txt']; only = [a for a in sys.argv[2:] if not a.endswith('.txt')]
items = [tuple(int(x) for x in re.match(r'\((\d+),(\d+)\): d = (-?\d+)', l).groups()) for l in open(scope[0]) if 'rank of E1 in [1, 1]' in l]
if only: items = [it for it in items if '%d,%d' % it[:2] in only]
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 2) as ex:
    for line in ex.map(one, items): print(line, flush=True)
