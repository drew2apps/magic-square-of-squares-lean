#!/usr/bin/env python3
"""Magma's two-cover descent against the first descent, on subsets of three forms.

magma_twocover_raw.txt: verbatim requests to the public Magma calculator and its answers.  Each request runs
TwoCoverDescent on the genus 2 curve Z^2 = f_a f_b f_c for three of the five forms of a ratio, and counts the elements of
the fake 2-Selmer set at which each of the three norms is a square (the elements that a solution could give, since each form
is a square at a solution).  subsets.txt, subsets_34.txt: the number of tuples that the first descent leaves for the same three
forms, up to a common rational factor (the two-cover descent works modulo Q*).  The two numbers must be equal.

Every curve that was asked for is listed, read from the input of each request (the ratio and the list of triples), also
when the calculator printed nothing for it: a request that ran out of time or memory is reported as not finished, with the
number of attempts.  The last part lists, for the ratios of pairs34.txt, whether Magma itself found an empty set for some
triple; a ratio for which no request was finished is named.  Magma is not run here.
usage: twocover_check.py [raw file [subsets file ...]]"""
import sys, re
NAMES = ['F0', 'F8', 'G4', 'GL', 'GU']
raw = open(sys.argv[1] if len(sys.argv) > 1 else 'magma_twocover_raw.txt').read()
ours = {}
for l in [x for f in (sys.argv[2:] or ['subsets.txt', 'subsets_34.txt']) for x in open(f)]:
    pair, rest = l.split(':', 1); m, n = (int(v) for v in pair.strip('() ').split(','))
    for item in rest.split(): k, v = item.split('='); ours[(m, n, k)] = int(v.split('/')[-1])            # tuples up to a common rational factor
asked = {}; seen = {}; versions = set(); blocks = raw.split('=== ')[1:]
for block in blocks:
    inp, resp = block.split('--- response\n', 1); versions.update(re.findall(r'<version>(.*?)</version>', resp))
    mn = re.search(r'm := (\d+); n := (\d+);', inp); ts = re.search(r'for T in \[(.*?)\] do', inp)
    assert mn and ts, 'request without a ratio or a list of triples'
    m, n = int(mn.group(1)), int(mn.group(2))
    for t in re.findall(r'\[([\d,]+)\]', ts.group(1)):
        key = (m, n, '.'.join(NAMES[int(j) - 1] for j in t.split(','))); asked[key] = asked.get(key, 0) + 1
    lines = re.findall(r'<line>(.*?)</line>', resp, re.S); cur = None
    for l in lines:
        mm = re.match(r'ratio (\d+) (\d+) forms \[ (.*?) \]', l)
        if mm:
            cur = (int(mm.group(1)), int(mm.group(2)), '.'.join(x.strip() for x in mm.group(3).split(',')))
            assert cur in asked, 'an answer for a curve that was not asked for'; res = {}
        elif cur and l.startswith('fake 2-Selmer set:'): res['sel'] = int(l.split(':')[1])
        elif cur and l.startswith('algebra modulus is f up to a constant:'): res['alg'] = l.split(':')[1].strip()
        elif cur and l.startswith('elements at which each of these forms has a square norm:'):
            res['cnt'] = int(l.split(':')[1])
            if cur in seen: assert seen[cur] == res, 'two finished runs of one curve disagree'
            seen[cur] = dict(res)
bad = 0; done = 0; unfinished = 0
for key in sorted(asked):
    d = seen.get(key)
    if d is None:
        unfinished += 1; print('(%d,%d) %s: Magma did not finish (%d request%s, no result)' % (key + (asked[key], 's' if asked[key] > 1 else ''))); continue
    o = ours.get(key); ok = o == d['cnt'] and d.get('alg') == 'true'; done += 1; bad += not ok
    print('(%d,%d) %s: fake 2-Selmer set of Magma %d, of these with the three norms square %d; tuples left by the first descent, up to a rational factor: %s  %s' % (key + (d['sel'], d['cnt'], o, 'equal' if ok else 'DIFFERENT')))
print('%d requests, %d curves asked for, %d finished, %d not finished; of the finished ones %d different; Magma %s' % (len(blocks), len(asked), done, unfinished, bad, ', '.join(sorted(versions))))
try: pairs = [tuple(int(v) for v in l.split(',')) for l in open('pairs34.txt').read().split()]
except FileNotFoundError: pairs = []
if pairs:
    empty = [p for p in pairs if any(k[:2] == p and d['sel'] == 0 for k, d in seen.items())]
    nosq = [p for p in pairs if p not in empty and any(k[:2] == p and d['cnt'] == 0 for k, d in seen.items())]
    never = [p for p in pairs if not any(k[:2] == p for k in seen)]
    other = [p for p in pairs if p not in empty and p not in nosq and p not in never]
    f = lambda L: ' '.join('%d,%d' % p for p in L) if L else 'none'
    print('the %d ratios of pairs34.txt: Magma finds an empty fake 2-Selmer set for some triple: %d' % (len(pairs), len(empty)))
    print('  a non-empty set, no element with the three norms square: %d: %s' % (len(nosq), f(nosq)))
    print('  no finished request: %d: %s' % (len(never), f(never)))
    print('  finished, some element left: %d: %s' % (len(other), f(other)))
assert bad == 0 and done > 0
