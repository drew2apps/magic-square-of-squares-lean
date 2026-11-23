#!/usr/bin/env python3
"""Magma's two-cover descent against the first descent, on subsets of three forms.

magma_twocover_raw.txt: verbatim requests to the public Magma calculator and its answers.  Each request runs
TwoCoverDescent on the genus 2 curve Z^2 = f_a f_b f_c for three of the five forms of a ratio, and counts the elements of
the fake 2-Selmer set at which each of the three norms is a square (the elements that a solution could give, since each form
is a square at a solution).  subsets.txt, subsets_34.txt: the number of tuples that the first descent leaves for the same three
forms, up to a common rational factor (the two-cover descent works modulo Q*).  The two numbers must be equal.  Magma is not run here.  usage: twocover_check.py [raw file [subsets file]]"""
import sys, re
raw = open(sys.argv[1] if len(sys.argv) > 1 else 'magma_twocover_raw.txt').read()
ours = {}
for l in [x for f in (sys.argv[2:] or ['subsets.txt', 'subsets_34.txt']) for x in open(f)]:
    pair, rest = l.split(':', 1); m, n = (int(v) for v in pair.strip('() ').split(','))
    for item in rest.split(): k, v = item.split('='); ours[(m, n, k)] = int(v.split('/')[-1])            # tuples up to a common rational factor
seen = {}; versions = set()
for block in raw.split('=== ')[1:]:
    resp = block.split('--- response\n', 1)[1]; versions.update(re.findall(r'<version>(.*?)</version>', resp))
    lines = re.findall(r'<line>(.*?)</line>', resp, re.S)
    cur = None
    for i, l in enumerate(lines):
        mm = re.match(r'ratio (\d+) (\d+) forms \[ (.*?) \]', l)
        if mm: cur = (int(mm.group(1)), int(mm.group(2)), '.'.join(x.strip() for x in mm.group(3).split(','))); seen.setdefault(cur, {})
        elif cur and l.startswith('fake 2-Selmer set:'): seen[cur]['sel'] = int(l.split(':')[1])
        elif cur and l.startswith('algebra modulus is f up to a constant:'): seen[cur]['alg'] = l.split(':')[1].strip()
        elif cur and l.startswith('elements at which each of these forms has a square norm:'): seen[cur]['cnt'] = int(l.split(':')[1])
bad = 0; done = 0
for key in sorted(seen):
    d = seen[key]
    if 'cnt' not in d: print('(%d,%d) %s: Magma did not finish' % key); continue
    o = ours.get(key); ok = o == d['cnt'] and d.get('alg') == 'true'; done += 1; bad += not ok
    print('(%d,%d) %s: fake 2-Selmer set of Magma %d, of these with the three norms square %d; tuples left by the first descent, up to a rational factor: %s  %s' % (key + (d['sel'], d['cnt'], o, 'equal' if ok else 'DIFFERENT')))
print('%d curves compared, %d different; Magma %s' % (done, bad, ', '.join(sorted(versions))))
assert bad == 0 and done > 0
