"""z15: the pilot's 201 ratios (as t10 drew them) and an independent fixed-seed sample of the other 736.
Written by an AI model (reviewer z15), not read by a person."""
import random, re, sys, json, os
G2 = '../roadmap-g/g2/'
def read_list(path):
    out = []
    for l in open(path):
        l = l.strip()
        if not l or l[0] == '#' or 'left' not in l: continue
        m, n = map(int, l.split()[0].split(',')); out.append((m, n))
    return out
L997 = read_list(G2 + 'g2_997_all.txt'); L15 = read_list(G2 + 'g2_15000_need_good.txt')
pilot15 = sorted(random.Random(1).sample(L15, 100))
pilot = set(L997) | set(pilot15)
# check against the pilot's own output files
seen = set()
for f in ('t10_997_a.txt', 't10_997_b.txt', 't10_15000_s100_a.txt', 't10_15000_s100_b.txt'):
    for l in open('../conductor-t10/' + f):
        mm = re.match(r'ratio \((\d+),(\d+)\)', l)
        if mm: seen.add((int(mm.group(1)), int(mm.group(2))))
print('lists:', len(L997), len(L15), 'pilot set', len(pilot), 'pilot files', len(seen), 'equal', seen == pilot)
rest = sorted(set(L997 + L15) - pilot)
print('non-pilot ratios', len(rest))
S = sorted(random.Random(20261127).sample(rest, 80))
open('sample80_seed20261127.txt', 'w').write('\n'.join('%d,%d' % r for r in S) + '\n')
print('sample of 80 written; first', S[:3])
