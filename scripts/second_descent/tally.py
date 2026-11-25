#!/usr/bin/env python3
"""The state of the 123 ratios after second_descent.py and lvalue_close.py.  usage: tally.py second_123.jsonl lvalue_open.jsonl"""
import sys, json
def partner(m, n):
    a, b = m + n, m - n
    return (a//2, b//2) if a % 2 == 0 else (a, b)
R = {tuple(json.loads(l)['pair']): json.loads(l) for l in open(sys.argv[1])}
Lv = {tuple(json.loads(l)['pair']): json.loads(l) for l in open(sys.argv[2])}
d2 = set(p for p in R if R[p]['status'] == 'closed'); dl = set(p for p in Lv if Lv[p]['status'].startswith('closed'))
assert not (d2 & dl) and set(Lv) == set(p for p in R if R[p]['status'] != 'closed')
closed = d2 | dl; full = set(p for p in R if p in closed or partner(*p) in closed)
OPEN24 = [(91,60),(151,31),(195,28),(209,120),(223,167),(329,89),(420,341),(440,279),(572,555),(608,105),(609,200),(627,364),(680,111),(713,503),(719,161),(741,580),(760,39),(761,79),(791,569),(799,721),(809,409),(861,620),(945,248),(991,263)]
print('%d ratios with a degenerate solution' % len(R))
print('closed with a rank bound 0 from descent (second_descent.py): %d' % len(d2))
print('closed with rank 0 from an L-value (lvalue_close.py): %d: %s' % (len(dl), ' '.join('%d,%d' % p for p in sorted(dl))))
print('closed, directly or through the partner: %d' % len(full))
print('not closed: %d: %s' % (len(R) - len(full), ' '.join('%d,%d' % p for p in sorted(R) if p not in full)))
print('of the 24 ratios that were open in the census at release 0.5.0: closed %d (by descent bounds alone %d), not closed: %s' % (
    sum(p in full for p in OPEN24), sum(p in d2 or partner(*p) in d2 for p in OPEN24), ' '.join('%d,%d' % p for p in OPEN24 if p not in full)))
