#!/usr/bin/env python3
"""That the test `oracle` of tuple_descent.py can fail (EXPERIMENT): two deliberately wrong class rules, and the right ones.
usage: tuple_mutation.py"""
import tuple_descent as td
forms, sup, known = td.magic(627, 364)
orig = td.key2
def bad(z):
    x, y = z; w = 0
    while (x + y) % 2 == 0: x, y = (x + y)//2, (y - x)//2; w += 1
    return w, (w % 2, min(((x*a - y*b) % 4, (x*b + y*a) % 4) for a, b in td.SQ8))
td.key2 = bad; td.D2 = {bad((d, 0))[1]: None for d in (1, 2, 3, 6)}
print('wrong rule (units of Z_2[i] compared modulo 4, not 8), prime 2: %d values tested, disagreements: %d' % td.oracle(forms, primes=((2, 6),)))
td.key2 = orig; td.D2 = {orig((d, 0))[1]: None for d in (1, 2, 3, 6)}
ol = td.leg; td.leg = lambda a, p: 1 - ol(a, p)
print('wrong rule (residue symbol reversed), primes 3, 5, 7: %d values tested, disagreements: %d' % td.oracle(forms, primes=((3, 4), (5, 3), (7, 3))))
td.leg = ol
print('the rules of the program: %d values tested, disagreements: %d' % td.oracle(forms))
