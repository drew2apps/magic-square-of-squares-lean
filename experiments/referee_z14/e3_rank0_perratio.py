"""z14: explore3's rank0.py bound() (its own code: the function text of e3_rank0.py, exec'd), per ratio, on the three
ratio sets of rank0_z14.py, for a line-by-line comparison with rank0_z14.jsonl.
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I e3_rank0_perratio.py | grep "^{" > e3_rank0_perratio.jsonl"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
src = open(os.path.join(HERE, 'e3_rank0.py')).read()
head = src.split('def is_sq')[0]          # imports, rank2(), bound()
exec(compile(head.replace("_os.path.dirname(_os.path.abspath(__file__))", repr(HERE)), 'e3_rank0.py', 'exec'))
src2 = open(os.path.join(HERE, 'rank0_z14.py')).read()
gen = src2[src2.index('def is_sq'):src2.index('def rank2')]
exec(gen)
for name, prs in (('nd5', ratios(5, 3000, 150, 'nd')), ('sq5', ratios(5, 3000, 80, 'sq')), ('nd6', ratios(6, 3000, 150, 'nd'))):
    for m, n in prs:
        A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
        try: dA, bA = bound((0, A * A, C * C), 7); dB, bB = bound((0, B * B, C * C), 7); err = None
        except Exception as ex: dA = bA = dB = bB = None; err = repr(ex)[:80]
        print(json.dumps(dict(set=name, m=m, n=n, dA=dA, bA=bA, dB=dB, bB=bB, err=err), sort_keys=True))
