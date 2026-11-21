#!/usr/bin/env python3
"""Evidence that the three kinds of local obstruction are the only ones (odd primes).

For each modulus M = p^k below, every projective class (a:b) modulo M is tested directly and
compared with the prediction: obstructed iff (p = 5 mod 8 and p divides a^2 + b^2 to an odd
power) or (p in {7, 23, 47} and (a:b) is one of the listed ratios). Classes where a^2 + b^2
is 0 modulo M are reported as undecided, because the power is then not visible modulo M.
"""
def val(x, p):
    k = 0
    while x % p == 0: x //= p; k += 1
    return k
def proj(M, p):
    return [(1, v) for v in range(M)] + [(u, 1) for u in range(0, M, p)]
def predicted(a, b, M, p):
    if p % 8 == 5:
        s = (a*a + b*b) % M
        return None if s == 0 else val(s, p) % 2 == 1
    c = {7: 4, 23: 4, 47: 25}.get(p)
    return c is not None and ((a*a - c*b*b) * (c*a*a - b*b)) % p == 0
MODULI = [(9, 3), (27, 3), (81, 3), (25, 5), (125, 5), (49, 7), (343, 7), (121, 11), (169, 13),
          (289, 17), (361, 19), (529, 23), (841, 29), (961, 31), (1369, 37), (2209, 47),
          (2809, 53)] + [(p, p) for p in (53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101)]
for (M, p) in MODULI:
    sq = [False] * M
    for x in range(M): sq[x * x % M] = True
    reps = proj(M, p); bad = obs = undecided = 0
    for (a, b) in reps:
        o = not any(sq[(a*a*V*V + b*b*U*U) % M] and sq[(2*(a*a+b*b)*(U*U+V*V)) % M]
                    and sq[(b*b*V*V + a*a*U*U) % M] for (U, V) in reps)
        obs += o
        pr = predicted(a, b, M, p)
        if pr is None: undecided += 1
        elif pr != o: bad += 1
    print(f'mod {M}: classes {len(reps)}, obstructed {obs}, disagreements {bad}, undecided {undecided}')
