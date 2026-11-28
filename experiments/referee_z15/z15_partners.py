"""z15: partner ratios (m, n) -> (m + n, m - n): E_08 of one and E_LU of the other have the same minimal model (PARI),
on the session's closed-form models (checked isomorphic to t10's in z15_models.txt).  Written by an AI model (reviewer
z15), not read by a person."""
import cypari2
P = cypari2.Pari()
def cubic(e): return [0, -(e[0] + e[1] + e[2]), 0, e[0] * e[1] + e[0] * e[2] + e[1] * e[2], -e[0] * e[1] * e[2]]
def mm(e): return [int(x) for x in P.ellminimalmodel(P.ellinit(cubic(e)))[:5]]
def curves(m, n):
    A, B, C = m * m - n * n, 2 * m * n, m * m + n * n
    return {'08': mm((0, A * A, C * C)), 'LU': mm((0, B * B, C * C)), 'mix': mm((0, 2 * (C * C - A * B), 2 * (C * C + A * B)))}
print('z15_partners.py: written by an AI model (reviewer z15), not read by a person.')
for a, b in (((3320, 2303), (5623, 1017)), ((9113, 3200), (12313, 5913))):
    x, y = curves(*a), curves(*b)
    print(a, b, 'E_08 ~ E_LU\':', x['08'] == y['LU'], '; E_LU ~ E_08\':', x['LU'] == y['08'], '; E_mix ~ E_mix\':', x['mix'] == y['mix'])
