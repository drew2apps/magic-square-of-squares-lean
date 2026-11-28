#!/usr/bin/env python3
"""Make the task list (one line per coordinate j and distinct class value b_j) and identify the torsion classes with PARI."""
import json, sys, cypari2
pari = cypari2.Pari(); pari.allocatemem(10**9)
e = [(0, 0), (273990441250, 0), (-102929464375, 118967496648)]
sel = json.load(open('selmer.json'))['selmer']
tasks = {}
for idx, tr in enumerate(sel):
    for j in range(3):
        b = tuple(tr[j]); tasks.setdefault((j, b), []).append(idx)
with open('tasks.txt', 'w') as f:
    for (j, b), idxs in sorted(tasks.items(), key=lambda kv: (kv[0][1][0]**2 + kv[0][1][1]**2, kv[0][0])):
        k, l = [x for x in range(3) if x != j]
        D2 = (e[j][0] - e[k][0], e[j][1] - e[k][1]); D3 = (e[j][0] - e[l][0], e[j][1] - e[l][1])
        f.write('j%d_c%s %d %d %d %d %d %d\n' % (j + 1, '.'.join(map(str, idxs)), b[0], b[1], D2[0], D2[1], D3[0], D3[1]))
print(len(tasks), 'tasks')
