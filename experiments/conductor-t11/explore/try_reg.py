import json, sys
import cypari2; pari = cypari2.Pari(); pari.allocatemem(2*10**9)
sys.path.insert(0, '.')
J = json.load(open('explore_rank_all.json'))
def core(n):
    return int(pari.core(n))
def info(D):
    D = core(D)
    disc = D if D % 4 == 1 else 4 * D
    try:
        q = pari.quadclassunit(disc)
        return '%d: h=%s R=%s' % (D, q[0], ('%.1f' % float(q[3])) if D > 0 else '-')
    except Exception as ex: return '%d: err %s' % (D, str(ex)[:40])
for lab, r in J.items():
    b = int(r['b']); a = int(r['a']); bp = a*a - 4*b
    b0 = core(b)
    print(lab, 'b0 =', b0)
    for side, op in (('E', r['E_open']), ("E'", r['Ep_open'])):
        for cos in op:
            for y in cos[:4]:
                if side == 'E': print('  E ', y, '|', info(y), '|', info(y * b0))
                else: print("  E'", y, '|', info(y), '|', info(-y))
