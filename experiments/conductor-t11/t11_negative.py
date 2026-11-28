"""t11_theta.py --negative: forged inputs that must be refused (and the true ones accepted).
Written by an AI model (build agent t11), not read by a person.  Experiment, not counted."""
import os, sys, json
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t11_lib as L
RI, R2, RF = L.RI, L.R2, L.RF
HINTS = os.path.join(HERE, 't11_hints.json')
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)
quiet = lambda s='': None

def main(write):
    say('t11_theta.py --negative: forged inputs that must be refused. Written by an AI model (build agent t11), not read '
        'by a person.')
    n = ok = 0
    def refuse(label, fn):
        nonlocal n, ok
        n += 1
        try: fn()
        except AssertionError as ex:
            ok += 1; say('  refused: %s (%s)' % (label, str(ex.args[0] if ex.args else type(ex).__name__)[:100])); return
        say('  NOT REFUSED: %s' % label)
    def accept(label, fn):
        nonlocal n, ok
        n += 1
        fn(); ok += 1; say('  accepted (positive control): %s' % label)
    lab = 'E1(10512,145)-'
    mi5 = L.BYLAB[lab][6]
    bad = list(mi5); bad[3] += 1
    refuse('a forged minimal model (a4 + 1) against the census curve',
           lambda: (_ for _ in ()).throw(AssertionError('not the same curve')) if not L.same_curve(L.census_model(lab)[0], bad) else None)
    C = L.curve(lab); Rd = RI.first_descent(C['a'], C['b'])
    b0 = RI.sqfree_part(C['b'])
    C3, _ = L.level3_lift(C, Rd, 'E', 2, say=quiet)
    H = json.load(open(HINTS)); name = 'E1(10512,145)-#2'; key = name + '|E|V'
    R = L.analyse3(name, 'E', C3, 2, hints=H)
    accept('the hinted conic point over Q(sqrt 2) on its cone (rank_fisher5 step (7))', lambda: R['pushout3'] is not None or (_ for _ in ()).throw(AssertionError('no pushout')))
    Hb = dict(H); h = dict(H[key]); h['Z'] = [str(Fr(h['Z'][0]) + 1), h['Z'][1]]; Hb[key] = h
    refuse('the conic point moved off its cone (Z -> Z + 1)', lambda: L.analyse3(name, 'E', C3, 2, hints=Hb))
    zs = [11878571404153441]
    refuse('a paired class 3 that is not in the Selmer group (not in the local image at some bad prime)',
           lambda: L.theta3(C, Rd, C3, R, 'V', [3], say=quiet))
    eta = 1534898            # in the E' side's S_1: a class of the wrong side for Theta_3 on the E side
    refuse('a class %d of the other side (E\') paired through the E-side C3: refused by the Selmer check' % eta,
           lambda: L.theta3(C, Rd, C3, R, 'V', [eta], say=quiet))
    refuse('the same class with the Selmer check switched off: the local term must depend on the point',
           lambda: L.theta3(C, Rd, C3, R, 'V', [eta], say=quiet, check_selmer=False))
    import sympy as sp
    import rank_fisher_theta4 as T
    orig = L.setup3
    def forged(R_, pair):
        D = orig(R_, pair)
        return dict(D, f0=sp.Poly(sp.expand(D['f0'].as_expr() + T.c * T.e), *T.VS))
    def run_forged():
        L.setup3 = forged
        try: L.theta3(C, Rd, C3, R, 'V', zs, say=quiet)
        finally: L.setup3 = orig
    refuse('f0 replaced by f0 + c e (not a pushout form): the local term depends on the point', run_forged)
    refuse('a level-3 twist by 3 (not in this side\'s S_1)', lambda: L.recertify_xi(C3, Rd, 'E', 3, say=quiet))
    C3f = dict(C3, xi=str(int(C3['xi']) * 10007))
    refuse('a forged xi (times 10007) with the stored local points', lambda: RF.check_C3_point(
        R2.DLam(int(C3['lam']), [int(x) for x in C3['U']], [int(x) for x in C3['V']]), int(C3f['xi']), 10007,
        tuple(int(x) for x in C3['local_points']['2'])) if '10007' in C3['local_points'] else
        RF.certify_C3(dict(D=R2.DLam(int(C3['lam']), [int(x) for x in C3['U']], [int(x) for x in C3['V']]), xi=int(C3f['xi']),
                           places=[('inf' if v == 'inf' else int(v)) for v in C3['local_points']],
                           cert={('inf' if v == 'inf' else int(v)): tuple(Fr(x) if v == 'inf' else int(x) for x in P) for v, P in C3['local_points'].items()})))
    say('negative tests: %d of %d as expected' % (ok, n))
    assert ok == n
    if write: open(os.path.join(HERE, 't11_theta_negative.txt'), 'w').write('\n'.join(OUT) + '\n')
    else:
        f = os.path.join(HERE, 't11_theta_negative.txt')
        if os.path.exists(f): print('# matches' if open(f).read() == '\n'.join(OUT) + '\n' else '# DIFFERS')

if __name__ == '__main__':
    main('--write' in sys.argv)
