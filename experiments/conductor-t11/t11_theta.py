#!/usr/bin/env python3
"""t11: Fisher's higher descents (Theta_2, Theta_3, Theta_4; T. Fisher, arXiv:1509.03234, Sects. 2, 3, 7) on the
nine curves of REPORT_I5 (the rows "conditional on a rank being 1"), in exact arithmetic, with the repository's
code (rank_indep, rank_second, rank_fisher, rank_fisher5, rank_fisher_theta4, rank_fisher_theta4_gen).

    python3 t11_theta.py [label ...] [--write]      the checked computation (hints from t11_hints.json; no finder)
    python3 t11_theta.py --negative [--write]       forged inputs that must be refused
    (outputs only with --write: t11_theta_<tag>.txt per curve, t11_theta_negative.txt; without it the script prints
    and compares with the committed file: '# matches' / '# DIFFERS')
    finders (PARI, not certifying; they wrote t11_hints.json): explore/try_l3.py, explore/try_find.py,
    explore/try_find2.py, run from this folder (see NOTES.md)

Written by an AI model (build agent t11), not read by a person.  Experiment, not counted.

WHAT IS DONE FOR EACH CURVE (see NOTES.md for the reading of Fisher's paper that this rests on):
 (0) the model is rebuilt from (m, n, d) with the census's code and compared with I5's minimal model (t11_lib);
     E: y^2 = x(x^2 + a x + b), a^2 - 4b = -square; a point of infinite order (n G != O, n <= 12, on E or E').
 (1) rank_indep.rank_bound (exact first and second descent): S_1, S'_1 and the open classes (the classes of the
     first descent not excluded by the second descent and not images of known points).  rank_indep's "E side" (d | b)
     is Fisher's S'_1, its "E' side" (d' | b') is Fisher's S_1.
 (2) Theta_2 (S_2 x S'_2, Thm 2.1): <D_lam, eta> by rank_second's Hilbert-symbol formula for checked level-2 lifts
     D_lam (local points at every place), both ways (a lift of each side) when both sides have open classes.
 (3) Theta_3 (alternating on S_3, or on S'_3): for a class y of one side, a level-3 lift C3 (rank_fisher.py step (4),
     checked locally soluble at R and every prime of S3), C3' and ONE pushout form f on it (rank_fisher5 steps (1)-(7):
     a conic point over a quadratic field, re-checked exactly; Lemma 7.1), and Theta_3(y, z) = sum_v (f(P_v), z)_v
     (Sect. 3) for z of the SAME side, at R and every prime of S (t11_lib.theta3; the term must not depend on the
     sampled point P_v, which is checked; outside S a good-reduction argument).  Either pushout of the two cone pairs
     may be used: f / f_dag = c h^2 (rank_fisher5 (8)), and sum_v (c, z)_v = 1.
 (4) Theta_4 (S_4 x S'_4): rank_fisher_theta4_gen.py, unchanged, with its curve source redirected (t11_lib
     .theta4_gen_module) and its twist sampler choose_Xi_lin replaced by one that closes the sampled classes with the
     local image (t11_lib.choose_Xi_lin_images; a sampler only, C4 is then checked locally soluble).
 (5) the bookkeeping: Thm 2.1 (kernels), Lemma 2.2 (rank <= dim S_{2n-1} + dim S'_{2n} - 2, and with the primes
     swapped).
"""
import os, sys, json, io, time, contextlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import t11_lib as L
RI, R2, RF = L.RI, L.R2, L.RF

WRITE = '--write' in sys.argv
HINTS = os.path.join(HERE, 't11_hints.json')
OUT = []
def say(s=''):
    print(s, flush=True); OUT.append(s)

# ---------------------------------------------------------------- the plans (what was found to work; see NOTES.md)
# Each plan: the side and classes lifted, the pair (U or V) whose pushout is used, the classes paired.
PLANS = {}

def theta2(C, Rd, side, dds, etas, say=say, nlam=2):
    """Theta_2(dd, eta) = <D_lam, eta> (rank_second.DLam.pairing: local points of D_lam at every place of S, Hilbert
    symbols) for checked level-2 lifts D_lam of each dd (side `side`) and classes eta of the other side; the value
    must not depend on the lift (nlam lifts each)."""
    AA, BB = (C['a'], C['b']) if side == 'E' else (C['ap'], C['bp'])
    M = {}
    for dd in dds:
        P0 = RI.find_conic_point(dd, AA, BB // dd)
        surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
        R2.check_conic_forms(dd, AA, BB // dd, info['U'], info['V'], P0[2])
        rows = []
        for lam in sorted(surv, key=abs):
            try:
                D = R2.DLam(lam, info['U'], info['V'], bad=info['primes'])
                pv = {eta: D.pairing(eta)[0] for eta in etas}
            except AssertionError as ex:
                if 'no local point' not in str(ex) and 'no conic point' not in str(ex): raise
                continue
            rows.append((lam, pv))
            if len(rows) >= nlam: break
        assert rows and all(r[1] == rows[0][1] for r in rows), 'Theta_2 depends on the lift'
        M[dd] = rows[0][1]
        say('  (2) Theta_2(%d, .) (%s side, lifts lam = %s): %s' % (dd, side, [r[0] for r in rows],
                                                                   ' '.join('%d:%+d' % kv for kv in rows[0][1].items())))
    return M

def plan_10512m(C, Rd, rec, dims, say=say):
    b0 = RI.sqfree_part(C['b'])
    zs = [2, b0, 11878571404153441, 23757142808306882]
    say('  type A: two open classes on the E side (S\'_1 = S\'_2 = <b0, 2, 11878571404153441>), the E\' side all known: '
        'Theta_2 = 0 (its left kernel contains S_2 = images of the rational points), S\'_3 = S\'_2; Theta\'_3 is '
        'alternating on S\'_3 with b0 in its kernel, so it is fixed by Theta\'_3(2, 11878571404153441)')
    r1 = theta3_row(C, Rd, 'E', 2, 'V', zs, 'E1(10512,145)-#2', say=say)
    say('  control: the coset partner 2 b0 (another C3, its U pair over Q(sqrt 2)); Theta\'_3(2 b0, .) = Theta\'_3(2, .)')
    r2 = theta3_row(C, Rd, 'E', RI.sqfree_part(2 * b0), 'U', zs, 'E1(10512,145)-#2b0', say=say)
    assert r1 == r2, 'the two lifts disagree'
    assert r1[2] == 1 and r1[b0] == 1, 'CONTROL FAILED: Theta_3(y, y) or the pairing with the image of T is not 0'
    th = r1[11878571404153441]
    if th == -1: return 'Theta\'_3 != 0: S\'_4 = known, rank <= 2 + 1 - 2 = 1'
    return ('Theta\'_3(2, 11878571404153441) = 0 (both lifts; controls 0): Theta\'_3 = 0, S\'_4 = S\'_3 (dim 3), '
            'S_4 = S_3 (dim 2, all known, so Theta_4 = 0 too); the 2-isogeny descents through level 4 give exactly 3; '
            'not closed (the next pairing, Theta\'_5, is beyond Fisher Sect. 7)')
PLANS['E1(10512,145)-'] = dict(run=plan_10512m)

def plan_8249(C, Rd, rec, dims, say=say):
    say('  type C: two open classes on each side.  Theta_2 first (rational: no quadratic field)')
    M = theta2(C, Rd, 'E', [13278577, 755496193], [113, 2991753233, -1, 108415202], say=say)
    Mt = theta2(C, Rd, "E'", [113, 2991753233], [13278577, 755496193, RI.sqfree_part(C['b'])], say=say)
    for y in (13278577, 755496193):
        for eta in (113, 2991753233): assert M[y][eta] == Mt[eta][y], 'Theta_2 computed from the two sides disagrees'
    assert M[13278577][-1] == M[13278577][108415202] == M[755496193][-1] == M[755496193][108415202] == 1
    assert (M[13278577][113], M[13278577][2991753233], M[755496193][113], M[755496193][2991753233]) == (1, -1, 1, 1)
    say('  Theta_2 on (open E side) x (open E\' side) = [[0, 1], [0, 0]] over F_2 (rank 1; both ways): S_3 = <known, 113> '
        '(dim 3), S\'_3 = <known, 755496193> (dim 2); Theta_3, Theta\'_3 alternating with one open class: 0; '
        '4-descent bound (Lemma 2.2, n = 2) 3 + 2 - 2 = 3 (from 5)')
    say('  Theta_4(113, 755496193): rank_fisher_theta4_gen on the E\' side level-3 lift of 113 (xi twisted by 11969: the '
        'first xi has a non-zero level-4 obstruction on the excluded classes of S_1, as Thm 2.1 allows)')
    C3 = level3(C, Rd, "E'", 113, mult=11969, say=say)
    G = L.theta4_gen_module(C, 'E8249x', "E'", 'E', 113, 755496193, C3, (3, 2), HINTS)
    G.choose_Xi_lin = L.choose_Xi_lin_images(G, Rd)
    G.say = lambda s='': (say(s), G.OUT.append(s))[0]
    th = G.main('E8249x')
    if th == -1: return 'bound 5 -> 3 (Theta_2) -> 1 (Theta_4 != 0)'
    return 'bound 5 -> 3 (Theta_2, rank 1, exact through rank_second); Theta_4(113, 755496193) = 0: stays 3; not closed'
PLANS['E1(8249,3720)'] = dict(run=plan_8249)

def plan_22057(C, Rd, rec, dims, say=say):
    zs = [-1, 14, 274, 339, 108146]
    say('  type A\' on the E\' side: S_1 = S_2 (no class excluded) of dim 5 = <-1> + 4 open (basis 14, 274, 339, 108146); '
        'the E side all known, so Theta_2 = 0 and S_3 = S_2; Theta_3 is alternating on S_3 with -1 in its kernel')
    r14 = theta3_row(C, Rd, "E'", 14, 'U', zs, 'S22057#14', say=say)
    r274 = theta3_row(C, Rd, "E'", -274, 'V', zs, 'S22057#-274', say=say)
    for r, y in ((r14, 14), (r274, 274)):
        assert r[-1] == 1 and r[y] == 1, 'CONTROL FAILED: Theta_3(y, y) or the pairing with the image of T\' is not 0'
    assert r14[274] == r274[14], 'CONTROL FAILED: Theta_3 is not symmetric'
    assert r14[339] == -1 and all(r274[z] == 1 for z in zs), 'not the values recorded in NOTES.md'
    return ('Theta_3(14, 339) != 0 and Theta_3(274, .) = 0: Theta_3 has rank exactly 2 on the 4-dim open part (274 in '
            'its kernel), dim S_4 = 3, rank <= dim S\'_3 + dim S_4 - 2 = 2 + 3 - 2 = 3 (from 5); Theta_4 = 0 (S\'_4 all '
            'known); not closed')
PLANS['Sig+(22057,2351)'] = dict(run=plan_22057)

def plan_4440(C, Rd, rec, dims, say=say):
    b0 = RI.sqfree_part(C['b'])
    zs = [2, b0, 970682689, 6846322457, 119601394217, 1620825056929]
    say('  type A\' on the E side: S\'_1 = S\'_2 of dim 5 = <b0> + 4 open; the E\' side all known (Theta_2 = 0); '
        'Theta\'_3 alternating on S\'_3; the only open class with a small field is 2 (Q(sqrt 2); the next is '
        'Q(sqrt 1941365378), regulator about 895)')
    r2 = theta3_row(C, Rd, 'E', 2, 'V', zs, 'E1(4440,2231)#2', say=say)
    assert r2[2] == 1 and r2[b0] == 1, 'CONTROL FAILED'
    assert all(v == 1 for v in r2.values()), 'not the value recorded in NOTES.md'
    return ('Theta\'_3(2, .) = 0: 2 is in the kernel, dim S\'_4 >= 2, so these descents cannot give less than '
            '2 + 2 - 2 = 2 (at best 3); bound stays 5 here (the other rows need fields with large regulators); not closed')
PLANS['E1(4440,2231)'] = dict(run=plan_4440)

def field_info(D):
    """class number and regulator of Q(sqrt D) (PARI quadclassunit, under GRH; for information only)."""
    D = RI.sqfree_part(D); disc = D if D % 4 == 1 else 4 * D
    q = L.pari.quadclassunit(disc)
    return 'Q(sqrt %d): h = %s%s' % (D, q[0], (', R = %.1f' % float(q[3])) if D > 0 else '')

def plan_fields(kind, side, classes, why, lift=None):
    def run(C, Rd, rec, dims, say=say):
        b0 = RI.sqfree_part(C['b'])
        say('  %s; the cone fields of a level-3 lift of y are Q(sqrt y) and Q(sqrt %s) (Fisher Sect. 5; rank_fisher5 (4\'))'
            % (kind, 'y b0' if side == 'E' else '-y'))
        for y in classes:
            w = RI.sqfree_part(y * b0) if side == 'E' else -y
            say('    y = %d: %s; %s   (PARI quadclassunit, GRH, for information only)' % (y, field_info(y), field_info(w)))
        if lift:
            C3 = level3(C, Rd, side, lift, say=say)
        return why
    return run

PLANS['E1(3540,2581)'] = dict(run=plan_fields(
    'type A: two open classes on the E side, the E\' side all known; Theta\'_3 needs one pushout over Q(sqrt w), w one '
    'of the six classes below', 'E', [234494823676754, 522222777812402, 1890333437563521337],
    'NOT ATTEMPTED (field too large): every cone field is real quadratic with regulator >= 396790, where the conic '
    'points have coefficients of more than 10^5 digits; bound stays 3'))
PLANS['D1(12007,1799)'] = dict(run=plan_fields(
    'type A: two open classes on the E side, the E\' side all known', 'E', [354384937, 3909914385316993, 4678765003011961],
    'FAILED (no conic point found): the level-3 lift of 354384937 exists (above); its V pair is over Q(sqrt 354384937) '
    '(regulator about 18900): after the square reduction the conic coefficients have 5115 and 15443 digits; the lattice '
    'finder found no point (20 trials), the norm finder stopped in its lattice reduction (precision); every other field '
    'is larger; bound stays 3', lift=354384937))
def plan_12204(C, Rd, rec, dims, say=say):
    run0 = plan_fields('type B: one open class on each side (E\' side: 959; E side: 2); Theta_4 needs both cone fields of '
                       'one side; the E side fields include Q(sqrt 7015612900073813401545347960788802), so the E\' side',
                       "E'", [959], None)
    run0(C, Rd, rec, dims, say=say)
    H = json.load(open(HINTS))
    for dd, skip, name in ((959, 0, 'S12204'), (-959, 0, 'S12204m')):
        C3 = level3(C, Rd, "E'", dd, say=say)
        R = L.analyse3(name, "E'", C3, dd, hints=H)
        assert R['pushout'] is not None
        say('  lift of %d (lam = %s, xi = %s): the U pair over Q(sqrt %d) has a checked conic point (rank_fisher5 (6)); '
            'the V pair over Q(sqrt %d): no point found' % (dd, C3['lam'], C3['xi'], RI.sqfree_part(int(R['dU'])),
                                                            RI.sqfree_part(int(R['dV']))))
    return ('UNFINISHED (finder): for each lift tried the U-pair conic point is found in seconds, but the V-pair conic '
            'has coefficient norms with cofactors of 289 and 364 digits (after the known primes and p < 1e5) that were '
            'not factored, so its square parts cannot be reduced and neither finder got a point in 20 minutes; '
            'Theta_4 not computed; bound stays 3')
PLANS['Sig+(12204,9853)'] = dict(run=plan_12204)
PLANS['E1(10512,145)+'] = dict(run=plan_fields(
    'type B: one open class on each side (E\' side 23868431349, E side 11878571404153441)', "E'", [23868431349],
    'UNFINISHED (time): the E\' side lift of 23868431349 exists (above); the finder over Q(sqrt -23868431349) (class '
    'group [4064, 4, 2, 2, 2], reduced norms of 134 and 76 digits) did not finish in 20 minutes; the E side fields have '
    'regulators >= 1.4e8; bound stays 3', lift=23868431349))
PLANS['E1(6480,4961)'] = dict(run=plan_fields(
    'type B: one open class on each side; Theta_4 needs both cone fields of one side', "E'", [34129302, 20942341989],
    'NOT ATTEMPTED (time): the E\' side lift of 34129302 exists (above); its fields Q(sqrt +-34129302) include a real '
    'one with regulator about 490; the E side fields have regulators >= 1.8e7; bound stays 3', lift=34129302))

def tagof(lab):
    return lab.replace('(', '_').replace(')', '').replace(',', '_').replace('+', 'p').replace('-', 'm')

def selmer(C, say=say):
    r, rec = RI.rank_bound(C['a'], C['b'], C['pts'], C['ptsp'], say=lambda s: None)
    kE = rec['E']['log2']; kEp = rec["E'"]['log2']
    hE = len(rec['E']['known']).bit_length() - 1; hEp = len(rec["E'"]['known']).bit_length() - 1
    say('  (1) rank_indep: |S\'_1| = %d (E side, d | b), |S_1| = %d (E\' side, d\' | b\'); after the second descent '
        'dim S\'_2 <= %d, dim S_2 <= %d; known images: dim %d (E side), %d (E\' side); exact bound %d + %d - 2 = %d'
        % (len(rec['E']['selmer']), len(rec["E'"]['selmer']), kE, kEp, hE, hEp, kE, kEp, r))
    say('      open cosets: E side %s; E\' side %s' % ([c[0] for c in rec['E']['open']], [c[0] for c in rec["E'"]['open']]))
    return r, rec, (kE, kEp, hE, hEp)

def level3(C, Rd, side, dd, mult=None, say=say):
    C3, pv = L.level3_lift(C, Rd, side, dd, say=say)
    if mult is not None:
        C3 = L.recertify_xi(C3, Rd, side, mult, say=say)
    return C3

def theta3_row(C, Rd, side, y, pair, zs, name, mult=None, say=say):
    """Theta_3(y, z) for z in zs (same side), through the level-3 lift of y and the pushout of `pair`."""
    C3 = level3(C, Rd, side, y, mult=mult, say=say)
    H = json.load(open(HINTS))
    R = L.analyse3(name, side, C3, y, hints=H)
    po = R['pushout'] if pair == 'U' else R['pushout3']
    assert po is not None, 'no checked conic point for the %s pair of %s (hint %s)' % (pair, name, name)
    say('  C3\' of y = %d (%s side; lam = %s, xi = %s): rank_fisher5 steps (1)-(%s) checked; the %s pair over Q(sqrt %d) '
        'has a checked conic point (tangent plane l, N(l) rational: a pushout form, Lemma 7.1)'
        % (y, side, C3['lam'], C3['xi'], '6' if pair == 'U' else '7', pair,
           RF.RI.sqfree_part(int(R['dU'] if pair == 'U' else R['dV']))))
    out, S, cl = L.theta3(C, Rd, C3, R, pair, zs, say=say, side=side)
    return out

def run_curve(lab, say=say):
    plan = PLANS[lab]
    say('%s: %s' % (lab, L.BYLAB[lab][5]))
    C = L.curve(lab)
    say('  (0) model rebuilt from (m, n, d) by the census code = I5\'s minimal model %s; E: y^2 = x(x^2 + a x + b), '
        'a = %d, b = %d (a^2 - 4b = -square)' % (C['mi'], C['a'], C['b']))
    G = [P for P in C['pts'] if L.infinite_order(C['a'], C['b'], P)]
    Gp = [P for P in C['ptsp'] if L.infinite_order(C['ap'], C['bp'], P)]
    assert G or Gp, 'no point of infinite order'
    say('      a point of infinite order on %s (PARI ellrank as a hint; on the curve, n G != O for n <= 12): rank >= 1'
        % ('E' if G else "E'"))
    r, rec, dims = selmer(C, say=say)
    Rd = RI.first_descent(C['a'], C['b'])
    res = plan['run'](C, Rd, rec, dims, say=say)
    say('  RESULT %s: %s' % (lab, res))
    return res

def main(labs):
    for lab in labs:
        OUT.clear()
        t0 = time.time()
        say('t11_theta.py: Fisher\'s higher descents on %s (arXiv:1509.03234 Sects. 2, 3, 7), exact. Written by an AI '
            'model (build agent t11), not read by a person. Experiment, not counted.' % lab)
        run_curve(lab)
        f = os.path.join(HERE, 't11_theta_%s.txt' % tagof(lab))
        s = '\n'.join(OUT) + '\n'
        if WRITE: open(f, 'w').write(s)
        elif os.path.exists(f): print('# matches' if open(f).read() == s else '# DIFFERS')
        print('# %s: %.0f s' % (lab, time.time() - t0), file=sys.stderr)

if __name__ == '__main__':
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--negative' in sys.argv:
        import t11_negative
        t11_negative.main(WRITE)
    else:
        main(args or [lab for lab in PLANS])
