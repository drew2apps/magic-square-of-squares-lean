"""t11: shared code for the Theta_3 / Theta_4 runs on the nine I5 curves (experiments/conductor-t11).
Written by an AI model (build agent t11), not read by a person.  Experiment, not counted.

The models are rebuilt from (m, n, d) with the census's own code, as REPORT_I5 did (census-wide/sd_pari.py `family`,
`weierstrass` for E1; census2000/dm_pari.py `all_curves` for D1, Sig+); PARI gives the minimal model (a finder: it is
checked by c4, c6 up to u^4, u^6) and it is compared with I5's printed minimal model; then rank_indep's
shift_to_two_torsion gives E: y^2 = x(x^2 + a x + b)."""
import os, sys, io, json, contextlib, random
from fractions import Fraction as Fr
from math import isqrt, gcd
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
G1 = os.path.join(REPO, 'scripts', 'g1_full')
for p_ in (os.path.join(REPO, 'experiments', 'census-wide'), os.path.join(REPO, 'experiments', 'census2000'), G1):
    if p_ not in sys.path: sys.path.insert(0, p_)
_argv = sys.argv; sys.argv = [sys.argv[0]]
with contextlib.redirect_stdout(io.StringIO()):
    import sd_pari as SD
    import dm_pari as DM
    import rank_indep as RI
    import rank_second as R2
    import rank_fisher as RF
sys.argv = _argv
pari = SD.pari

# (label, source, m, n, d, rows resting on its rank (REPORT_I5's table), I5's printed minimal model)
CURVES = [
 ('E1(3540,2581)', 'E1', 3540, 2581, -11740078, '(3540,2581) [CENSUS5000], (6121,959) [CENSUS7000, partner]',
  [0, -1, 0, 16899659070115573331857545392554393999753919, -4577250886777635699050726927293136713005349082509501848957193919]),
 ('E1(4440,2231)', 'E1', 4440, 2231, -13342, '(4440,2231) [CENSUS5000], (6671,2209) [CENSUS7000, partner]',
  [0, -1, 0, 39795335239688639860836320471793848959, -84072399232288381781467536481424339431003081071435104895]),
 ('E1(6480,4961)', 'E1', 6480, 4961, 709342, '(6480,4961) [CENSUS7000]',
  [0, 1, 0, 9211122592674457236286980493932726649216639, 1318503358214771909347173196865869794685175648543198847781303935]),
 ('E1(8249,3720)', 'E1', 8249, 3720, -108415202, '(8249,3720) [CENSUS9000]',
  [0, 1, 0, 263016164015259849023980430859354403964456553599, -65439936015720822957240852006058045665149706287739221875757099862115201]),
 ('E1(10512,145)+', 'E1', 10512, 145, 220962238, '(10512,145), (10657,10367) [CENSUS11000]',
  [0, 1, 0, -2419208314318668146199775264884775233139683805441, 1463244133933325813781975262027627717377710145889707390761192388296636159]),
 ('E1(10512,145)-', 'E1', 10512, 145, -220962238, '(10512,145), (10657,10367) [CENSUS11000]',
  [0, -1, 0, -2419208314318668146199775264884775233139683805441, -1463244133933325813781975262027627717377710145889707390761192388296636159]),
 ('D1(12007,1799)', 'D1', 12007, 1799, None, '(6903,5104) [CENSUS7000], (12007,1799) [G3]',
  [0, 0, 0, -2100795130954721348682439136479163983865084, 1304206922206538200817834175771821647493983141228332733879641600]),
 ('Sig+(12204,9853)', 'Sig+', 12204, 9853, None, '(12204,9853) [G3 open]',
  [0, -1, 0, -86753333016712673790059436059407475720840059034305, -3978647120759912663046198465625130039355274066092469254836088590418933674303]),
 ('Sig+(22057,2351)', 'Sig+', 22057, 2351, None, '(12204,9853) [G3 open] through its partner (22057,2351)',
  [1, 0, 0, -3273186841081526246623448342297018808717045938715, 5262479891378770656266792116117772613394291382315357263651247762914134225]),
]
BYLAB = {c[0]: c for c in CURVES}

def c46(ai):
    a1, a2, a3, a4, a6 = [Fr(x) for x in ai]
    b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
    b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
    return b2 * b2 - 24 * b4, -b2 ** 3 + 36 * b2 * b4 - 216 * b6

def same_curve(ai, bi):
    """c4, c6 of the two models agree up to u^4, u^6 with u^2 a rational square (exact)."""
    (c4, c6), (d4, d6) = c46(ai), c46(bi)
    if c4 == 0 or c6 == 0: return False
    u2 = (c6 / d6) / (c4 / d4)
    return u2 * u2 == c4 / d4 and u2 ** 3 == c6 / d6 and u2 > 0 and isqrt(u2.numerator) ** 2 == u2.numerator \
        and isqrt(u2.denominator) ** 2 == u2.denominator

def census_model(lab):
    _, src, m, n, d, rows, mi5 = BYLAB[lab]
    if src == 'E1':
        Q, r = SD.family(m, n)['+13']; ain, A = SD.weierstrass(d * Q, r)
    else:
        ain = [Fr(int(v)) for v in [c[0] for c in DM.all_curves(m, n) if c[4] == src][0]]
    Em = pari.ellminimalmodel(pari.ellinit([pari(x.numerator) / x.denominator for x in ain]))   # a finder
    mi = [int(v) for v in Em[:5]]
    assert same_curve(ain, mi), 'the minimal model is not the census curve'
    assert mi == mi5, "the rebuilt minimal model differs from I5's"
    return ain, mi

def curve(lab):
    """E: y^2 = x(x^2 + a x + b) (the census curve, one rational 2-torsion point), with a^2 - 4b = -square; a point G
    of E(Q) of infinite order (PARI ellrank on E or on E' as a hint, checked on the curve and n G != O for n <= 12,
    Mazur); returns dict."""
    ain, mi = census_model(lab)
    a1, a2, a3, a4, a6 = mi
    if a1 == 0 and a3 == 0: A2, A4, A6 = a2, a4, a6
    else:
        b2, b4, b6 = a1 * a1 + 4 * a2, 2 * a4 + a1 * a3, a3 * a3 + 4 * a6
        A2, A4, A6 = b2, 8 * b4, 16 * b6            # y^2 = x^3 + b2 x^2 + 8 b4 x + 16 b6, isomorphic over Q
    assert same_curve(mi, [0, A2, 0, A4, A6])
    e, a, b = RI.shift_to_two_torsion(A2, A4, A6)
    ap, bp = -2 * a, a * a - 4 * b
    assert bp < 0 and isqrt(-bp) ** 2 == -bp, 'a^2 - 4b is not minus a square'
    pts = RI.find_points([0, a, 0, b, 0]); ptsp = RI.find_points([0, ap, 0, bp, 0])
    for P in pts: assert RI.on_curve([0, a, 0, b, 0], P)
    for P in ptsp: assert RI.on_curve([0, ap, 0, bp, 0], P)
    return dict(lab=lab, a=a, b=b, ap=ap, bp=bp, pts=pts, ptsp=ptsp, mi=mi)

def ell_add(a, b, P, Q):
    """y^2 = x^3 + a x^2 + b x, affine points or None (O)."""
    if P is None: return Q
    if Q is None: return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and y1 == -y2: return None
    lam = (3 * x1 * x1 + 2 * a * x1 + b) / (2 * y1) if (x1, y1) == (x2, y2) else (y2 - y1) / (x2 - x1)
    x3 = lam * lam - a - x1 - x2
    return (x3, lam * (x1 - x3) - y1)

def infinite_order(a, b, P):
    """n P != O for n <= 12 (exact; Mazur: a rational torsion point has order <= 12)."""
    P = (Fr(P[0]), Fr(P[1])); Q = None
    for n in range(1, 13):
        Q = ell_add(a, b, Q, P)
        if Q is None: return False
    return True

EQ2, EQ3 = RF.EQ2, RF.EQ3

def dl_json(D):
    return dict(equations=EQ2, lam=str(D.lam), U=[str(c_) for c_ in D.U], V=[str(c_) for c_ in D.V],
                l_X=[str(c_) for c_ in D.X['l']], l_Y=[str(c_) for c_ in D.Y['l']],
                p_X=[str(c_) for c_ in D.X['P']], p_Y=[str(c_) for c_ in D.Y['P']])

def other_basis(C, side, R):
    """an F_2 basis of the other side's first Selmer group (classes certified locally soluble), torsion class first."""
    tors_other = RI.sqfree_part(C['bp'] if side == 'E' else C['b'])
    grp = R['Sp'] if side == 'E' else R['S']
    return RF.rref_ints([tors_other] + sorted(grp, key=abs))

def level3_lift(C, R, side, dd, lam_cap=12, say=print, skip=0):
    """rank_fisher.py steps (2)-(4) for the class dd of the given side (copied from rank_fisher.analyse, made
    general: any number of open classes): the level-2 lifts D_lam of dd not excluded by rank_indep's second descent,
    their pairings with a basis of the other side's S_1 (rank_second's Hilbert-symbol formula; Theta_2 and the
    level-3 obstruction), and the level-3 Selmer set S(D_lam, f) (rank_fisher.selmer_set_level3, sampled twice);
    the first D_lam with a certified C3 is returned (rank_fisher.certify_C3: a checked local point at R and at
    every prime of S3)."""
    a, b, ap, bp = C['a'], C['b'], C['ap'], C['bp']
    AA, BB = (a, b) if side == 'E' else (ap, bp)
    for p in R['bad']:
        assert R2.certified(dd, AA, BB, p, R['locS'] if side == 'E' else R['locSp']), ('class not locally soluble', dd, p)
    assert RI.real_ok(dd, AA, BB)
    c = BB // dd
    P0 = RI.find_conic_point(dd, AA, c)
    surv, info = RI.second_descent(dd, AA, BB, P0, stop_at_first=False)
    U, V = info['U'], info['V']
    R2.check_conic_forms(dd, AA, c, U, V, P0[2])
    Bo = other_basis(C, side, R)
    say('  %s side, class %d: lifts D_lam not excluded by the second descent: %d (tried: up to %d); basis of the other '
        'side\'s S_1: %s' % (side, dd, len(surv), lam_cap, Bo))
    tried = 0
    for lam in sorted(surv, key=abs):
        if tried >= lam_cap: break
        try: D = R2.DLam(lam, U, V, bad=info['primes'])
        except AssertionError as ex:
            if 'no conic point' not in str(ex): raise
            continue
        pv = {}
        try:
            for eta in Bo: pv[eta] = D.pairing(eta)[0]
        except AssertionError as ex:
            if 'no local point found' not in str(ex) and 'not constant' not in str(ex): raise
            continue
        tried += 1
        assert not D.special_rational
        res = RF.selmer_set_level3(D, R, side)
        res2 = RF.selmer_set_level3(D, R, side, seed=101)
        assert (res2['nonempty'] is False) == (res['nonempty'] is False) and res2['dim'] == res['dim'], \
            ('the level-3 Selmer set depends on the sample', lam)
        pred = all(x == 1 for x in pv.values())
        agree = (res['nonempty'] is not False) == pred if res['nonempty'] is not None else None
        assert agree is not False, 'the level-3 Selmer set disagrees with the pairings'
        st = {True: 'non-empty (dim %s), xi = %s, C3 checked locally soluble at R and %d primes' % (res['dim'], res.get('xi'), len(res['places']) - 1),
              False: 'empty', None: 'not decided', 'theory': 'non-empty, no sampled point for one xi'}[res['nonempty']]
        say('    lam = %d: pairings with the other side\'s S_1 %s; level-3 Selmer set %s (%s)'
            % (lam, ' '.join('%d:%+d' % kv for kv in pv.items()), st, 'agrees' if agree else 'not compared'))
        if res['nonempty'] is True and skip > 0:
            skip -= 1; continue
        if res['nonempty'] is True:
            RF.certify_C3(res)
            r3 = res
            return dict(dl_json(r3['D']), xi=str(r3['xi']), equations=EQ3,
                        local_points={str(v): [str(x_) for x_ in P] for v, P in r3['cert'].items()}), pv
    raise AssertionError('no level-3 lift found for %d (%s side)' % (dd, side))

def analyse3(name, side, C3, dd, hints=None):
    """rank_fisher5.analyse (steps (1)-(8), every identity checked exactly) on a level-3 curve, with the conic
    points of `hints` ({'<name>|<side>': .., '<name>|<side>|V': ..}) put into rank_fisher5.HINTS."""
    import rank_fisher5 as R5
    for k in ('%s|%s' % (name, side), '%s|%s|V' % (name, side)):
        R5.HINTS.pop(k, None)
        if hints and k in hints: R5.HINTS[k] = hints[k]
    return R5.analyse(name, side, C3, dd=dd)

def find_conic(R, PAIR, key, verbose=True, big=None):
    """FINDER (PARI; not certifying): a point of the cone conic Gamma_1 (PAIR 'U', over Q(sqrt disc U)) or Gamma_3
    ('V', over Q(sqrt disc V)) of rank_fisher5.analyse's output R.  The body of rank_fisher_theta4_find.find, taking
    R instead of a curve index of rank_fisher_curves.json (otherwise unchanged); the caller re-checks the point
    exactly (rank_fisher5.analyse)."""
    import time
    import sympy as sp
    import rank_fisher5 as F
    import rank_fisher_theta4_find as TF
    import conic_lattice_find as CL
    import cypari2
    pari_ = TF.pari
    t = TF.t
    if PAIR == 'U': h0, Bq, v1, v2, dd = R['h0'], R['Bq'], F.e, F.g, R['dU']
    else: h0, Bq, v1, v2, dd = R['coneV'][3], R['coneV'][4], F.c, F.d, R['dV']
    D0, m = TF.sqfree_split(dd)
    K = TF.QF(D0)
    t0 = time.time()
    def topari(expr):
        P_ = sp.Poly(sp.expand(expr), t)
        return K.el(P_.coeff_monomial(1), '(%s)*%d' % (P_.coeff_monomial(t), m))
    P = sp.Poly(Bq, v1, v2)
    b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (v1 ** 2, v1 * v2, v2 ** 2)]
    h0p = topari(h0)
    assert b20 != 0
    A1 = -4 * b20 * h0p; B1 = b11 ** 2 - 4 * b20 * b02
    try:
        A2, sa = K.sqred(A1); B2, sb = K.sqred(B1)
    except cypari2.PariError:
        A2, sa = K.sqred_rational(A1); B2, sb = K.sqred_rational(B1)
    nA, nB = abs(int(pari_.norm(A2))), abs(int(pari_.norm(B2)))
    if verbose: print('  finder %s %s: L = Q(sqrt %d), class group %s; N(A2) %d digits, N(B2) %d digits (%.1fs)'
                      % (key, PAIR, D0, K.bnf.bnf_get_cyc(), len(str(nA)), len(str(nB)), time.time() - t0), flush=True)
    if big: pt = CL.find_point(D0, A2, B2, tries=big, nvec=20000, Rmul=40, npair=2500, verbose=verbose)
    else: pt = CL.find_point(D0, A2, B2, tries=20, nvec=3000, Rmul=8, npair=400, verbose=verbose)
    if pt is None: return None
    E_, Zr, gr = pt
    E_ = pari_.Mod(pari_.lift(E_), K.pol); Zr = pari_.Mod(pari_.lift(Zr), K.pol); gr = pari_.Mod(pari_.lift(gr), K.pol)
    assert E_ ** 2 - A2 * Zr ** 2 == B2 * gr ** 2
    Z = Zr / sa; gg = gr / sb
    assert E_ ** 2 - A1 * Z ** 2 == B1 * gg ** 2
    ee = (E_ - b11 * gg) / (2 * b20)
    assert h0p * Z ** 2 + b20 * ee ** 2 + b11 * ee * gg + b02 * gg ** 2 == 0
    def out(z):
        z = z.lift(); return [str(pari_.polcoef(z, 0, 'y')), str(pari_.polcoef(z, 1, 'y'))]
    if PAIR == 'U': hint = dict(D0=D0, m=m, Z=out(Z), e=out(ee), g=out(gg))
    else: hint = dict(D0=D0, m=m, Z=out(Z), c=out(ee), d=out(gg))
    if verbose: print('  point found (%.1fs)' % (time.time() - t0), flush=True)
    return hint

def find_conic_norm(R, PAIR, target=10 ** 11, log=print, seeds=range(1, 12), rnf_time=60):
    """FINDER (PARI under GRH; not certifying): rank_fisher_theta4_norm.find_hint's body (class-group reduction, a
    Legendre-type descent over L, rnfisnorm, quaternion reduction of the point), taking rank_fisher5.analyse's
    output R instead of the json's (945,248) curve; otherwise unchanged.  The caller re-checks the point exactly."""
    import time
    import sympy as sp
    import rank_fisher5 as R5
    import rank_fisher_theta4_norm as NM
    from cysignals.alarm import alarm, AlarmInterrupt, cancel_alarm
    pari_ = NM.pari; t_ = NM.t_
    t0 = time.time()
    if PAIR == 'U': h0, Bq, v1, v2, dd = R['h0'], R['Bq'], R5.e, R5.g, R['dU']
    else: h0, Bq, v1, v2, dd = R['coneV'][3], R['coneV'][4], R5.c, R5.d, R['dV']
    D0, m = NM.sqfree_split(dd)
    def topari(expr):
        P_ = sp.Poly(sp.expand(expr), t_)
        return pari_('Mod((%s) + (%s)*%d*y, y^2 - (%d))' % (P_.coeff_monomial(1), P_.coeff_monomial(t_), m, D0))
    P = sp.Poly(Bq, v1, v2)
    b20, b11, b02 = [topari(P.coeff_monomial(mm)) for mm in (v1 ** 2, v1 * v2, v2 ** 2)]
    h0p = topari(h0)
    Cn = dict(D0=D0, m=m, b20=b20, b11=b11, b02=b02, h0=h0p, A1=-4 * b20 * h0p, B1=b11 ** 2 - 4 * b20 * b02)
    K = NM.Fld(Cn['D0'])
    log('  norm finder %s: L = Q(sqrt %d), class group %s (GRH)' % (PAIR, D0, K.bnf.bnf_get_cyc()))
    A1, B1 = Cn['A1'], Cn['B1']
    known = set()
    for z in (A1, B1):
        zl = pari_.lift(z); den = pari_.denominator(pari_.content(zl))
        known |= set(NM.factor_int(pari_.norm(pari_.Mod(zl * den * den, K.pol))))
    known = tuple(sorted(known))
    sol = None
    for seed in seeds:
        A, B, steps, nA, nB = NM.descend(K, A1, B1, target=target, known=known, seed=seed)
        log('  %s, descent seed %d: N(J_A), N(J_B) of %d, %d digits after %d steps' % (PAIR, seed, len(str(nA)), len(str(nB)), len(steps) - 1))
        if nA * nB > target: continue
        swap = nB < nA
        al, be = (B, A) if swap else (A, B)
        alarm(rnf_time)
        try:
            r = pari_.rnfisnorm(pari_.rnfisnorminit(K.pol, pari_('x^2') - pari_.lift(al), 2), pari_.lift(be), 0)
        except AlarmInterrupt:
            log('  %s: rnfisnorm over %d s, next seed' % (PAIR, rnf_time)); continue
        finally:
            cancel_alarm()
        if r[1] != 1: continue
        a = pari_.lift(r[0])
        a0 = K.el(pari_.polcoef(a, 0, 'x')); a1 = K.el(pari_.polcoef(a, 1, 'x'))
        P0 = (a0, K.el(1), a1) if swap else (a0, a1, K.el(1))
        if NM.digits(P0) > 20000: continue
        sol = P0; break
    if sol is None: return None
    assert sol[0] ** 2 - A * sol[1] ** 2 - B * sol[2] ** 2 == 0
    X, Y, Z = NM.unwind(K, steps[1:], NM.small_points(K, A, B, sol)[0])
    sA, sB = steps[0]['sA'], steps[0]['sB']
    X, Y, Z = NM.small_points(K, A1 / sA ** 2, B1 / sB ** 2, (X, Y, Z))[0]
    X, Y, Z = NM.reduce_point(K, (X, Y / sA, Z / sB))
    assert X ** 2 - A1 * Y ** 2 - B1 * Z ** 2 == 0
    ee = (X - Cn['b11'] * Z) / (2 * Cn['b20'])
    assert Cn['h0'] * Y ** 2 + Cn['b20'] * ee ** 2 + Cn['b11'] * ee * Z + Cn['b02'] * Z ** 2 == 0
    out = lambda z: [str(pari_.polcoef(pari_.lift(z), 0, 'y')), str(pari_.polcoef(pari_.lift(z), 1, 'y'))]
    log('  point found (%.1fs)' % (time.time() - t0))
    if PAIR == 'U': return dict(D0=Cn['D0'], m=Cn['m'], Z=out(Y), e=out(ee), g=out(Z))
    return dict(D0=Cn['D0'], m=Cn['m'], Z=out(Y), c=out(ee), d=out(Z))

# ---------------------------------------------------------------- Theta_3 on one side (type A)
def setup3(R, pair):
    """the rational data of one pushout form on C3' (rank_fisher5 step (6) for pair 'U', step (7) for 'V'):
    QA, QB integral; f = N(l_1) (or f_dag = N(l_3)), f0 = f / s0 primitive integral."""
    import sympy as sp
    import rank_fisher_theta4 as T
    VS = T.VS
    QA, QB = sp.expand(R['QA']), sp.expand(R['QB'])
    for Q_ in (QA, QB):
        assert all(sp.Rational(k).q == 1 for k in sp.Poly(Q_, *VS).coeffs())
    po = R['pushout'] if pair == 'U' else R['pushout3']
    assert po is not None, 'no pushout form for the pair %s' % pair
    f = sp.expand(po['f'])
    s0, f0 = T.poly_int(f, VS)
    return dict(QA=QA, QB=QB, f=f, f0=f0, s0=Fr(int(s0.p), int(s0.q)), pair=pair,
                field=T.sqfree(int(R['dU'] if pair == 'U' else R['dV'])))

def theta3(C, Rd, C3, R, pair, ys, say=print, nseeds=(1, 2), side='E', check_selmer=True):
    """Theta_3(y, z) = sum_v (f(P_v), z)_v (Fisher Sect. 3) on the level-3 lift C3' of y, f the pushout of the
    given pair (Lemma 7.1, checked by rank_fisher5.analyse), P_v in C3'(Q_v), for each z in ys (E-side classes,
    certified in S'_1).  Places: S = {R} u {p < 47} u B0 u {f0 mod p in span(QA, QB)} u primes(ys), B0 = {2} u bad
    primes of E, E' u primes of the pencil discriminant (factored completely).  At p outside S, C3' has good reduction,
    f0 mod p is a non-zero form on it (8 zeros at most), a point of C3'(F_p) off them lifts (p >= 47) and there f0 and
    z are units: the term is 1 (written argument).  Returns ({z: product}, S, classes)."""
    import rank_fisher_theta4 as T
    import rank_fisher_theta4_gen as G
    D = setup3(R, pair)
    AA, BB, loc = (C['a'], C['b'], Rd['locS']) if side == 'E' else (C['ap'], C['bp'], Rd['locSp'])
    for z in (ys if check_selmer else []):
        for p in Rd['bad']:
            assert R2.certified(z, AA, BB, p, loc), '%d is not in the local image at %d (%s side)' % (z, p, side)
        assert RI.real_ok(z, AA, BB)
    nd, ps = G.pencil_primes(D, known=sorted(set([2] + Rd['bad'])))
    B0 = sorted(set([2] + Rd['bad'] + ps))
    fbad = G.f0_degenerate_primes(D)
    small = [p for p in range(2, 47) if pari.isprime(p)]
    S = sorted(set(small + B0 + fbad + [q for z in ys for q in T.primes_of(z)]))
    say('  pushout f of the %s pair (over Q(sqrt %d)), f = s0 f0 with f0 primitive; pencil discriminant %d digits, '
        'primes %s; f0 mod p in span(QA, QB) only for p in %s' % (pair, D['field'], nd, ps, fbad))
    say('  S = {R} u %s' % S)
    M = T.Model(D['QA'], D['QB'])
    terms, classes = {}, {}
    for p in S:
        cls = set()
        for sd in list(nseeds) + [3, 4, 5, 6, 7, 8]:
            for P in T.c3_points(M, C3, p, 100, sd):
                v = T.ev_poly(D['f0'], P)
                if v == 0: continue
                try: cls.add(T.padic_class(v, p))
                except AssertionError: continue
            if len(cls) >= 2 and sd >= nseeds[-1]: break
        assert cls, 'no point of C3\'(Q_%d) found' % p
        terms[p] = {}
        for z in ys:
            vals = {T.hilbert(k, T.int_class(z, p), p) for k in cls}
            assert len(vals) == 1, 'the term (f0(P), %d)_%d depends on the point P: %s (f is not a pushout for this pairing)' % (z, p, vals)
            terms[p][z] = vals.pop()
        classes[p] = sorted(cls)
    Rs = T.real_signs(M, D, C=C3)
    assert Rs, 'no real point of C3\' found'
    assert all(z > 0 for z in ys) or len(set(Rs)) == 1, 'the real term is not determined'
    rt = {z: (-1 if (z < 0 and set(Rs) == {-1}) else 1) for z in ys}
    p4 = next((p for p in S if len(classes[p]) > 1), None)
    assert p4 is not None, '(I4) not shown: f0 might be a constant times a square'
    say('  R: real points of C3\' (interval enclosures), signs of f0 there %s; (I4) at p = %d f0 takes the classes %s'
        % (sorted(Rs), p4, classes[p4]))
    out = {}
    for z in ys:
        prod = rt[z]
        for p in S: prod *= terms[p][z]
        out[z] = prod
        nontriv = [p for p in S if terms[p][z] == -1] + (['R'] if rt[z] == -1 else [])
        say('  sum_v (f(P_v), %d)_v: the term is -1 exactly at %s; product = %+d' % (z, nontriv if nontriv else 'no place', prod))
    for p in S:
        say('    p = %d: classes of f0 on C3\'(Q_p) %s; terms %s' % (p, classes[p], ' '.join('%d:%+d' % (z, terms[p][z]) for z in ys)))
    return out, S, classes

# ---------------------------------------------------------------- Theta_4 (type B), through rank_fisher_theta4_gen
def theta4_gen_module(C, name, side, other, dd_lift, y_other, C3, dims, hints_path):
    """rank_fisher_theta4_gen.py, unchanged, with its curve source redirected (module attributes replaced, no file
    edited): curve_index -> this curve's record (the json format of rank_fisher_curves.json: theta4_reduction,
    sides[side].level3_curve, sides[...].dd), base_curve -> (a, b, G) of this curve, HINTS_GEN -> t11_hints.json,
    DIMS[name] -> (dim S_4, dim S'_4) of this curve (used only in the printed conclusion)."""
    import rank_fisher_theta4_gen as G
    Gpt = next(P for P in C['pts'] if P[0] != 0)
    cv = dict(name=name, theta4_reduction=[side, other],
              sides={side: dict(dd=str(dd_lift), level3_curve=C3), other: dict(dd=str(y_other))})
    G.curve_index = lambda nm: (0, cv)
    G.base_curve = lambda cvname: (C['a'], C['b'], (Fr(Gpt[0]), Fr(Gpt[1])))
    G.HINTS_GEN = hints_path
    G.DIMS[name] = dims
    G.SIDE = None
    G.OUT.clear()
    return G

def choose_Xi_lin_images(G, RD):
    """a replacement for rank_fisher_theta4_gen.choose_Xi_lin (a sampler that chooses Xi; nothing rests on it, the
    level-4 curve C4 is then checked locally soluble at every place of S with found points): the same F_2 linear
    algebra, but at each place the sampled classes of f0 on C3'(Q_p) are first closed up with the local image
    Im(delta_p) of the twisting side (Fisher Sect. 3: they form a coset of it): rank_indep's certified local image
    (RD['locS'], the first descent of the working curve) at its bad primes, the unit classes at an odd prime of good
    reduction.  Change from the original: only that closure (the original closes the sampled classes up to the
    coset they span, which is too small where the sampler finds one class only, as at p = 1332819889 for E1(8249,3720))."""
    T = G.T
    orig = G.choose_Xi_lin
    def choose(M, D, C, places, say=G.say):
        K = {p: T.f0_classes(M, D, C, p) for p in places}
        for p in places:
            sd = 2
            while not K[p] and sd <= 8:
                K[p] = T.f0_classes(M, D, C, p, seed=sd); sd += 1
        Rs = T.real_signs(M, D, C=C)
        if not all(K[p] for p in places) or not Rs: return [], K, Rs
        for p in places:
            if p in RD['locS']:
                poss, cert = RD['locS'][p]
                assert poss == cert, 'local image at %d not decided' % p
                im = [T.int_class(dl, p) for dl in RI.class_reps(p) if RI.qp_class(dl, p) in cert]
            else:
                assert p != 2
                im = [(0, 0), (0, 1)]
            k0 = T.cls_vec(next(iter(K[p])), p)
            add = {T.vec_cls(tuple((x + y) % 2 for x, y in zip(k0, T.cls_vec(i, p))), p): None for i in im}
            K2 = dict(K[p]); K2.update(add); K[p] = T.coset_closure(K2, p)
            assert len(K[p]) == len(im), ('the sampled classes of f0 at %d are not one coset of the local image' % p, len(K[p]), len(im))
        gens = [-1] + list(places)
        n = len(gens)
        rows, rhs = [], []
        for p in places:
            ks = [T.cls_vec(k, p) for k in K[p]]
            k0 = ks[0]
            span = []
            for k in ks[1:]:
                dv = [(x + y) % 2 for x, y in zip(k, k0)]
                for b in span:
                    lead = next(i for i, x in enumerate(b) if x)
                    if dv[lead]: dv = [(x + y) % 2 for x, y in zip(dv, b)]
                if any(dv): span.append(dv)
            dim = len(k0)
            ann = G.gf2_solve(span, [0] * len(span), dim) if span else [[1 if i == j else 0 for i in range(dim)] for j in range(dim)]
            basis = []
            for a in ann:
                v = list(a)
                for b in basis:
                    lead = next(i for i, x in enumerate(b) if x)
                    if v[lead]: v = [(x + y) % 2 for x, y in zip(v, b)]
                if any(v): basis.append(v)
            gv = [T.cls_vec(T.int_class(q_, p), p) for q_ in gens]
            for a in basis:
                rows.append([sum(ai * gi for ai, gi in zip(a, gv[j])) % 2 for j in range(n)])
                rhs.append(sum(ai * ki for ai, ki in zip(a, k0)) % 2)
        rs = set(Rs)
        if rs == {1}: rows.append([1] + [0] * (n - 1)); rhs.append(0)
        elif rs == {-1}: rows.append([1] + [0] * (n - 1)); rhs.append(1)
        sols = []
        for e_ in G.gf2_solve(rows, rhs, n):
            Xi = 1
            for j, q_ in enumerate(gens):
                if e_[j]: Xi *= q_
            sols.append(Xi)
        return sorted(sols, key=lambda z: (abs(z), z)), K, Rs
    choose.orig = orig
    return choose

def recertify_xi(C3, R, side, s, say=print, rounds=4):
    """another member xi' = xi s (s in this side's S_1, the twisting group of the level-3 Selmer set) of the level-3
    Selmer set of the same D_lam: local points P_v of D_lam with xi' f(P_v) a non-zero square at every place of S3
    (rank_fisher's samplers; each point checked exactly by rank_fisher.check_C3_point), returned as a new C3 record
    (xi', its local points).  The places S3 and the argument outside S3 are rank_fisher.py's (step (4))."""
    lam = int(C3['lam']); U = [int(x) for x in C3['U']]; V = [int(x) for x in C3['V']]
    AA_S1 = R['Sp'] if side == "E'" else R['S']
    assert RI.sqfree_part(s) in AA_S1, 's is not in this side\'s S_1'
    xi2 = RI.sqfree_part(int(C3['xi']) * s)
    # the D_lam of the record (its bad primes from rank_indep's second descent are not stored: the places of S3
    # are those of the record's local points, which rank_fisher.certify_C3 used)
    D = R2.DLam(lam, U, V)
    assert [str(x) for x in D.X['l']] == C3['l_X'] and [str(x) for x in D.Y['l']] == C3['l_Y'], 'the record is not this D_lam'
    places = [('inf' if v == 'inf' else int(v)) for v in C3['local_points']]
    assert all(p in places for p in RI.factor_int(abs(xi2))), 'xi\' has a prime outside S3'
    rng = random.Random(abs(xi2) % 99991)
    cert = {}
    for v in places:
        smp = RF.sample(D, v, rng)
        for rd in range(rounds):
            cv = R2.sq_cls(xi2, v)
            hit = [P for (c_, P) in smp if c_ == cv]
            if hit: cert[v] = hit[0]; break
            smp += RF.sample_tree(D, v, rng) if (rd == 0 and v != 'inf') else RF.sample(D, v, rng, 3 * RF.NSAMP)
        assert v in cert, 'no sampled point of D_lam with xi\' f(P) a square at %s' % v
        assert RF.check_C3_point(D, xi2, v, cert[v])
    say('  level-3 twist xi\' = %d = xi * %d (mod squares): C3 checked locally soluble at R and %d primes' % (xi2, s, len(places) - 1))
    return dict(C3, xi=str(xi2), local_points={str(v): [str(x_) for x_ in P] for v, P in cert.items()})
