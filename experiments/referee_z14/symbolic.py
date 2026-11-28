"""z14 tasks 2 and 3: the algebra of EXPLORATION_3 §1-§2 and of EXPLORATION_4 §5 (tangent_norm.py), by sympy.
Written by an AI model (reviewer z14); not read by a person.   usage: python3 -I symbolic.py | diff - symbolic.txt"""
from sympy import symbols, expand, factor, simplify, cancel, I, Poly, degree, sqrt, together, numer, denom, gcd, cyclotomic_poly, QQ, roots, minimal_polynomial, discriminant, resultant, zeta
m, n, r, t, x, T, Z, d = symbols('m n r t x T Z d')
A, B, C = m**2 - n**2, 2*m*n, m**2 + n**2
Rp, Rm = C**2 + A*B, C**2 - A*B
W = expand(Rp * Rm)
def say(*a): print(*a)
def disc_cubic(e1, e2, e3): return 16 * ((e1 - e2) * (e1 - e3) * (e2 - e3))**2
def c4_cubic(e1, e2, e3):
    # y^2 = (x - e1)(x - e2)(x - e3) = x^3 + a2 x^2 + a4 x + a6: c4 = 16 (a2^2 - 3 a4)
    a2 = -(e1 + e2 + e3); a4 = e1*e2 + e1*e3 + e2*e3
    return 16 * (a2**2 - 3*a4)
def j_cubic(e):
    return cancel(c4_cubic(*e)**3 / disc_cubic(*e))
say('# z14 symbolic checks (sympy). Written by an AI model (reviewer z14); not read by a person.')
say('W = R+ R- =', W, '; equals m^8 + 14 m^4 n^4 + n^8:', expand(W - (m**8 + 14*m**4*n**4 + n**8)) == 0)
say('W = A^4 - A^2 C^2 + C^4 = C^4 Phi_12(A/C):', expand(W - (A**4 - A**2*C**2 + C**4)) == 0, '; Phi_12 =', cyclotomic_poly(12, x))
Wr = Poly(expand(W.subs({m: r, n: 1})), r)
for nm, R_ in (('R+', Rp), ('R-', Rm)):
    Rr = Poly(expand(R_.subs({m: r, n: 1})), r)
    say('%s(r,1) =' % nm, Rr.as_expr(), '; irreducible over Q:', Rr.is_irreducible, '; linear factors over Q(i, sqrt 3) = Q(zeta_12):', len(factor(Rr.as_expr(), extension=[I, sqrt(3)]).as_ordered_factors()))
say('  so (Q(zeta_12) abelian, the root field of R+-) an odd prime p | R+-(m,n), p not dividing mn, splits completely in Q(zeta_12):')
say('  p = 1 mod 12 (not mod 24).  R+ - R- = 2AB =', factor(2*A*B), '= 0 mod 8 when mn is even, so W = R+^2 = 1 mod 8')
say('  (odd-odd ratios: the same after removing 4^2); hence the Jacobi symbols (-1/W0), (2/W0) of the odd-exponent part W0')
say('  of W are 1 for every ratio, and the controls S*J(-1), S*J(2) of task23.py cannot fail')
# §1: the mixed curve and the conic
e_mix = (-2*C**2, -2*A*B, 2*A*B)
say('mixed curve e2 - e1 = 2R-:', expand(e_mix[1] - e_mix[0] - 2*Rm) == 0, '; e3 - e1 = 2R+:', expand(e_mix[2] - e_mix[0] - 2*Rp) == 0, '; e3 - e2 = 4AB:', expand(e_mix[2] - e_mix[1] - 4*A*B) == 0)
Aa, Bb, Cc = symbols('A B C')
s_, D_ = d*T + 4*Aa*Bb*Z, d*T - 4*Aa*Bb*Z
ident = expand((Cc**2 + Aa*Bb) * s_**2 - (Cc**2 - Aa*Bb) * D_**2 - 2*Aa*Bb*d*(d*T**2 + 8*Cc**2*T*Z + 16*Aa**2*Bb**2/d*Z**2))
say('conic identity R+ s^2 - R- D^2 = 2AB d (d T^2 + 8C^2 T Z + 16A^2B^2/d Z^2) (A, B, C free):', ident == 0)
say('2AB R+ = R+^2 - W (so 2AB R+ is a norm from Q(sqrt W) and the conic is soluble iff d is):', expand(2*A*B*Rp - (Rp**2 - W)) == 0)
# the covering of the class (1, d, d): x - e1 = u^2 on E_mix gives the 2-isogenous curve's coefficients
a2_, b2_ = -2*(Rp + Rm), 4*W
say('E_mix shifted, x(x^2 + a x + b): a = -4C^2:', expand(a2_ + 4*C**2) == 0, '; b = 4W; isogenous a\' = -2a = 8C^2, b\' = a^2 - 4b = 16A^2B^2:', expand(a2_**2 - 4*b2_ - 16*A**2*B**2) == 0)
say('  -> w^2 = d t^4 + 8C^2 t^2 + 16A^2B^2/d is the quartic of the class d for the 2-isogenous curve Y^2 = X(X^2 + 8C^2 X + 16A^2B^2)')
# degenerate-ratio identities
a1, b1, a2, b2 = m*(m + n), n*(m - n), m*(m - n), n*(m + n)
say('alpha1 beta1 = alpha2 beta2 = AB/2:', expand(a1*b1 - A*B/2) == 0 and expand(a2*b2 - A*B/2) == 0,
    '; R+ = alpha1^2 + beta1^2:', expand(a1**2 + b1**2 - Rp) == 0, '; R- = alpha2^2 + beta2^2:', expand(a2**2 + b2**2 - Rm) == 0,
    '; m^2 + n^2 = alpha1 - beta1 = alpha2 + beta2:', expand(a1 - b1 - C) == 0 and expand(a2 + b2 - C) == 0)
# EXPLORATION_4 §5: tangent form as a norm
u1, u3, q1, q3, qT, c = symbols('u1 u3 q1 q3 qT c')
L = q1*u1 - d*q3*u3 - c*qT*T
Nw = (qT*u1 - T*q1)**2 - d*(qT*u3 - T*q3)**2
diff = expand(-2*qT*T*L - Nw)
diff = expand(diff.subs(u1**2, c*T**2 + d*u3**2).subs(q1**2, c*qT**2 + d*q3**2))
say('tangent lemma: -2 qT T L(P) = N(qT(u1 + sqrt d u3) - T(q1 + sqrt d q3)) on u1^2 - d u3^2 = c T^2 through Q:', diff == 0)
say('  tangent form of cassels.py on H_2 (a = (1,d,d)): (a1 q1, -a3 q3, -(e3 - e1) qT) = (q1, -d q3, -2R+ qT), i.e. c = 2R+; on H_3 c = 2R-')
# §2: the modular picture
e08 = (0, A**2, C**2); eLU = (0, B**2, C**2)
D08 = factor(disc_cubic(*e08)); DLU = factor(disc_cubic(*eLU)); Dmix = factor(disc_cubic(*e_mix))
say('Delta(E_08) =', D08)
say('Delta(E_LU) =', DLU)
say('Delta(E_mix) =', Dmix)
for nm, e, Dl in (('E_08', e08, D08), ('E_LU', eLU, DLU), ('E_mix', e_mix, Dmix)):
    Dr = Poly(expand(disc_cubic(*e).subs({m: r, n: 1})), r)
    c4r = expand(c4_cubic(*e).subs({m: r, n: 1}))
    a = [expand(-(e[0] + e[1] + e[2])), expand(e[0]*e[1] + e[0]*e[2] + e[1]*e[2]), expand(-e[0]*e[1]*e[2])]
    degs = [Poly(a_, m, n).total_degree() for a_ in a]
    # fibre types: at each root of Delta (and at r = infinity), ord of Delta and whether c4 vanishes there
    fac = factor(Dr.as_expr(), extension=[I])
    rows = []
    for f_, k_ in fac.as_powers_dict().items():
        if Poly(f_, r).degree() == 0: continue
        rts = roots(Poly(f_, r))
        for rt, mult in rts.items():
            c4v = simplify(c4r.subs(r, rt))
            rows.append((str(rt), k_ * mult, 'c4 != 0' if c4v != 0 else 'c4 = 0'))
        if not rts: rows.append(('roots of ' + str(f_), k_, 'c4 != 0' if resultant(Poly(f_, r), Poly(c4r, r)) != 0 else 'c4 = 0 somewhere'))
    ordinf = 24 - Dr.degree()
    rows.append(('infinity', ordinf, 'c4 != 0' if Poly(c4r, r).degree() == 8 else 'c4 = 0'))
    say('%s: (deg a2, deg a4, deg a6) = %s (K3 range: <= (4, 8, 12)); fibres (place, ord Delta, c4): %s; sum of ord Delta = %d'
        % (nm, degs, sorted(rows), sum(x[1] for x in rows)))
j08 = j_cubic(e08); jmix = j_cubic(e_mix)
say('j(E_08) = 256 (A^4 - A^2 C^2 + C^4)^3 / (A^4 B^4 C^4):', cancel(j08 - 256*W**3/(A**4*B**4*C**4)) == 0, '; so j(E_08) = 0 exactly at the 8 zeros of W, deg_r j(E_08) =', max(Poly(numer(together(j08.subs({m: r, n: 1}))), r).degree(), Poly(denom(together(j08.subs({m: r, n: 1}))), r).degree()))
say('j(E_LU) = j(E_08):', cancel(j_cubic(eLU) - j08) == 0)
jm = cancel(jmix)
say('j(E_mix) =', factor(jm), '-> at the zeros of W: j(E_mix) = infinity (I2 fibres), not 0; j(E_mix) = 0 where', factor(numer(together(jm))))
# twist: E_08 twisted by -1 is y^2 = x(x + A^2)(x + C^2); translate x -> x - C^2
tw = expand((x - C**2) * (x - C**2 + A**2) * x - x * (x - B**2) * (x - C**2))
say('E_08 twisted by -1 (y^2 = x(x+A^2)(x+C^2)), x -> x - C^2, is y^2 = x(x - B^2)(x - C^2) (the model of E_LU up to 2^2):', tw == 0)
# 4-torsion of E_08 over Q(i)(r): halves of the three 2-torsion points
f08 = lambda X: X * (X - A**2) * (X - C**2)
halves = {'(0,0)': [A*C, -A*C], '(A^2,0)': [A**2 + I*A*B, A**2 - I*A*B], '(C^2,0)': [C**2 + C*B, C**2 - C*B]}
ys = {'(0,0)': [I*A*C*(C - A), I*A*C*(A + C)], '(A^2,0)': [None, None], '(C^2,0)': [C*B*(C + B), I*C*B*(C - B)]}
for k_, xs in halves.items():
    out = []
    for X in xs:
        y2 = factor(expand(f08(X)))
        out.append(str(y2))
    say('halves of %s on E_08: x = %s, y^2 = %s' % (k_, [str(factor(X)) for X in xs], out))
yA = factor(expand(f08(A**2 + I*A*B)))
cand = expand((A**2 + I*A*B) * (I*A*B) * (A**2 + I*A*B - C**2))
say('  y^2 at x = A^2 + iAB is the square of 2i mn(m^2 - n^2)(m + in)^2:', expand(cand - (2*I*m*n*(m**2 - n**2)*(m + I*n)**2)**2) == 0)
say('  (0,0): y = i A C (C - A) checks:', expand((I*A*C*(C - A))**2 - f08(A*C)) == 0)
# the four corner forms vanish on graphs of automorphisms
F0 = m**2 + n**2 * t**2; F8 = n**2 + m**2 * t**2
FL = (m + n)**2 * (t - 1)**2 + (m - n)**2 * (t + 1)**2; FU = (m - n)**2 * (t - 1)**2 + (m + n)**2 * (t + 1)**2
phi = lambda z: (z - 1) / (z + 1); phinv = lambda z: (1 + z) / (1 - z)
for nm, F, g in (('F_0', F0, lambda z: I * z), ('F_8', F8, lambda z: I / z), ('F_L', FL, lambda z: phinv(I * phi(z))), ('F_U', FU, lambda z: phinv(I / phi(z)))):
    ok = all(simplify(F.subs({m: r, n: 1}).subs(t, gg)) == 0 for gg in (g(r), g(r).subs(I, -I)))
    say('%s(r, t) vanishes on t = g(r) and on its complex conjugate, g = %s:' % (nm, simplify(g(r))), ok)
say('Delta(E_mix) at r = +-i (C = 0): ', simplify(expand(disc_cubic(*e_mix)).subs({m: I, n: 1})), '(non-zero: the cusps +-i carry smooth fibres of E_mix)')
