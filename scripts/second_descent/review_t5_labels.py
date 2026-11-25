#!/usr/bin/env python3
"""Hostile review of branch sess_t5 (slice v5, run 3): is SECOND_DESCENT.md section 6 right that the four quartics with
2-Selmer bound 2 over Q(i) for (627,364) are q+ L1 L3, q- L1 L3, q+ L2 L4, q- L2 L4, "curves defined over Q"?

    python3 review_t5_labels.py

Written by an AI model (Claude); not read by a person.  sympy only (no Sage, no PARI); reads kdist_627_364.txt.

1. The labels of kclose.factors (copied here as data, not imported: kclose needs Sage):
     L1 = m - n s, L2 = m + n s, L3 = n + m s, L4 = n - m s, L5 = 1 - s, L6 = 1 + s,
     M1 = p s - pb, M2 = pb s + p, M3 = p s + pb, M4 = pb s - p,   p = m + i n, pb = m - i n,
   against SECOND_DESCENT.md section 2:  q+ = S s^2 - 2 D s + S,  q- = S s^2 + 2 D s + S  (S = m^2 + n^2, D = m^2 - n^2).
   Which products Mi Mj are q+, q-, and which are rational in s, or in t (s = -i t)?
2. For each of the 210 quartics (products of four of the ten factors) whether it is a constant times a polynomial with
   rational coefficients in s, or in t; and for (627,364) the j-invariant of its Jacobian (from the invariants I, J of the
   quartic; j = 6912 I^3 / (4 I^3 - J^2), independent of the twist), exactly in Q(i).
3. The bound-2 quartics of kdist_627_364.txt and the bounds of q+ L1 L3 = ?, etc., read from that file.
Negative test: the same j test applied to q+ L1 L3 must give j in Q (it is a curve over Q in s).
"""
import os, sys, itertools
import sympy as sp

here = os.path.dirname(os.path.abspath(__file__))
I = sp.I
s, t, m, n = sp.symbols('s t m n')

def factors(m, n):
    p = m + I*n; pb = m - I*n
    return {'L1': m - n*s, 'L2': m + n*s, 'L3': n + m*s, 'L4': n - m*s, 'L5': 1 - s, 'L6': 1 + s,
            'M1': p*s - pb, 'M2': pb*s + p, 'M3': p*s + pb, 'M4': pb*s - p}

def say(*a): print(*a, flush=True)

def rational_up_to_constant(expr, var):
    """Is expr = c * g(var) with g in Q[var] (c a constant of Q(i))?"""
    P = sp.Poly(sp.expand(expr), var)
    co = P.all_coeffs()
    lead = next(c for c in co if c != 0)
    for c in co:
        q = sp.expand(sp.radsimp(sp.expand(c) / sp.expand(lead)))
        re_, im_ = q.as_real_imag()
        if im_ != 0: return False
    return True

def j_of_quartic(Q, var):
    P = sp.Poly(sp.expand(Q), var); co = P.all_coeffs()
    while len(co) < 5: co = [0] + co
    a, b, c, d, e = co
    Iv = sp.expand(12*a*e - 3*b*d + c*c)
    Jv = sp.expand(72*a*c*e + 9*b*c*d - 27*a*d*d - 27*e*b*b - 2*c**3)
    den = sp.expand(4*Iv**3 - Jv**2)
    assert den != 0
    j = sp.expand(sp.radsimp(6912*Iv**3 / den))
    jr, ji = j.as_real_imag()
    return sp.Rational(jr), sp.Rational(ji)

def main():
    S = m*m + n*n; D = m*m - n*n
    F = factors(m, n)
    qp = S*s*s - 2*D*s + S; qm = S*s*s + 2*D*s + S
    GL = 2*S*t*t - 8*m*n*t + 2*S
    say("1. Products of two M's (general m, n), with s = -i t for the t-column:")
    names_M = ['M1', 'M2', 'M3', 'M4']
    for a_, b_ in itertools.combinations(names_M, 2):
        pr = sp.expand(F[a_]*F[b_])
        tags = []
        if sp.expand(pr - qp) == 0: tags.append('= q+')
        if sp.expand(pr - qm) == 0: tags.append('= q-')
        prt = sp.expand(pr.subs(s, -I*t))
        if sp.expand(prt + GL/2) == 0: tags.append('= -GL/2 in t')
        if sp.expand(prt - GL/2) == 0: tags.append('= GL/2 in t')
        say(f"   {a_}{b_} = {sp.factor(pr)}   {' '.join(tags)}")
    assert sp.expand(F['M1']*F['M4'] - qp) == 0 and sp.expand(F['M2']*F['M3'] - qm) == 0
    assert sp.expand((F['M1']*F['M2']).subs(s, -I*t) + GL/2) == 0
    say("   => with kclose's labels q+ = M1 M4 and q- = M2 M3 (identities in Q(i)[m, n, s]); M1 M2 = -GL/2 at s = -i t.")
    L1L3t = sp.expand((F['L1']*F['L3']).subs(s, -I*t))
    say(f"   L1 L3 at s = -i t: {L1L3t}")

    mm, nn = 627, 364
    Fn = factors(sp.Integer(mm), sp.Integer(nn))
    names = sorted(Fn)
    say(f"2. ({mm},{nn}): the 210 quartics, rational (up to a constant) in s / in t, and j of the bound-2 ones")
    # kdist table
    tab = {}
    for line in open(os.path.join(here, 'kdist_627_364.txt')):
        w = line.split()
        if len(w) >= 2 and w[0][0] in 'LM' and w[1].startswith('['):
            tab[w[0]] = eval(' '.join(w[1:]))
    assert len(tab) == 210, len(tab)
    two = sorted(k for k, v in tab.items() if max(v) == 2)
    say(f"   kdist_627_364.txt: quartics with larger bound 2: {two}")
    claimed = {'q+L1L3': 'L1L3M1M4', 'q-L1L3': 'L1L3M2M3', 'q+L2L4': 'L2L4M1M4', 'q-L2L4': 'L2L4M2M3'}
    for k, v in claimed.items():
        say(f"   {k} = {v} (kclose labels, sorted): bound {tab[v]}")
    cnt_s = cnt_t = 0
    rat = {}
    for J in itertools.combinations(names, 4):
        Q = sp.expand(sp.prod([Fn[k] for k in J]))
        rs = rational_up_to_constant(Q, s)
        rt = rational_up_to_constant(sp.expand(Q.subs(s, -I*t)), t)
        rat[''.join(J)] = (rs, rt)
        cnt_s += rs; cnt_t += rt
    say(f"   rational in s: {cnt_s} of 210; rational in t: {cnt_t} of 210")
    for k in two:
        rs, rt = rat[k]
        Q = sp.expand(sp.prod([Fn[x] for x in (k[0:2], k[2:4], k[4:6], k[6:8])]))
        jr, ji = j_of_quartic(Q, s)
        say(f"   {k}: rational in s: {rs}, in t: {rt}; j = {float(jr):.6g} + {float(ji):.6g} i  (Im j = 0: {ji == 0})")
        assert not rs and not rt and ji != 0
    for k, v in claimed.items():
        Q = sp.expand(sp.prod([Fn[x] for x in (v[0:2], v[2:4], v[4:6], v[6:8])]))
        jr, ji = j_of_quartic(Q, s)
        say(f"   {k} = {v}: rational in s: {rat[v][0]}; j in Q: {ji == 0}  (negative control: must be True)")
        assert rat[v][0] and ji == 0
    # are the four bound-2 quartics conjugate / related?
    js = {}
    for k in two:
        Q = sp.expand(sp.prod([Fn[x] for x in (k[0:2], k[2:4], k[4:6], k[6:8])]))
        js[k] = j_of_quartic(Q, s)
    say("   j-invariants of the four bound-2 quartics: " + ", ".join(f"{k}: {'conj of ' if False else ''}{js[k][0]}{'+' if js[k][1] >= 0 else '-'}{abs(js[k][1])}i"[:0] + f"{k}: Re j = {float(js[k][0]):.6g}, Im j = {float(js[k][1]):.6g}" for k in two))
    say("CONCLUSION: SECOND_DESCENT.md section 6 identifies L1L3M1M2, L1L3M3M4, L2L4M1M2, L2L4M3M4 with q+L1L3, q-L1L3,")
    say("  q+L2L4, q-L2L4; with kclose's labels q+ = M1M4, q- = M2M3, so the identification is wrong; the four bound-2")
    say("  quartics are rational neither in s nor in t and have j not in Q (not defined over Q in any coordinate);")
    say("  q+L1L3 = L1L3M1M4 etc. have bound 11 in the same table.")

if __name__ == '__main__':
    main()
