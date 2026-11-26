#!/usr/bin/env python3
"""Exact cross-check of LeanRank/CensusClosure.lean, independent of PARI and of the generator's code.

Written by an AI model; not read by a person.

For every pair section of the Lean file it parses p, the reduced curve, the count N, the list of
points, the list of affine points, the tau cases and the statements, and checks with integers and
Fractions only:
  1. p is an odd prime (trial division) not dividing 16(ab(a - b))^2, a = m^4, b = n^4;
  2. the reduced curve is (0, (a + b) mod p, 0, ab mod p, 0);
  3. #E~(F_p) = N, counted by brute force (O and the affine solutions);
  4. every listed point lies on y^2 = x(x + a)(x + b); the N - 1 affine points are distinct;
  5. the N points are closed under the chord-and-tangent law and under negation, so they form a
     subgroup of E(Q) of order N (consistent with the Lean bound #E(Q) <= #E~(F_p) = N);
  6. the affine list of `affine_m_n` is the point list of `pts_m_n`, in order;
  7. each tau case matches its abscissa: x < 0 'neg', x = 0 'zero', x = m^2n^2 'one', otherwise
     x/(m^2n^2) = P/Q in lowest terms, N = PQ, k = isqrt(N), k^2 < N < (k + 1)^2, x Q^2 = m^2n^2 N;
  8. the statements carry the right numbers (curve (-a) (-b), m^2 + n^2 tau^2, the ratio, c = m^2n^2);
  9. the pair is an exception of the census (exceptions1000 of the main library), and the file has
     a `#print axioms` line for each headline theorem;
 10. each partner theorem ratio_ne_M_N uses ratio_ne_partner with (M + N)·n = (M − N)·m for a checked
     pair (m, n), and (M, N) is coprime and in the census;
 11. (search, a sanity check, not a proof) for t = u/v with |u|, |v| <= 300, gcd = 1, the integers
     m^2 v^2 + n^2 u^2 and n^2 v^2 + m^2 u^2 are both squares only at t in {0, 1, -1} -- as tau_m_n says.
Then the negative test: 12 mutated copies of the parsed data (dropped point, duplicated point, point
off the curve, wrong N, bad prime, prime 2, wrong reduced coefficient, wrong tau case, square
passed as non-square, wrong bracket, wrong ratio, a pair outside the census) and 4 mutated partner
records (wrong relation, unchecked source pair, wrong statement, partner outside the census), and 4
edited copies of the Lean file (count 16 -> 15, a point moved off the curve, a tau case deleted, the
ratio of a theorem changed) must each be refused.

usage: python3 -I census_closure_check.py ../LeanRank/CensusClosure.lean ../../MagicSquareOfSquares/G1Census1000.lean
"""
import copy
import re
import sys
from fractions import Fraction as Fr
from math import gcd, isqrt


def is_prime(p):
    if p < 2:
        return False
    d = 2
    while d * d <= p:
        if p % d == 0:
            return False
        d += 1
    return True


def is_sq_int(x):
    return x >= 0 and isqrt(x) ** 2 == x


def is_sq(q):
    q = Fr(q)
    return q >= 0 and is_sq_int(q.numerator) and is_sq_int(q.denominator)


def num(s):
    s = s.strip()
    while s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()
    return int(s)


def parse(path):
    src = open(path).read()
    heads = list(re.finditer(r"/-! ### `\((\d+), (\d+)\)`: `Y² = X\(X \+ (\d+)\)\(X \+ (\d+)\)`, (\d+) points, "
                             r"reduction modulo (\d+) -/", src))
    secs = []
    for j, h in enumerate(heads):
        end = heads[j + 1].start() if j + 1 < len(heads) else src.index("end LeanRank.CensusClosure")
        body = src[h.start():end]
        m, n, a, b, N, p = map(int, h.groups())
        tag = f"{m}_{n}"
        ep = re.search(rf"def Ep_{tag} : WeierstrassCurve \(ZMod (\d+)\) := ⟨0, (\d+), 0, (\d+), 0⟩", body)
        cardp = re.search(rf"theorem cardp_{tag} : Fintype.card Ep_{tag}.toAffine.Point = (\d+) :=", body)
        ptsdef = body[body.index(f"noncomputable def pts_{tag}"):body.index(f"theorem pts_{tag}_injective")]
        finN = int(re.search(r": Fin (\d+) →", ptsdef).group(1))
        pts = [(num(x), num(y)) for x, y in
               re.findall(r"\.some (\(-?\d+\)|\d+) (\(-?\d+\)|\d+) \(ns \(by norm_num\)\)", ptsdef)]
        assert ptsdef.split("![", 1)[1].lstrip().startswith("0,"), "first entry must be O"
        aff = body[body.index(f"theorem affine_{tag}"):]
        affl = aff[aff.index("([") + 2:aff.index("] : List (ℚ × ℚ))")]
        affpts = [(num(x), num(y)) for x, y in re.findall(r"\((\(-?\d+\)|\d+), (\(-?\d+\)|\d+)\)", affl)]
        tau = body[body.index(f"theorem tau_{tag}"):body.index(f"theorem ratio_ne_{tag}")]
        th = re.search(r"\(h0 : r0 \^ 2 = (\d+) \+ (\d+) \* τ \^ 2\)\s+\(h8 : r8 \^ 2 = (\d+) \+ (\d+) \* τ \^ 2\)",
                       tau)
        cc = re.search(r"have hE : \((\d+) \* τ \* r0 \* r8\) \^ 2 =\s+\((\d+) \* τ \^ 2\) \* \(\((\d+) \* τ \^ 2\) "
                       r"\+ (\d+)\) \* \(\((\d+) \* τ \^ 2\) \+ (\d+)\)", tau)
        cases = []
        for line in tau.split("\n"):
            line = line.strip()
            if not line.startswith("· exact"):
                continue
            if "tau_zero" in line:
                cases.append(("zero",))
            elif "tau_one" in line:
                cases.append(("one",))
            elif "tau_neg" in line:
                cases.append(("neg",))
            else:
                g = re.search(r"\(N := (\d+)\) \(k := (\d+)\) \(Q := (\d+)\)", line)
                cases.append(("nonsq", int(g.group(1)), int(g.group(2)), int(g.group(3))))
        ratio = re.search(rf"theorem ratio_ne_{tag}[\s\S]*?: \(b \+ d\) \* (\d+) ≠ \(h \+ f\) \* (\d+) :=", body)
        curves = set(re.findall(r"curve \((-\d+)\) \((-\d+)\)", body))
        secs.append(dict(m=m, n=n, a=a, b=b, N=N, p=p, ep=tuple(map(int, ep.groups())),
                         cardp=int(cardp.group(1)), finN=finN, pts=pts, aff=affpts,
                         tauh=tuple(map(int, th.groups())), hE=tuple(map(int, cc.groups())), cases=cases,
                         ratio=tuple(map(int, ratio.groups())), curves=curves))
    axioms = set(re.findall(r"#print axioms LeanRank\.CensusClosure\.(\S+)", src))
    partners = []
    psrc = src[src.index("/-! ### Partner ratios -/"):] if "/-! ### Partner ratios -/" in src else ""
    for g in re.finditer(r"theorem ratio_ne_(\d+)_(\d+) \{a b c d e f g h i T : ℤ\}[\s\S]*?"
                         r": \(b \+ d\) \* (\d+) ≠ \(h \+ f\) \* (\d+) := by\n"
                         r"  have := ratio_ne_partner \(m := (\d+)\) \(n := (\d+)\) \(M := (\d+)\) \(N := (\d+)\)"
                         r"[\s\S]*?have := ratio_ne_(\d+)_(\d+) ", psrc):
        partners.append(tuple(int(x) for x in g.groups()))
    return secs, axioms, partners


def src_partner_count(path):
    return len(re.findall(r"ratio_ne_partner \(m :=", open(path).read()))


def add(P, Q, A, B):
    if P is None:
        return Q
    if Q is None:
        return P
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and y1 == -y2:
        return None
    lam = (y2 - y1) / (x2 - x1) if x1 != x2 else (3 * x1 * x1 + 2 * A * x1 + B) / (2 * y1)
    x3 = lam * lam - A - x1 - x2
    return (x3, lam * (x1 - x3) - y1)


def check_section(s, census):
    m, n, a, b, N, p = s["m"], s["n"], s["a"], s["b"], s["N"], s["p"]
    c = m * m * n * n
    assert (m, n) in census, "pair not in the census exceptions"
    assert a == m ** 4 and b == n ** 4
    assert s["curves"] == {(str(-a), str(-b))}, s["curves"]
    # 1, 2
    assert is_prime(p) and p != 2, "p must be an odd prime"
    assert (16 * (a * b * (a - b)) ** 2) % p != 0, "bad reduction"
    assert s["ep"] == (p, (a + b) % p, (a * b) % p), "reduced curve"
    # 3
    sq = [0] * p
    for y in range(p):
        sq[y * y % p] += 1
    cnt = 1 + sum(sq[(x ** 3 + (a + b) * x * x + a * b * x) % p] for x in range(p))
    assert cnt == N == s["cardp"] == s["finN"], ("count", cnt, N, s["cardp"], s["finN"])
    # 4
    pts = s["pts"]
    assert len(pts) == N - 1 and len(set(pts)) == N - 1, "points not distinct or wrong number"
    for (x, y) in pts:
        assert y * y == x * (x + a) * (x + b), ("off the curve", x, y)
    # 5
    G = {None} | {(Fr(x), Fr(y)) for (x, y) in pts}
    for P in G:
        if P is not None:
            assert (P[0], -P[1]) in G, "not closed under negation"
        for Q in G:
            assert add(P, Q, a + b, a * b) in G, "not closed under addition"
    # 6
    assert s["aff"] == pts, "affine list differs from the point list"
    # 7
    assert len(s["cases"]) == len(pts)
    for (x, _), cs in zip(pts, s["cases"]):
        if x < 0:
            assert cs == ("neg",), (x, cs)
        elif x == 0:
            assert cs == ("zero",), (x, cs)
        elif x == c:
            assert cs == ("one",), (x, cs)
        else:
            w = Fr(x, c)
            assert not is_sq(w)
            assert cs[0] == "nonsq", (x, cs)
            _, Nn, k, Q = cs
            assert k * k < Nn < (k + 1) * (k + 1), "bracket"
            assert Fr(x) * Q * Q == c * Nn, "x Q^2 = c N"
    # 8
    assert s["tauh"] == (m * m, n * n, n * n, m * m), s["tauh"]
    assert s["hE"] == (c, c, c, a, c, b), s["hE"]
    assert s["ratio"] == (n, m), s["ratio"]
    return True


def check_partner(t, done, census):
    """t = (pm, pn, n', m', m'', n'', M, N, m0, n0) parsed from a partner theorem."""
    pm, pn, n1, m1, m2, n2, M, N, m0, n0 = t
    assert (n1, m1) == (pn, pm) and (m2, n2) == (pm, pn), "statement and ratio_ne_partner disagree"
    assert (M, N) == (m0, n0), "the partner argument is not the theorem it uses"
    assert (m0, n0) in done, "the partner's own theorem is not a checked section"
    assert pm > pn > 0 and gcd(pm, pn) == 1
    assert (pm + pn) * N == (pm - pn) * M, "not the partner ratio"
    assert (pm, pn) in census, "partner pair not in the census exceptions"
    return True


def search(m, n, B=300):
    hits = set()
    for v in range(1, B + 1):
        for u in range(-B, B + 1):
            if gcd(u, v) != 1:
                continue
            if is_sq_int(m * m * v * v + n * n * u * u) and is_sq_int(n * n * v * v + m * m * u * u):
                hits.add(Fr(u, v))
    return hits


def refused(s, census):
    try:
        check_section(s, census)
    except AssertionError as e:
        return str(e) or "assertion"
    return None


def main():
    secs, axioms, partners = parse(sys.argv[1])
    src = open(sys.argv[2]).read()
    line = src[src.index("def exceptions1000"):].split("\n", 1)[0]
    census = {(int(x), int(y)) for x, y in re.findall(r"\((\d+), (\d+)\)", line)}
    assert len(census) == 425
    assert secs, "no sections parsed"
    for th in ["mss_forms", "mss_point", "ratio_ne_of_tau", "ratio_ne_partner"]:
        assert th in axioms, th
    for s in secs:
        tag = f"{s['m']}_{s['n']}"
        for th in ["points_", "card_", "affine_", "tau_", "ratio_ne_"]:
            assert th + tag in axioms, th + tag
        check_section(s, census)
        hits = search(s["m"], s["n"])
        assert hits <= {Fr(0), Fr(1), Fr(-1)}, hits
        print(f"({s['m']},{s['n']}): p = {s['p']}, #E~(F_p) = {s['N']} = #list, subgroup closed, "
              f"{len(s['cases'])} tau cases ok, search |u|,|v| <= 300 finds t in {sorted(hits)}")
    done = {(s["m"], s["n"]) for s in secs}
    assert len(partners) == src_partner_count(sys.argv[1]), "unparsed partner theorem"
    for t in partners:
        check_partner(t, done, census)
        assert f"ratio_ne_{t[0]}_{t[1]}" in axioms
        print(f"({t[0]},{t[1]}): partner of ({t[8]},{t[9]}), (m + n) N = (m - n) M, in the census")
    # negative test
    s0 = next(s for s in secs if (s["m"], s["n"]) == (4, 3))
    muts = []

    def mut(f):
        t = copy.deepcopy(s0)
        f(t)
        muts.append(t)
    mut(lambda t: t["pts"].pop())                                             # dropped point
    mut(lambda t: t["pts"].__setitem__(-1, t["pts"][-2]))                     # duplicated point
    mut(lambda t: t["pts"].__setitem__(5, (t["pts"][5][0], t["pts"][5][1] + 1)))   # off the curve
    mut(lambda t: t.update(N=t["N"] + 1, cardp=t["cardp"] + 1, finN=t["finN"] + 1))  # wrong N
    mut(lambda t: t.update(p=5, ep=(5, (256 + 81) % 5, (256 * 81) % 5)))    # bad prime
    mut(lambda t: t.update(p=2, ep=(2, 1, 0)))                               # prime 2
    mut(lambda t: t.update(ep=(t["ep"][0], t["ep"][1], (t["ep"][2] + 1) % t["p"])))  # wrong coefficient
    mut(lambda t: t["cases"].__setitem__(1, ("zero",)))                      # wrong tau case
    mut(lambda t: t["cases"].__setitem__(13, ("nonsq", 4, 1, 1)))            # N = 4 is a square
    mut(lambda t: t["cases"].__setitem__(3, ("nonsq", 6, 1, 6)))             # wrong bracket
    mut(lambda t: t.update(ratio=(4, 3)))                                    # wrong ratio
    mut(lambda t: t.update(m=5, n=4, a=625, b=256))                          # not in the census
    names = ["dropped point", "duplicated point", "point off the curve", "wrong N", "bad prime 5",
             "prime 2", "wrong reduced coefficient", "wrong tau case", "square as non-square",
             "wrong bracket", "wrong ratio", "pair outside the census"]
    pmuts = []
    if partners:
        t0 = partners[0]
        pmuts.append(("partner relation (18 : 7 for 12 : 5)", (18, 7, 7, 18, 18, 7) + t0[6:]))
        pmuts.append(("partner of a pair not checked", t0[:6] + (5, 2, 5, 2)))
        pmuts.append(("statement ratio differs", (t0[0], t0[1], t0[2] + 1) + t0[3:]))
        pmuts.append(("partner outside the census (census without it)", t0))
    for name, t in pmuts:
        cen = census - {(t[0], t[1])} if "census without" in name else census
        try:
            check_partner(t, done, cen)
            why = None
        except AssertionError as e:
            why = str(e) or "assertion"
        print(f"  mutation {name!r}: " + (f"refused ({why})" if why else "ACCEPTED"))
        assert why, name
    bad = []
    for name, t in zip(names, muts):
        why = refused(t, census)
        print(f"  mutation {name!r}: " + (f"refused ({why})" if why else "ACCEPTED"))
        if not why:
            bad.append(name)
    assert not bad, f"mutations accepted: {bad}"
    # file-level mutations: the parser and the checks together must refuse an edited Lean file
    import os
    import tempfile
    text = open(sys.argv[1]).read()
    fmuts = [("point count 16 -> 15 in cardp_4_3",
              "theorem cardp_4_3 : Fintype.card Ep_4_3.toAffine.Point = 16",
              "theorem cardp_4_3 : Fintype.card Ep_4_3.toAffine.Point = 15"),
             ("point (24, 840) -> (24, 841) in pts_4_3", ".some 24 840 (ns", ".some 24 841 (ns"),
             ("one tau case of (4, 3) deleted", "  · exact Or.inl (tau_zero (by norm_num) h)\n", ""),
             ("ratio of ratio_ne_12_5 changed", "(b + d) * 5 ≠ (h + f) * 12 :=", "(b + d) * 5 ≠ (h + f) * 13 :=")]
    nf = 0
    for name, a, b in fmuts:
        assert text.count(a) >= 1, name
        with tempfile.NamedTemporaryFile("w", suffix=".lean", delete=False) as fh:
            fh.write(text.replace(a, b, 1))
            tmp = fh.name
        try:
            ss, _, pp = parse(tmp)
            for x in ss:
                check_section(x, census)
            dd = {(x["m"], x["n"]) for x in ss}
            for t in pp:
                check_partner(t, dd, census)
            assert len(ss) == len(secs), "a section disappeared"
            why = None
        except (AssertionError, AttributeError, ValueError) as e:
            why = str(e) or type(e).__name__
        finally:
            os.unlink(tmp)
        print(f"  file mutation {name!r}: " + (f"refused ({why})" if why else "ACCEPTED"))
        assert why, name
        nf += 1
    print(f"negative test: {len(muts) + len(pmuts) + nf} of {len(muts) + len(pmuts) + nf} mutated inputs refused")
    print(f"ALL OK: {len(secs)} pairs, {len(partners)} partner ratios")


if __name__ == "__main__":
    main()
