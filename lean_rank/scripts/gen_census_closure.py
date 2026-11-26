#!/usr/bin/env python3
"""Generate LeanRank/CensusClosure.lean: for census pairs (m, n) whose closing curve is w·f0·f8
(G1_FULL_SYSTEM.md, Result 3), the explicit list of all rational points of
E : Y^2 = X(X + m^4)(X + n^4) = curve (-m^4) (-n^4), and from it, in Lean,
  tau_<m>_<n>      : m^2 + n^2 t^2 and n^2 + m^2 t^2 both squares of rationals  =>  t in {0, 1, -1},
  ratio_ne_<m>_<n> : a 3x3 magic square of nine integer squares with b^2 != d^2, b^2 != f^2 has
                     (b + d) * n != (h + f) * m.

Written by an AI model; not read by a person.  Untrusted generator, pure Python (no PARI): it finds
the 2-power torsion points by repeated halving (x(Q) = x0 + r1 r2 + r1 r3 + r2 r3, r_i^2 = x0 - e_i)
and the smallest odd prime p of good reduction with #E~(F_p) equal to their number, by brute
force.  (For this family E(Q)_tors contains Z/2 x Z/4, so by Mazur it is Z/2 x Z/4 or Z/2 x Z/8 and
halving finds all of it; Lean does not rely on this.)  Lean re-checks everything (the point count
over ZMod p by `decide +kernel`, every point on the curve by `norm_num`, their distinctness by
`decide +kernel`, every non-square by a bracket k^2 < N < (k+1)^2).

Rank 0 is taken from LeanRank/RankZeroCensus.lean (run 3) for the pairs it contains; for other
pairs the rank-0 certificate is produced by the functions of the run-3 generator
gen_rank_zero_certs.py (imported, not edited), which refuses a pair whose local test leaves
more than 4 classes.

usage: python3 -I gen_census_closure.py OUT.lean m,n [m,n ...]
"""
import importlib.util
import itertools
import os
import sys
from fractions import Fraction
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("gen_rank_zero_certs",
                                               os.path.join(HERE, "gen_rank_zero_certs.py"))
G3 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(G3)

# pairs whose rank-0 theorem is in LeanRank/RankZeroCensus.lean (namespace LeanRank.RankZeroCensus)
RUN3 = {(4, 3), (7, 1), (12, 5), (24, 7)}


def is_sq(q):
    q = Fraction(q)
    if q < 0:
        return False
    a, b = q.numerator, q.denominator
    return isqrt(a) ** 2 == a and isqrt(b) ** 2 == b


def sqrt_q(q):
    q = Fraction(q)
    if q < 0:
        return None
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    return Fraction(a, b) if a * a == q.numerator and b * b == q.denominator else None


def ell_add(P, Q, A, B):
    """Chord and tangent on y^2 = x^3 + A x^2 + B x (None is O), exact."""
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


def torsion_points(m, n):
    """The 2-power torsion of Y^2 = X(X + m^4)(X + n^4), by repeated halving; affine points as integers."""
    a, b = m ** 4, n ** 4
    A, B = a + b, a * b
    S = {None, (Fraction(0), Fraction(0)), (Fraction(-a), Fraction(0)), (Fraction(-b), Fraction(0))}
    frontier = [P for P in S if P is not None]
    while frontier:
        new = []
        for (x0, y0) in frontier:
            rs = [sqrt_q(x0 - e) for e in (0, -a, -b)]
            if any(r is None for r in rs):
                continue
            for sg in itertools.product((1, -1), repeat=3):
                r1, r2, r3 = (si * ri for si, ri in zip(sg, rs))
                x = x0 + r1 * r2 + r1 * r3 + r2 * r3
                y = sqrt_q(x * (x + a) * (x + b))
                if y is None:
                    continue
                for Q in ((x, y), (x, -y)):
                    if Q not in S and ell_add(Q, Q, A, B) == (x0, y0):
                        S.add(Q)
                        new.append(Q)
        frontier = new
    for P in S:
        for Q in S:
            assert ell_add(P, Q, A, B) in S                     # a group
    out = []
    for P in S:
        if P is None:
            continue
        x, y = P
        assert x.denominator == 1 and y.denominator == 1        # Nagell-Lutz; Lean does not need it
        out.append((int(x), int(y)))
    return len(S), out


def good_prime(m, n, N):
    """Smallest odd prime p not dividing 16 (ab(a - b))^2 with #E~(F_p) = N (brute force)."""
    a, b = m ** 4, n ** 4
    disc = 16 * (a * b * (a - b)) ** 2
    for p in range(3, 400):
        if any(p % d == 0 for d in range(2, isqrt(p) + 1)) or disc % p == 0:
            continue
        if count_mod_p((a + b) % p, (a * b) % p, p) == N:
            return p
    return None


def count_mod_p(A2, A4, p):
    """#E(F_p) for y^2 = x^3 + A2 x^2 + A4 x, by brute force (independent of PARI)."""
    sq = [0] * p
    for y in range(p):
        sq[y * y % p] += 1
    return 1 + sum(sq[(x ** 3 + A2 * x * x + A4 * x) % p] for x in range(p))


def tau_case(m, n, x):
    """How the abscissa x excludes t with m^2 n^2 t^2 = x: returns a Lean term."""
    c = m * m * n * n
    if x < 0:
        return "neg", None
    if x == 0:
        return "zero", None
    if x == c:
        return "one", None
    w = Fraction(x, c)
    assert not is_sq(w), (m, n, x)
    P, Q = w.numerator, w.denominator
    N = P * Q
    k = isqrt(N)
    assert k * k < N < (k + 1) * (k + 1)
    assert Fraction(x) * Q * Q == c * N
    return "nonsq", (N, k, Q)


def lq(v):
    return f"({v})" if v < 0 else f"{v}"


def section(m, n):
    tag = f"{m}_{n}"
    a, b = m ** 4, n ** 4
    c = m * m * n * n
    N, pts = torsion_points(m, n)
    p = good_prime(m, n, N)
    assert p is not None, (m, n)
    A2, A4 = (a + b) % p, (a * b) % p
    assert count_mod_p(A2, A4, p) == N
    for (x, y) in pts:
        assert y * y == x * (x + a) * (x + b)
    pts = sorted(pts, key=lambda t: (t[1] != 0, abs(t[0]), t[0], -t[1]))
    xs = []
    for (x, _) in pts:
        if x not in xs:
            xs.append(x)
    finrank = (f"LeanRank.RankZeroCensus.finrank_{tag}" if (m, n) in RUN3 else f"finrank_{tag}")
    L = []
    L.append(f"/-! ### `({m}, {n})`: `Y² = X(X + {a})(X + {b})`, {N} points, reduction modulo {p} -/\n")
    L.append(f"instance : (curve ({-a}) ({-b})).IsElliptic :=\n"
             f"  LeanRank.RankZeroGen.isElliptic (by decide) (by decide) (by decide)\n")
    L.append(f"instance : Fact (Nat.Prime {p}) := ⟨by norm_num⟩\n")
    L.append(f"instance : (redCurve (intPrime {p}) (W₀ {a} {b})).IsElliptic :=\n"
             f"  redCurve_isElliptic {p} {a} {b} (by norm_num)\n")
    L.append(f"/-- `Y² = X³ + {a + b}X² + {a * b}X` modulo `{p}`. -/\n"
             f"def Ep_{tag} : WeierstrassCurve (ZMod {p}) := ⟨0, {A2}, 0, {A4}, 0⟩\n")
    L.append(f"theorem base_{tag} :\n"
             f"    (((W₀ {a} {b}).toAffine ⁄ (ZMod {p})) : WeierstrassCurve _).toAffine =\n"
             f"      Ep_{tag}.toAffine := by\n  ext <;> decide +kernel\n")
    L.append(f"/-- The reduction modulo `{p}` has `{N}` points. -/\n"
             f"theorem cardp_{tag} : Fintype.card Ep_{tag}.toAffine.Point = {N} := by decide +kernel\n")
    L.append(f"theorem fin_card_le_{tag} :\n"
             f"    Finite (curve ({-a}) ({-b})).Point ∧ Nat.card (curve ({-a}) ({-b})).Point ≤ {N} :=\n"
             f"  card_le {p} {a} {b} (by norm_num) {finrank} Ep_{tag} base_{tag} cardp_{tag}\n")
    L.append(f"instance : Finite (curve ({-a}) ({-b})).Point := fin_card_le_{tag}.1\n")
    L.append(f"theorem card_le_{tag} : Nat.card (curve ({-a}) ({-b})).Point ≤ {N} := fin_card_le_{tag}.2\n")
    entries = ["0"] + [f".some {lq(x)} {lq(y)} (ns (by norm_num))" for (x, y) in pts]
    body = ",\n    ".join(entries)
    L.append(f"/-- The {N} points of `Y² = X(X + {a})(X + {b})`. -/\n"
             f"noncomputable def pts_{tag} : Fin {N} → (curve ({-a}) ({-b})).Point :=\n"
             f"  ![{body}]\n")
    L.append(f"theorem pts_{tag}_injective : Function.Injective pts_{tag} := by\n"
             f"  have : Function.Injective (key ∘ pts_{tag}) := by decide +kernel\n"
             f"  exact this.of_comp\n")
    L.append(f"/-- **`E(ℚ)` is exactly the {N} points `pts_{tag}`**, for "
             f"`E : Y² = X(X + {a})(X + {b})`. -/\n"
             f"theorem points_{tag} (P : (curve ({-a}) ({-b})).Point) : ∃ i, pts_{tag} i = P :=\n"
             f"  (pts_{tag}_injective.bijective_of_nat_card_le (by rw [Nat.card_fin]; exact card_le_{tag})).2 P\n")
    L.append(f"theorem card_{tag} : Nat.card (curve ({-a}) ({-b})).Point = {N} :=\n"
             f"  le_antisymm card_le_{tag} (by\n"
             f"    have := Nat.card_le_card_of_injective _ pts_{tag}_injective\n"
             f"    rwa [Nat.card_fin] at this)\n")
    plist = ", ".join(f"({lq(x)}, {lq(y)})" for (x, y) in pts)
    L.append(f"/-- **The affine rational points of `Y² = X(X + {a})(X + {b})`.** -/\n"
             f"theorem affine_{tag} {{x y : ℚ}} (h : y ^ 2 = x * (x + {a}) * (x + {b})) :\n"
             f"    (x, y) ∈ ([{plist}] : List (ℚ × ℚ)) := by\n"
             f"  obtain ⟨i, hi⟩ := points_{tag} (.some x y (ns (by push_cast; linear_combination h)))\n"
             f"  have hk := congrArg key hi\n"
             f"  fin_cases i <;> simp [pts_{tag}, key] at hk <;> obtain ⟨rfl, rfl⟩ := hk <;> simp\n")
    cases = []
    for (x, y) in pts:
        kind, data = tau_case(m, n, x)
        if kind == "zero":
            cases.append("exact Or.inl (tau_zero (by norm_num) h)")
        elif kind == "one":
            cases.append("exact Or.inr (tau_one (by norm_num) h)")
        elif kind == "neg":
            cases.append("exact (tau_neg (by norm_num) h (by norm_num)).elim")
        else:
            Nn, k, Q = data
            cases.append(f"exact (tau_nonsq (by norm_num) h (N := {Nn}) (k := {k}) (Q := {Q})"
                         f" (by norm_num) (by norm_num) (by norm_num)).elim")
    pat = " | ".join(["⟨h, -⟩"] * len(pts))
    L.append(f"/-- **`{m * m} + {n * n}τ²` and `{n * n} + {m * m}τ²` are both squares of rationals only for "
             f"`τ ∈ {{0, 1, −1}}`** (curve `w·f0·f8`, `w = τ²`, of the census pair `({m}, {n})`). -/\n"
             f"theorem tau_{tag} {{τ r0 r8 : ℚ}} (h0 : r0 ^ 2 = {m * m} + {n * n} * τ ^ 2)\n"
             f"    (h8 : r8 ^ 2 = {n * n} + {m * m} * τ ^ 2) : τ = 0 ∨ τ = 1 ∨ τ = -1 := by\n"
             f"  have hE : ({c} * τ * r0 * r8) ^ 2 =\n"
             f"      ({c} * τ ^ 2) * (({c} * τ ^ 2) + {a}) * (({c} * τ ^ 2) + {b}) := by\n"
             f"    have e : ({c} * τ * r0 * r8) ^ 2 = {c} ^ 2 * τ ^ 2 * r0 ^ 2 * r8 ^ 2 := by ring\n"
             f"    rw [e, h0, h8]; ring\n"
             f"  have hx := affine_{tag} hE\n"
             f"  simp only [List.mem_cons, List.not_mem_nil, or_false, Prod.mk.injEq] at hx\n"
             f"  rcases hx with {pat}\n" +
             "".join(f"  · {s}\n" for s in cases))
    L.append(f"/-- **No magic square of squares at the ratio `{m} : {n}`.** In a 3×3 magic square of the\n"
             f"squares of nine integers `a b c / d e f / g h i` (any signs, magic sum `T`) with\n"
             f"`b² ≠ d²` and `b² ≠ f²`, the edge roots do not satisfy `(b + d) : (h + f) = {m} : {n}`. -/\n"
             f"theorem ratio_ne_{tag} {{a b c d e f g h i T : ℤ}}\n"
             f"    (r₁ : a ^ 2 + b ^ 2 + c ^ 2 = T) (r₂ : d ^ 2 + e ^ 2 + f ^ 2 = T)\n"
             f"    (r₃ : g ^ 2 + h ^ 2 + i ^ 2 = T) (c₁ : a ^ 2 + d ^ 2 + g ^ 2 = T)\n"
             f"    (c₂ : b ^ 2 + e ^ 2 + h ^ 2 = T) (c₃ : c ^ 2 + f ^ 2 + i ^ 2 = T)\n"
             f"    (dg : a ^ 2 + e ^ 2 + i ^ 2 = T) (ad : c ^ 2 + e ^ 2 + g ^ 2 = T)\n"
             f"    (hbd : b ^ 2 ≠ d ^ 2) (hbf : b ^ 2 ≠ f ^ 2) : (b + d) * {n} ≠ (h + f) * {m} := by\n"
             f"  have := ratio_ne_of_tau (m := {m}) (n := {n}) (by norm_num) (by norm_num)\n"
             f"    (fun τ r0 r8 h0 h8 => tau_{tag} (τ := τ) (r0 := r0) (r8 := r8)\n"
             f"      (by push_cast at h0; linear_combination h0)\n"
             f"      (by push_cast at h8; linear_combination h8))\n"
             f"    r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf\n"
             f"  exact_mod_cast this\n")
    return "\n".join(L) + "\n", dict(m=m, n=n, N=N, p=p, pts=pts)


def partner(m, n):
    """The partner ratio (m + n) : (m - n), in lowest terms (G1_FULL_SYSTEM.md, "The system")."""
    from math import gcd
    M, N = m + n, m - n
    g = gcd(M, N)
    return M // g, N // g


def partner_section(m, n):
    """ratio_ne for the partner P of (m, n), from ratio_ne_m_n through ratio_ne_partner."""
    pm, pn = partner(m, n)
    assert (pm + pn) * n == (pm - pn) * m and pm > pn > 0
    tag, ptag = f"{m}_{n}", f"{pm}_{pn}"
    return (f"/-- **No magic square of squares at the ratio `{pm} : {pn}`**, the partner of `{m} : {n}`\n"
            f"(`ratio_ne_partner`); hypotheses as in `ratio_ne_{tag}`. -/\n"
            f"theorem ratio_ne_{ptag} {{a b c d e f g h i T : ℤ}}\n"
            f"    (r₁ : a ^ 2 + b ^ 2 + c ^ 2 = T) (r₂ : d ^ 2 + e ^ 2 + f ^ 2 = T)\n"
            f"    (r₃ : g ^ 2 + h ^ 2 + i ^ 2 = T) (c₁ : a ^ 2 + d ^ 2 + g ^ 2 = T)\n"
            f"    (c₂ : b ^ 2 + e ^ 2 + h ^ 2 = T) (c₃ : c ^ 2 + f ^ 2 + i ^ 2 = T)\n"
            f"    (dg : a ^ 2 + e ^ 2 + i ^ 2 = T) (ad : c ^ 2 + e ^ 2 + g ^ 2 = T)\n"
            f"    (hbd : b ^ 2 ≠ d ^ 2) (hbf : b ^ 2 ≠ f ^ 2) : (b + d) * {pn} ≠ (h + f) * {pm} := by\n"
            f"  have := ratio_ne_partner (m := {pm}) (n := {pn}) (M := {m}) (N := {n}) (by norm_num) (by norm_num)\n"
            f"    (fun a b c d e f g h i T r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf => by\n"
            f"      have := ratio_ne_{tag} r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf\n"
            f"      exact_mod_cast this)\n"
            f"    r₁ r₂ r₃ c₁ c₂ c₃ dg ad hbd hbf\n"
            f"  exact_mod_cast this\n\n"), (pm, pn)


def rank_part(pairs):
    """Rank-0 certificates for the pairs not in RankZeroCensus.lean, by the run-3 generator's own
    functions (cert, emit), moved into this file's namespace."""
    certs = []
    for (m, n) in pairs:
        c = G3.cert(m, n)
        if c is None:
            return None, (m, n)
        certs.append(c)
    if not certs:
        return "", None
    txt = G3.emit(certs)
    start = txt.index("namespace LeanRank.RankZeroCensus\n") + len("namespace LeanRank.RankZeroCensus\n")
    end = txt.index("end LeanRank.RankZeroCensus\n")
    body = txt[start:end].strip("\n")
    return ("/-! ### Rank-0 certificates (functions of `gen_rank_zero_certs.py`, run 3) -/\n\n"
            + body + "\n\n"), None


def main():
    path = sys.argv[1]
    head = open(os.path.join(HERE, "census_closure_head.lean.in")).read()
    pairs = [tuple(map(int, s.split(","))) for s in sys.argv[2:]]
    assert len(set(pairs)) == len(pairs), "duplicate pair"
    new = [pq for pq in pairs if pq not in RUN3]
    rk, bad = rank_part(new)
    if rk is None:
        print(f"{bad}: REFUSED (rank-0 local test leaves more than 4 classes)", file=sys.stderr)
        return 1
    out = [head, rk]
    infos = []
    for (m, n) in pairs:
        s, info = section(m, n)
        out.append(s)
        infos.append(info)
        print(f"({m},{n}): {info['N']} points, p = {info['p']}", file=sys.stderr)
    partners = []
    for (m, n) in pairs:
        P = partner(m, n)
        if P not in pairs and P not in [q for q, _ in partners]:
            partners.append((P, (m, n)))
    if partners:
        out.append("/-! ### Partner ratios -/\n\n")
        for (P, (m, n)) in partners:
            txt, _ = partner_section(m, n)
            out.append(txt)
            print(f"({P[0]},{P[1]}): partner of ({m},{n})", file=sys.stderr)
    out.append("end LeanRank.CensusClosure\n\n")
    for th in ["mss_forms", "mss_point", "ratio_ne_of_tau", "ratio_ne_partner"]:
        out.append(f"#print axioms LeanRank.CensusClosure.{th}\n")
    for info in infos:
        tag = f"{info['m']}_{info['n']}"
        ths = (["finrank_"] if (info['m'], info['n']) not in RUN3 else []) + \
            ["points_", "card_", "affine_", "tau_", "ratio_ne_"]
        for th in ths:
            out.append(f"#print axioms LeanRank.CensusClosure.{th}{tag}\n")
    for (P, _) in partners:
        out.append(f"#print axioms LeanRank.CensusClosure.ratio_ne_{P[0]}_{P[1]}\n")
    open(path, "w").write("".join(out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
