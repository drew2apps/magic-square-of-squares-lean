#!/usr/bin/env python3
"""Witnesses for the positive controls of ClassHensel.lean (exact integer arithmetic).

For a prime p, a ratio m : n, a point (u, v) and a class d, find roots w_j of the four norms modulo
p^K with p^3 not dividing w_j, and Gaussian roots (a, b) with d*Q_f = (a + ib)^2 modulo p^K for the
eight products, then check every congruence exactly and print Lean source. Negative test: a mutated
witness must be refused by check().
Written by an AI model; not read by a person.
"""
import sys

def P(m, n, u, v):
    return [(m*v, n*u), (n*v, m*u), ((m+n)*(u-v), (m-n)*(u+v)), ((m-n)*(u-v), (m+n)*(u+v))]

def gm(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def cj(a): return (a[0], -a[1])
def cjIf(c, a): return cj(a) if c else a
FLAGS = [(a, b, c) for a in (False, True) for b in (False, True) for c in (False, True)]
def Q(m, n, u, v, f):
    p0, p8, pl, pu = P(m, n, u, v)
    return gm(gm(gm(p0, cjIf(f[0], p8)), cjIf(f[1], pl)), cjIf(f[2], pu))
def nrm(a): return a[0]**2 + a[1]**2
def val(p, x):
    if x == 0: return 10**9
    e = 0
    while x % p == 0: x //= p; e += 1
    return e

def int_root(x, M, p, K):
    """t with t^2 = x mod p^K (= M), by content, a root mod p (search) and Newton."""
    x %= M
    if x == 0:
        return 0
    e = val(p, x)
    assert e % 2 == 0 and e < K, ("odd valuation", e)
    A = x // p**e
    c = next((c for c in range(1, p) if (c*c - A) % p == 0), None)
    assert c is not None, "unit part not a square mod p"
    t = c
    for j in range(1, K):
        Pj = p**j
        s = (A - t*t) // Pj
        assert (A - t*t) % Pj == 0
        t = (t + Pj * (pow(2*t, -1, p) * s)) % (Pj*p)
    assert (t*t - A) % M == 0
    return (p**(e//2) * t) % M

def gauss_root(p, K, X, Y):
    """(a, b) with X = a^2 - b^2, Y = 2ab mod p^K."""
    M = p**K
    iota = int_root(-1, M, p, K) if pow(p - 1, (p - 1)//2, p) == 1 else None
    if iota is not None:   # split
        z1, z2 = (X + iota*Y) % M, (X - iota*Y) % M
        x1, x2 = int_root(z1, M, p, K), int_root(z2, M, p, K)
        h = pow(2, -1, M)
        a = h*(x1 + x2) % M
        b = (-iota*h*(x1 - x2)) % M
        return a, b
    # inert: content, unit part, root mod p, Newton
    c = min(val(p, X), val(p, Y))
    assert c % 2 == 0 and 2*(c//2) < K
    r = c // 2
    A, B = X // p**(2*r), Y // p**(2*r)
    base = None
    for c1 in range(p):
        for c2 in range(p):
            if (c1, c2) != (0, 0) and (c1*c1 - c2*c2 - A) % p == 0 and (2*c1*c2 - B) % p == 0:
                base = (c1, c2); break
        if base: break
    assert base is not None, "unit part not a square mod p"
    x, y = base
    for j in range(1, K):        # lift from p^j to p^(j+1)
        Pj = p**j
        t1 = (A - (x*x - y*y)) // Pj
        t2 = (B - 2*x*y) // Pj
        assert (A - (x*x - y*y)) % Pj == 0 and (B - 2*x*y) % Pj == 0
        alpha = pow(2*(x*x + y*y), -1, p)
        x, y = x + Pj*(alpha*(t1*x + t2*y)), y + Pj*(alpha*(t2*x - t1*y))
        x, y = x % (Pj*p), y % (Pj*p)
    return (p**r * x) % M, (p**r * y) % M

def witnesses(p, K, m, n, u, v, d):
    M = p**K
    ws = []
    for q in P(m, n, u, v):
        N = nrm(q)
        e = val(p, N)
        assert e % 2 == 0 and e <= 4, ("norm valuation", e)
        w = int_root(N, M, p, K)
        assert w % p**3 != 0
        ws.append(w)
    gs = [gauss_root(p, K, d*Q(m, n, u, v, f)[0], d*Q(m, n, u, v, f)[1]) for f in FLAGS]
    return ws, gs

def check(p, K, m, n, u, v, d, ws, gs):
    M = p**K
    ok = True
    for q, w in zip(P(m, n, u, v), ws):
        ok &= (nrm(q) - w*w) % M == 0 and w % p**3 != 0
    for f, (a, b) in zip(FLAGS, gs):
        X, Y = Q(m, n, u, v, f)
        ok &= (d*X - (a*a - b*b)) % M == 0 and (d*Y - 2*a*b) % M == 0
    ok &= d % (p*p) != 0
    return ok

def lean(name, p, K, m, n, u, v, d, ws, gs):
    b = lambda x: "false" if not x else "true"
    lines = [f"theorem {name} : ClassAt ((({p} : ℕ) : ℤ) ^ {K}) ((({p} : ℕ) : ℤ) ^ 3) {d} {m} {n} {u} {v} := by",
             f"  refine ⟨{ws[0]}, {ws[1]}, {ws[2]}, {ws[3]}, by decide +kernel, by decide +kernel, ?_⟩",
             "  rintro ⟨_ | _, _ | _, _ | _⟩"]
    for (a, bb) in gs:
        lines.append(f"  · exact ⟨{a}, {bb}, by decide +kernel, by decide +kernel⟩")
    return "\n".join(lines)

if __name__ == "__main__":
    K = 18
    cases = [("classAt_seven", 7, 49, 1, 7, 1, 7*3), ("classAt_seventeen", 17, 289, 1, 17, 1, 17*2)]
    if len(sys.argv) > 1 and sys.argv[1] == "search":
        for (name, p, m, n, u, v, _) in cases:
            for d1 in range(1, 3*p):
                if d1 % p == 0: continue
                for d in (d1, p*d1):
                    try:
                        ws, gs = witnesses(p, K, m, n, u, v, d)
                        print(name, "d =", d, "ok", check(p, K, m, n, u, v, d, ws, gs))
                    except AssertionError as e:
                        pass
        sys.exit()
    for (name, p, m, n, u, v, d) in cases:
        ws, gs = witnesses(p, K, m, n, u, v, d)
        assert check(p, K, m, n, u, v, d, ws, gs)
        # negative tests: a mutated Gaussian root and a norm root times p^3 are refused by the congruences;
        # d times p^2 is refused by the hypothesis p^2 not dividing d (the congruences alone would not refuse it)
        bad = list(gs); bad[3] = (bad[3][0] + 1, bad[3][1])
        assert not check(p, K, m, n, u, v, d, ws, bad)
        badw = list(ws); badw[1] = ws[1] * p**3
        assert not check(p, K, m, n, u, v, d, badw, gs)
        assert not check(p, K, m, n, u, v, d * p * p, ws, gs)
        print(lean(name, p, K, m, n, u, v, d, ws, gs))
        print()
    print("-- all checks passed")
