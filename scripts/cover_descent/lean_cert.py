#!/usr/bin/env python3
"""Writes MagicSquareOfSquares/CoverDescentTable.lean: certificates for the descent of CoverDescent.lean.

For each pair (m, n) of G1.exceptions1000 (read from G1Census1000.lean) with m^2 + n^2 neither a square nor twice a
square, the candidates are the products d of 2 and of the primes = +-1 mod 8 dividing m n (m^2 - n^2) (Theorem E:
mod_eight_of_dv and exists_class_prod8 in CoverDescent.lean; `keep8` there).  For a modulus M = p^k from a fixed list
(the prime powers up to 128), the program finds the candidates that are excluded modulo M by the check `kills` of
CoverDescent.lean, which it mirrors: at every point (1, v), and (u, 1) with p | u, of the projective line modulo M, the
four forms are not all squares modulo M, or d times one of the eight products of the four linear forms (with
conjugates) is not the square of a Gaussian integer modulo M.  A ratio gets an entry when every candidate is excluded by
some modulus.  The Lean kernel then checks every entry again; this program is not trusted.
usage (from the root of the repository): python3 scripts/cover_descent/lean_cert.py [--check]"""
import sys, re, itertools
from math import isqrt
FLAGS = list(itertools.product((0, 1), repeat=3))
def gmul(p, q): return (p[0]*q[0] - p[1]*q[1], p[0]*q[1] + p[1]*q[0])
def cj(p): return (p[0], -p[1])
def forms(m, n, u, v): return [(m*v, n*u), (n*v, m*u), ((m+n)*(u-v), (m-n)*(u+v)), ((m-n)*(u-v), (m+n)*(u+v))]
def Q(P, f):
    r = gmul(P[0], cj(P[1]) if f[0] else P[1]); r = gmul(r, cj(P[2]) if f[1] else P[2]); return gmul(r, cj(P[3]) if f[2] else P[3])
MASK = {}
def masks(M):
    if M not in MASK:
        MASK[M] = (set(w*w % M for w in range(M)), set(((a*a - b*b) % M, (2*a*b) % M) for a in range(M) for b in range(M)))
    return MASK[M]
def survivors(m, n, M, p):
    sq, g = masks(M); out = []
    for u, v in [(1, v) for v in range(M)] + [(u, 1) for u in range(0, M, p)]:
        P = forms(m, n, u, v)
        if all((z[0]*z[0] + z[1]*z[1]) % M in sq for z in P): out.append([Q(P, f) for f in FLAGS])
    return out
def killed(d, surv, M):
    g = masks(M)[1]
    return not any(all(((d*X) % M, (d*Y) % M) in g for X, Y in pr) for pr in surv)
def factor(x):
    out = set(); q = 2
    while q*q <= x:
        while x % q == 0: out.add(q); x //= q
        q += 1
    if x > 1: out.add(x)
    return out
CAP = 128
PR = [q for q in range(5, CAP) if all(q % r for r in range(2, int(q**0.5) + 1))]
MODS = [(2, 3), (3, 2)] + sorted([(q, 1) for q in PR] + [(2, 4), (5, 2), (3, 3), (2, 5), (7, 2), (2, 6), (3, 4), (11, 2), (5, 3), (2, 7)], key=lambda t: t[0]**t[1])
def one(mn):
    m, n = mn
    P = sorted(factor(2*m*n*(m*m - n*n))); P8 = [q for q in P if q == 2 or q % 8 in (1, 7)]; cands = [1]
    for q in P8: cands = cands + [c*q for c in cands]
    left = set(cands); cert = []
    for p, k in MODS:
        M = p**k; surv = survivors(m, n, M, p)
        kl = sorted(d for d in left if killed(d, surv, M))
        if kl: cert.append((p, k, kl)); left -= set(kl)
        if not left: break
    for p in EPR:                                           # lines by Euler's criterion (k = 0 in the certificate)
        if not left: break
        io = iota(p); kl = sorted(d for d in left if killed_e(m, n, p, io, d))
        if kl: cert.append((p, 0, kl)); left -= set(kl)
    return dict(m=m, n=n, primes=P, primes8=P8, cands=len(cands), cert=cert, left=sorted(left))
ECAP = 400
EPR = [q for q in range(CAP, ECAP) if q % 4 == 1 and all(q % r for r in range(2, int(q**0.5) + 1))]
def eul(x, p): return pow(x % p, (p - 1)//2, p) <= 1      # `eul` in CoverDescent.lean
def iota(p): return next(j for j in range(p) if j*j % p == p - 1)    # `iotaOf`
def killed_e(m, n, p, io, d):                               # `killsE` for one class
    for u, v in [(1, v) for v in range(p)] + [(0, 1)]:
        P = forms(m, n, u, v)
        if not all(eul(z[0]*z[0] + z[1]*z[1], p) for z in P): continue
        if all(eul(d*(X + io*Y), p) and eul(d*(X + (p - io)*Y), p) for X, Y in (Q(P, f) for f in FLAGS)): return False
    return True
def four(mn):                                               # four classes (`entryOk4`): the lines (p, k) that remove tuples
    m, n = mn; P = sorted(factor(2*m*n*(m*m - n*n))); keep = lambda q: q == 2 or q % 8 in (1, 7)
    A = prods([q for q in P if keep(q) and (2*m*n) % q == 0]); B = prods([q for q in P if keep(q) and (2*(m + n)*(m - n)) % q == 0])
    C = [A, A, B, B]; left = set(itertools.product(*C)); lines = []
    for p, k in MODS:
        M = p**k; sq, g = masks(M); allowed = set()
        for u in range(M):
            for v in range(M):
                if u % p == 0 and v % p == 0: continue
                Pf = forms(m, n, u, v)
                if not all((z[0]*z[0] + z[1]*z[1]) % M in sq for z in Pf): continue
                S = [[d for d in C[j] if ((d*Pf[j][0]) % M, (d*Pf[j][1]) % M) in g] for j in range(4)]
                if all(S): allowed |= set(itertools.product(*S))
        if len(left & allowed) < len(left): lines.append((p, k)); left &= allowed
        if not left: break
    return dict(m=m, n=n, primes=P, lines=lines, left=len(left), tuples=len(A)**2*len(B)**2)
SBOUND = 100000
SPR = [q for q in range(3, SBOUND, 2) if all(q % r for r in range(3, int(q**0.5) + 1, 2))]
def struct(mn):                                             # `structOk`: primes of (m^2+n^2)(m^8+14m^4n^4+n^8) that leave one candidate
    m, n = mn; P = sorted(factor(2*m*n*(m*m - n*n))); C = [d for d in prods([q for q in P if q == 2 or q % 8 in (1, 7)]) if d % 2 == 0]
    N = (m*m + n*n)*(m**8 + 14*m**4*n**4 + n**8); qs = []; left = set(C)
    for q in SPR:
        if len(left) <= 1: break
        if N % q == 0:
            new = set(d for d in left if eul(d, q))
            if len(new) < len(left): qs.append(q); left = new
    return dict(m=m, n=n, primes=P, d0=min(left), qs=qs, left=len(left))
def smask(M):
    S = 0
    for r in range(M): S |= 1 << (r*r % M)
    return S
def gmask(M):
    G = 0
    for a in range(M):
        for b in range(M): G |= 1 << (M*((a*a - b*b) % M) + (2*a*b) % M)
    return G
def prods(ps):                                              # the order of `prods` in CoverDescent.lean
    if not ps: return [1]
    r = prods(ps[1:]); return r + [ps[0]*x for x in r]
def cmask(ps, ds):
    pos = {d: i for i, d in enumerate(prods(ps))}; c = 0
    for d in ds: c |= 1 << pos[d]
    return c
HEAD = 'import MagicSquareOfSquares.CoverDescent\nimport MagicSquareOfSquares.G1Census1000\n\n/-!\n# Certificates for the descent of `CoverDescent.lean`\n\nGenerated by `scripts/cover_descent/lean_cert.py`; do not edit by hand.\n\n* `masks`: for each modulus used, the squares and the Gaussian squares modulo it, as bit masks;\n  `mask8`, …: each mask is the one computed by `sMask`, `gMask` (evaluated by the kernel).\n* `table`: for each ratio `m : n`, the primes of `2mn(m − n)(m + n)` and a certificate: lines\n  `(p, k, c)` saying that modulo `p^k` the candidates at the bits of `c` are excluded (for\n  `k = 0`: modulo the prime `p`, by the criterion of Euler, `killsE`). The\n  candidates are the products of those of the primes that are 2 or `≡ ±1 (mod 8)` (Theorem E,\n  `exists_class_prod8`). `table_ok`: every entry passes `entryOk8` (evaluated by the kernel, in\n  pieces `chunk0`, …).\n* `table4`: certificates with one class for each form (`entryOk4`), for the ratios that have no\n  certificate with one class.\n* `structTable`, `class_of_mem_structTable`: for 150 ratios with a degenerate solution, among them\n  `627 : 364` and `991 : 263` (class `521266`), the class of every solution, from the theorems on\n  the class alone (`structOk`).\n* `no_fullSquares_of_mem_ratios`: for the ratios of the two tables, no rational `t` makes the\n  five free cells of the G1 family squares.\n* `IsMSS.ratio_ne_of_coverDescent`, `IsMSS.ratio_ne_of_le_1000_descent`: the same for a magic\n  square of squares; with the census by residue classes, for coprime `1 ≤ n < m ≤ 1000` the ratio\n  is one of the pairs of `exceptions`, and `m² + n²` is a square or twice a square\n  (`exceptions_degenerate`, `IsMSS.ratio_norm_sq_or_two_sq_of_le_1000`).\n\nNo `native_decide`: every check is `decide +kernel`.\n-/\n\nset_option Elab.async false\n\nnamespace MagicSquareOfSquares.CoverDescent\n\n'
TAIL = "\n/-- **No solution of the five-form system, from a certificate.** -/\ntheorem no_fullSquares_of_entry {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks)\n    {e : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)} (he : entryOk8 masks e = true) (t : ℚ) :\n    ¬G1.FullSquares e.1 e.2.1 t := by\n  intro h\n  obtain ⟨u, V, hV, hcop, h0, _, h8, hL, hU⟩ := G1.exists_coprime_squares_full h\n  have hv : (V : ℤ) ≠ 0 := by exact_mod_cast hV.ne'\n  have huv : (u : ℤ) + V ≠ 0 := by\n    have : (0 : ℤ) < (u : ℤ) + V := by positivity\n    exact this.ne'\n  refine no_point_of_entry8 hmasks he (Nat.isCoprime_iff_coprime.mpr hcop) hv huv ?_ ?_ ?_ ?_\n  · obtain ⟨r, hr⟩ := h0\n    exact ⟨r, by dsimp only [nrm, P0]; linear_combination hr⟩\n  · obtain ⟨r, hr⟩ := h8\n    exact ⟨r, by dsimp only [nrm, P8]; linear_combination hr⟩\n  · obtain ⟨r, hr⟩ := hL\n    exact ⟨r, by dsimp only [nrm, PL]; linear_combination hr⟩\n  · obtain ⟨r, hr⟩ := hU\n    exact ⟨r, by dsimp only [nrm, PU]; linear_combination hr⟩\n\n/-- **No solution of the five-form system, from a certificate with four classes.** -/\ntheorem no_fullSquares_of_entry4 {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks)\n    {e : ℕ × ℕ × List ℕ × List (ℕ × ℕ)} (he : entryOk4 masks e = true) (t : ℚ) :\n    ¬G1.FullSquares e.1 e.2.1 t := by\n  intro h\n  obtain ⟨u, V, hV, hcop, h0, _, h8, hL, hU⟩ := G1.exists_coprime_squares_full h\n  have hv : (V : ℤ) ≠ 0 := by exact_mod_cast hV.ne'\n  have huv : (u : ℤ) + V ≠ 0 := by\n    have : (0 : ℤ) < (u : ℤ) + V := by positivity\n    exact this.ne'\n  refine no_point_of_entry4 hmasks he (Nat.isCoprime_iff_coprime.mpr hcop) hv huv ?_ ?_ ?_ ?_\n  · obtain ⟨r, hr⟩ := h0\n    exact ⟨r, by dsimp only [nrm, P0]; linear_combination hr⟩\n  · obtain ⟨r, hr⟩ := h8\n    exact ⟨r, by dsimp only [nrm, P8]; linear_combination hr⟩\n  · obtain ⟨r, hr⟩ := hL\n    exact ⟨r, by dsimp only [nrm, PL]; linear_combination hr⟩\n  · obtain ⟨r, hr⟩ := hU\n    exact ⟨r, by dsimp only [nrm, PU]; linear_combination hr⟩\n\n/-- **For the ratios of the two tables, no rational `t` makes the five free cells of the G1 family\nsquares.** Among them are `820 : 799` and `889 : 320`. -/\ntheorem no_fullSquares_of_mem_ratios {m n : ℕ} (h : (m, n) ∈ ratios) (t : ℚ) :\n    ¬G1.FullSquares m n t := by\n  rw [← ratios_eq, List.mem_append, List.mem_map, List.mem_map] at h\n  rcases h with ⟨e, he, heq⟩ | ⟨e, he, heq⟩\n  · have h1 : entryOk8 masks e = true := (List.all_eq_true.mp table_ok) e he\n    have h2 := no_fullSquares_of_entry masks_ok h1 t\n    rw [Prod.mk.injEq] at heq\n    rw [← heq.1, ← heq.2]\n    exact h2\n  · have h1 : entryOk4 masks e = true := (List.all_eq_true.mp table4_ok) e he\n    have h2 := no_fullSquares_of_entry4 masks_ok h1 t\n    rw [Prod.mk.injEq] at heq\n    rw [← heq.1, ← heq.2]\n    exact h2\n\n/-- The pairs of `G1.exceptions1000` that are not in the table. -/\ndef exceptions : List (ℕ × ℕ) := G1.exceptions1000.filter fun x => decide (x ∉ ratios)\n\ntheorem exceptions_length : exceptions.length = EXCEPTIONS_LEFT := by decide +kernel\n\n/-- **Every pair that is left has `m² + n²` a square or twice a square.** -/\ntheorem exceptions_degenerate : ∀ x ∈ exceptions,\n    Nat.sqrt (x.1 ^ 2 + x.2 ^ 2) ^ 2 = x.1 ^ 2 + x.2 ^ 2 ∨\n      2 * Nat.sqrt ((x.1 ^ 2 + x.2 ^ 2) / 2) ^ 2 = x.1 ^ 2 + x.2 ^ 2 := by decide +kernel\n\nend MagicSquareOfSquares.CoverDescent\n\nnamespace MagicSquareOfSquares\n\n/-- **In a 3×3 magic square of squares the ratio `(b + d) : (h + f)` of the edge roots is not one\nof the ratios of the table of `CoverDescent`.** -/\ntheorem IsMSS.ratio_ne_of_coverDescent {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)\n    {m n : ℕ} (hn0 : 0 < n) (hmem : (m, n) ∈ CoverDescent.ratios) (hne : h + f ≠ 0) (hbd : b ≠ d) :\n    (b + d) * n ≠ (h + f) * m :=\n  fun hmn => CoverDescent.no_fullSquares_of_mem_ratios hmem _ (hs.fullSquares hn0 hmn hne hbd)\n\n/-- **In a 3×3 magic square of squares the ratio `(b + d) : (h + f)` of the edge roots is not\n`m : n` for any coprime `1 ≤ n < m ≤ 1000` outside the pairs of `CoverDescent.exceptions`.** -/\ntheorem IsMSS.ratio_ne_of_le_1000_descent {a b c d e f g h i T : ℤ}\n    (hs : IsMSS a b c d e f g h i T) {m n : ℕ} (hn0 : 0 < n) (hnm : n < m) (hm : m ≤ 1000)\n    (hcop : Nat.Coprime m n) (hex : (m, n) ∉ CoverDescent.exceptions) (hne : h + f ≠ 0)\n    (hbd : b ≠ d) : (b + d) * n ≠ (h + f) * m := by\n  by_cases hr : (m, n) ∈ CoverDescent.ratios\n  · exact hs.ratio_ne_of_coverDescent hn0 hr hne hbd\n  · refine hs.ratio_ne_of_le_1000 hn0 hnm hm hcop ?_ hne hbd\n    intro hmem\n    apply hex\n    unfold CoverDescent.exceptions\n    rw [List.mem_filter]\n    exact ⟨hmem, decide_eq_true hr⟩\n\n/-- **In a 3×3 magic square of squares whose ratio `(b + d) : (h + f)` is `m : n` with coprime\n`1 ≤ n < m ≤ 1000`, the number `m² + n²` is a square or twice a square.** -/\ntheorem IsMSS.ratio_norm_sq_or_two_sq_of_le_1000 {a b c d e f g h i T : ℤ}\n    (hs : IsMSS a b c d e f g h i T) {m n : ℕ} (hn0 : 0 < n) (hnm : n < m) (hm : m ≤ 1000)\n    (hcop : Nat.Coprime m n) (hne : h + f ≠ 0) (hbd : b ≠ d) (hr : (b + d) * n = (h + f) * m) :\n    (∃ k : ℕ, m ^ 2 + n ^ 2 = k ^ 2) ∨ ∃ k : ℕ, m ^ 2 + n ^ 2 = 2 * k ^ 2 := by\n  by_contra hcon\n  refine hs.ratio_ne_of_le_1000_descent hn0 hnm hm hcop ?_ hne hbd hr\n  intro hmem\n  have h1 := CoverDescent.exceptions_degenerate (m, n) hmem\n  dsimp only at h1\n  rcases h1 with h1 | h1\n  · exact hcon (Or.inl ⟨_, h1.symm⟩)\n  · exact hcon (Or.inr ⟨_, h1.symm⟩)\n\nend MagicSquareOfSquares\n"
if __name__ == '__main__':
    src = open('MagicSquareOfSquares/G1Census1000.lean').read()
    body = src[src.index('def exceptions1000'):]; body = body[:body.index(']')]
    pairs = [tuple(int(x) for x in t) for t in re.findall(r'\((\d+), (\d+)\)', body)]
    def triv(m, n):
        S = m*m + n*n; return isqrt(S)**2 == S or (S % 2 == 0 and isqrt(S//2)**2 == S//2)
    res = [one(p) for p in pairs if not triv(*p)]
    ok = [r for r in res if not r['left']]
    res4 = [four((r['m'], r['n'])) for r in res if r['left']]; ok4 = [r for r in res4 if not r['left']]
    mods = sorted(set(p**k for r in ok for p, k, _ in r['cert'] if k) | set(p**k for r in ok4 for p, k in r['lines']))
    out = []
    for M in mods:
        out += ['def S%d : ℕ := 0x%x' % (M, smask(M)), 'def G%d : ℕ := 0x%x' % (M, gmask(M)),
                'theorem mask%d : S%d = sMask %d ∧ G%d = gMask %d := by decide +kernel' % (M, M, M, M, M), '']
    out += ['/-- The masks of the moduli used: `(M, squares modulo M, Gaussian squares modulo M)`. -/',
            'def masks : List (ℕ × ℕ × ℕ) := [' + ', '.join('(%d, S%d, G%d)' % (M, M, M) for M in mods) + ']', '',
            'theorem masks_ok : MasksOk masks := by\n  intro e he\n  simp only [masks, List.mem_cons, List.not_mem_nil, or_false] at he\n  rcases he with ' + ' | '.join('rfl' for M in mods) + '\n' + ''.join('  · exact mask%d\n' % M for M in mods), '']
    CH = 8; chunks = [ok[i:i + CH] for i in range(0, len(ok), CH)]
    for j, ch in enumerate(chunks):
        out.append('def chunk%d : List (ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)) := [\n' % j + ',\n'.join(
            '  (%d, %d, %s, [%s])' % (r['m'], r['n'], r['primes'], ', '.join('(%d, %d, 0x%x)' % (p, k, cmask(r['primes8'], ds)) for p, k, ds in r['cert'])) for r in ch) + ']')
        out.append(''); out.append('theorem chunk%d_ok : chunk%d.all (entryOk8 masks) = true := by decide +kernel' % (j, j)); out.append('')
    out.append('/-- The table: `(m, n, primes, certificate)` for %d ratios. -/' % len(ok))
    out.append('def table : List (ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)) :=\n  ' + ' ++ ('.join('chunk%d' % j for j in range(len(chunks))) + ')'*(len(chunks) - 1)); out.append('')
    out.append('theorem table_ok : table.all (entryOk8 masks) = true := by\n  unfold table\n' + ''.join('  rw [List.all_append, chunk%d_ok, Bool.true_and]\n' % j for j in range(len(chunks) - 1)) + '  exact chunk%d_ok' % (len(chunks) - 1)); out.append('')
    out.append('/-- Certificates with four classes (`entryOk4`), for the ratios that have none with one class:\n`(m, n, primes, lines)`; a line `(p, k)` removes the tuples that have no point modulo `p^k`. -/')
    out.append('def table4 : List (ℕ × ℕ × List ℕ × List (ℕ × ℕ)) := [\n' + ',\n'.join('  (%d, %d, %s, [%s])' % (r['m'], r['n'], r['primes'], ', '.join('(%d, %d)' % l for l in r['lines'])) for r in ok4) + ']'); out.append('')
    out.append('theorem table4_ok : table4.all (entryOk4 masks) = true := by decide +kernel'); out.append('')
    st = [r for r in (struct(p) for p in pairs if triv(*p)) if r['left'] == 1]
    out.append('/-- Ratios with a degenerate solution for which the theorems on the class leave one candidate:\n`(m, n, primes, class, primes of (m² + n²)(m⁸ + 14m⁴n⁴ + n⁸))`. -/')
    SCH = 10; schunks = [st[i:i + SCH] for i in range(0, len(st), SCH)]
    for j, ch in enumerate(schunks):
        out.append('def structChunk%d : List (ℕ × ℕ × List ℕ × ℕ × List ℕ) := [\n' % j + ',\n'.join('  (%d, %d, %s, %d, %s)' % (r['m'], r['n'], r['primes'], r['d0'], r['qs']) for r in ch) + ']'); out.append('')
        out.append('theorem structChunk%d_ok : structChunk%d.all structOk = true := by decide +kernel' % (j, j)); out.append('')
    out.append('/-- The table of the classes: %d ratios. -/\ndef structTable : List (ℕ × ℕ × List ℕ × ℕ × List ℕ) :=\n  ' % len(st) + ' ++ ('.join('structChunk%d' % j for j in range(len(schunks))) + ')'*(len(schunks) - 1)); out.append('')
    out.append('theorem structTable_ok : structTable.all structOk = true := by\n  unfold structTable\n' + ''.join('  rw [List.all_append, structChunk%d_ok, Bool.true_and]\n' % j for j in range(len(schunks) - 1)) + '  exact structChunk%d_ok' % (len(schunks) - 1)); out.append('')
    out.append('/-- **For the %d ratios of `structTable`, among them `627 : 364` and `991 : 263` with the class\n`521266`, the class of every solution on the four forms is the one of the table**: that number times\neach of the eight products is a Gaussian square. -/\ntheorem class_of_mem_structTable {e : ℕ × ℕ × List ℕ × ℕ × List ℕ} (he : e ∈ structTable) {u v : ℤ}\n    (hcop : IsCoprime u v)\n    (h0 : ∃ w, nrm (P0 e.1 e.2.1 u v) = w ^ 2) (h8 : ∃ w, nrm (P8 e.1 e.2.1 u v) = w ^ 2)\n    (hL : ∃ w, nrm (PL e.1 e.2.1 u v) = w ^ 2) (hU : ∃ w, nrm (PU e.1 e.2.1 u v) = w ^ 2) :\n    ∀ f : Bool × Bool × Bool, GSqP (e.2.2.2.1 : ℤ) (Q e.1 e.2.1 u v f) :=\n  class_of_structOk ((List.all_eq_true.mp structTable_ok) e he) hcop h0 h8 hL hU' % len(st)); out.append('')
    out.append('theorem structTable_length : structTable.length = %d := by decide +kernel' % len(st)); out.append('')
    out.append('theorem structTable_last_pair : (627, 364, 521266) ∈ structTable.map (fun e => (e.1, e.2.1, e.2.2.2.1)) ∧\n    (991, 263, 521266) ∈ structTable.map (fun e => (e.1, e.2.1, e.2.2.2.1)) := by decide +kernel'); out.append('')
    out.append('/-- The ratios of the two tables. -/\ndef ratios : List (ℕ × ℕ) := [' + ', '.join('(%d, %d)' % (r['m'], r['n']) for r in ok + ok4) + ']'); out.append('')
    out.append('theorem ratios_eq :\n    table.map (fun e => (e.1, e.2.1)) ++ table4.map (fun e => (e.1, e.2.1)) = ratios := by decide +kernel'); out.append('')
    text = HEAD + '\n'.join(out) + TAIL.replace('EXCEPTIONS_LEFT', str(len(pairs) - len(ok) - len(ok4)))
    if '--check' in sys.argv:
        same = open('MagicSquareOfSquares/CoverDescentTable.lean').read() == text
        print('CoverDescentTable.lean is the output of this program: %s' % same); sys.exit(0 if same else 1)
    open('MagicSquareOfSquares/CoverDescentTable.lean', 'w').write(text)
    print('pairs of exceptions1000: %d; not of the two exceptional types: %d; with a certificate: %d; moduli: %s' % (len(pairs), len(res), len(ok), mods))
    print('primes of the lines by Euler\'s criterion: %s' % sorted(set(p for r in ok for p, k, _ in r['cert'] if not k)))
    print('without a certificate with one class (moduli up to %d, Euler lines at primes up to %d): %s' % (CAP, ECAP, [((r['m'], r['n']), r['left']) for r in res if r['left']]))
    print('ratios with a degenerate solution whose class is fixed by E, F, G, O (primes below %d): %d; largest prime used: %d' % (SBOUND, len(st), max(q for r in st for q in r['qs'])))
    print('with four classes: %s' % [((r['m'], r['n']), r['tuples'], r['lines'], r['left']) for r in res4])
