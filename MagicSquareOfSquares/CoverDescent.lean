import Mathlib

/-!
# A descent on the four forms, with a finite check (`LAST_PAIR.md`, sections 11 and 15)

For a ratio `m : n` and coprime integers `u, v`, the four forms of the G1 family that are not
`G4` are norms of pairs `x + i y`:

`F0 = N(mv + i·nu)`, `F8 = N(nv + i·mu)`, `GL = N((m+n)(u−v) + i·(m−n)(u+v))`,
`GU = N((m−n)(u−v) + i·(m+n)(u+v))`.

Everything is done in integer coordinates; no Gaussian integers, no unique factorisation.

* `GSq d x y`: `d (x + i y)` is the square of some `a + i b`.
* `exists_gsq_of_norm`: a pair of square norm satisfies `GSq k x y` for some `k > 0`
  (take `a + i b = x + w + i y`).
* `GSq.reduce`: a square factor of the multiplier can be removed (`e ∣ a` and `e ∣ b`).
* `GSq.support`: an odd prime dividing a multiplier without square factors divides `x` and `y`.
* `exists_class_prod`: **if the four forms are squares, some `d > 0` without square factors,
  whose prime factors divide `2mn(m − n)(m + n)`, makes `d` times each of the eight products
  `P0 · P8^± · PL^± · PU^±` a Gaussian square.**
* `kills`: a computable check modulo `M = p^k`, over the points of the projective line, with bit
  masks of the squares and of the Gaussian squares modulo `M`.
* `not_gsq_of_kills`: if the check succeeds for a list of `d`, none of them is the class of a
  solution.
* `no_point_of_entry`: if every product of the primes of `2mn(m − n)(m + n)` is excluded by some
  modulus, there is no solution at all.
* `exists_class_struct`: **for coprime `m, n` the class `d` has three more properties, for every
  ratio:** its odd prime factors are `≡ ±1 (mod 8)` (`mod_eight_of_dv`), it is even
  (`two_dvd_class`), it is a square modulo every odd prime dividing `m² + n²`
  (`isSquare_class_mod`), and modulo every odd prime dividing `m⁸ + 14m⁴n⁴ + n⁸`
  (`isSquare_class_mod_octic`).
* `killsE`, `not_gsq_of_killsE`: a second check, modulo an odd prime `p`, by the criterion of
  Euler in the two residue fields above `p`; no mask, so `p` may be large.
* `entryOk8`, `no_point_of_entry8`: the certificates of the table, with the candidates that the
  first of the three properties leaves and with lines of both kinds.
* `structOk`, `class_of_structOk`: for a ratio where those properties leave one candidate, the class
  of every solution is that candidate; no local computation at the points.
* `exists_classes4`, `recs4`, `entryOk4`, `no_point_of_entry4`: the same with one class for each
  of the four forms, for the ratios where one class is not enough.

The certificates and the statements about magic squares are in `CoverDescentTable.lean`.
-/

namespace MagicSquareOfSquares.CoverDescent

/-! ### Gaussian squares, in integer coordinates -/

/-- `d (x + i y) = (a + i b)²` for some integers `a, b`. -/
def GSq (d x y : ℤ) : Prop := ∃ a b : ℤ, d * x = a ^ 2 - b ^ 2 ∧ d * y = 2 * a * b

theorem GSq.conj {d x y : ℤ} (h : GSq d x y) : GSq d x (-y) := by
  obtain ⟨a, b, h1, h2⟩ := h
  exact ⟨a, -b, by linear_combination h1, by linear_combination -h2⟩

theorem GSq.neg_left {d x y : ℤ} (h : GSq d x y) : GSq (-d) x y := by
  obtain ⟨a, b, h1, h2⟩ := h
  exact ⟨b, -a, by linear_combination -h1, by linear_combination -h2⟩

theorem GSq.mul {d₁ d₂ x₁ y₁ x₂ y₂ : ℤ} (h₁ : GSq d₁ x₁ y₁) (h₂ : GSq d₂ x₂ y₂) :
    GSq (d₁ * d₂) (x₁ * x₂ - y₁ * y₂) (x₁ * y₂ + y₁ * x₂) := by
  obtain ⟨a₁, b₁, e₁, f₁⟩ := h₁
  obtain ⟨a₂, b₂, e₂, f₂⟩ := h₂
  refine ⟨a₁ * a₂ - b₁ * b₂, a₁ * b₂ + b₁ * a₂, ?_, ?_⟩
  · linear_combination (d₂ * x₂) * e₁ + (a₁ ^ 2 - b₁ ^ 2) * e₂ - (d₂ * y₂) * f₁ - (2 * a₁ * b₁) * f₂
  · linear_combination (d₂ * y₂) * e₁ + (a₁ ^ 2 - b₁ ^ 2) * f₂ + (d₂ * x₂) * f₁ + (2 * a₁ * b₁) * e₂

/-- A point of square norm is a nonzero integer times a Gaussian square. -/
theorem exists_gsq_of_norm {x y w : ℤ} (h : x ^ 2 + y ^ 2 = w ^ 2) (hne : x ≠ 0 ∨ y ≠ 0) :
    ∃ k : ℕ, 0 < k ∧ GSq k x y := by
  have key : ∃ r : ℤ, r ≠ 0 ∧ GSq r x y := by
    by_cases hw : x + w = 0
    · refine ⟨2 * (x - w), ?_, x - w, y, ?_, ?_⟩
      · intro h0
        have hx : x = 0 := by linarith
        have hw0 : w = 0 := by linarith
        rw [hx, hw0] at h
        have hy : y = 0 := by
          have : y ^ 2 = 0 := by linarith
          exact pow_eq_zero_iff (two_ne_zero) |>.mp this
        rcases hne with h1 | h1
        · exact h1 hx
        · exact h1 hy
      · linear_combination h
      · ring
    · refine ⟨2 * (x + w), ?_, x + w, y, ?_, ?_⟩
      · intro h0
        apply hw
        linarith
      · linear_combination h
      · ring
  obtain ⟨r, hr, hg⟩ := key
  obtain ⟨k, rfl | rfl⟩ := Int.eq_nat_or_neg r
  · exact ⟨k, Nat.pos_of_ne_zero (by rintro rfl; exact hr rfl), hg⟩
  · refine ⟨k, Nat.pos_of_ne_zero (by rintro rfl; exact hr (by simp)), ?_⟩
    have := hg.neg_left
    rwa [neg_neg] at this

theorem dvd_both {q : ℤ} (hq : Prime q) {a b : ℤ} (hab : q ∣ a * b) (h1 : q ∣ a ^ 2 - b ^ 2) :
    q ∣ a ∧ q ∣ b := by
  rcases hq.dvd_or_dvd hab with h | h
  · refine ⟨h, ?_⟩
    have h2 : q ∣ a ^ 2 - (a ^ 2 - b ^ 2) := dvd_sub (dvd_pow h two_ne_zero) h1
    rw [sub_sub_cancel] at h2
    exact hq.dvd_of_dvd_pow h2
  · refine ⟨?_, h⟩
    have h2 : q ∣ a ^ 2 - b ^ 2 + b ^ 2 := dvd_add h1 (dvd_pow h two_ne_zero)
    rw [sub_add_cancel] at h2
    exact hq.dvd_of_dvd_pow h2

theorem dvd_of_sq_dvd {q : ℕ} (hq : q.Prime) {a b : ℤ} (h1 : (q : ℤ) ^ 2 ∣ a ^ 2 - b ^ 2)
    (h2 : (q : ℤ) ^ 2 ∣ 2 * a * b) : (q : ℤ) ∣ a ∧ (q : ℤ) ∣ b := by
  have hqi : Prime (q : ℤ) := Nat.prime_iff_prime_int.mp hq
  have h1' : (q : ℤ) ∣ a ^ 2 - b ^ 2 := dvd_trans (dvd_pow_self _ two_ne_zero) h1
  have hab : (q : ℤ) ∣ a * b := by
    obtain ⟨c, hc⟩ := h2
    have h3 : (q : ℤ) ∣ 2 * (a * b) := ⟨q * c, by linear_combination hc⟩
    rcases hqi.dvd_or_dvd h3 with h | h
    · have hq2' : q ∣ 2 := by exact_mod_cast h
      have hq2 : q = 2 := (Nat.prime_dvd_prime_iff_eq hq Nat.prime_two).mp hq2'
      subst hq2
      refine ⟨c, ?_⟩
      push_cast at hc ⊢
      linarith
    · exact h
  exact dvd_both hqi hab h1'

/-- A square factor of the multiplier can be removed. -/
theorem GSq.reduce : ∀ e : ℕ, 0 < e → ∀ d x y : ℤ, GSq (d * (e : ℤ) ^ 2) x y → GSq d x y := by
  intro e
  induction e using Nat.strong_induction_on with
  | _ e ih =>
    intro he d x y h
    by_cases h1 : e = 1
    · subst h1
      simpa using h
    obtain ⟨q, hq, e', rfl⟩ := Nat.exists_prime_and_dvd h1
    obtain ⟨a, b, ha, hb⟩ := h
    have hq0 : (q : ℤ) ≠ 0 := by exact_mod_cast hq.ne_zero
    push_cast at ha hb
    have d1 : (q : ℤ) ^ 2 ∣ a ^ 2 - b ^ 2 := ⟨d * (e' : ℤ) ^ 2 * x, by linear_combination -ha⟩
    have d2 : (q : ℤ) ^ 2 ∣ 2 * a * b := ⟨d * (e' : ℤ) ^ 2 * y, by linear_combination -hb⟩
    obtain ⟨⟨a', rfl⟩, ⟨b', rfl⟩⟩ := dvd_of_sq_dvd hq d1 d2
    have he' : 0 < e' := Nat.pos_of_ne_zero (by rintro rfl; simp at he)
    have hlt : e' < q * e' := lt_mul_of_one_lt_left he' hq.one_lt
    apply ih e' hlt he' d x y
    refine ⟨a', b', ?_, ?_⟩
    · apply mul_left_cancel₀ (pow_ne_zero 2 hq0)
      linear_combination ha
    · apply mul_left_cancel₀ (pow_ne_zero 2 hq0)
      linear_combination hb

/-- No square of a prime divides `d`. -/
def NoSq (d : ℕ) : Prop := ∀ q : ℕ, q.Prime → ¬q ^ 2 ∣ d

theorem exists_nosq : ∀ n : ℕ, 0 < n → ∃ d e : ℕ, 0 < d ∧ 0 < e ∧ n = d * e ^ 2 ∧ NoSq d := by
  intro n
  induction n using Nat.strong_induction_on with
  | _ n ih =>
    intro hn
    by_cases h : ∃ q : ℕ, q.Prime ∧ q ^ 2 ∣ n
    swap
    · exact ⟨n, 1, hn, one_pos, by ring, fun q hq hdvd => h ⟨q, hq, hdvd⟩⟩
    · obtain ⟨q, hq, n', rfl⟩ := h
      have hn' : 0 < n' := Nat.pos_of_ne_zero (by rintro rfl; simp at hn)
      have hq2 : 1 < q ^ 2 := Nat.one_lt_pow two_ne_zero hq.one_lt
      have hlt : n' < q ^ 2 * n' := lt_mul_of_one_lt_left hn' hq2
      obtain ⟨d, e, hd, he, hde, hns⟩ := ih n' hlt hn'
      exact ⟨d, q * e, hd, Nat.mul_pos hq.pos he, by rw [hde]; ring, hns⟩

/-- An odd prime that divides a multiplier without square factors divides both coordinates. -/
theorem GSq.support {d : ℕ} {x y : ℤ} (h : GSq d x y) {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hpd : p ∣ d) (hsq : ¬p ^ 2 ∣ d) : (p : ℤ) ∣ x ∧ (p : ℤ) ∣ y := by
  obtain ⟨a, b, ha, hb⟩ := h
  obtain ⟨d', rfl⟩ := hpd
  have hpi : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hp0 : (p : ℤ) ≠ 0 := by exact_mod_cast hp.ne_zero
  have hpd' : ¬(p : ℤ) ∣ (d' : ℤ) := by
    intro hh
    have h3 : p ∣ d' := by exact_mod_cast hh
    obtain ⟨c, rfl⟩ := h3
    exact hsq ⟨c, by ring⟩
  have hp2' : ¬(p : ℤ) ∣ 2 := by
    intro hh
    have h3 : p ∣ 2 := by exact_mod_cast hh
    exact hp2 ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp h3)
  push_cast at ha hb
  have e1 : (p : ℤ) ∣ a ^ 2 - b ^ 2 := ⟨d' * x, by linear_combination -ha⟩
  have e2 : (p : ℤ) ∣ 2 * (a * b) := ⟨d' * y, by linear_combination -hb⟩
  have hab : (p : ℤ) ∣ a * b := (hpi.dvd_or_dvd e2).resolve_left hp2'
  obtain ⟨⟨a', rfl⟩, ⟨b', rfl⟩⟩ := dvd_both hpi hab e1
  constructor
  · have h3 : (p : ℤ) ∣ d' * x := ⟨a' ^ 2 - b' ^ 2, by
      apply mul_left_cancel₀ hp0
      linear_combination ha⟩
    exact (hpi.dvd_or_dvd h3).resolve_left hpd'
  · have h3 : (p : ℤ) ∣ d' * y := ⟨2 * a' * b', by
      apply mul_left_cancel₀ hp0
      linear_combination hb⟩
    exact (hpi.dvd_or_dvd h3).resolve_left hpd'

/-- The class of a point of square norm: a positive multiplier without square factors whose odd
prime factors divide both coordinates. -/
theorem exists_class {x y w : ℤ} (h : x ^ 2 + y ^ 2 = w ^ 2) (hne : x ≠ 0 ∨ y ≠ 0) :
    ∃ d : ℕ, 0 < d ∧ GSq d x y ∧
      ∀ p : ℕ, p.Prime → p ∣ d → p = 2 ∨ ((p : ℤ) ∣ x ∧ (p : ℤ) ∣ y) := by
  obtain ⟨k, hk, hg⟩ := exists_gsq_of_norm h hne
  obtain ⟨d, e, hd, he, hde, hns⟩ := exists_nosq k hk
  have hg' : GSq (d : ℤ) x y := by
    apply GSq.reduce e he
    have e1 : ((k : ℕ) : ℤ) = (d : ℤ) * (e : ℤ) ^ 2 := by exact_mod_cast hde
    rw [← e1]
    exact hg
  refine ⟨d, hd, hg', ?_⟩
  intro p hp hpd
  by_cases hp2 : p = 2
  · exact Or.inl hp2
  · exact Or.inr (hg'.support hp hp2 hpd (hns p hp))

/-! ### The four linear forms and two products -/

/-- Multiplication of `x + i y` and `x' + i y'`, on pairs. -/
def gm (p q : ℤ × ℤ) : ℤ × ℤ := (p.1 * q.1 - p.2 * q.2, p.1 * q.2 + p.2 * q.1)

/-- Conjugation, on pairs. -/
def cj (p : ℤ × ℤ) : ℤ × ℤ := (p.1, -p.2)

/-- The norm, on pairs. -/
def nrm (p : ℤ × ℤ) : ℤ := p.1 ^ 2 + p.2 ^ 2

def P0 (m n u v : ℤ) : ℤ × ℤ := (m * v, n * u)
def P8 (m n u v : ℤ) : ℤ × ℤ := (n * v, m * u)
def PL (m n u v : ℤ) : ℤ × ℤ := ((m + n) * (u - v), (m - n) * (u + v))
def PU (m n u v : ℤ) : ℤ × ℤ := ((m - n) * (u - v), (m + n) * (u + v))

/-- Conjugation, or not. -/
def cjIf (c : Bool) (p : ℤ × ℤ) : ℤ × ℤ := if c then cj p else p

/-- The eight products `P0 · P8^± · PL^± · PU^±`; the flags say which factors are conjugated. -/
def Q (m n u v : ℤ) (f : Bool × Bool × Bool) : ℤ × ℤ :=
  gm (gm (gm (P0 m n u v) (cjIf f.1 (P8 m n u v))) (cjIf f.2.1 (PL m n u v))) (cjIf f.2.2 (PU m n u v))

/-- `GSq` on a pair. -/
def GSqP (d : ℤ) (p : ℤ × ℤ) : Prop := GSq d p.1 p.2

theorem GSqP.gm {d₁ d₂ : ℤ} {p q : ℤ × ℤ} (h₁ : GSqP d₁ p) (h₂ : GSqP d₂ q) :
    GSqP (d₁ * d₂) (gm p q) := GSq.mul h₁ h₂

theorem GSqP.cj {d : ℤ} {p : ℤ × ℤ} (h : GSqP d p) : GSqP d (cj p) := GSq.conj h

theorem GSqP.cjIf {d : ℤ} {p : ℤ × ℤ} (h : GSqP d p) (c : Bool) : GSqP d (cjIf c p) := by
  cases c
  · exact h
  · exact GSqP.cj h

theorem dvd_bad_of_dvd_mn {m n : ℤ} {q : ℤ} (h : q ∣ m * n) : q ∣ 2 * m * n * (m - n) * (m + n) := by
  obtain ⟨c, hc⟩ := h
  exact ⟨c * (2 * (m - n) * (m + n)), by linear_combination (2 * (m - n) * (m + n)) * hc⟩

theorem dvd_bad_of_dvd_diff {m n : ℤ} {q : ℤ} (h : q ∣ 2 * (m + n) * (m - n)) :
    q ∣ 2 * m * n * (m - n) * (m + n) := by
  obtain ⟨c, hc⟩ := h
  exact ⟨c * (m * n), by linear_combination (m * n) * hc⟩

/-- **The class of a solution.** If the four forms are squares at coprime `u, v`, some positive
`d` without square factors, whose prime factors divide `2mn(m − n)(m + n)`, makes `d` times each of
the eight products the square of a Gaussian integer. -/
theorem exists_class_prod {m n u v : ℤ} (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (n0 : (P0 m n u v).1 ≠ 0 ∨ (P0 m n u v).2 ≠ 0) (n8 : (P8 m n u v).1 ≠ 0 ∨ (P8 m n u v).2 ≠ 0)
    (nL : (PL m n u v).1 ≠ 0 ∨ (PL m n u v).2 ≠ 0) (nU : (PU m n u v).1 ≠ 0 ∨ (PU m n u v).2 ≠ 0) :
    ∃ d : ℕ, 0 < d ∧ NoSq d ∧
      (∀ p : ℕ, p.Prime → p ∣ d → (p : ℤ) ∣ 2 * m * n * (m - n) * (m + n)) ∧
      ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f) := by
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  obtain ⟨d0, hd0, g0, s0⟩ := exists_class hw0 n0
  obtain ⟨d8, hd8, g8, s8⟩ := exists_class hw8 n8
  obtain ⟨dL, hdL, gL, sL⟩ := exists_class hwL nL
  obtain ⟨dU, hdU, gU, sU⟩ := exists_class hwU nU
  obtain ⟨a, b, hab⟩ := hcop
  have hN : 0 < d0 * d8 * dL * dU := Nat.mul_pos (Nat.mul_pos (Nat.mul_pos hd0 hd8) hdL) hdU
  obtain ⟨d, e, hd, he, hde, hns⟩ := exists_nosq _ hN
  have e1 : (d0 : ℤ) * d8 * dL * dU = (d : ℤ) * (e : ℤ) ^ 2 := by exact_mod_cast hde
  have two_dvd : (2 : ℤ) ∣ 2 * m * n * (m - n) * (m + n) := ⟨m * n * (m - n) * (m + n), by ring⟩
  refine ⟨d, hd, hns, ?_, ?_⟩
  · intro p hp hpd
    have hpN : p ∣ d0 * d8 * dL * dU := by
      rw [hde]
      exact Dvd.dvd.mul_right hpd _
    rcases (Nat.Prime.dvd_mul hp).mp hpN with h3 | h3
    · rcases (Nat.Prime.dvd_mul hp).mp h3 with h4 | h4
      · rcases (Nat.Prime.dvd_mul hp).mp h4 with h5 | h5
        · rcases s0 p hp h5 with rfl | ⟨hx, hy⟩
          · exact_mod_cast two_dvd
          · apply dvd_bad_of_dvd_mn
            have e2 : m * n = b * n * (m * v) + a * m * (n * u) := by
              linear_combination (-(m * n)) * hab
            rw [e2]
            exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
        · rcases s8 p hp h5 with rfl | ⟨hx, hy⟩
          · exact_mod_cast two_dvd
          · apply dvd_bad_of_dvd_mn
            have e2 : m * n = b * m * (n * v) + a * n * (m * u) := by
              linear_combination (-(m * n)) * hab
            rw [e2]
            exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
      · rcases sL p hp h4 with rfl | ⟨hx, hy⟩
        · exact_mod_cast two_dvd
        · apply dvd_bad_of_dvd_diff
          have e2 : 2 * (m + n) * (m - n) =
              (a - b) * (m - n) * ((m + n) * (u - v)) + (a + b) * (m + n) * ((m - n) * (u + v)) := by
            linear_combination (-(2 * (m + n) * (m - n))) * hab
          rw [e2]
          exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
    · rcases sU p hp h3 with rfl | ⟨hx, hy⟩
      · exact_mod_cast two_dvd
      · apply dvd_bad_of_dvd_diff
        have e2 : 2 * (m + n) * (m - n) =
            (a - b) * (m + n) * ((m - n) * (u - v)) + (a + b) * (m - n) * ((m + n) * (u + v)) := by
          linear_combination (-(2 * (m + n) * (m - n))) * hab
        rw [e2]
        exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
  · intro f
    apply GSq.reduce e he
    rw [← e1]
    exact GSqP.gm (GSqP.gm (GSqP.gm g0 (GSqP.cjIf g8 f.1)) (GSqP.cjIf gL f.2.1)) (GSqP.cjIf gU f.2.2)

/-! ### The primes of the class: Theorem E -/

/-- `p` divides both coordinates of a pair. -/
def Dv (p : ℕ) (q : ℤ × ℤ) : Prop := (p : ℤ) ∣ q.1 ∧ (p : ℤ) ∣ q.2

/-- The class of a solution, with the reason for each of its primes: a prime dividing it is 2 or
divides both coordinates of one of the four pairs. -/
theorem exists_class_core {m n u v : ℤ}
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (n0 : (P0 m n u v).1 ≠ 0 ∨ (P0 m n u v).2 ≠ 0) (n8 : (P8 m n u v).1 ≠ 0 ∨ (P8 m n u v).2 ≠ 0)
    (nL : (PL m n u v).1 ≠ 0 ∨ (PL m n u v).2 ≠ 0) (nU : (PU m n u v).1 ≠ 0 ∨ (PU m n u v).2 ≠ 0) :
    ∃ d : ℕ, 0 < d ∧ NoSq d ∧
      (∀ p : ℕ, p.Prime → p ∣ d → p = 2 ∨ Dv p (P0 m n u v) ∨ Dv p (P8 m n u v) ∨
        Dv p (PL m n u v) ∨ Dv p (PU m n u v)) ∧
      ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f) := by
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  obtain ⟨d0, hd0, g0, s0⟩ := exists_class hw0 n0
  obtain ⟨d8, hd8, g8, s8⟩ := exists_class hw8 n8
  obtain ⟨dL, hdL, gL, sL⟩ := exists_class hwL nL
  obtain ⟨dU, hdU, gU, sU⟩ := exists_class hwU nU
  have hN : 0 < d0 * d8 * dL * dU := Nat.mul_pos (Nat.mul_pos (Nat.mul_pos hd0 hd8) hdL) hdU
  obtain ⟨d, e, hd, he, hde, hns⟩ := exists_nosq _ hN
  have e1 : (d0 : ℤ) * d8 * dL * dU = (d : ℤ) * (e : ℤ) ^ 2 := by exact_mod_cast hde
  refine ⟨d, hd, hns, ?_, ?_⟩
  · intro p hp hpd
    have hpN : p ∣ d0 * d8 * dL * dU := by
      rw [hde]
      exact Dvd.dvd.mul_right hpd _
    rcases (Nat.Prime.dvd_mul hp).mp hpN with h3 | h3
    · rcases (Nat.Prime.dvd_mul hp).mp h3 with h4 | h4
      · rcases (Nat.Prime.dvd_mul hp).mp h4 with h5 | h5
        · rcases s0 p hp h5 with h | h
          · exact Or.inl h
          · exact Or.inr (Or.inl h)
        · rcases s8 p hp h5 with h | h
          · exact Or.inl h
          · exact Or.inr (Or.inr (Or.inl h))
      · rcases sL p hp h4 with h | h
        · exact Or.inl h
        · exact Or.inr (Or.inr (Or.inr (Or.inl h)))
    · rcases sU p hp h3 with h | h
      · exact Or.inl h
      · exact Or.inr (Or.inr (Or.inr (Or.inr h)))
  · intro f
    apply GSq.reduce e he
    rw [← e1]
    exact GSqP.gm (GSqP.gm (GSqP.gm g0 (GSqP.cjIf g8 f.1)) (GSqP.cjIf gL f.2.1)) (GSqP.cjIf gU f.2.2)

/-- If twice a nonzero square is a square modulo an odd prime `p`, then `p ≡ ±1 (mod 8)`. -/
theorem mod_eight_of_two_sq {p : ℕ} [Fact p.Prime] (hp2 : p ≠ 2) {a w : ZMod p} (ha : a ≠ 0)
    (h : w ^ 2 = 2 * a ^ 2) : p % 8 = 1 ∨ p % 8 = 7 := by
  apply (ZMod.exists_sq_eq_two_iff hp2).mp
  refine ⟨w / a, ?_⟩
  field_simp
  linear_combination -h

/-- **Theorem E.** Let `m, n` be coprime and `u, v` coprime, with `F0` and `GL` squares. An odd
prime that divides both coordinates of one of the four pairs is `≡ ±1 (mod 8)`: at such a point one
of the two forms is congruent to twice a nonzero square. -/
theorem mod_eight_of_dv {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (hL : ∃ w, nrm (PL m n u v) = w ^ 2)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (h : Dv p (P0 m n u v) ∨ Dv p (P8 m n u v) ∨ Dv p (PL m n u v) ∨ Dv p (PU m n u v)) :
    p % 8 = 1 ∨ p % 8 = 7 := by
  have : Fact p.Prime := ⟨hp⟩
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨wL, hwL⟩ := hL
  have z : ∀ x : ℤ, (p : ℤ) ∣ x → ((x : ℤ) : ZMod p) = 0 := fun x hx =>
    (ZMod.intCast_zmod_eq_zero_iff_dvd x p).mpr hx
  have e0 : ((m : ZMod p) * v) ^ 2 + ((n : ZMod p) * u) ^ 2 = (w0 : ZMod p) ^ 2 := by
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hw0
    simp only [nrm, P0] at h1
    push_cast at h1
    exact h1
  have eL : (((m : ZMod p) + n) * (u - v)) ^ 2 + (((m : ZMod p) - n) * (u + v)) ^ 2 =
      (wL : ZMod p) ^ 2 := by
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hwL
    simp only [nrm, PL] at h1
    push_cast at h1
    exact h1
  have hmn' : ¬((m : ZMod p) = 0 ∧ (n : ZMod p) = 0) := by
    rintro ⟨hm, hn⟩
    obtain ⟨a, b, hab⟩ := hmn
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hm, hn, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have hcop' : ¬((u : ZMod p) = 0 ∧ (v : ZMod p) = 0) := by
    rintro ⟨hu, hv⟩
    obtain ⟨a, b, hab⟩ := hcop
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hu, hv, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have two_ne : (2 : ZMod p) ≠ 0 := by
    intro h2
    have h3 : ((2 : ℕ) : ZMod p) = 0 := by exact_mod_cast h2
    have h4 : p ∣ 2 := (ZMod.natCast_eq_zero_iff 2 p).mp h3
    have h5 : p ≤ 2 := Nat.le_of_dvd (by norm_num) h4
    have h6 := hp.two_le
    omega
  have half : ∀ x y : ZMod p, x + y = 0 → x - y = 0 → x = 0 ∧ y = 0 := by
    intro x y h1 h2
    have hx : (2 : ZMod p) * x = 0 := by linear_combination h1 + h2
    have hx' : x = 0 := (mul_eq_zero.mp hx).resolve_left two_ne
    refine ⟨hx', ?_⟩
    rw [hx', zero_add] at h1
    exact h1
  rcases h with ⟨hx, hy⟩ | ⟨hx, hy⟩ | ⟨hx, hy⟩ | ⟨hx, hy⟩
  · -- (m v, n u)
    have h1 := z _ hx
    have h2 := z _ hy
    simp only [P0] at h1 h2
    push_cast at h1 h2
    rcases mul_eq_zero.mp h1 with hm | hv
    · rcases mul_eq_zero.mp h2 with hn | hu
      · exact absurd ⟨hm, hn⟩ hmn'
      · have hn : (n : ZMod p) ≠ 0 := fun hn => hmn' ⟨hm, hn⟩
        have hv : (v : ZMod p) ≠ 0 := fun hv => hcop' ⟨hu, hv⟩
        refine mod_eight_of_two_sq (w := (wL : ZMod p)) hp2 (mul_ne_zero hn hv) ?_
        rw [← eL, hm, hu]
        ring
    · rcases mul_eq_zero.mp h2 with hn | hu
      · have hm : (m : ZMod p) ≠ 0 := fun hm => hmn' ⟨hm, hn⟩
        have hu : (u : ZMod p) ≠ 0 := fun hu => hcop' ⟨hu, hv⟩
        refine mod_eight_of_two_sq (w := (wL : ZMod p)) hp2 (mul_ne_zero hm hu) ?_
        rw [← eL, hn, hv]
        ring
      · exact absurd ⟨hu, hv⟩ hcop'
  · -- (n v, m u)
    have h1 := z _ hx
    have h2 := z _ hy
    simp only [P8] at h1 h2
    push_cast at h1 h2
    rcases mul_eq_zero.mp h1 with hn | hv
    · rcases mul_eq_zero.mp h2 with hm | hu
      · exact absurd ⟨hm, hn⟩ hmn'
      · have hm : (m : ZMod p) ≠ 0 := fun hm => hmn' ⟨hm, hn⟩
        have hv : (v : ZMod p) ≠ 0 := fun hv => hcop' ⟨hu, hv⟩
        refine mod_eight_of_two_sq (w := (wL : ZMod p)) hp2 (mul_ne_zero hm hv) ?_
        rw [← eL, hn, hu]
        ring
    · rcases mul_eq_zero.mp h2 with hm | hu
      · have hn : (n : ZMod p) ≠ 0 := fun hn => hmn' ⟨hm, hn⟩
        have hu : (u : ZMod p) ≠ 0 := fun hu => hcop' ⟨hu, hv⟩
        refine mod_eight_of_two_sq (w := (wL : ZMod p)) hp2 (mul_ne_zero hn hu) ?_
        rw [← eL, hm, hv]
        ring
      · exact absurd ⟨hu, hv⟩ hcop'
  · -- ((m + n)(u − v), (m − n)(u + v))
    have h1 := z _ hx
    have h2 := z _ hy
    simp only [PL] at h1 h2
    push_cast at h1 h2
    rcases mul_eq_zero.mp h1 with hA | hA
    · rcases mul_eq_zero.mp h2 with hB | hB
      · exact absurd (half _ _ hA hB) hmn'
      · have hn : (n : ZMod p) = -m := by linear_combination hA
        have hv : (v : ZMod p) = -u := by linear_combination hB
        have hm : (m : ZMod p) ≠ 0 := fun hm => hmn' ⟨hm, by rw [hn, hm, neg_zero]⟩
        have hu : (u : ZMod p) ≠ 0 := fun hu => hcop' ⟨hu, by rw [hv, hu, neg_zero]⟩
        refine mod_eight_of_two_sq (w := (w0 : ZMod p)) hp2 (mul_ne_zero hm hu) ?_
        rw [← e0, hn, hv]
        ring
    · rcases mul_eq_zero.mp h2 with hB | hB
      · have hn : (n : ZMod p) = m := by linear_combination -hB
        have hv : (v : ZMod p) = u := by linear_combination -hA
        have hm : (m : ZMod p) ≠ 0 := fun hm => hmn' ⟨hm, by rw [hn, hm]⟩
        have hu : (u : ZMod p) ≠ 0 := fun hu => hcop' ⟨hu, by rw [hv, hu]⟩
        refine mod_eight_of_two_sq (w := (w0 : ZMod p)) hp2 (mul_ne_zero hm hu) ?_
        rw [← e0, hn, hv]
        ring
      · exact absurd (half _ _ hB hA) hcop'
  · -- ((m − n)(u − v), (m + n)(u + v))
    have h1 := z _ hx
    have h2 := z _ hy
    simp only [PU] at h1 h2
    push_cast at h1 h2
    rcases mul_eq_zero.mp h1 with hA | hA
    · rcases mul_eq_zero.mp h2 with hB | hB
      · exact absurd (half _ _ hB hA) hmn'
      · have hn : (n : ZMod p) = m := by linear_combination -hA
        have hv : (v : ZMod p) = -u := by linear_combination hB
        have hm : (m : ZMod p) ≠ 0 := fun hm => hmn' ⟨hm, by rw [hn, hm]⟩
        have hu : (u : ZMod p) ≠ 0 := fun hu => hcop' ⟨hu, by rw [hv, hu, neg_zero]⟩
        refine mod_eight_of_two_sq (w := (w0 : ZMod p)) hp2 (mul_ne_zero hm hu) ?_
        rw [← e0, hn, hv]
        ring
    · rcases mul_eq_zero.mp h2 with hB | hB
      · have hn : (n : ZMod p) = -m := by linear_combination hB
        have hv : (v : ZMod p) = u := by linear_combination -hA
        have hm : (m : ZMod p) ≠ 0 := fun hm => hmn' ⟨hm, by rw [hn, hm, neg_zero]⟩
        have hu : (u : ZMod p) ≠ 0 := fun hu => hcop' ⟨hu, by rw [hv, hu]⟩
        refine mod_eight_of_two_sq (w := (w0 : ZMod p)) hp2 (mul_ne_zero hm hu) ?_
        rw [← e0, hn, hv]
        ring
      · exact absurd (half _ _ hB hA) hcop'

/-- **The class of a solution, with Theorem E.** For coprime `m, n`: the class `d` of
`exists_class_prod`, and in addition every prime dividing `d` is 2 or `≡ ±1 (mod 8)`. -/
theorem exists_class_prod8 {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (n0 : (P0 m n u v).1 ≠ 0 ∨ (P0 m n u v).2 ≠ 0) (n8 : (P8 m n u v).1 ≠ 0 ∨ (P8 m n u v).2 ≠ 0)
    (nL : (PL m n u v).1 ≠ 0 ∨ (PL m n u v).2 ≠ 0) (nU : (PU m n u v).1 ≠ 0 ∨ (PU m n u v).2 ≠ 0) :
    ∃ d : ℕ, 0 < d ∧ NoSq d ∧
      (∀ p : ℕ, p.Prime → p ∣ d → (p : ℤ) ∣ 2 * m * n * (m - n) * (m + n)) ∧
      (∀ p : ℕ, p.Prime → p ∣ d → p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧
      ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f) := by
  obtain ⟨d, hd, hns, hcore, g⟩ := exists_class_core h0 h8 hL hU n0 n8 nL nU
  obtain ⟨a, b, hab⟩ := hcop
  have two_dvd : (2 : ℤ) ∣ 2 * m * n * (m - n) * (m + n) := ⟨m * n * (m - n) * (m + n), by ring⟩
  refine ⟨d, hd, hns, ?_, ?_, g⟩
  · intro p hp hpd
    rcases hcore p hp hpd with rfl | ⟨hx, hy⟩ | ⟨hx, hy⟩ | ⟨hx, hy⟩ | ⟨hx, hy⟩
    · exact_mod_cast two_dvd
    · apply dvd_bad_of_dvd_mn
      have e2 : m * n = b * n * (m * v) + a * m * (n * u) := by
        linear_combination (-(m * n)) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
    · apply dvd_bad_of_dvd_mn
      have e2 : m * n = b * m * (n * v) + a * n * (m * u) := by
        linear_combination (-(m * n)) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
    · apply dvd_bad_of_dvd_diff
      have e2 : 2 * (m + n) * (m - n) =
          (a - b) * (m - n) * ((m + n) * (u - v)) + (a + b) * (m + n) * ((m - n) * (u + v)) := by
        linear_combination (-(2 * (m + n) * (m - n))) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
    · apply dvd_bad_of_dvd_diff
      have e2 : 2 * (m + n) * (m - n) =
          (a - b) * (m + n) * ((m - n) * (u - v)) + (a + b) * (m - n) * ((m + n) * (u + v)) := by
        linear_combination (-(2 * (m + n) * (m - n))) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
  · intro p hp hpd
    by_cases hp2 : p = 2
    · exact Or.inl hp2
    · exact Or.inr (mod_eight_of_dv hmn ⟨a, b, hab⟩ h0 hL hp hp2
        ((hcore p hp hpd).resolve_left hp2))

/-! ### The finite check -/

theorem testBit_foldl_or (f : ℕ → ℕ) (l : List ℕ) (acc i : ℕ)
    (h : acc.testBit i = true ∨ ∃ r ∈ l, f r = i) :
    (l.foldl (fun acc r => acc ||| 1 <<< f r) acc).testBit i = true := by
  induction l generalizing acc with
  | nil => simpa using h
  | cons r l ih =>
    rw [List.foldl_cons]
    apply ih
    rcases h with h | ⟨ρ, hρ, hi⟩
    · left
      rw [Nat.testBit_or, h, Bool.true_or]
    · rcases List.mem_cons.mp hρ with hr | hρ'
      · left
        subst hi
        subst hr
        rw [Nat.testBit_or, Nat.one_shiftLeft, Nat.testBit_two_pow_self, Bool.or_true]
      · right
        exact ⟨ρ, hρ', hi⟩

/-- The residue of an integer modulo `M`, as a natural number. -/
def res (M : ℕ) (z : ℤ) : ℕ := (z % (M : ℤ)).toNat

theorem res_lt {M : ℕ} (hM : 0 < M) (z : ℤ) : res M z < M := by
  have hM' : (0 : ℤ) < M := by exact_mod_cast hM
  have h1 := Int.emod_lt_of_pos z hM'
  have h2 := Int.emod_nonneg z hM'.ne'
  unfold res
  omega

theorem res_cast {M : ℕ} (hM : 0 < M) (z : ℤ) : ((res M z : ℕ) : ℤ) ≡ z [ZMOD M] := by
  have hM' : (0 : ℤ) < M := by exact_mod_cast hM
  unfold res
  rw [Int.toNat_of_nonneg (Int.emod_nonneg z hM'.ne')]
  exact Int.mod_modEq z M

theorem res_congr {M : ℕ} {z z' : ℤ} (h : z ≡ z' [ZMOD M]) : res M z = res M z' := by
  unfold res
  rw [h]

/-- The squares modulo `M`, as a bit mask. -/
def sMask (M : ℕ) : ℕ :=
  (List.range M).foldl (fun (acc r : ℕ) => acc ||| 1 <<< res M ((r : ℤ) ^ 2)) 0

/-- The position of the pair `(x, y)` modulo `M` in the mask of Gaussian squares. -/
def gIdx (M : ℕ) (x y : ℤ) : ℕ := M * res M x + res M y

/-- The squares `(a + i b)²` modulo `M`, as a bit mask: bit `M·x + y` for `(x, y) = (a² − b², 2ab)`. -/
def gMask (M : ℕ) : ℕ :=
  (List.range (M * M)).foldl
    (fun (acc k : ℕ) => acc ||| 1 <<< gIdx M (((k / M : ℕ) : ℤ) ^ 2 - ((k % M : ℕ) : ℤ) ^ 2)
      (2 * ((k / M : ℕ) : ℤ) * ((k % M : ℕ) : ℤ))) 0

def sqOk (M S : ℕ) (z : ℤ) : Bool := S.testBit (res M z)

theorem sqOk_of {M : ℕ} (hM : 0 < M) {z w : ℤ} (h : z ≡ w ^ 2 [ZMOD M]) :
    sqOk M (sMask M) z = true := by
  unfold sqOk sMask
  apply testBit_foldl_or (fun r => res M ((r : ℤ) ^ 2))
  right
  refine ⟨res M w, List.mem_range.mpr (res_lt hM w), res_congr ?_⟩
  exact ((res_cast hM w).pow 2).trans h.symm

theorem gOk_of {M : ℕ} (hM : 0 < M) {X Y a b : ℤ} (h1 : X ≡ a ^ 2 - b ^ 2 [ZMOD M])
    (h2 : Y ≡ 2 * a * b [ZMOD M]) : (gMask M).testBit (gIdx M X Y) = true := by
  unfold gMask
  apply testBit_foldl_or (fun k => gIdx M (((k / M : ℕ) : ℤ) ^ 2 - ((k % M : ℕ) : ℤ) ^ 2)
      (2 * ((k / M : ℕ) : ℤ) * ((k % M : ℕ) : ℤ)))
  right
  have hA := res_lt hM a
  have hB := res_lt hM b
  refine ⟨M * res M a + res M b, List.mem_range.mpr ?_, ?_⟩
  · calc M * res M a + res M b < M * res M a + M := by omega
      _ = M * (res M a + 1) := by ring
      _ ≤ M * M := Nat.mul_le_mul_left M hA
  · have hdiv : (M * res M a + res M b) / M = res M a := by
      rw [Nat.mul_add_div hM, Nat.div_eq_of_lt hB, add_zero]
    have hmod : (M * res M a + res M b) % M = res M b := by
      rw [Nat.mul_add_mod, Nat.mod_eq_of_lt hB]
    simp only [hdiv, hmod]
    unfold gIdx
    have ca := res_cast hM a
    have cb := res_cast hM b
    rw [res_congr (((ca.pow 2).sub (cb.pow 2)).trans h1.symm),
      res_congr (((ca.mul_left 2).mul cb).trans h2.symm)]

def formsOk (m n : ℤ) (M S : ℕ) (u v : ℤ) : Bool :=
  sqOk M S (nrm (P0 m n u v)) && sqOk M S (nrm (P8 m n u v)) &&
    sqOk M S (nrm (PL m n u v)) && sqOk M S (nrm (PU m n u v))

def gOk (M G d : ℕ) (p : ℤ × ℤ) : Bool := G.testBit (gIdx M ((d : ℤ) * p.1) ((d : ℤ) * p.2))

/-- The eight choices of conjugates. -/
def flags : List (Bool × Bool × Bool) :=
  [(false, false, false), (false, false, true), (false, true, false), (false, true, true),
   (true, false, false), (true, false, true), (true, true, false), (true, true, true)]

theorem mem_flags (f : Bool × Bool × Bool) : f ∈ flags := by
  obtain ⟨a, b, c⟩ := f
  cases a <;> cases b <;> cases c <;> simp [flags]

theorem res_eq_of {M k : ℕ} {z : ℤ} (hk : k < M) (h : (k : ℤ) ≡ z [ZMOD M]) : res M z = k := by
  rw [← res_congr h]
  unfold res
  have h1 : (k : ℤ) % (M : ℤ) = k := Int.emod_eq_of_lt (by positivity) (by exact_mod_cast hk)
  rw [h1]
  rfl

theorem res_mul {M : ℕ} (hM : 0 < M) (d : ℕ) (z : ℤ) : res M ((d : ℤ) * z) = d * res M z % M := by
  apply res_eq_of (Nat.mod_lt _ hM)
  have h1 : ((d * res M z % M : ℕ) : ℤ) ≡ ((d * res M z : ℕ) : ℤ) [ZMOD M] := by
    push_cast
    exact Int.mod_modEq _ _
  refine h1.trans ?_
  push_cast
  exact (res_cast hM z).mul_left _

/-- The residues of a pair. -/
def resQ (M : ℕ) (q : ℤ × ℤ) : ℕ × ℕ := (res M q.1, res M q.2)

/-- `gOk` on residues, in natural numbers. -/
def gOkN (M G d : ℕ) (q : ℕ × ℕ) : Bool := G.testBit (M * (d * q.1 % M) + d * q.2 % M)

theorem gOkN_eq {M : ℕ} (hM : 0 < M) (G d : ℕ) (q : ℤ × ℤ) : gOkN M G d (resQ M q) = gOk M G d q := by
  unfold gOkN gOk gIdx resQ
  rw [res_mul hM, res_mul hM]

/-- For every `d` of the list, `d` times one of the products is not a Gaussian square. -/
def dsOk (M G : ℕ) (ds : List ℕ) (qs : List (ℕ × ℕ)) : Bool :=
  ds.all fun d => !(qs.all fun q => gOkN M G d q)

/-- At the point `(u, v)` modulo `M`, no `d` of the list survives: the four forms are not all
squares, or for each `d` one of the eight products times `d` is not a Gaussian square. -/
def pointOk (m n : ℤ) (M S G : ℕ) (ds : List ℕ) (u v : ℤ) : Bool :=
  !formsOk m n M S u v || dsOk M G ds (flags.map fun f => resQ M (Q m n u v f))

/-- `kills m n M p S G ds`: the check `pointOk` at every point of the projective line modulo `M`,
a power of `p`: the points `(1, v)`, and the points `(u, 1)` with `p ∣ u`. -/
def kills (m n : ℤ) (M p S G : ℕ) (ds : List ℕ) : Bool :=
  ((List.range M).all fun v => pointOk m n M S G ds 1 v) &&
    ((List.range M).all fun u => !(u % p == 0) || pointOk m n M S G ds u 1)

theorem nrm_congr {M : ℕ} {m n u v u' v' : ℤ} (hu : u ≡ u' [ZMOD M]) (hv : v ≡ v' [ZMOD M]) :
    nrm (P0 m n u v) ≡ nrm (P0 m n u' v') [ZMOD M] ∧ nrm (P8 m n u v) ≡ nrm (P8 m n u' v') [ZMOD M] ∧
    nrm (PL m n u v) ≡ nrm (PL m n u' v') [ZMOD M] ∧ nrm (PU m n u v) ≡ nrm (PU m n u' v') [ZMOD M] := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> dsimp only [nrm, P0, P8, PL, PU] <;> gcongr

theorem Q_congr {M : ℕ} {m n u v u' v' : ℤ} (hu : u ≡ u' [ZMOD M]) (hv : v ≡ v' [ZMOD M])
    (f : Bool × Bool × Bool) :
    (Q m n u v f).1 ≡ (Q m n u' v' f).1 [ZMOD M] ∧ (Q m n u v f).2 ≡ (Q m n u' v' f).2 [ZMOD M] := by
  obtain ⟨a, b, c⟩ := f
  cases a <;> cases b <;> cases c <;> constructor <;>
    dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] <;> simp only [if_true, if_false, Bool.false_eq_true] <;>
    gcongr

theorem gOk_of_GSqP {M : ℕ} (hM : 0 < M) {d : ℕ} {p p' : ℤ × ℤ} (h : GSqP d p)
    (h1 : p'.1 ≡ p.1 [ZMOD M]) (h2 : p'.2 ≡ p.2 [ZMOD M]) : gOk M (gMask M) d p' = true := by
  obtain ⟨a, b, ha, hb⟩ := h
  unfold gOk
  apply gOk_of hM (a := a) (b := b)
  · rw [← ha]
    exact h1.mul_left _
  · rw [← hb]
    exact h2.mul_left _

/-- A solution, read modulo `M`, is a point where `pointOk` fails for its class `d`. -/
theorem pointOk_false {m n : ℤ} {M : ℕ} {ds : List ℕ} (hM : 0 < M) {d : ℕ} (hd : d ∈ ds) {u v : ℤ}
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) {u' v' : ℤ}
    (cu : u' ≡ u [ZMOD M]) (cv : v' ≡ v [ZMOD M]) :
    pointOk m n M (sMask M) (gMask M) ds u' v' = false := by
  obtain ⟨c0, c8, cL, cU⟩ := nrm_congr (m := m) (n := n) cu cv
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  have f : formsOk m n M (sMask M) u' v' = true := by
    unfold formsOk
    rw [sqOk_of hM (w := w0) (hw0 ▸ c0), sqOk_of hM (w := w8) (hw8 ▸ c8),
      sqOk_of hM (w := wL) (hwL ▸ cL), sqOk_of hM (w := wU) (hwU ▸ cU)]
    rfl
  have o : ((flags.map fun f => resQ M (Q m n u' v' f)).all fun q => gOkN M (gMask M) d q) = true := by
    rw [List.all_eq_true]
    intro q hq
    rw [List.mem_map] at hq
    obtain ⟨f, _, rfl⟩ := hq
    obtain ⟨q1, q2⟩ := Q_congr (m := m) (n := n) cu cv f
    rw [gOkN_eq hM]
    exact gOk_of_GSqP hM (g f) q1 q2
  unfold pointOk
  rw [f, Bool.not_true, Bool.false_or]
  rw [Bool.eq_false_iff]
  intro hall
  unfold dsOk at hall
  rw [List.all_eq_true] at hall
  have := hall d hd
  rw [o] at this
  simp at this

theorem scale_sq {m n u v : ℤ} (l : ℤ) (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2)
    (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2) (hL : ∃ w, nrm (PL m n u v) = w ^ 2)
    (hU : ∃ w, nrm (PU m n u v) = w ^ 2) :
    (∃ w, nrm (P0 m n (l * u) (l * v)) = w ^ 2) ∧ (∃ w, nrm (P8 m n (l * u) (l * v)) = w ^ 2) ∧
    (∃ w, nrm (PL m n (l * u) (l * v)) = w ^ 2) ∧ (∃ w, nrm (PU m n (l * u) (l * v)) = w ^ 2) := by
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  dsimp only [nrm, P0, P8, PL, PU] at *
  exact ⟨⟨l * w0, by linear_combination l ^ 2 * hw0⟩, ⟨l * w8, by linear_combination l ^ 2 * hw8⟩,
    ⟨l * wL, by linear_combination l ^ 2 * hwL⟩, ⟨l * wU, by linear_combination l ^ 2 * hwU⟩⟩

theorem scale_gsq {m n u v : ℤ} {d : ℤ} (l : ℤ) (f : Bool × Bool × Bool) (h : GSqP d (Q m n u v f)) :
    GSqP d (Q m n (l * u) (l * v) f) := by
  obtain ⟨a, b, ha, hb⟩ := h
  refine ⟨l ^ 2 * a, l ^ 2 * b, ?_, ?_⟩
  · have e : (Q m n (l * u) (l * v) f).1 = l ^ 4 * (Q m n u v f).1 := by
      obtain ⟨x, y, z⟩ := f
      cases x <;> cases y <;> cases z <;> dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] <;>
        simp only [if_true, if_false, Bool.false_eq_true] <;> ring
    rw [e]
    linear_combination l ^ 4 * ha
  · have e : (Q m n (l * u) (l * v) f).2 = l ^ 4 * (Q m n u v f).2 := by
      obtain ⟨x, y, z⟩ := f
      cases x <;> cases y <;> cases z <;> dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] <;>
        simp only [if_true, if_false, Bool.false_eq_true] <;> ring
    rw [e]
    linear_combination l ^ 4 * hb

theorem exists_inv {p k : ℕ} (hp : p.Prime) {x : ℤ} (hx : ¬(p : ℤ) ∣ x) :
    ∃ l : ℤ, l * x ≡ 1 [ZMOD ((p ^ k : ℕ) : ℤ)] := by
  have h1 : IsCoprime x (p : ℤ) := by
    rw [Int.isCoprime_iff_gcd_eq_one]
    have h2 : ¬p ∣ x.natAbs := fun h => hx (Int.natCast_dvd.mpr h)
    have h3 : Nat.Coprime p x.natAbs := (Nat.Prime.coprime_iff_not_dvd hp).mpr h2
    have h4 : Int.gcd x (p : ℤ) = Nat.gcd x.natAbs p := rfl
    rw [h4]
    exact Nat.Coprime.symm h3
  obtain ⟨a, b, hab⟩ := (h1.pow_right (n := k))
  refine ⟨a, ?_⟩
  rw [Int.modEq_iff_dvd]
  refine ⟨b, ?_⟩
  push_cast
  linear_combination -hab

theorem res_one {M : ℕ} (hM : 1 < M) {x : ℤ} (h : x ≡ 1 [ZMOD M]) : res M x = 1 := by
  rw [res_congr h]
  unfold res
  have h1 : (1 : ℤ) % (M : ℤ) = 1 := Int.emod_eq_of_lt (by norm_num) (by exact_mod_cast hM)
  rw [h1]
  rfl

/-- If the check succeeds, no `d` of the list is the class of a solution. -/
theorem not_gsq_of_kills {m n : ℤ} {M p k S G : ℕ} {ds : List ℕ} (hp : p.Prime) (hpM : M = p ^ k)
    (hk0 : 0 < k) (hS : S = sMask M) (hG : G = gMask M) (hk : kills m n M p S G ds = true)
    {d : ℕ} (hd : d ∈ ds) {u v : ℤ} (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : False := by
  subst hS hG
  have hM1 : 1 < M := by
    rw [hpM]
    exact Nat.one_lt_pow hk0.ne' hp.one_lt
  have hM : 0 < M := by omega
  have hpi : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  unfold kills at hk
  rw [Bool.and_eq_true, List.all_eq_true, List.all_eq_true] at hk
  obtain ⟨hk1, hk2⟩ := hk
  by_cases hu : (p : ℤ) ∣ u
  · -- then `p` does not divide `v`
    have hv : ¬(p : ℤ) ∣ v := fun hv => hpi.not_isUnit (hcop.isUnit_of_dvd' hu hv)
    obtain ⟨l, hl⟩ := exists_inv (k := k) hp hv
    rw [← hpM] at hl
    obtain ⟨s0, s8, sL, sU⟩ := scale_sq l h0 h8 hL hU
    have key := pointOk_false (ds := ds) hM hd s0 s8 sL sU (fun f => scale_gsq l f (g f))
      (res_cast hM (l * u)) hl.symm
    have h3 := hk2 (res M (l * u)) (List.mem_range.mpr (res_lt hM _))
    have h4 : res M (l * u) % p = 0 := by
      have h5 : (p : ℤ) ∣ ((res M (l * u) : ℕ) : ℤ) := by
        have h6 : (M : ℤ) ∣ l * u - (res M (l * u) : ℕ) := Int.ModEq.dvd (res_cast hM (l * u))
        have h7 : (p : ℤ) ∣ (M : ℤ) := by
          rw [hpM]
          push_cast
          exact dvd_pow_self _ hk0.ne'
        have h8 := dvd_sub (Dvd.dvd.mul_left hu l) (dvd_trans h7 h6)
        rwa [sub_sub_cancel] at h8
      exact Nat.mod_eq_zero_of_dvd (by exact_mod_cast h5)
    rw [h4] at h3
    simp only [beq_self_eq_true, Bool.not_true, Bool.false_or] at h3
    have key' : pointOk m n M (sMask M) (gMask M) ds ((res M (l * u) : ℕ) : ℤ) 1 = false := key
    rw [key'] at h3
    exact Bool.false_ne_true h3
  · obtain ⟨l, hl⟩ := exists_inv (k := k) hp hu
    rw [← hpM] at hl
    obtain ⟨s0, s8, sL, sU⟩ := scale_sq l h0 h8 hL hU
    have key := pointOk_false (ds := ds) hM hd s0 s8 sL sU (fun f => scale_gsq l f (g f))
      hl.symm (res_cast hM (l * v))
    have h3 := hk1 (res M (l * v)) (List.mem_range.mpr (res_lt hM _))
    rw [key] at h3
    exact Bool.false_ne_true h3

/-! ### The class is even -/

theorem nrm_congr4 {M : ℕ} {m n u v m' n' u' v' : ℤ} (hm : m ≡ m' [ZMOD M]) (hn : n ≡ n' [ZMOD M])
    (hu : u ≡ u' [ZMOD M]) (hv : v ≡ v' [ZMOD M]) :
    nrm (P0 m n u v) ≡ nrm (P0 m' n' u' v') [ZMOD M] ∧ nrm (P8 m n u v) ≡ nrm (P8 m' n' u' v') [ZMOD M] ∧
    nrm (PL m n u v) ≡ nrm (PL m' n' u' v') [ZMOD M] ∧ nrm (PU m n u v) ≡ nrm (PU m' n' u' v') [ZMOD M] := by
  refine ⟨?_, ?_, ?_, ?_⟩ <;> dsimp only [nrm, P0, P8, PL, PU] <;> gcongr

theorem Q_congr4 {M : ℕ} {m n u v m' n' u' v' : ℤ} (hm : m ≡ m' [ZMOD M]) (hn : n ≡ n' [ZMOD M])
    (hu : u ≡ u' [ZMOD M]) (hv : v ≡ v' [ZMOD M]) (f : Bool × Bool × Bool) :
    (Q m n u v f).1 ≡ (Q m' n' u' v' f).1 [ZMOD M] ∧ (Q m n u v f).2 ≡ (Q m' n' u' v' f).2 [ZMOD M] := by
  obtain ⟨a, b, c⟩ := f
  cases a <;> cases b <;> cases c <;> constructor <;>
    dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] <;> simp only [if_true, if_false, Bool.false_eq_true] <;>
    gcongr

theorem gOkN_mod (M G d : ℕ) (q : ℕ × ℕ) : gOkN M G (d % M) q = gOkN M G d q := by
  unfold gOkN
  rw [Nat.mod_mul_mod, Nat.mod_mul_mod]

/-- The check behind `two_dvd_class`: modulo 8, at every `(m, n, u, v)` with `m, n` not both even
and `u, v` not both even, the four forms are not all squares, or no odd `d` passes. -/
def evenCheck : Bool :=
  (List.range 8).all fun m => (List.range 8).all fun n => (List.range 8).all fun u =>
    (List.range 8).all fun v =>
      (m % 2 == 0 && n % 2 == 0) || (u % 2 == 0 && v % 2 == 0) ||
        !formsOk m n 8 (sMask 8) u v ||
        dsOk 8 (gMask 8) [1, 3, 5, 7] (flags.map fun f => resQ 8 (Q m n u v f))

theorem evenCheck_eq : evenCheck = true := by decide +kernel

/-- **The class is even.** For coprime `m, n` and coprime `u, v` at which the four forms are
squares, a `d` such that `d` times each of the eight products is a Gaussian square is even. No
condition on the ratio. -/
theorem two_dvd_class {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    {d : ℕ} (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : 2 ∣ d := by
  by_contra hd
  have hM : 0 < 8 := by norm_num
  have cm := res_cast hM m
  have cn := res_cast hM n
  have cu := res_cast hM u
  have cv := res_cast hM v
  have h := evenCheck_eq
  unfold evenCheck at h
  rw [List.all_eq_true] at h
  have h1 := h (res 8 m) (List.mem_range.mpr (res_lt hM m))
  rw [List.all_eq_true] at h1
  have h2 := h1 (res 8 n) (List.mem_range.mpr (res_lt hM n))
  rw [List.all_eq_true] at h2
  have h3 := h2 (res 8 u) (List.mem_range.mpr (res_lt hM u))
  rw [List.all_eq_true] at h3
  have h4 := h3 (res 8 v) (List.mem_range.mpr (res_lt hM v))
  obtain ⟨c0, c8, cL, cU⟩ := nrm_congr4 cm cn cu cv
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  have f : formsOk (res 8 m) (res 8 n) 8 (sMask 8) (res 8 u) (res 8 v) = true := by
    unfold formsOk
    rw [sqOk_of hM (w := w0) (hw0 ▸ c0), sqOk_of hM (w := w8) (hw8 ▸ c8),
      sqOk_of hM (w := wL) (hwL ▸ cL), sqOk_of hM (w := wU) (hwU ▸ cU)]
    rfl
  have hmem : d % 8 ∈ [1, 3, 5, 7] := by
    have : d % 8 = 1 ∨ d % 8 = 3 ∨ d % 8 = 5 ∨ d % 8 = 7 := by omega
    rcases this with h | h | h | h <;> simp [h]
  have o : ((flags.map fun f => resQ 8 (Q (res 8 m) (res 8 n) (res 8 u) (res 8 v) f)).all
      fun q => gOkN 8 (gMask 8) (d % 8) q) = true := by
    rw [List.all_eq_true]
    intro q hq
    rw [List.mem_map] at hq
    obtain ⟨f, _, rfl⟩ := hq
    obtain ⟨q1, q2⟩ := Q_congr4 cm cn cu cv f
    rw [gOkN_mod, gOkN_eq hM]
    exact gOk_of_GSqP hM (g f) q1 q2
  have nd : dsOk 8 (gMask 8) [1, 3, 5, 7]
      (flags.map fun f => resQ 8 (Q (res 8 m) (res 8 n) (res 8 u) (res 8 v) f)) = false := by
    rw [Bool.eq_false_iff]
    intro hall
    unfold dsOk at hall
    rw [List.all_eq_true] at hall
    have := hall (d % 8) hmem
    rw [o] at this
    simp at this
  have e8 : ∀ {x : ℤ}, ((res 8 x : ℕ) : ℤ) ≡ x [ZMOD ((8 : ℕ) : ℤ)] → res 8 x % 2 = 0 → (2 : ℤ) ∣ x := by
    intro x hx h2
    have := Int.ModEq.dvd hx
    omega
  have bA : (res 8 m % 2 == 0 && res 8 n % 2 == 0) = false := by
    rw [Bool.eq_false_iff]
    intro h
    simp only [Bool.and_eq_true, beq_iff_eq] at h
    exact Int.prime_two.not_isUnit (hmn.isUnit_of_dvd' (e8 cm h.1) (e8 cn h.2))
  have bB : (res 8 u % 2 == 0 && res 8 v % 2 == 0) = false := by
    rw [Bool.eq_false_iff]
    intro h
    simp only [Bool.and_eq_true, beq_iff_eq] at h
    exact Int.prime_two.not_isUnit (hcop.isUnit_of_dvd' (e8 cu h.1) (e8 cv h.2))
  rw [bA, bB, f, nd] at h4
  exact absurd h4 (by decide)

/-! ### The class modulo the primes of `m² + n²` -/

theorem isSquare_of_mul_sq {F : Type*} [Field F] {d s t : F} (hs : s ≠ 0) (h : d * s ^ 2 = t ^ 2) :
    IsSquare d :=
  ⟨t / s, by field_simp; linear_combination h⟩

/-- **The class is a square modulo every odd prime dividing `m² + n²`.** For coprime `m, n` and
coprime `u, v`, a `d` such that `d` times each of the eight products is a Gaussian square is a
square modulo `p`, for every odd prime `p` dividing `m² + n²`. No condition on the ratio, and the
four forms need not be squares. -/
theorem isSquare_class_mod {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) (hdvd : (p : ℤ) ∣ m ^ 2 + n ^ 2)
    {d : ℕ} (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : IsSquare (d : ZMod p) := by
  have : Fact p.Prime := ⟨hp⟩
  have hS : (m : ZMod p) ^ 2 + (n : ZMod p) ^ 2 = 0 := by
    have h1 := (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).mpr hdvd
    push_cast at h1
    exact h1
  have hmn' : ¬((m : ZMod p) = 0 ∧ (n : ZMod p) = 0) := by
    rintro ⟨hm, hn⟩
    obtain ⟨a, b, hab⟩ := hmn
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hm, hn, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have hcop' : ¬((u : ZMod p) = 0 ∧ (v : ZMod p) = 0) := by
    rintro ⟨hu, hv⟩
    obtain ⟨a, b, hab⟩ := hcop
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hu, hv, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have hn : (n : ZMod p) ≠ 0 := by
    intro h0
    apply hmn'
    refine ⟨?_, h0⟩
    rw [h0] at hS
    have : (m : ZMod p) ^ 2 = 0 := by linear_combination hS
    exact pow_eq_zero_iff (two_ne_zero) |>.mp this
  have two_ne : (2 : ZMod p) ≠ 0 := by
    intro h2
    have h3 : ((2 : ℕ) : ZMod p) = 0 := by exact_mod_cast h2
    rw [ZMod.natCast_eq_zero_iff] at h3
    have := (Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp h3
    exact hp2 this
  obtain ⟨ι, hι, hm⟩ : ∃ ι : ZMod p, ι ^ 2 = -1 ∧ (m : ZMod p) = ι * n := by
    refine ⟨m / n, ?_, (div_mul_cancel₀ (m : ZMod p) hn).symm⟩
    rw [div_pow, div_eq_iff (pow_ne_zero 2 hn)]
    linear_combination hS
  have r8 : (-8 * ι : ZMod p) = (2 * (1 - ι)) ^ 2 := by linear_combination (-4) * hι
  have hι1 : (1 - ι : ZMod p) ≠ 0 := by
    intro h1
    have h2 : ι = 1 := by linear_combination -h1
    rw [h2] at hι
    apply two_ne
    linear_combination hι
  have hpm : (u : ZMod p) + v ≠ 0 ∨ (u : ZMod p) - v ≠ 0 := by
    by_contra hc
    rw [not_or, not_not, not_not] at hc
    obtain ⟨c1, c2⟩ := hc
    apply hcop'
    have hu2 : (2 : ZMod p) * u = 0 := by linear_combination c1 + c2
    have hv2 : (2 : ZMod p) * v = 0 := by linear_combination c1 - c2
    exact ⟨(mul_eq_zero.mp hu2).resolve_left two_ne, (mul_eq_zero.mp hv2).resolve_left two_ne⟩
  have huv : (u : ZMod p) ≠ 0 ∨ (v : ZMod p) ≠ 0 := by
    by_contra hc
    rw [not_or, not_not, not_not] at hc
    exact hcop' hc
  -- the two products used, read in `ZMod p`
  obtain ⟨a1, b1, ha1, hb1⟩ := g (true, false, true)
  obtain ⟨a2, b2, ha2, hb2⟩ := g (true, true, false)
  have ha1' := congrArg (Int.cast : ℤ → ZMod p) ha1
  have hb1' := congrArg (Int.cast : ℤ → ZMod p) hb1
  have ha2' := congrArg (Int.cast : ℤ → ZMod p) ha2
  have hb2' := congrArg (Int.cast : ℤ → ZMod p) hb2
  dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] at ha1' hb1' ha2' hb2'
  simp only [if_true, if_false, Bool.false_eq_true] at ha1' hb1' ha2' hb2'
  push_cast at ha1' hb1' ha2' hb2'
  rw [hm] at ha1' hb1' ha2' hb2'
  set D : ZMod p := (d : ZMod p) with hD
  set N : ZMod p := (n : ZMod p)
  set U : ZMod p := (u : ZMod p)
  set V : ZMod p := (v : ZMod p)
  rcases hpm with hp1 | hp1 <;> rcases huv with hq1 | hq1
  · -- u + v and u
    refine isSquare_of_mul_sq (s := 2 * (1 - ι) * N ^ 2 * (U + V) * U) (t := (a2 : ZMod p) + ι * b2)
      (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne hι1) (pow_ne_zero 2 hn)) hp1) hq1) ?_
    have key : D * (-8 * ι * N ^ 4 * (U + V) ^ 2 * U ^ 2) = ((a2 : ZMod p) + ι * b2) ^ 2 := by
      linear_combination ha2' + ι * hb2' + (-(D * (6*N^4*U^4*ι - 2*N^4*U^3*V*ι^3 + 10*N^4*U^3*V*ι + 4*N^4*U^2*V^2*ι - 2*N^4*U*V^3*ι^3 + 2*N^4*U*V^3*ι - 2*N^4*V^4*ι)) - (b2 : ZMod p) ^ 2) * hι
    linear_combination key - (D * N ^ 4 * (U + V) ^ 2 * U ^ 2) * r8
  · -- u + v and v
    refine isSquare_of_mul_sq (s := 2 * (1 - ι) * N ^ 2 * (U + V) * V) (t := (a1 : ZMod p) + ι * b1)
      (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne hι1) (pow_ne_zero 2 hn)) hp1) hq1) ?_
    have key : D * (-8 * ι * N ^ 4 * (U + V) ^ 2 * V ^ 2) = ((a1 : ZMod p) + ι * b1) ^ 2 := by
      linear_combination ha1' + ι * hb1' + (-(D * (-2*N^4*U^4*ι - 2*N^4*U^3*V*ι^3 + 2*N^4*U^3*V*ι + 4*N^4*U^2*V^2*ι - 2*N^4*U*V^3*ι^3 + 10*N^4*U*V^3*ι + 6*N^4*V^4*ι)) - (b1 : ZMod p) ^ 2) * hι
    linear_combination key - (D * N ^ 4 * (U + V) ^ 2 * V ^ 2) * r8
  · -- u - v and u
    refine isSquare_of_mul_sq (s := 2 * (1 - ι) * N ^ 2 * (U - V) * U) (t := (a1 : ZMod p) - ι * b1)
      (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne hι1) (pow_ne_zero 2 hn)) hp1) hq1) ?_
    have key : D * (-8 * ι * N ^ 4 * (U - V) ^ 2 * U ^ 2) = ((a1 : ZMod p) - ι * b1) ^ 2 := by
      linear_combination ha1' - ι * hb1' + (-(D * (6*N^4*U^4*ι + 2*N^4*U^3*V*ι^3 - 10*N^4*U^3*V*ι + 4*N^4*U^2*V^2*ι + 2*N^4*U*V^3*ι^3 - 2*N^4*U*V^3*ι - 2*N^4*V^4*ι)) - (b1 : ZMod p) ^ 2) * hι
    linear_combination key - (D * N ^ 4 * (U - V) ^ 2 * U ^ 2) * r8
  · -- u - v and v
    refine isSquare_of_mul_sq (s := 2 * (1 - ι) * N ^ 2 * (U - V) * V) (t := (a2 : ZMod p) - ι * b2)
      (mul_ne_zero (mul_ne_zero (mul_ne_zero (mul_ne_zero two_ne hι1) (pow_ne_zero 2 hn)) hp1) hq1) ?_
    have key : D * (-8 * ι * N ^ 4 * (U - V) ^ 2 * V ^ 2) = ((a2 : ZMod p) - ι * b2) ^ 2 := by
      linear_combination ha2' - ι * hb2' + (-(D * (-2*N^4*U^4*ι + 2*N^4*U^3*V*ι^3 - 2*N^4*U^3*V*ι + 4*N^4*U^2*V^2*ι + 2*N^4*U*V^3*ι^3 - 10*N^4*U*V^3*ι + 6*N^4*V^4*ι)) - (b2 : ZMod p) ^ 2) * hι
    linear_combination key - (D * N ^ 4 * (U - V) ^ 2 * V ^ 2) * r8

/-! ### The class modulo the primes of `m⁸ + 14m⁴n⁴ + n⁸` -/

/-! `m⁸ + 14m⁴n⁴ + n⁸ = R₁·R₂` with `R₂ = m⁴ + 2m³n + 2m²n² − 2mn³ + n⁴` and `R₁` the same with
`−n`. With `A = m² + n²` and `J = m² + 2mn − n²`: `J² + A² = 2R₂` and `R₂ = (m² − n²)² + 2mnJ`.
Modulo a prime `p` of `R₂` the four linear forms are proportional in pairs in each of the two
residue fields above `p`, and one of the eight products, `X + iY`, satisfies
`d·((m² − n²)·G₀·G₈)² = ((Aa + Jb)·mn·A)²` and `d·(J(m² − n²)·H₀·H₈)² = (A(Aa − Jb)·mn·A)²`
modulo `p`, where `G₀ = A·mv + J·nu`, `G₈ = A·nv + J·mu` and `H₀`, `H₈` are the same with `−J`.
`G₀G₈` and `H₀H₈` do not both vanish modulo `p` when `u, v` are coprime. -/

theorem isSquare_class_mod_R2 {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hdvd : (p : ℤ) ∣ m ^ 4 + 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 - 2 * m * n ^ 3 + n ^ 4)
    {d : ℕ} (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : IsSquare (d : ZMod p) := by
  have : Fact p.Prime := ⟨hp⟩
  have hmn' : ¬((m : ZMod p) = 0 ∧ (n : ZMod p) = 0) := by
    rintro ⟨hm, hn⟩
    obtain ⟨a, b, hab⟩ := hmn
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hm, hn, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have hcop' : ¬((u : ZMod p) = 0 ∧ (v : ZMod p) = 0) := by
    rintro ⟨hu, hv⟩
    obtain ⟨a, b, hab⟩ := hcop
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hu, hv, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have two_ne : (2 : ZMod p) ≠ 0 := by
    intro h2
    have h3 : ((2 : ℕ) : ZMod p) = 0 := by exact_mod_cast h2
    rw [ZMod.natCast_eq_zero_iff] at h3
    exact hp2 ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp h3)
  have hR := (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).mpr hdvd
  push_cast at hR
  obtain ⟨a, b, ha, hb⟩ := g (false, false, false)
  have ha' := congrArg (Int.cast : ℤ → ZMod p) ha
  have hb' := congrArg (Int.cast : ℤ → ZMod p) hb
  dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] at ha' hb'
  simp only [if_false, Bool.false_eq_true] at ha' hb'
  push_cast at ha' hb'
  set D : ZMod p := (d : ZMod p)
  set M : ZMod p := (m : ZMod p)
  set N : ZMod p := (n : ZMod p)
  set U : ZMod p := (u : ZMod p)
  set V : ZMod p := (v : ZMod p)
  have hN : N ≠ 0 := by
    intro h
    apply hmn'
    refine ⟨?_, h⟩
    have h4 : M ^ 4 = 0 := by
      rw [h] at hR
      linear_combination hR
    exact pow_eq_zero_iff (by norm_num) |>.mp h4
  have hM : M ≠ 0 := by
    intro h
    apply hmn'
    refine ⟨h, ?_⟩
    have h4 : N ^ 4 = 0 := by
      rw [h] at hR
      linear_combination hR
    exact pow_eq_zero_iff (by norm_num) |>.mp h4
  have hE : M ^ 2 - N ^ 2 ≠ 0 := by
    intro h
    have h4 : (2 * M * N) ^ 2 = 0 := by linear_combination hR - (M ^ 2 - N ^ 2 + 2 * M * N) * h
    exact mul_ne_zero (mul_ne_zero two_ne hM) hN (pow_eq_zero_iff (by norm_num) |>.mp h4)
  have hA : M ^ 2 + N ^ 2 ≠ 0 := by
    intro h
    have hJ2 : (M ^ 2 + 2 * M * N - N ^ 2) ^ 2 = 0 := by linear_combination 2 * hR - (M ^ 2 + N ^ 2) * h
    have hJ := pow_eq_zero_iff (by norm_num) |>.mp hJ2
    apply hE
    have h6 : 2 * N * (M - N) = 0 := by linear_combination hJ - h
    have h7 := (mul_eq_zero.mp h6).resolve_left (mul_ne_zero two_ne hN)
    have h8' : M = N := by linear_combination h7
    rw [h8']
    ring
  have hJ : M ^ 2 + 2 * M * N - N ^ 2 ≠ 0 := by
    intro h
    have h2 : (M ^ 2 + N ^ 2) ^ 2 = 0 := by linear_combination 2 * hR - (M ^ 2 + 2 * M * N - N ^ 2) * h
    exact hA (pow_eq_zero_iff (by norm_num) |>.mp h2)
  by_cases hz : ((M ^ 2 + N ^ 2) * (M * V) + (M ^ 2 + 2 * M * N - N ^ 2) * (N * U)) * ((M ^ 2 + N ^ 2) * (N * V) + (M ^ 2 + 2 * M * N - N ^ 2) * (M * U)) = 0
  · -- the point is at a common zero of two of the forms: the other residue field
    have hAAJ : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 + 2 * M * N - N ^ 2) ≠ 0 := mul_ne_zero (pow_ne_zero 2 hA) hJ
    have hH : ((M ^ 2 + N ^ 2) * (M * V) - (M ^ 2 + 2 * M * N - N ^ 2) * (N * U)) * ((M ^ 2 + N ^ 2) * (N * V) - (M ^ 2 + 2 * M * N - N ^ 2) * (M * U)) ≠ 0 := by
      intro hH
      apply hcop'
      rcases mul_eq_zero.mp hz with h1 | h1 <;> rcases mul_eq_zero.mp hH with h2 | h2
      · have eU : (M ^ 2 + 2 * M * N - N ^ 2) * N * U * 2 = 0 := by linear_combination h1 - h2
        have eV : (M ^ 2 + N ^ 2) * M * V * 2 = 0 := by linear_combination h1 + h2
        exact ⟨by simpa [hJ, hN, two_ne] using eU, by simpa [hA, hM, two_ne] using eV⟩
      · have eU : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 + 2 * M * N - N ^ 2) * U = 0 := by linear_combination ((M ^ 2 + N ^ 2) * N) * h1 - ((M ^ 2 + N ^ 2) * M) * h2
        have eV : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 + 2 * M * N - N ^ 2) * V = 0 := by linear_combination ((M ^ 2 + 2 * M * N - N ^ 2) * M) * h1 + ((M ^ 2 + 2 * M * N - N ^ 2) * N) * h2
        exact ⟨(mul_eq_zero.mp eU).resolve_left hAAJ, (mul_eq_zero.mp eV).resolve_left hAAJ⟩
      · have eU : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 + 2 * M * N - N ^ 2) * U = 0 := by linear_combination ((M ^ 2 + N ^ 2) * M) * h1 - ((M ^ 2 + N ^ 2) * N) * h2
        have eV : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 + 2 * M * N - N ^ 2) * V = 0 := by linear_combination ((M ^ 2 + 2 * M * N - N ^ 2) * N) * h1 + ((M ^ 2 + 2 * M * N - N ^ 2) * M) * h2
        exact ⟨(mul_eq_zero.mp eU).resolve_left hAAJ, (mul_eq_zero.mp eV).resolve_left hAAJ⟩
      · have eU : (M ^ 2 + 2 * M * N - N ^ 2) * M * U * 2 = 0 := by linear_combination h1 - h2
        have eV : (M ^ 2 + N ^ 2) * N * V * 2 = 0 := by linear_combination h1 + h2
        exact ⟨by simpa [hJ, hM, two_ne] using eU, by simpa [hA, hN, two_ne] using eV⟩
    have key : D * ((M ^ 2 + 2 * M * N - N ^ 2) * (M ^ 2 - N ^ 2) * (((M ^ 2 + N ^ 2) * (M * V) - (M ^ 2 + 2 * M * N - N ^ 2) * (N * U)) *
        ((M ^ 2 + N ^ 2) * (N * V) - (M ^ 2 + 2 * M * N - N ^ 2) * (M * U)))) ^ 2 =
        ((M ^ 2 + N ^ 2) * ((M ^ 2 + N ^ 2) * (a : ZMod p) - (M ^ 2 + 2 * M * N - N ^ 2) * (b : ZMod p)) * M * N * (M ^ 2 + N ^ 2)) ^ 2 := by
      linear_combination ((M ^ 2 + N ^ 2) ^ 6 * M ^ 2 * N ^ 2) * ha' -
        ((M ^ 2 + N ^ 2) ^ 5 * M ^ 2 * N ^ 2 * (M ^ 2 + 2 * M * N - N ^ 2)) * hb' +
        (
          D*M^16*U^2*V^2 - 2*D*M^15*N*U^3*V + 6*D*M^15*N*U^2*V^2 - 2*D*M^15*N*U*V^3 +
          D*M^14*N^2*U^4 - 14*D*M^14*N^2*U^3*V + 6*D*M^14*N^2*U^2*V^2 - 10*D*M^14*N^2*U*V^3 +
          D*M^14*N^2*V^4 + 8*D*M^13*N^3*U^4 - 42*D*M^13*N^3*U^3*V + 14*D*M^13*N^3*U^2*V^2 +
          6*D*M^13*N^3*U*V^3 + 30*D*M^12*N^4*U^4 + 44*D*M^12*N^4*U^3*V + 16*D*M^12*N^4*U^2*V^2
          - 28*D*M^12*N^4*U*V^3 - 2*D*M^12*N^4*V^4 - 8*D*M^11*N^5*U^4 + 14*D*M^11*N^5*U^3*V -
          58*D*M^11*N^5*U^2*V^2 + 30*D*M^11*N^5*U*V^3 - 145*D*M^10*N^6*U^4 +
          46*D*M^10*N^6*U^3*V + 26*D*M^10*N^6*U^2*V^2 - 22*D*M^10*N^6*U*V^3 - 17*D*M^10*N^6*V^4
          - 16*D*M^9*N^7*U^4 + 54*D*M^9*N^7*U^3*V - 66*D*M^9*N^7*U^2*V^2 + 22*D*M^9*N^7*U*V^3 +
          164*D*M^8*N^8*U^4 - 24*D*M^8*N^8*U^3*V + 30*D*M^8*N^8*U^2*V^2 - 8*D*M^8*N^8*U*V^3 -
          28*D*M^8*N^8*V^4 + 16*D*M^7*N^9*U^4 - 54*D*M^7*N^9*U^3*V + 66*D*M^7*N^9*U^2*V^2 -
          22*D*M^7*N^9*U*V^3 - 145*D*M^6*N^10*U^4 + 46*D*M^6*N^10*U^3*V + 26*D*M^6*N^10*U^2*V^2
          - 22*D*M^6*N^10*U*V^3 - 17*D*M^6*N^10*V^4 + 8*D*M^5*N^11*U^4 - 14*D*M^5*N^11*U^3*V +
          58*D*M^5*N^11*U^2*V^2 - 30*D*M^5*N^11*U*V^3 + 30*D*M^4*N^12*U^4 + 44*D*M^4*N^12*U^3*V
          + 16*D*M^4*N^12*U^2*V^2 - 28*D*M^4*N^12*U*V^3 - 2*D*M^4*N^12*V^4 - 8*D*M^3*N^13*U^4 +
          42*D*M^3*N^13*U^3*V - 14*D*M^3*N^13*U^2*V^2 - 6*D*M^3*N^13*U*V^3 + D*M^2*N^14*U^4 -
          14*D*M^2*N^14*U^3*V + 6*D*M^2*N^14*U^2*V^2 - 10*D*M^2*N^14*U*V^3 + D*M^2*N^14*V^4 +
          2*D*M*N^15*U^3*V - 6*D*M*N^15*U^2*V^2 + 2*D*M*N^15*U*V^3 + D*N^16*U^2*V^2 -
          2*M^10*N^2*(b : ZMod p)^2 - 8*M^8*N^4*(b : ZMod p)^2 - 12*M^6*N^6*(b : ZMod p)^2 -
          8*M^4*N^8*(b : ZMod p)^2 - 2*M^2*N^10*(b : ZMod p)^2
        ) * hR
    exact isSquare_of_mul_sq (mul_ne_zero (mul_ne_zero hJ hE) hH) key
  · have key : D * ((M ^ 2 - N ^ 2) * (((M ^ 2 + N ^ 2) * (M * V) + (M ^ 2 + 2 * M * N - N ^ 2) * (N * U)) *
        ((M ^ 2 + N ^ 2) * (N * V) + (M ^ 2 + 2 * M * N - N ^ 2) * (M * U)))) ^ 2 =
        (((M ^ 2 + N ^ 2) * (a : ZMod p) + (M ^ 2 + 2 * M * N - N ^ 2) * (b : ZMod p)) * M * N * (M ^ 2 + N ^ 2)) ^ 2 := by
      linear_combination (M ^ 2 * N ^ 2 * (M ^ 2 + N ^ 2) ^ 4) * ha' +
        (M ^ 2 * N ^ 2 * (M ^ 2 + N ^ 2) ^ 3 * (M ^ 2 + 2 * M * N - N ^ 2)) * hb' +
        (
          D*M^12*U^2*V^2 + 2*D*M^11*N*U^3*V + 2*D*M^11*N*U^2*V^2 + 2*D*M^11*N*U*V^3 +
          D*M^10*N^2*U^4 + 10*D*M^10*N^2*U^3*V + 4*D*M^10*N^2*U^2*V^2 - 2*D*M^10*N^2*U*V^3 +
          D*M^10*N^2*V^4 + 8*D*M^9*N^3*U^4 - 10*D*M^9*N^3*U^3*V + 6*D*M^9*N^3*U^2*V^2 +
          6*D*M^9*N^3*U*V^3 + 4*D*M^8*N^4*U^4 + 8*D*M^8*N^4*U^3*V - 9*D*M^8*N^4*U^2*V^2 -
          8*D*M^8*N^4*U*V^3 + 4*D*M^8*N^4*V^4 - 24*D*M^7*N^5*U^4 - 12*D*M^7*N^5*U^3*V +
          4*D*M^7*N^5*U^2*V^2 + 4*D*M^7*N^5*U*V^3 + 6*D*M^6*N^6*U^4 - 4*D*M^6*N^6*U^3*V -
          24*D*M^6*N^6*U^2*V^2 - 12*D*M^6*N^6*U*V^3 + 6*D*M^6*N^6*V^4 + 24*D*M^5*N^7*U^4 +
          12*D*M^5*N^7*U^3*V - 4*D*M^5*N^7*U^2*V^2 - 4*D*M^5*N^7*U*V^3 + 4*D*M^4*N^8*U^4 +
          8*D*M^4*N^8*U^3*V - 9*D*M^4*N^8*U^2*V^2 - 8*D*M^4*N^8*U*V^3 + 4*D*M^4*N^8*V^4 -
          8*D*M^3*N^9*U^4 + 10*D*M^3*N^9*U^3*V - 6*D*M^3*N^9*U^2*V^2 - 6*D*M^3*N^9*U*V^3 +
          D*M^2*N^10*U^4 + 10*D*M^2*N^10*U^3*V + 4*D*M^2*N^10*U^2*V^2 - 2*D*M^2*N^10*U*V^3 +
          D*M^2*N^10*V^4 - 2*D*M*N^11*U^3*V - 2*D*M*N^11*U^2*V^2 - 2*D*M*N^11*U*V^3 +
          D*N^12*U^2*V^2 - 2*M^6*N^2*(b : ZMod p)^2 - 4*M^4*N^4*(b : ZMod p)^2 - 2*M^2*N^6*(b :
          ZMod p)^2
        ) * hR
    exact isSquare_of_mul_sq (mul_ne_zero hE hz) key

theorem isSquare_class_mod_R1 {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hdvd : (p : ℤ) ∣ m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4)
    {d : ℕ} (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : IsSquare (d : ZMod p) := by
  have : Fact p.Prime := ⟨hp⟩
  have hmn' : ¬((m : ZMod p) = 0 ∧ (n : ZMod p) = 0) := by
    rintro ⟨hm, hn⟩
    obtain ⟨a, b, hab⟩ := hmn
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hm, hn, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have hcop' : ¬((u : ZMod p) = 0 ∧ (v : ZMod p) = 0) := by
    rintro ⟨hu, hv⟩
    obtain ⟨a, b, hab⟩ := hcop
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hab
    push_cast at h1
    rw [hu, hv, mul_zero, mul_zero, add_zero] at h1
    exact zero_ne_one h1
  have two_ne : (2 : ZMod p) ≠ 0 := by
    intro h2
    have h3 : ((2 : ℕ) : ZMod p) = 0 := by exact_mod_cast h2
    rw [ZMod.natCast_eq_zero_iff] at h3
    exact hp2 ((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp h3)
  have hR := (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).mpr hdvd
  push_cast at hR
  obtain ⟨a, b, ha, hb⟩ := g (false, true, true)
  have ha' := congrArg (Int.cast : ℤ → ZMod p) ha
  have hb' := congrArg (Int.cast : ℤ → ZMod p) hb
  dsimp only [Q, gm, cj, cjIf, P0, P8, PL, PU] at ha' hb'
  simp only [if_true, if_false, Bool.false_eq_true] at ha' hb'
  push_cast at ha' hb'
  set D : ZMod p := (d : ZMod p)
  set M : ZMod p := (m : ZMod p)
  set N : ZMod p := (n : ZMod p)
  set U : ZMod p := (u : ZMod p)
  set V : ZMod p := (v : ZMod p)
  have hN : N ≠ 0 := by
    intro h
    apply hmn'
    refine ⟨?_, h⟩
    have h4 : M ^ 4 = 0 := by
      rw [h] at hR
      linear_combination hR
    exact pow_eq_zero_iff (by norm_num) |>.mp h4
  have hM : M ≠ 0 := by
    intro h
    apply hmn'
    refine ⟨h, ?_⟩
    have h4 : N ^ 4 = 0 := by
      rw [h] at hR
      linear_combination hR
    exact pow_eq_zero_iff (by norm_num) |>.mp h4
  have hE : M ^ 2 - N ^ 2 ≠ 0 := by
    intro h
    have h4 : (2 * M * N) ^ 2 = 0 := by linear_combination hR - (M ^ 2 - N ^ 2 - 2 * M * N) * h
    exact mul_ne_zero (mul_ne_zero two_ne hM) hN (pow_eq_zero_iff (by norm_num) |>.mp h4)
  have hA : M ^ 2 + N ^ 2 ≠ 0 := by
    intro h
    have hJ2 : (M ^ 2 - 2 * M * N - N ^ 2) ^ 2 = 0 := by linear_combination 2 * hR - (M ^ 2 + N ^ 2) * h
    have hJ := pow_eq_zero_iff (by norm_num) |>.mp hJ2
    apply hE
    have h6 : 2 * N * (M + N) = 0 := by linear_combination h - hJ
    have h7 := (mul_eq_zero.mp h6).resolve_left (mul_ne_zero two_ne hN)
    have h8' : M = -N := by linear_combination h7
    rw [h8']
    ring
  have hJ : M ^ 2 - 2 * M * N - N ^ 2 ≠ 0 := by
    intro h
    have h2 : (M ^ 2 + N ^ 2) ^ 2 = 0 := by linear_combination 2 * hR - (M ^ 2 - 2 * M * N - N ^ 2) * h
    exact hA (pow_eq_zero_iff (by norm_num) |>.mp h2)
  by_cases hz : ((M ^ 2 + N ^ 2) * (M * V) + (M ^ 2 - 2 * M * N - N ^ 2) * (N * U)) * ((M ^ 2 + N ^ 2) * (N * V) + (M ^ 2 - 2 * M * N - N ^ 2) * (M * U)) = 0
  · -- the point is at a common zero of two of the forms: the other residue field
    have hAAJ : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 - 2 * M * N - N ^ 2) ≠ 0 := mul_ne_zero (pow_ne_zero 2 hA) hJ
    have hH : ((M ^ 2 + N ^ 2) * (M * V) - (M ^ 2 - 2 * M * N - N ^ 2) * (N * U)) * ((M ^ 2 + N ^ 2) * (N * V) - (M ^ 2 - 2 * M * N - N ^ 2) * (M * U)) ≠ 0 := by
      intro hH
      apply hcop'
      rcases mul_eq_zero.mp hz with h1 | h1 <;> rcases mul_eq_zero.mp hH with h2 | h2
      · have eU : (M ^ 2 - 2 * M * N - N ^ 2) * N * U * 2 = 0 := by linear_combination h1 - h2
        have eV : (M ^ 2 + N ^ 2) * M * V * 2 = 0 := by linear_combination h1 + h2
        exact ⟨by simpa [hJ, hN, two_ne] using eU, by simpa [hA, hM, two_ne] using eV⟩
      · have eU : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 - 2 * M * N - N ^ 2) * U = 0 := by linear_combination ((M ^ 2 + N ^ 2) * N) * h1 - ((M ^ 2 + N ^ 2) * M) * h2
        have eV : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 - 2 * M * N - N ^ 2) * V = 0 := by linear_combination ((M ^ 2 - 2 * M * N - N ^ 2) * M) * h1 + ((M ^ 2 - 2 * M * N - N ^ 2) * N) * h2
        exact ⟨(mul_eq_zero.mp eU).resolve_left hAAJ, (mul_eq_zero.mp eV).resolve_left hAAJ⟩
      · have eU : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 - 2 * M * N - N ^ 2) * U = 0 := by linear_combination ((M ^ 2 + N ^ 2) * M) * h1 - ((M ^ 2 + N ^ 2) * N) * h2
        have eV : (M ^ 2 + N ^ 2) ^ 2 * (M ^ 2 - 2 * M * N - N ^ 2) * V = 0 := by linear_combination ((M ^ 2 - 2 * M * N - N ^ 2) * N) * h1 + ((M ^ 2 - 2 * M * N - N ^ 2) * M) * h2
        exact ⟨(mul_eq_zero.mp eU).resolve_left hAAJ, (mul_eq_zero.mp eV).resolve_left hAAJ⟩
      · have eU : (M ^ 2 - 2 * M * N - N ^ 2) * M * U * 2 = 0 := by linear_combination h1 - h2
        have eV : (M ^ 2 + N ^ 2) * N * V * 2 = 0 := by linear_combination h1 + h2
        exact ⟨by simpa [hJ, hM, two_ne] using eU, by simpa [hA, hN, two_ne] using eV⟩
    have key : D * ((M ^ 2 - 2 * M * N - N ^ 2) * (M ^ 2 - N ^ 2) * (((M ^ 2 + N ^ 2) * (M * V) - (M ^ 2 - 2 * M * N - N ^ 2) * (N * U)) *
        ((M ^ 2 + N ^ 2) * (N * V) - (M ^ 2 - 2 * M * N - N ^ 2) * (M * U)))) ^ 2 =
        ((M ^ 2 + N ^ 2) * ((M ^ 2 + N ^ 2) * (a : ZMod p) - (M ^ 2 - 2 * M * N - N ^ 2) * (b : ZMod p)) * M * N * (M ^ 2 + N ^ 2)) ^ 2 := by
      linear_combination ((M ^ 2 + N ^ 2) ^ 6 * M ^ 2 * N ^ 2) * ha' -
        ((M ^ 2 + N ^ 2) ^ 5 * M ^ 2 * N ^ 2 * (M ^ 2 - 2 * M * N - N ^ 2)) * hb' +
        (
          D*M^16*U^2*V^2 - 2*D*M^15*N*U^3*V - 6*D*M^15*N*U^2*V^2 - 2*D*M^15*N*U*V^3 +
          D*M^14*N^2*U^4 + 14*D*M^14*N^2*U^3*V + 6*D*M^14*N^2*U^2*V^2 + 10*D*M^14*N^2*U*V^3 +
          D*M^14*N^2*V^4 - 8*D*M^13*N^3*U^4 - 42*D*M^13*N^3*U^3*V - 14*D*M^13*N^3*U^2*V^2 +
          6*D*M^13*N^3*U*V^3 + 30*D*M^12*N^4*U^4 - 44*D*M^12*N^4*U^3*V + 16*D*M^12*N^4*U^2*V^2
          + 28*D*M^12*N^4*U*V^3 - 2*D*M^12*N^4*V^4 + 8*D*M^11*N^5*U^4 + 14*D*M^11*N^5*U^3*V +
          58*D*M^11*N^5*U^2*V^2 + 30*D*M^11*N^5*U*V^3 - 145*D*M^10*N^6*U^4 -
          46*D*M^10*N^6*U^3*V + 26*D*M^10*N^6*U^2*V^2 + 22*D*M^10*N^6*U*V^3 - 17*D*M^10*N^6*V^4
          + 16*D*M^9*N^7*U^4 + 54*D*M^9*N^7*U^3*V + 66*D*M^9*N^7*U^2*V^2 + 22*D*M^9*N^7*U*V^3 +
          164*D*M^8*N^8*U^4 + 24*D*M^8*N^8*U^3*V + 30*D*M^8*N^8*U^2*V^2 + 8*D*M^8*N^8*U*V^3 -
          28*D*M^8*N^8*V^4 - 16*D*M^7*N^9*U^4 - 54*D*M^7*N^9*U^3*V - 66*D*M^7*N^9*U^2*V^2 -
          22*D*M^7*N^9*U*V^3 - 145*D*M^6*N^10*U^4 - 46*D*M^6*N^10*U^3*V + 26*D*M^6*N^10*U^2*V^2
          + 22*D*M^6*N^10*U*V^3 - 17*D*M^6*N^10*V^4 - 8*D*M^5*N^11*U^4 - 14*D*M^5*N^11*U^3*V -
          58*D*M^5*N^11*U^2*V^2 - 30*D*M^5*N^11*U*V^3 + 30*D*M^4*N^12*U^4 - 44*D*M^4*N^12*U^3*V
          + 16*D*M^4*N^12*U^2*V^2 + 28*D*M^4*N^12*U*V^3 - 2*D*M^4*N^12*V^4 + 8*D*M^3*N^13*U^4 +
          42*D*M^3*N^13*U^3*V + 14*D*M^3*N^13*U^2*V^2 - 6*D*M^3*N^13*U*V^3 + D*M^2*N^14*U^4 +
          14*D*M^2*N^14*U^3*V + 6*D*M^2*N^14*U^2*V^2 + 10*D*M^2*N^14*U*V^3 + D*M^2*N^14*V^4 +
          2*D*M*N^15*U^3*V + 6*D*M*N^15*U^2*V^2 + 2*D*M*N^15*U*V^3 + D*N^16*U^2*V^2 -
          2*M^10*N^2*(b : ZMod p)^2 - 8*M^8*N^4*(b : ZMod p)^2 - 12*M^6*N^6*(b : ZMod p)^2 -
          8*M^4*N^8*(b : ZMod p)^2 - 2*M^2*N^10*(b : ZMod p)^2
        ) * hR
    exact isSquare_of_mul_sq (mul_ne_zero (mul_ne_zero hJ hE) hH) key
  · have key : D * ((M ^ 2 - N ^ 2) * (((M ^ 2 + N ^ 2) * (M * V) + (M ^ 2 - 2 * M * N - N ^ 2) * (N * U)) *
        ((M ^ 2 + N ^ 2) * (N * V) + (M ^ 2 - 2 * M * N - N ^ 2) * (M * U)))) ^ 2 =
        (((M ^ 2 + N ^ 2) * (a : ZMod p) + (M ^ 2 - 2 * M * N - N ^ 2) * (b : ZMod p)) * M * N * (M ^ 2 + N ^ 2)) ^ 2 := by
      linear_combination (M ^ 2 * N ^ 2 * (M ^ 2 + N ^ 2) ^ 4) * ha' +
        (M ^ 2 * N ^ 2 * (M ^ 2 + N ^ 2) ^ 3 * (M ^ 2 - 2 * M * N - N ^ 2)) * hb' +
        (
          D*M^12*U^2*V^2 + 2*D*M^11*N*U^3*V - 2*D*M^11*N*U^2*V^2 + 2*D*M^11*N*U*V^3 +
          D*M^10*N^2*U^4 - 10*D*M^10*N^2*U^3*V + 4*D*M^10*N^2*U^2*V^2 + 2*D*M^10*N^2*U*V^3 +
          D*M^10*N^2*V^4 - 8*D*M^9*N^3*U^4 - 10*D*M^9*N^3*U^3*V - 6*D*M^9*N^3*U^2*V^2 +
          6*D*M^9*N^3*U*V^3 + 4*D*M^8*N^4*U^4 - 8*D*M^8*N^4*U^3*V - 9*D*M^8*N^4*U^2*V^2 +
          8*D*M^8*N^4*U*V^3 + 4*D*M^8*N^4*V^4 + 24*D*M^7*N^5*U^4 - 12*D*M^7*N^5*U^3*V -
          4*D*M^7*N^5*U^2*V^2 + 4*D*M^7*N^5*U*V^3 + 6*D*M^6*N^6*U^4 + 4*D*M^6*N^6*U^3*V -
          24*D*M^6*N^6*U^2*V^2 + 12*D*M^6*N^6*U*V^3 + 6*D*M^6*N^6*V^4 - 24*D*M^5*N^7*U^4 +
          12*D*M^5*N^7*U^3*V + 4*D*M^5*N^7*U^2*V^2 - 4*D*M^5*N^7*U*V^3 + 4*D*M^4*N^8*U^4 -
          8*D*M^4*N^8*U^3*V - 9*D*M^4*N^8*U^2*V^2 + 8*D*M^4*N^8*U*V^3 + 4*D*M^4*N^8*V^4 +
          8*D*M^3*N^9*U^4 + 10*D*M^3*N^9*U^3*V + 6*D*M^3*N^9*U^2*V^2 - 6*D*M^3*N^9*U*V^3 +
          D*M^2*N^10*U^4 - 10*D*M^2*N^10*U^3*V + 4*D*M^2*N^10*U^2*V^2 + 2*D*M^2*N^10*U*V^3 +
          D*M^2*N^10*V^4 - 2*D*M*N^11*U^3*V + 2*D*M*N^11*U^2*V^2 - 2*D*M*N^11*U*V^3 +
          D*N^12*U^2*V^2 - 2*M^6*N^2*(b : ZMod p)^2 - 4*M^4*N^4*(b : ZMod p)^2 - 2*M^2*N^6*(b :
          ZMod p)^2
        ) * hR
    exact isSquare_of_mul_sq (mul_ne_zero hE hz) key

/-- **The class modulo the primes of the octic.** For coprime `m, n` and coprime `u, v`, a `d` such
that `d` times each of the eight products is a Gaussian square is a square modulo every odd prime
dividing `m⁸ + 14m⁴n⁴ + n⁸`. No condition on the ratio, and the four forms need not be squares. -/
theorem isSquare_class_mod_octic {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) (hdvd : (p : ℤ) ∣ m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8)
    {d : ℕ} (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : IsSquare (d : ZMod p) := by
  have e : m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8 =
      (m ^ 4 - 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 + 2 * m * n ^ 3 + n ^ 4) *
      (m ^ 4 + 2 * m ^ 3 * n + 2 * m ^ 2 * n ^ 2 - 2 * m * n ^ 3 + n ^ 4) := by ring
  rw [e] at hdvd
  rcases (Nat.prime_iff_prime_int.mp hp).dvd_or_dvd hdvd with h | h
  · exact isSquare_class_mod_R1 hmn hcop hp hp2 h g
  · exact isSquare_class_mod_R2 hmn hcop hp hp2 h g

/-! ### Lines by Euler's criterion, for primes too large for a mask -/

/-- A necessary condition for `x` to be a square modulo the odd prime `p`:
`x^((p − 1)/2)` is `0` or `1` modulo `p`. -/
def eul (p x : ℕ) : Bool := decide (x ^ ((p - 1) / 2) % p ≤ 1)

theorem eul_of_sq {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {x : ℕ} {c : ZMod p}
    (h : (x : ZMod p) = c ^ 2) : eul p x = true := by
  have : Fact p.Prime := ⟨hp⟩
  have hodd : 2 * ((p - 1) / 2) = p - 1 := by
    have := hp.eq_two_or_odd
    have := hp.two_le
    omega
  have key : ((x ^ ((p - 1) / 2) : ℕ) : ZMod p) = c ^ (p - 1) := by
    push_cast
    rw [h, ← pow_mul, hodd]
  unfold eul
  rw [decide_eq_true_eq]
  by_cases hc : c = 0
  · have h0 : ((x ^ ((p - 1) / 2) : ℕ) : ZMod p) = ((0 : ℕ) : ZMod p) := by
      rw [key, hc, zero_pow (by have := hp.two_le; omega)]
      simp
    have h1 := (ZMod.natCast_eq_natCast_iff _ _ _).mp h0
    unfold Nat.ModEq at h1
    rw [h1]
    simp
  · have h0 : ((x ^ ((p - 1) / 2) : ℕ) : ZMod p) = ((1 : ℕ) : ZMod p) := by
      rw [key, ZMod.pow_card_sub_one_eq_one hc]
      simp
    have h1 := (ZMod.natCast_eq_natCast_iff _ _ _).mp h0
    unfold Nat.ModEq at h1
    rw [h1, Nat.mod_eq_of_lt hp.one_lt]

theorem res_zmod {p : ℕ} (hp0 : 0 < p) (z : ℤ) : ((res p z : ℕ) : ZMod p) = (z : ZMod p) := by
  have h := (ZMod.intCast_eq_intCast_iff _ _ _).mpr (res_cast hp0 z)
  rw [Int.cast_natCast] at h
  exact h

def formsOkE (m n : ℤ) (p : ℕ) (u v : ℤ) : Bool :=
  eul p (res p (nrm (P0 m n u v))) && eul p (res p (nrm (P8 m n u v))) &&
    eul p (res p (nrm (PL m n u v))) && eul p (res p (nrm (PU m n u v)))

/-- `d` times the pair `q` (residues modulo `p`), read in the two residue fields above `p`:
`io` is a square root of `−1` modulo `p`. -/
def gOkE (p io d : ℕ) (q : ℕ × ℕ) : Bool :=
  eul p (d * (q.1 + io * q.2) % p) && eul p (d * (q.1 + (p - io) * q.2) % p)

def dsOkE (p io : ℕ) (ds : List ℕ) (qs : List (ℕ × ℕ)) : Bool :=
  ds.all fun d => !(qs.all fun q => gOkE p io d q)

def pointOkE (m n : ℤ) (p io : ℕ) (ds : List ℕ) (u v : ℤ) : Bool :=
  !formsOkE m n p u v || dsOkE p io ds (flags.map fun f => resQ p (Q m n u v f))

/-- The check `pointOkE` at every point of the projective line modulo `p`. -/
def killsE (m n : ℤ) (p io : ℕ) (ds : List ℕ) : Bool :=
  ((List.range p).all fun v => pointOkE m n p io ds 1 v) && pointOkE m n p io ds 0 1

theorem iota_sq {p io : ℕ} (hp : 1 ≤ p) (h : io * io % p = p - 1) : (io : ZMod p) * io = -1 := by
  have h1 : ((io * io % p : ℕ) : ZMod p) = ((p - 1 : ℕ) : ZMod p) := by rw [h]
  rw [ZMod.natCast_mod, Nat.cast_sub hp, ZMod.natCast_self] at h1
  push_cast at h1
  rw [h1]
  ring

theorem eul_comp {p : ℕ} (hp : p.Prime) (hp2 : p ≠ 2) {j : ℕ} (hj : (j : ZMod p) * j = -1) {d : ℕ}
    {q : ℤ × ℤ} (h : GSqP d q) : eul p (d * (res p q.1 + j * res p q.2) % p) = true := by
  obtain ⟨a, b, ha, hb⟩ := h
  have ha' : (d : ZMod p) * (q.1 : ZMod p) = (a : ZMod p) ^ 2 - (b : ZMod p) ^ 2 := by
    have h1 := congrArg (Int.cast : ℤ → ZMod p) ha
    push_cast at h1
    exact h1
  have hb' : (d : ZMod p) * (q.2 : ZMod p) = 2 * (a : ZMod p) * (b : ZMod p) := by
    have h1 := congrArg (Int.cast : ℤ → ZMod p) hb
    push_cast at h1
    exact h1
  apply eul_of_sq hp hp2 (c := (a : ZMod p) + j * b)
  rw [ZMod.natCast_mod]
  push_cast
  rw [res_zmod hp.pos, res_zmod hp.pos]
  linear_combination ha' + (j : ZMod p) * hb' - (b : ZMod p) ^ 2 * hj

theorem pointOkE_false {m n : ℤ} {p io : ℕ} {ds : List ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hio : io * io % p = p - 1) (hlt : io < p) {d : ℕ} (hd : d ∈ ds) {u v : ℤ}
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) {u' v' : ℤ}
    (cu : u' ≡ u [ZMOD p]) (cv : v' ≡ v [ZMOD p]) : pointOkE m n p io ds u' v' = false := by
  obtain ⟨c0, c8, cL, cU⟩ := nrm_congr (m := m) (n := n) cu cv
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  have sq : ∀ {z z' w : ℤ}, z' ≡ z [ZMOD p] → z = w ^ 2 → eul p (res p z') = true := by
    intro z z' w hc hz
    apply eul_of_sq hp hp2 (c := (w : ZMod p))
    rw [res_zmod hp.pos, (ZMod.intCast_eq_intCast_iff _ _ _).mpr hc, hz]
    exact Int.cast_pow w 2
  have f : formsOkE m n p u' v' = true := by
    unfold formsOkE
    rw [sq c0 hw0, sq c8 hw8, sq cL hwL, sq cU hwU]
    rfl
  have hj1 := iota_sq (p := p) hp.one_lt.le hio
  have hj2 : ((p - io : ℕ) : ZMod p) * ((p - io : ℕ) : ZMod p) = -1 := by
    rw [Nat.cast_sub hlt.le, ZMod.natCast_self, zero_sub, neg_mul_neg]
    exact hj1
  have o : ((flags.map fun f => resQ p (Q m n u' v' f)).all fun q => gOkE p io d q) = true := by
    rw [List.all_eq_true]
    intro q hq
    rw [List.mem_map] at hq
    obtain ⟨f, _, rfl⟩ := hq
    obtain ⟨q1, q2⟩ := Q_congr (m := m) (n := n) cu cv f
    have e1 : res p (Q m n u' v' f).1 = res p (Q m n u v f).1 := res_congr q1
    have e2 : res p (Q m n u' v' f).2 = res p (Q m n u v f).2 := res_congr q2
    unfold gOkE resQ
    dsimp only
    rw [e1, e2, eul_comp hp hp2 hj1 (g f), eul_comp hp hp2 hj2 (g f)]
    rfl
  unfold pointOkE
  rw [f, Bool.not_true, Bool.false_or, Bool.eq_false_iff]
  intro hall
  unfold dsOkE at hall
  rw [List.all_eq_true] at hall
  have := hall d hd
  rw [o] at this
  simp at this

/-- If the check succeeds, no `d` of the list is the class of a solution. -/
theorem not_gsq_of_killsE {m n : ℤ} {p io : ℕ} {ds : List ℕ} (hp : p.Prime) (hp2 : p ≠ 2)
    (hio : io * io % p = p - 1) (hlt : io < p) (hk : killsE m n p io ds = true)
    {d : ℕ} (hd : d ∈ ds) {u v : ℤ} (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (g : ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f)) : False := by
  have hM : 0 < p := hp.pos
  have hpi : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  unfold killsE at hk
  rw [Bool.and_eq_true, List.all_eq_true] at hk
  obtain ⟨hk1, hk2⟩ := hk
  by_cases hu : (p : ℤ) ∣ u
  · have hv : ¬(p : ℤ) ∣ v := fun hv => hpi.not_isUnit (hcop.isUnit_of_dvd' hu hv)
    obtain ⟨l, hl⟩ := exists_inv (k := 1) hp hv
    rw [pow_one] at hl
    obtain ⟨s0, s8, sL, sU⟩ := scale_sq l h0 h8 hL hU
    have cu : (0 : ℤ) ≡ l * u [ZMOD p] := by
      rw [Int.modEq_iff_dvd, sub_zero]
      exact Dvd.dvd.mul_left hu l
    have key := pointOkE_false (ds := ds) hp hp2 hio hlt hd s0 s8 sL sU
      (fun f => scale_gsq l f (g f)) cu hl.symm
    rw [key] at hk2
    exact Bool.false_ne_true hk2
  · obtain ⟨l, hl⟩ := exists_inv (k := 1) hp hu
    rw [pow_one] at hl
    obtain ⟨s0, s8, sL, sU⟩ := scale_sq l h0 h8 hL hU
    have key := pointOkE_false (ds := ds) hp hp2 hio hlt hd s0 s8 sL sU
      (fun f => scale_gsq l f (g f)) hl.symm (res_cast hM (l * v))
    have h3 := hk1 (res p (l * v)) (List.mem_range.mpr (res_lt hM _))
    rw [key] at h3
    exact Bool.false_ne_true h3

/-- A square root of `−1` modulo `p`, by search. -/
def iotaOf (p : ℕ) : Option ℕ := (List.range p).find? fun j => j * j % p == p - 1

theorem iotaOf_spec {p io : ℕ} (h : iotaOf p = some io) : io * io % p = p - 1 ∧ io < p := by
  unfold iotaOf at h
  have h1 := List.find?_some h
  have h2 := List.mem_of_find?_eq_some h
  exact ⟨by simpa using h1, List.mem_range.mp h2⟩

/-! ### The structure of the class, for every ratio -/

/-- **The class of a solution, with what is known about it for every ratio.** For coprime `m, n`
and coprime `u, v` at which the four forms are squares there is a `d > 0` without square factors
such that `d` times each of the eight products is a Gaussian square, and:
every prime factor of `d` divides `2mn(m − n)(m + n)`; every odd prime factor of `d` is
`≡ ±1 (mod 8)`; `d` is even; `d` is a square modulo every odd prime dividing `m² + n²` or
`m⁸ + 14m⁴n⁴ + n⁸`. -/
theorem exists_class_struct {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (n0 : (P0 m n u v).1 ≠ 0 ∨ (P0 m n u v).2 ≠ 0) (n8 : (P8 m n u v).1 ≠ 0 ∨ (P8 m n u v).2 ≠ 0)
    (nL : (PL m n u v).1 ≠ 0 ∨ (PL m n u v).2 ≠ 0) (nU : (PU m n u v).1 ≠ 0 ∨ (PU m n u v).2 ≠ 0) :
    ∃ d : ℕ, 0 < d ∧ NoSq d ∧
      (∀ p : ℕ, p.Prime → p ∣ d → (p : ℤ) ∣ 2 * m * n * (m - n) * (m + n)) ∧
      (∀ p : ℕ, p.Prime → p ∣ d → p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧
      2 ∣ d ∧
      (∀ p : ℕ, p.Prime → p ≠ 2 → (p : ℤ) ∣ m ^ 2 + n ^ 2 → IsSquare (d : ZMod p)) ∧
      (∀ p : ℕ, p.Prime → p ≠ 2 → (p : ℤ) ∣ m ^ 8 + 14 * m ^ 4 * n ^ 4 + n ^ 8 →
        IsSquare (d : ZMod p)) ∧
      ∀ f : Bool × Bool × Bool, GSqP d (Q m n u v f) := by
  obtain ⟨d, hd, hns, hsup, hE, g⟩ := exists_class_prod8 hmn hcop h0 h8 hL hU n0 n8 nL nU
  exact ⟨d, hd, hns, hsup, hE, two_dvd_class hmn hcop h0 h8 hL hU g,
    fun p hp hp2 hdvd => isSquare_class_mod hmn hcop hp hp2 hdvd g,
    fun p hp hp2 hdvd => isSquare_class_mod_octic hmn hcop hp hp2 hdvd g, g⟩

/-! ### The candidates and the certificate -/

/-- All products of sublists of `ps`. -/
def prods : List ℕ → List ℕ
  | [] => [1]
  | p :: ps => prods ps ++ (prods ps).map (p * ·)

theorem mem_prods (ps : List ℕ) (hps : ∀ p ∈ ps, p.Prime) :
    ∀ d : ℕ, 0 < d → NoSq d → (∀ q : ℕ, q.Prime → q ∣ d → q ∈ ps) → d ∈ prods ps := by
  induction ps with
  | nil =>
    intro d _ _ hsup
    have h1 : d = 1 := by
      by_contra h1
      obtain ⟨q, hq, hqd⟩ := Nat.exists_prime_and_dvd h1
      simpa using hsup q hq hqd
    simp [prods, h1]
  | cons p ps ih =>
    intro d hd hns hsup
    have hp : p.Prime := hps p (by simp)
    have hps' : ∀ q ∈ ps, q.Prime := fun q hq => hps q (by simp [hq])
    unfold prods
    rw [List.mem_append, List.mem_map]
    by_cases hpd : p ∣ d
    · right
      obtain ⟨d', rfl⟩ := hpd
      have hd' : 0 < d' := Nat.pos_of_ne_zero (by rintro rfl; simp at hd)
      have hpd' : ¬p ∣ d' := by
        rintro ⟨c, rfl⟩
        exact hns p hp ⟨c, by ring⟩
      refine ⟨d', ih hps' d' hd' ?_ ?_, rfl⟩
      · intro q hq hdvd
        exact hns q hq (Dvd.dvd.mul_left hdvd p)
      · intro q hq hqd
        rcases List.mem_cons.mp (hsup q hq (Dvd.dvd.mul_left hqd p)) with rfl | h
        · exact absurd hqd hpd'
        · exact h
    · left
      refine ih hps' d hd hns ?_
      intro q hq hqd
      rcases List.mem_cons.mp (hsup q hq hqd) with rfl | h
      · exact absurd hqd hpd
      · exact h

theorem mem_of_prime_dvd_prod {ps : List ℕ} (hps : ∀ p ∈ ps, p.Prime) {q : ℕ} (hq : q.Prime)
    (h : q ∣ ps.prod) : q ∈ ps := by
  induction ps with
  | nil =>
    rw [List.prod_nil, Nat.dvd_one] at h
    exact absurd h hq.one_lt.ne'
  | cons p ps ih =>
    rw [List.prod_cons] at h
    rcases (Nat.Prime.dvd_mul hq).mp h with h1 | h1
    · have e := (Nat.prime_dvd_prime_iff_eq hq (hps p (by simp))).mp h1
      simp [e]
    · exact List.mem_cons_of_mem _ (ih (fun p hp => hps p (by simp [hp])) h1)

/-- A primality test by trial division, for the kernel. -/
def isPrime (p : ℕ) : Bool := decide (2 ≤ p) && (List.range p).all fun k => decide (k < 2) || !(p % k == 0)

theorem prime_of_isPrime {p : ℕ} (h : isPrime p = true) : p.Prime := by
  unfold isPrime at h
  rw [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at h
  obtain ⟨h2, hall⟩ := h
  rw [Nat.prime_def_lt]
  refine ⟨h2, ?_⟩
  intro k hk hdvd
  have h3 := hall k (List.mem_range.mpr hk)
  rw [Bool.or_eq_true, decide_eq_true_eq] at h3
  rcases h3 with h3 | h3
  · have hk0 : k ≠ 0 := by
      rintro rfl
      rw [Nat.zero_dvd] at hdvd
      omega
    omega
  · rw [Nat.mod_eq_zero_of_dvd hdvd] at h3
    simp at h3

/-- The masks of a modulus, in a list of triples `(M, S, G)`. -/
def findMask (M : ℕ) : List (ℕ × ℕ × ℕ) → Option (ℕ × ℕ)
  | [] => none
  | e :: l => if e.1 = M then some e.2 else findMask M l

theorem mem_of_findMask {M S G : ℕ} {l : List (ℕ × ℕ × ℕ)} (h : findMask M l = some (S, G)) :
    (M, S, G) ∈ l := by
  induction l with
  | nil => simp [findMask] at h
  | cons e l ih =>
    unfold findMask at h
    split at h
    · rename_i he
      obtain ⟨M', S', G'⟩ := e
      simp only [Option.some.injEq] at h
      simp only at he
      subst he
      have h1 : (S', G') = (S, G) := h
      rw [h1]
      exact List.mem_cons_self
    · exact List.mem_cons_of_mem _ (ih h)

/-- The masks in the table are the masks of their moduli. -/
def MasksOk (masks : List (ℕ × ℕ × ℕ)) : Prop :=
  ∀ e ∈ masks, e.2.1 = sMask e.1 ∧ e.2.2 = gMask e.1

/-- The elements of `l` at the positions of the bits of `k`. -/
def pick : List ℕ → ℕ → List ℕ
  | [], _ => []
  | d :: l, k => if k % 2 = 1 then d :: pick l (k / 2) else pick l (k / 2)

/-- Every position of `l` is a bit of one of the numbers `ks`. -/
def coveredBy : List ℕ → List ℕ → Bool
  | [], _ => true
  | _ :: l, ks => (ks.any fun k => decide (k % 2 = 1)) && coveredBy l (ks.map (· / 2))

theorem mem_pick_of_coveredBy : ∀ (l ks : List ℕ), coveredBy l ks = true →
    ∀ d ∈ l, ∃ k ∈ ks, d ∈ pick l k := by
  intro l
  induction l with
  | nil =>
    intro ks _ d hd
    simp at hd
  | cons x l ih =>
    intro ks h d hd
    unfold coveredBy at h
    rw [Bool.and_eq_true, List.any_eq_true] at h
    obtain ⟨⟨k, hk, hodd⟩, hrest⟩ := h
    rw [decide_eq_true_eq] at hodd
    rcases List.mem_cons.mp hd with rfl | hd'
    · refine ⟨k, hk, ?_⟩
      unfold pick
      rw [if_pos hodd]
      exact List.mem_cons_self
    · obtain ⟨k', hk', hmem⟩ := ih (ks.map (· / 2)) hrest d hd'
      rw [List.mem_map] at hk'
      obtain ⟨k0, hk0, rfl⟩ := hk'
      refine ⟨k0, hk0, ?_⟩
      unfold pick
      split
      · exact List.mem_cons_of_mem _ hmem
      · exact hmem

/-- One line of a certificate: `(p, k, c)`. Modulo `M = p^k` the classes at the positions of the
bits of `c` in the list of candidates `cands` are excluded. -/
def lineOk (masks : List (ℕ × ℕ × ℕ)) (m n : ℤ) (cands : List ℕ) (c : ℕ × ℕ × ℕ) : Bool :=
  match findMask (c.1 ^ c.2.1) masks with
  | some (S, G) =>
    isPrime c.1 && decide (0 < c.2.1) && kills m n (c.1 ^ c.2.1) c.1 S G (pick cands c.2.2)
  | none => false

/-- An entry of the table: `(m, n, ps, cert)`. `ps` lists primes; every prime factor of
`2mn(m − n)(m + n)` is among them; every product of some of them is excluded by a line. -/
def entryOk (masks : List (ℕ × ℕ × ℕ)) (e : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)) : Bool :=
  decide (0 < e.2.1) && decide (e.2.1 < e.1) && e.2.2.1.all isPrime &&
    decide (2 * e.1 * e.2.1 * (e.1 - e.2.1) * (e.1 + e.2.1) ∣ e.2.2.1.prod ^ 64) &&
    e.2.2.2.all (lineOk masks (e.1 : ℤ) (e.2.1 : ℤ) (prods e.2.2.1)) &&
    coveredBy (prods e.2.2.1) (e.2.2.2.map fun c => c.2.2)

/-- **No point, from a certificate.** If `entryOk` holds for `(m, n, ps, cert)`, no coprime `u, v`
with `v ≠ 0`, `u + v ≠ 0` make the four forms squares. -/
theorem no_point_of_entry {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks)
    {e : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)} (he : entryOk masks e = true) {u v : ℤ}
    (hcop : IsCoprime u v) (hv : v ≠ 0) (huv : u + v ≠ 0)
    (h0 : ∃ w, nrm (P0 e.1 e.2.1 u v) = w ^ 2) (h8 : ∃ w, nrm (P8 e.1 e.2.1 u v) = w ^ 2)
    (hL : ∃ w, nrm (PL e.1 e.2.1 u v) = w ^ 2) (hU : ∃ w, nrm (PU e.1 e.2.1 u v) = w ^ 2) : False := by
  obtain ⟨m, n, ps, cert⟩ := e
  unfold entryOk at he
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at he
  obtain ⟨⟨⟨⟨⟨hn, hnm⟩, hps⟩, hdvd⟩, hlines⟩, hcov⟩ := he
  dsimp only at h0 h8 hL hU
  have hps' : ∀ p ∈ ps, p.Prime := fun p hp => prime_of_isPrime (hps p hp)
  have hm0 : (m : ℤ) ≠ 0 := by
    have : 0 < m := by omega
    exact_mod_cast this.ne'
  have hn0 : (n : ℤ) ≠ 0 := by exact_mod_cast hn.ne'
  have hmn0 : (m : ℤ) - n ≠ 0 := by
    have : (n : ℤ) < m := by exact_mod_cast hnm
    intro h
    omega
  have hbad : ∀ q : ℕ, q.Prime → (q : ℤ) ∣ 2 * (m : ℤ) * n * (m - n) * (m + n) → q ∈ ps := by
    intro q hq hdiv
    have e1 : (2 * (m : ℤ) * n * (m - n) * (m + n)) = ((2 * m * n * (m - n) * (m + n) : ℕ) : ℤ) := by
      push_cast [Nat.cast_sub hnm.le]
      ring
    rw [e1] at hdiv
    have h1 : q ∣ 2 * m * n * (m - n) * (m + n) := Int.natCast_dvd_natCast.mp hdiv
    exact mem_of_prime_dvd_prod hps' hq (hq.dvd_of_dvd_pow (dvd_trans h1 hdvd))
  obtain ⟨d, hd, hns, hsup, g⟩ := exists_class_prod hcop h0 h8 hL hU
    (Or.inl (by dsimp only [P0]; exact mul_ne_zero hm0 hv))
    (Or.inl (by dsimp only [P8]; exact mul_ne_zero hn0 hv))
    (Or.inr (by dsimp only [PL]; exact mul_ne_zero hmn0 huv))
    (Or.inr (by
      dsimp only [PU]
      refine mul_ne_zero ?_ huv
      have : (0 : ℤ) < (m : ℤ) + n := by positivity
      exact this.ne'))
  have hmem : d ∈ prods ps := mem_prods ps hps' d hd hns (fun q hq hqd => hbad q hq (hsup q hq hqd))
  obtain ⟨k0, hk0, hdc⟩ := mem_pick_of_coveredBy _ _ hcov d hmem
  rw [List.mem_map] at hk0
  obtain ⟨c, hc, rfl⟩ := hk0
  have h2 := hlines c hc
  unfold lineOk at h2
  split at h2
  · rename_i S G hfind
    simp only [Bool.and_eq_true, decide_eq_true_eq] at h2
    obtain ⟨⟨hp, hk0⟩, hk⟩ := h2
    have h3 := hmasks _ (mem_of_findMask hfind)
    exact not_gsq_of_kills (prime_of_isPrime hp) rfl hk0 h3.1 h3.2 hk hdc hcop h0 h8 hL hU g
  · exact Bool.false_ne_true h2

/-! ### Certificates with Theorem E: fewer candidates; lines of two kinds -/

/-- One line of a certificate, of one of two kinds. `(p, k, c)` with `k > 0`: as `lineOk`, modulo
`p^k` with the masks. `(p, 0, c)`: modulo the odd prime `p`, by Euler's criterion in the two
residue fields above `p` (`killsE`); no mask is needed, so `p` may be large. -/
def lineOkX (masks : List (ℕ × ℕ × ℕ)) (m n : ℤ) (cands : List ℕ) (c : ℕ × ℕ × ℕ) : Bool :=
  if c.2.1 = 0 then
    match iotaOf c.1 with
    | some io => isPrime c.1 && decide (c.1 ≠ 2) && killsE m n c.1 io (pick cands c.2.2)
    | none => false
  else lineOk masks m n cands c


/-- The primes that can divide a class: 2 and the primes `≡ ±1 (mod 8)`. -/
def keep8 (p : ℕ) : Bool := p == 2 || p % 8 == 1 || p % 8 == 7

/-- An entry of the table, checked with the candidates of Theorem E: `(m, n, ps, cert)`, with
`m, n` coprime; `ps` lists primes, among them every prime factor of `2mn(m − n)(m + n)`; every
product of those primes of `ps` that are 2 or `≡ ±1 (mod 8)` is excluded by a line (`lineOkX`). -/
def entryOk8 (masks : List (ℕ × ℕ × ℕ)) (e : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)) : Bool :=
  decide (0 < e.2.1) && decide (e.2.1 < e.1) && decide (Nat.Coprime e.1 e.2.1) &&
    e.2.2.1.all isPrime &&
    decide (2 * e.1 * e.2.1 * (e.1 - e.2.1) * (e.1 + e.2.1) ∣ e.2.2.1.prod ^ 64) &&
    e.2.2.2.all (lineOkX masks (e.1 : ℤ) (e.2.1 : ℤ) (prods (e.2.2.1.filter keep8))) &&
    coveredBy (prods (e.2.2.1.filter keep8)) (e.2.2.2.map fun c => c.2.2)

/-- **No point, from a certificate with the candidates of Theorem E.** -/
theorem no_point_of_entry8 {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks)
    {e : ℕ × ℕ × List ℕ × List (ℕ × ℕ × ℕ)} (he : entryOk8 masks e = true) {u v : ℤ}
    (hcop : IsCoprime u v) (hv : v ≠ 0) (huv : u + v ≠ 0)
    (h0 : ∃ w, nrm (P0 e.1 e.2.1 u v) = w ^ 2) (h8 : ∃ w, nrm (P8 e.1 e.2.1 u v) = w ^ 2)
    (hL : ∃ w, nrm (PL e.1 e.2.1 u v) = w ^ 2) (hU : ∃ w, nrm (PU e.1 e.2.1 u v) = w ^ 2) : False := by
  obtain ⟨m, n, ps, cert⟩ := e
  unfold entryOk8 at he
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at he
  obtain ⟨⟨⟨⟨⟨⟨hn, hnm⟩, hmn⟩, hps⟩, hdvd⟩, hlines⟩, hcov⟩ := he
  dsimp only at h0 h8 hL hU
  have hps' : ∀ p ∈ ps, p.Prime := fun p hp => prime_of_isPrime (hps p hp)
  have hps8 : ∀ p ∈ ps.filter keep8, p.Prime := fun p hp => hps' p (List.mem_filter.mp hp).1
  have hm0 : (m : ℤ) ≠ 0 := by
    have : 0 < m := by omega
    exact_mod_cast this.ne'
  have hn0 : (n : ℤ) ≠ 0 := by exact_mod_cast hn.ne'
  have hmn0 : (m : ℤ) - n ≠ 0 := by
    have : (n : ℤ) < m := by exact_mod_cast hnm
    intro h
    omega
  have hbad : ∀ q : ℕ, q.Prime → (q : ℤ) ∣ 2 * (m : ℤ) * n * (m - n) * (m + n) → q ∈ ps := by
    intro q hq hdiv
    have e1 : (2 * (m : ℤ) * n * (m - n) * (m + n)) = ((2 * m * n * (m - n) * (m + n) : ℕ) : ℤ) := by
      push_cast [Nat.cast_sub hnm.le]
      ring
    rw [e1] at hdiv
    have h1 : q ∣ 2 * m * n * (m - n) * (m + n) := Int.natCast_dvd_natCast.mp hdiv
    exact mem_of_prime_dvd_prod hps' hq (hq.dvd_of_dvd_pow (dvd_trans h1 hdvd))
  obtain ⟨d, hd, hns, hsup, hE, g⟩ := exists_class_prod8 (Nat.isCoprime_iff_coprime.mpr hmn)
    hcop h0 h8 hL hU
    (Or.inl (by dsimp only [P0]; exact mul_ne_zero hm0 hv))
    (Or.inl (by dsimp only [P8]; exact mul_ne_zero hn0 hv))
    (Or.inr (by dsimp only [PL]; exact mul_ne_zero hmn0 huv))
    (Or.inr (by
      dsimp only [PU]
      refine mul_ne_zero ?_ huv
      have : (0 : ℤ) < (m : ℤ) + n := by positivity
      exact this.ne'))
  have hmem : d ∈ prods (ps.filter keep8) := by
    refine mem_prods _ hps8 d hd hns (fun q hq hqd => ?_)
    refine List.mem_filter.mpr ⟨hbad q hq (hsup q hq hqd), ?_⟩
    unfold keep8
    rcases hE q hq hqd with h | h | h
    · simp [h]
    · simp [h]
    · simp [h]
  obtain ⟨k0, hk0, hdc⟩ := mem_pick_of_coveredBy _ _ hcov d hmem
  rw [List.mem_map] at hk0
  obtain ⟨c, hc, rfl⟩ := hk0
  have h2 := hlines c hc
  unfold lineOkX at h2
  split at h2
  · split at h2
    · rename_i io hfind
      simp only [Bool.and_eq_true, decide_eq_true_eq] at h2
      obtain ⟨⟨hp, hp2⟩, hk⟩ := h2
      obtain ⟨hio, hlt⟩ := iotaOf_spec hfind
      exact not_gsq_of_killsE (prime_of_isPrime hp) hp2 hio hlt hk hdc hcop h0 h8 hL hU g
    · exact Bool.false_ne_true h2
  · unfold lineOk at h2
    split at h2
    · rename_i S G hfind
      simp only [Bool.and_eq_true, decide_eq_true_eq] at h2
      obtain ⟨⟨hp, hk0⟩, hk⟩ := h2
      have h3 := hmasks _ (mem_of_findMask hfind)
      exact not_gsq_of_kills (prime_of_isPrime hp) rfl hk0 h3.1 h3.2 hk hdc hcop h0 h8 hL hU g
    · exact Bool.false_ne_true h2

/-! ### The class from its structure alone -/

/-- An entry `(m, n, ps, d0, qs)`: `m, n` coprime; `ps` lists primes, among them every prime factor
of `2mn(m − n)(m + n)`; `qs` lists odd primes dividing `(m² + n²)(m⁸ + 14m⁴n⁴ + n⁸)`; and every
candidate (a product of the primes of `ps` that are 2 or `≡ ±1 (mod 8)`) other than `d0` is odd,
or is not a square modulo one of the primes of `qs` (by Euler's criterion). -/
def structOk (e : ℕ × ℕ × List ℕ × ℕ × List ℕ) : Bool :=
  decide (0 < e.2.1) && decide (e.2.1 < e.1) && decide (Nat.Coprime e.1 e.2.1) &&
    e.2.2.1.all isPrime &&
    decide (2 * e.1 * e.2.1 * (e.1 - e.2.1) * (e.1 + e.2.1) ∣ e.2.2.1.prod ^ 64) &&
    e.2.2.2.2.all (fun p => isPrime p && decide (p ≠ 2) &&
      decide (p ∣ (e.1 ^ 2 + e.2.1 ^ 2) * (e.1 ^ 8 + 14 * e.1 ^ 4 * e.2.1 ^ 4 + e.2.1 ^ 8))) &&
    (prods (e.2.2.1.filter keep8)).all (fun d => d == e.2.2.2.1 || d % 2 == 1 ||
      e.2.2.2.2.any fun p => !eul p (d % p))

/-- **The class of a solution, from the theorems on its structure alone.** If `structOk` holds
for `(m, n, ps, d0, qs)`, then at every solution on the four forms `d0` times each of the eight
products is a Gaussian square: the class is `d0`. No local computation at the points is used. -/
theorem class_of_structOk {e : ℕ × ℕ × List ℕ × ℕ × List ℕ} (he : structOk e = true) {u v : ℤ}
    (hcop : IsCoprime u v) (hv : v ≠ 0) (huv : u + v ≠ 0)
    (h0 : ∃ w, nrm (P0 e.1 e.2.1 u v) = w ^ 2) (h8 : ∃ w, nrm (P8 e.1 e.2.1 u v) = w ^ 2)
    (hL : ∃ w, nrm (PL e.1 e.2.1 u v) = w ^ 2) (hU : ∃ w, nrm (PU e.1 e.2.1 u v) = w ^ 2) :
    ∀ f : Bool × Bool × Bool, GSqP (e.2.2.2.1 : ℤ) (Q e.1 e.2.1 u v f) := by
  obtain ⟨m, n, ps, d0, qs⟩ := e
  unfold structOk at he
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at he
  obtain ⟨⟨⟨⟨⟨⟨hn, hnm⟩, hmn⟩, hps⟩, hdvd⟩, hqs⟩, hall⟩ := he
  dsimp only at h0 h8 hL hU ⊢
  have hps' : ∀ p ∈ ps, p.Prime := fun p hp => prime_of_isPrime (hps p hp)
  have hps8 : ∀ p ∈ ps.filter keep8, p.Prime := fun p hp => hps' p (List.mem_filter.mp hp).1
  have hm0 : (m : ℤ) ≠ 0 := by
    have : 0 < m := by omega
    exact_mod_cast this.ne'
  have hn0 : (n : ℤ) ≠ 0 := by exact_mod_cast hn.ne'
  have hmn0 : (m : ℤ) - n ≠ 0 := by
    have : (n : ℤ) < m := by exact_mod_cast hnm
    intro h
    omega
  have hbad : ∀ q : ℕ, q.Prime → (q : ℤ) ∣ 2 * (m : ℤ) * n * (m - n) * (m + n) → q ∈ ps := by
    intro q hq hdiv
    have e1 : (2 * (m : ℤ) * n * (m - n) * (m + n)) = ((2 * m * n * (m - n) * (m + n) : ℕ) : ℤ) := by
      push_cast [Nat.cast_sub hnm.le]
      ring
    rw [e1] at hdiv
    have h1 : q ∣ 2 * m * n * (m - n) * (m + n) := Int.natCast_dvd_natCast.mp hdiv
    exact mem_of_prime_dvd_prod hps' hq (hq.dvd_of_dvd_pow (dvd_trans h1 hdvd))
  obtain ⟨d, hd, hns, hsup, hE, h2, hG, hO, g⟩ := exists_class_struct
    (Nat.isCoprime_iff_coprime.mpr hmn) hcop h0 h8 hL hU
    (Or.inl (by dsimp only [P0]; exact mul_ne_zero hm0 hv))
    (Or.inl (by dsimp only [P8]; exact mul_ne_zero hn0 hv))
    (Or.inr (by dsimp only [PL]; exact mul_ne_zero hmn0 huv))
    (Or.inr (by
      dsimp only [PU]
      refine mul_ne_zero ?_ huv
      have : (0 : ℤ) < (m : ℤ) + n := by positivity
      exact this.ne'))
  have hmem : d ∈ prods (ps.filter keep8) := by
    refine mem_prods _ hps8 d hd hns (fun q hq hqd => ?_)
    refine List.mem_filter.mpr ⟨hbad q hq (hsup q hq hqd), ?_⟩
    unfold keep8
    rcases hE q hq hqd with h | h | h
    · simp [h]
    · simp [h]
    · simp [h]
  have h3 := hall d hmem
  rcases (Bool.or_eq_true _ _).mp h3 with h4 | h4
  · rcases (Bool.or_eq_true _ _).mp h4 with h5 | h5
    · have h6 : d = d0 := beq_iff_eq.mp h5
      rw [← h6]
      exact g
    · exfalso
      have h6 : d % 2 = 1 := beq_iff_eq.mp h5
      omega
  · exfalso
    obtain ⟨p, hp, hne⟩ := List.any_eq_true.mp h4
    obtain ⟨⟨hpp, hp2⟩, hdv⟩ := hqs p hp
    have hp' := prime_of_isPrime hpp
    have hdvZ : (p : ℤ) ∣ ((m : ℤ) ^ 2 + (n : ℤ) ^ 2) *
        ((m : ℤ) ^ 8 + 14 * (m : ℤ) ^ 4 * (n : ℤ) ^ 4 + (n : ℤ) ^ 8) := by
      have h7 := Int.natCast_dvd_natCast.mpr hdv
      push_cast at h7
      exact h7
    have hsq : IsSquare (d : ZMod p) := by
      rcases (Nat.prime_iff_prime_int.mp hp').dvd_or_dvd hdvZ with h7 | h7
      · exact hG p hp' hp2 h7
      · exact hO p hp' hp2 h7
    obtain ⟨r, hr⟩ := hsq
    have h8' : eul p (d % p) = true := by
      apply eul_of_sq hp' hp2 (c := r)
      rw [ZMod.natCast_mod, hr]
      ring
    rw [h8'] at hne
    exact Bool.false_ne_true hne

/-! ### Four classes, one for each form -/

/-- As `exists_class`, keeping the fact that `d` has no square factor. -/
theorem exists_class_nosq {x y w : ℤ} (h : x ^ 2 + y ^ 2 = w ^ 2) (hne : x ≠ 0 ∨ y ≠ 0) :
    ∃ d : ℕ, 0 < d ∧ NoSq d ∧ GSq d x y ∧
      ∀ p : ℕ, p.Prime → p ∣ d → p = 2 ∨ ((p : ℤ) ∣ x ∧ (p : ℤ) ∣ y) := by
  obtain ⟨k, hk, hg⟩ := exists_gsq_of_norm h hne
  obtain ⟨d, e, hd, he, hde, hns⟩ := exists_nosq k hk
  have hg' : GSq (d : ℤ) x y := by
    apply GSq.reduce e he
    have e1 : ((k : ℕ) : ℤ) = (d : ℤ) * (e : ℤ) ^ 2 := by exact_mod_cast hde
    rw [← e1]
    exact hg
  refine ⟨d, hd, hns, hg', ?_⟩
  intro p hp hpd
  by_cases hp2 : p = 2
  · exact Or.inl hp2
  · exact Or.inr (hg'.support hp hp2 hpd (hns p hp))

theorem two_cast_dvd (x : ℤ) : ((2 : ℕ) : ℤ) ∣ 2 * x := by
  rw [Nat.cast_ofNat]
  exact dvd_mul_right 2 x

/-- **Four classes.** For coprime `m, n` and coprime `u, v` at which the four forms are squares,
each of the four pairs has its own class: a positive integer without square factors that makes
the pair a Gaussian square, whose prime factors are 2 or `≡ ±1 (mod 8)` and divide `2mn` (for
`P0`, `P8`) or `2(m + n)(m − n)` (for `PL`, `PU`). -/
theorem exists_classes4 {m n u v : ℤ} (hmn : IsCoprime m n) (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    (n0 : (P0 m n u v).1 ≠ 0 ∨ (P0 m n u v).2 ≠ 0) (n8 : (P8 m n u v).1 ≠ 0 ∨ (P8 m n u v).2 ≠ 0)
    (nL : (PL m n u v).1 ≠ 0 ∨ (PL m n u v).2 ≠ 0) (nU : (PU m n u v).1 ≠ 0 ∨ (PU m n u v).2 ≠ 0) :
    ∃ d0 d8 dL dU : ℕ,
      (0 < d0 ∧ NoSq d0 ∧ GSqP d0 (P0 m n u v) ∧ ∀ p : ℕ, p.Prime → p ∣ d0 →
        (p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧ (p : ℤ) ∣ 2 * (m * n)) ∧
      (0 < d8 ∧ NoSq d8 ∧ GSqP d8 (P8 m n u v) ∧ ∀ p : ℕ, p.Prime → p ∣ d8 →
        (p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧ (p : ℤ) ∣ 2 * (m * n)) ∧
      (0 < dL ∧ NoSq dL ∧ GSqP dL (PL m n u v) ∧ ∀ p : ℕ, p.Prime → p ∣ dL →
        (p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧ (p : ℤ) ∣ 2 * ((m + n) * (m - n))) ∧
      (0 < dU ∧ NoSq dU ∧ GSqP dU (PU m n u v) ∧ ∀ p : ℕ, p.Prime → p ∣ dU →
        (p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧ (p : ℤ) ∣ 2 * ((m + n) * (m - n))) := by
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  obtain ⟨d0, hd0, hs0, g0, s0⟩ := exists_class_nosq hw0 n0
  obtain ⟨d8, hd8, hs8, g8, s8⟩ := exists_class_nosq hw8 n8
  obtain ⟨dL, hdL, hsL, gL, sL⟩ := exists_class_nosq hwL nL
  obtain ⟨dU, hdU, hsU, gU, sU⟩ := exists_class_nosq hwU nU
  obtain ⟨a, b, hab⟩ := hcop
  have E := fun (p : ℕ) (hp : p.Prime) (hp2 : p ≠ 2)
      (h : Dv p (P0 m n u v) ∨ Dv p (P8 m n u v) ∨ Dv p (PL m n u v) ∨ Dv p (PU m n u v)) =>
    mod_eight_of_dv hmn ⟨a, b, hab⟩ ⟨w0, hw0⟩ ⟨wL, hwL⟩ hp hp2 h
  refine ⟨d0, d8, dL, dU, ⟨hd0, hs0, g0, ?_⟩, ⟨hd8, hs8, g8, ?_⟩, ⟨hdL, hsL, gL, ?_⟩,
    ⟨hdU, hsU, gU, ?_⟩⟩
  · intro p hp hpd
    by_cases hp2 : p = 2
    · subst hp2
      exact ⟨Or.inl rfl, two_cast_dvd _⟩
    · have hD := (s0 p hp hpd).resolve_left hp2
      refine ⟨Or.inr (E p hp hp2 (Or.inl hD)), Dvd.dvd.mul_left ?_ 2⟩
      obtain ⟨hx, hy⟩ := hD
      have e2 : m * n = b * n * (m * v) + a * m * (n * u) := by
        linear_combination (-(m * n)) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
  · intro p hp hpd
    by_cases hp2 : p = 2
    · subst hp2
      exact ⟨Or.inl rfl, two_cast_dvd _⟩
    · have hD := (s8 p hp hpd).resolve_left hp2
      refine ⟨Or.inr (E p hp hp2 (Or.inr (Or.inl hD))), Dvd.dvd.mul_left ?_ 2⟩
      obtain ⟨hx, hy⟩ := hD
      have e2 : m * n = b * m * (n * v) + a * n * (m * u) := by
        linear_combination (-(m * n)) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
  · intro p hp hpd
    by_cases hp2 : p = 2
    · subst hp2
      exact ⟨Or.inl rfl, two_cast_dvd _⟩
    · have hD := (sL p hp hpd).resolve_left hp2
      refine ⟨Or.inr (E p hp hp2 (Or.inr (Or.inr (Or.inl hD)))), ?_⟩
      obtain ⟨hx, hy⟩ := hD
      have e2 : 2 * ((m + n) * (m - n)) =
          (a - b) * (m - n) * ((m + n) * (u - v)) + (a + b) * (m + n) * ((m - n) * (u + v)) := by
        linear_combination (-(2 * (m + n) * (m - n))) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)
  · intro p hp hpd
    by_cases hp2 : p = 2
    · subst hp2
      exact ⟨Or.inl rfl, two_cast_dvd _⟩
    · have hD := (sU p hp hpd).resolve_left hp2
      refine ⟨Or.inr (E p hp hp2 (Or.inr (Or.inr (Or.inr hD)))), ?_⟩
      obtain ⟨hx, hy⟩ := hD
      have e2 : 2 * ((m + n) * (m - n)) =
          (a - b) * (m + n) * ((m - n) * (u - v)) + (a + b) * (m - n) * ((m + n) * (u + v)) := by
        linear_combination (-(2 * (m + n) * (m - n))) * hab
      rw [e2]
      exact dvd_add (Dvd.dvd.mul_left hx _) (Dvd.dvd.mul_left hy _)

theorem P_congr {M : ℕ} {m n u v u' v' : ℤ} (hu : u ≡ u' [ZMOD M]) (hv : v ≡ v' [ZMOD M]) :
    ((P0 m n u v).1 ≡ (P0 m n u' v').1 [ZMOD M] ∧ (P0 m n u v).2 ≡ (P0 m n u' v').2 [ZMOD M]) ∧
    ((P8 m n u v).1 ≡ (P8 m n u' v').1 [ZMOD M] ∧ (P8 m n u v).2 ≡ (P8 m n u' v').2 [ZMOD M]) ∧
    ((PL m n u v).1 ≡ (PL m n u' v').1 [ZMOD M] ∧ (PL m n u v).2 ≡ (PL m n u' v').2 [ZMOD M]) ∧
    ((PU m n u v).1 ≡ (PU m n u' v').1 [ZMOD M] ∧ (PU m n u v).2 ≡ (PU m n u' v').2 [ZMOD M]) := by
  refine ⟨⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩, ⟨?_, ?_⟩⟩ <;> dsimp only [P0, P8, PL, PU] <;> gcongr

/-- The candidates for the classes of `P0`, `P8`. -/
def candA (m n : ℕ) (ps : List ℕ) : List ℕ :=
  prods (ps.filter fun p => keep8 p && decide (p ∣ 2 * (m * n)))

/-- The candidates for the classes of `PL`, `PU`. -/
def candB (m n : ℕ) (ps : List ℕ) : List ℕ :=
  prods (ps.filter fun p => keep8 p && decide (p ∣ 2 * ((m + n) * (m - n))))

/-- Four lists of classes, one for each form. -/
abbrev Box := List ℕ × List ℕ × List ℕ × List ℕ

/-- The tuple `t` has each of its four classes in the list of the box. -/
abbrev In4 (t : ℕ × ℕ × ℕ × ℕ) (b : Box) : Prop :=
  t.1 ∈ b.1 ∧ t.2.1 ∈ b.2.1 ∧ t.2.2.1 ∈ b.2.2.1 ∧ t.2.2.2 ∈ b.2.2.2

/-- The classes of the lists `A`, `A`, `B`, `B` that make the four pairs Gaussian squares at the
pair of residues `(u, v)` modulo `M`. -/
def rec4 (m n : ℤ) (M G : ℕ) (A B : List ℕ) (u v : ℤ) : Box :=
  (A.filter fun d => gOkN M G d (resQ M (P0 m n u v)),
   A.filter fun d => gOkN M G d (resQ M (P8 m n u v)),
   B.filter fun d => gOkN M G d (resQ M (PL m n u v)),
   B.filter fun d => gOkN M G d (resQ M (PU m n u v)))

/-- The same elements, without repetitions. -/
def dedupL : List Box → List Box
  | [] => []
  | x :: l => if (dedupL l).contains x then dedupL l else x :: dedupL l

theorem mem_dedupL {x : Box} : ∀ {l : List Box}, x ∈ l → x ∈ dedupL l := by
  intro l
  induction l with
  | nil =>
    intro h
    exact h
  | cons y l ih =>
    intro h
    rw [dedupL]
    rcases List.mem_cons.mp h with rfl | h'
    · split
      · rename_i hc
        simpa using hc
      · exact List.mem_cons_self
    · split
      · exact ih h'
      · exact List.mem_cons_of_mem _ (ih h')

/-- The boxes of all pairs of residues `(u, v)` modulo `M`, a power of `p`, not both divisible by
`p`, at which the four forms are squares modulo `M`; each box once. -/
def recs4 (m n : ℤ) (M p S G : ℕ) (A B : List ℕ) : List Box :=
  dedupL ((List.range M).flatMap fun u => (List.range M).filterMap fun v =>
    if !(u % p == 0 && v % p == 0) && formsOk m n M S u v then some (rec4 m n M G A B u v)
    else none)

/-- The classes of a solution are in one of the boxes, modulo every prime power. -/
theorem exists_rec4 {m n : ℤ} {M p k : ℕ} (hp : p.Prime) (hpM : M = p ^ k) (hk0 : 0 < k)
    {u v : ℤ} (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    {d0 d8 dL dU : ℕ} (g0 : GSqP d0 (P0 m n u v)) (g8 : GSqP d8 (P8 m n u v))
    (gL : GSqP dL (PL m n u v)) (gU : GSqP dU (PU m n u v))
    {A B : List ℕ} (hd0 : d0 ∈ A) (hd8 : d8 ∈ A) (hdL : dL ∈ B) (hdU : dU ∈ B) :
    ∃ r ∈ recs4 m n M p (sMask M) (gMask M) A B, In4 (d0, d8, dL, dU) r := by
  have hM : 0 < M := by
    rw [hpM]
    exact pow_pos hp.pos k
  have cu := res_cast hM u
  have cv := res_cast hM v
  obtain ⟨c0, c8, cL, cU⟩ := nrm_congr (m := m) (n := n) cu cv
  obtain ⟨⟨a0, b0⟩, ⟨a8, b8⟩, ⟨aL, bL⟩, ⟨aU, bU⟩⟩ := P_congr (m := m) (n := n) cu cv
  obtain ⟨w0, hw0⟩ := h0
  obtain ⟨w8, hw8⟩ := h8
  obtain ⟨wL, hwL⟩ := hL
  obtain ⟨wU, hwU⟩ := hU
  have f : formsOk m n M (sMask M) (res M u) (res M v) = true := by
    unfold formsOk
    rw [sqOk_of hM (w := w0) (hw0 ▸ c0), sqOk_of hM (w := w8) (hw8 ▸ c8),
      sqOk_of hM (w := wL) (hwL ▸ cL), sqOk_of hM (w := wU) (hwU ▸ cU)]
    rfl
  have hpi : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
  have hpM' : (p : ℤ) ∣ (M : ℤ) := by
    rw [hpM]
    push_cast
    exact dvd_pow_self _ hk0.ne'
  have e : ∀ {x : ℤ}, ((res M x : ℕ) : ℤ) ≡ x [ZMOD (M : ℤ)] → res M x % p = 0 → (p : ℤ) ∣ x := by
    intro x hx h2
    have h3 : (p : ℤ) ∣ ((res M x : ℕ) : ℤ) := by exact_mod_cast Nat.dvd_of_mod_eq_zero h2
    have h4 : (M : ℤ) ∣ x - (res M x : ℕ) := Int.ModEq.dvd hx
    have h5 := dvd_add (dvd_trans hpM' h4) h3
    rwa [sub_add_cancel] at h5
  have nb : (res M u % p == 0 && res M v % p == 0) = false := by
    rw [Bool.eq_false_iff]
    intro h
    simp only [Bool.and_eq_true, beq_iff_eq] at h
    exact hpi.not_isUnit (hcop.isUnit_of_dvd' (e cu h.1) (e cv h.2))
  have hc : (!(res M u % p == 0 && res M v % p == 0) &&
      formsOk m n M (sMask M) (res M u) (res M v)) = true := by
    rw [nb, f]
    rfl
  refine ⟨rec4 m n M (gMask M) A B (res M u) (res M v), ?_, ?_⟩
  · unfold recs4
    apply mem_dedupL
    rw [List.mem_flatMap]
    refine ⟨res M u, List.mem_range.mpr (res_lt hM u), ?_⟩
    rw [List.mem_filterMap]
    exact ⟨res M v, List.mem_range.mpr (res_lt hM v), if_pos hc⟩
  · refine ⟨List.mem_filter.mpr ⟨hd0, ?_⟩, List.mem_filter.mpr ⟨hd8, ?_⟩,
      List.mem_filter.mpr ⟨hdL, ?_⟩, List.mem_filter.mpr ⟨hdU, ?_⟩⟩
    · show gOkN M (gMask M) d0 _ = true
      rw [gOkN_eq hM]
      exact gOk_of_GSqP hM g0 a0 b0
    · show gOkN M (gMask M) d8 _ = true
      rw [gOkN_eq hM]
      exact gOk_of_GSqP hM g8 a8 b8
    · show gOkN M (gMask M) dL _ = true
      rw [gOkN_eq hM]
      exact gOk_of_GSqP hM gL aL bL
    · show gOkN M (gMask M) dU _ = true
      rw [gOkN_eq hM]
      exact gOk_of_GSqP hM gU aU bU

/-- The elements of `l` that are in `l'`. -/
def inter (l l' : List ℕ) : List ℕ := l.filter fun x => l'.contains x

theorem mem_inter {x : ℕ} {l l' : List ℕ} (h : x ∈ l) (h' : x ∈ l') : x ∈ inter l l' := by
  unfold inter
  rw [List.mem_filter]
  exact ⟨h, by simpa using h'⟩

theorem isEmpty_false_of_mem {x : ℕ} {l : List ℕ} (h : x ∈ l) : l.isEmpty = false := by
  cases l with
  | nil => simp at h
  | cons a t => rfl

def inter4 (b r : Box) : Box :=
  (inter b.1 r.1, inter b.2.1 r.2.1, inter b.2.2.1 r.2.2.1, inter b.2.2.2 r.2.2.2)

def nonempty4 (b : Box) : Bool :=
  !b.1.isEmpty && !b.2.1.isEmpty && !b.2.2.1.isEmpty && !b.2.2.2.isEmpty

/-- Intersect every box with every box of `rs`, and keep the results in which no list is empty. -/
def step4 (boxes rs : List Box) : List Box :=
  boxes.flatMap fun b => rs.filterMap fun r =>
    if nonempty4 (inter4 b r) then some (inter4 b r) else none

theorem exists_step4 {t : ℕ × ℕ × ℕ × ℕ} {boxes rs : List Box} (hb : ∃ b ∈ boxes, In4 t b)
    (hr : ∃ r ∈ rs, In4 t r) : ∃ c ∈ step4 boxes rs, In4 t c := by
  obtain ⟨b, hb, h1, h2, h3, h4⟩ := hb
  obtain ⟨r, hr, k1, k2, k3, k4⟩ := hr
  have hin : In4 t (inter4 b r) :=
    ⟨mem_inter h1 k1, mem_inter h2 k2, mem_inter h3 k3, mem_inter h4 k4⟩
  have hne : nonempty4 (inter4 b r) = true := by
    obtain ⟨i1, i2, i3, i4⟩ := hin
    unfold nonempty4
    rw [isEmpty_false_of_mem i1, isEmpty_false_of_mem i2, isEmpty_false_of_mem i3,
      isEmpty_false_of_mem i4]
    rfl
  refine ⟨inter4 b r, ?_, hin⟩
  unfold step4
  rw [List.mem_flatMap]
  refine ⟨b, hb, ?_⟩
  rw [List.mem_filterMap]
  exact ⟨r, hr, if_pos hne⟩

/-- One line `(p, k)` of a certificate with four classes: intersect the boxes with the boxes of
the points modulo `p^k`. A line that cannot be read changes nothing. -/
def line4 (masks : List (ℕ × ℕ × ℕ)) (m n : ℤ) (A B : List ℕ) (boxes : List Box) (c : ℕ × ℕ) :
    List Box :=
  match findMask (c.1 ^ c.2) masks with
  | some (S, G) =>
    if isPrime c.1 && decide (0 < c.2) then step4 boxes (recs4 m n (c.1 ^ c.2) c.1 S G A B)
    else boxes
  | none => boxes

theorem exists_line4 {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks) {m n : ℤ}
    {u v : ℤ} (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    {d0 d8 dL dU : ℕ} (g0 : GSqP d0 (P0 m n u v)) (g8 : GSqP d8 (P8 m n u v))
    (gL : GSqP dL (PL m n u v)) (gU : GSqP dU (PU m n u v))
    {A B : List ℕ} (hd0 : d0 ∈ A) (hd8 : d8 ∈ A) (hdL : dL ∈ B) (hdU : dU ∈ B)
    (boxes : List Box) (c : ℕ × ℕ) (hb : ∃ b ∈ boxes, In4 (d0, d8, dL, dU) b) :
    ∃ b ∈ line4 masks m n A B boxes c, In4 (d0, d8, dL, dU) b := by
  unfold line4
  split
  · rename_i S G hfind
    split
    · rename_i hc
      simp only [Bool.and_eq_true, decide_eq_true_eq] at hc
      obtain ⟨hS, hG⟩ : S = sMask (c.1 ^ c.2) ∧ G = gMask (c.1 ^ c.2) :=
        hmasks _ (mem_of_findMask hfind)
      subst hS hG
      exact exists_step4 hb (exists_rec4 (prime_of_isPrime hc.1) rfl hc.2 hcop h0 h8 hL hU
        g0 g8 gL gU hd0 hd8 hdL hdU)
    · exact hb
  · exact hb

theorem exists_foldl_line4 {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks) {m n : ℤ}
    {u v : ℤ} (hcop : IsCoprime u v)
    (h0 : ∃ w, nrm (P0 m n u v) = w ^ 2) (h8 : ∃ w, nrm (P8 m n u v) = w ^ 2)
    (hL : ∃ w, nrm (PL m n u v) = w ^ 2) (hU : ∃ w, nrm (PU m n u v) = w ^ 2)
    {d0 d8 dL dU : ℕ} (g0 : GSqP d0 (P0 m n u v)) (g8 : GSqP d8 (P8 m n u v))
    (gL : GSqP dL (PL m n u v)) (gU : GSqP dU (PU m n u v))
    {A B : List ℕ} (hd0 : d0 ∈ A) (hd8 : d8 ∈ A) (hdL : dL ∈ B) (hdU : dU ∈ B) :
    ∀ (cs : List (ℕ × ℕ)) (boxes : List Box), (∃ b ∈ boxes, In4 (d0, d8, dL, dU) b) →
      ∃ b ∈ cs.foldl (line4 masks m n A B) boxes, In4 (d0, d8, dL, dU) b := by
  intro cs
  induction cs with
  | nil =>
    intro boxes hb
    exact hb
  | cons c cs ih =>
    intro boxes hb
    rw [List.foldl_cons]
    exact ih _ (exists_line4 hmasks hcop h0 h8 hL hU g0 g8 gL gU hd0 hd8 hdL hdU boxes c hb)

/-- An entry with four classes: `(m, n, ps, lines)`, with `m, n` coprime; `ps` lists primes, among
them every prime factor of `2mn(m − n)(m + n)`; after the lines no box of candidates is left. -/
def entryOk4 (masks : List (ℕ × ℕ × ℕ)) (e : ℕ × ℕ × List ℕ × List (ℕ × ℕ)) : Bool :=
  decide (0 < e.2.1) && decide (e.2.1 < e.1) && decide (Nat.Coprime e.1 e.2.1) &&
    e.2.2.1.all isPrime &&
    decide (2 * e.1 * e.2.1 * (e.1 - e.2.1) * (e.1 + e.2.1) ∣ e.2.2.1.prod ^ 64) &&
    (e.2.2.2.foldl
      (line4 masks (e.1 : ℤ) (e.2.1 : ℤ) (candA e.1 e.2.1 e.2.2.1) (candB e.1 e.2.1 e.2.2.1))
      [(candA e.1 e.2.1 e.2.2.1, candA e.1 e.2.1 e.2.2.1, candB e.1 e.2.1 e.2.2.1,
        candB e.1 e.2.1 e.2.2.1)]).isEmpty

/-- **No point, from a certificate with four classes.** -/
theorem no_point_of_entry4 {masks : List (ℕ × ℕ × ℕ)} (hmasks : MasksOk masks)
    {e : ℕ × ℕ × List ℕ × List (ℕ × ℕ)} (he : entryOk4 masks e = true) {u v : ℤ}
    (hcop : IsCoprime u v) (hv : v ≠ 0) (huv : u + v ≠ 0)
    (h0 : ∃ w, nrm (P0 e.1 e.2.1 u v) = w ^ 2) (h8 : ∃ w, nrm (P8 e.1 e.2.1 u v) = w ^ 2)
    (hL : ∃ w, nrm (PL e.1 e.2.1 u v) = w ^ 2) (hU : ∃ w, nrm (PU e.1 e.2.1 u v) = w ^ 2) : False := by
  obtain ⟨m, n, ps, cert⟩ := e
  unfold entryOk4 at he
  simp only [Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at he
  obtain ⟨⟨⟨⟨⟨hn, hnm⟩, hmn⟩, hps⟩, hdvd⟩, hemp⟩ := he
  dsimp only at h0 h8 hL hU hemp
  have hps' : ∀ p ∈ ps, p.Prime := fun p hp => prime_of_isPrime (hps p hp)
  have hm0 : (m : ℤ) ≠ 0 := by
    have : 0 < m := by omega
    exact_mod_cast this.ne'
  have hn0 : (n : ℤ) ≠ 0 := by exact_mod_cast hn.ne'
  have hmn0 : (m : ℤ) - n ≠ 0 := by
    have : (n : ℤ) < m := by exact_mod_cast hnm
    intro h
    omega
  have hbad : ∀ q : ℕ, q.Prime → (q : ℤ) ∣ 2 * (m : ℤ) * n * (m - n) * (m + n) → q ∈ ps := by
    intro q hq hdiv
    have e1 : (2 * (m : ℤ) * n * (m - n) * (m + n)) = ((2 * m * n * (m - n) * (m + n) : ℕ) : ℤ) := by
      push_cast [Nat.cast_sub hnm.le]
      ring
    rw [e1] at hdiv
    have h1 : q ∣ 2 * m * n * (m - n) * (m + n) := Int.natCast_dvd_natCast.mp hdiv
    exact mem_of_prime_dvd_prod hps' hq (hq.dvd_of_dvd_pow (dvd_trans h1 hdvd))
  obtain ⟨d0, d8, dL, dU, ⟨hd0, hs0, g0, t0⟩, ⟨hd8, hs8, g8, t8⟩, ⟨hdL, hsL, gL, tL⟩,
      ⟨hdU, hsU, gU, tU⟩⟩ := exists_classes4 (Nat.isCoprime_iff_coprime.mpr hmn) hcop h0 h8 hL hU
    (Or.inl (by dsimp only [P0]; exact mul_ne_zero hm0 hv))
    (Or.inl (by dsimp only [P8]; exact mul_ne_zero hn0 hv))
    (Or.inr (by dsimp only [PL]; exact mul_ne_zero hmn0 huv))
    (Or.inr (by
      dsimp only [PU]
      refine mul_ne_zero ?_ huv
      have : (0 : ℤ) < (m : ℤ) + n := by positivity
      exact this.ne'))
  have hA : ∀ d : ℕ, 0 < d → NoSq d → (∀ p : ℕ, p.Prime → p ∣ d →
      (p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧ (p : ℤ) ∣ 2 * ((m : ℤ) * n)) → d ∈ candA m n ps := by
    intro d hd hns ht
    unfold candA
    refine mem_prods _ (fun p hp => hps' p (List.mem_filter.mp hp).1) d hd hns (fun q hq hqd => ?_)
    obtain ⟨h8', hdv⟩ := ht q hq hqd
    refine List.mem_filter.mpr ⟨hbad q hq ?_, ?_⟩
    · have e2 : 2 * (m : ℤ) * n * (m - n) * (m + n) = (2 * ((m : ℤ) * n)) * ((m - n) * (m + n)) := by
        ring
      rw [e2]
      exact Dvd.dvd.mul_right hdv _
    · have hnat : q ∣ 2 * (m * n) := by exact_mod_cast hdv
      rw [Bool.and_eq_true, decide_eq_true_eq]
      refine ⟨?_, hnat⟩
      unfold keep8
      rcases h8' with h | h | h <;> simp [h]
  have hB : ∀ d : ℕ, 0 < d → NoSq d → (∀ p : ℕ, p.Prime → p ∣ d →
      (p = 2 ∨ p % 8 = 1 ∨ p % 8 = 7) ∧ (p : ℤ) ∣ 2 * (((m : ℤ) + n) * (m - n))) →
      d ∈ candB m n ps := by
    intro d hd hns ht
    unfold candB
    refine mem_prods _ (fun p hp => hps' p (List.mem_filter.mp hp).1) d hd hns (fun q hq hqd => ?_)
    obtain ⟨h8', hdv⟩ := ht q hq hqd
    refine List.mem_filter.mpr ⟨hbad q hq ?_, ?_⟩
    · have e2 : 2 * (m : ℤ) * n * (m - n) * (m + n) = (2 * (((m : ℤ) + n) * (m - n))) * (m * n) := by
        ring
      rw [e2]
      exact Dvd.dvd.mul_right hdv _
    · have e3 : ((2 * ((m + n) * (m - n)) : ℕ) : ℤ) = 2 * (((m : ℤ) + n) * (m - n)) := by
        push_cast [Nat.cast_sub hnm.le]
        ring
      have hnat : q ∣ 2 * ((m + n) * (m - n)) := by
        rw [← e3] at hdv
        exact Int.natCast_dvd_natCast.mp hdv
      rw [Bool.and_eq_true, decide_eq_true_eq]
      refine ⟨?_, hnat⟩
      unfold keep8
      rcases h8' with h | h | h <;> simp [h]
  have hfin := exists_foldl_line4 hmasks hcop h0 h8 hL hU g0 g8 gL gU (hA d0 hd0 hs0 t0)
    (hA d8 hd8 hs8 t8) (hB dL hdL hsL tL) (hB dU hdU hsU tU) cert
    [(candA m n ps, candA m n ps, candB m n ps, candB m n ps)]
    ⟨(candA m n ps, candA m n ps, candB m n ps, candB m n ps), List.mem_singleton.mpr rfl,
      hA d0 hd0 hs0 t0, hA d8 hd8 hs8 t8, hB dL hdL hsL tL, hB dU hdU hsU tU⟩
  rw [List.isEmpty_iff] at hemp
  rw [hemp] at hfin
  obtain ⟨b, hb, _⟩ := hfin
  simp at hb

end MagicSquareOfSquares.CoverDescent
