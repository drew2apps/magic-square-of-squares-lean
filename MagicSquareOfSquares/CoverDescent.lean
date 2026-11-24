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

end MagicSquareOfSquares.CoverDescent
