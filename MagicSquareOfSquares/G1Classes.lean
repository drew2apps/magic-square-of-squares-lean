import MagicSquareOfSquares.G1

/-!
# The G1 family: whole remainder classes with no rational point

`G1.lean` gives a check that rules out one pair `(m, n)` at a time. Because that check only
depends on the remainders of `m` and `n`, a finite table settles every pair in a remainder
class, of any size. This file records the classes found so far.

## Main results

* `no_point_of_two_mod_four`: one of `m`, `n` is `2` modulo `4` and the other is odd
  (a table modulo `8`).
* `no_point_of_odd`: both odd, one is `±1` and the other `±3` modulo `8` (a table modulo `32`).
* `no_point_mod_seven`, `no_point_mod_twenty_three`: `m² ≡ 4n²` or `4m² ≡ n²` modulo `7` or
  `23`, with `m`, `n` not both divisible by that prime.
* `no_point_mod_forty_seven`: `m² ≡ 25n²` or `25m² ≡ n²` modulo `47`, likewise.

Every table is evaluated by the Lean kernel (`decide +kernel`); no `native_decide`.
-/

namespace MagicSquareOfSquares.G1

/-- `classTable M p S c`: for every pair of remainders `(a, b)` modulo `M` that satisfies `c`,
the check `obstructsProj a b M p S` succeeds. -/
def classTable (M p S : ℕ) (c : ℕ → ℕ → Bool) : Bool :=
  (List.range M).all fun a => (List.range M).all fun b => !c a b || obstructsProj a b M p S

private theorem of_classTable {m n M p S : ℕ} {c : ℕ → ℕ → Bool} (hp : p.Prime) (hM1 : 1 < M)
    (hpM : p ∣ M) (hcopM : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x M) (hS : S = sqMask M)
    (htab : classTable M p S c = true) (hc : c (m % M) (n % M) = true) (t : ℚ) :
    ¬DiagonalSquares m n t := by
  have hM : 0 < M := by omega
  unfold classTable at htab
  have h1 := List.all_eq_true.mp htab (m % M) (List.mem_range.mpr (Nat.mod_lt _ hM))
  have h2 := List.all_eq_true.mp h1 (n % M) (List.mem_range.mpr (Nat.mod_lt _ hM))
  rw [hc] at h2
  exact no_point_of_obstructsProj_mod hp hM1 hpM hcopM hS (by simpa using h2) t

/-! ## The prime 2 -/

/-- One of the two remainders is `2` modulo `4` and the other is odd. -/
def cond8 (a b : ℕ) : Bool :=
  decide ((a % 4 = 2 ∧ b % 2 = 1) ∨ (b % 4 = 2 ∧ a % 2 = 1))

/-- Both odd, one is `±1` and the other `±3` modulo `8`. -/
def cond32 (a b : ℕ) : Bool :=
  decide (((a % 8 = 1 ∨ a % 8 = 7) ∧ (b % 8 = 3 ∨ b % 8 = 5)) ∨
    ((a % 8 = 3 ∨ a % 8 = 5) ∧ (b % 8 = 1 ∨ b % 8 = 7)))

theorem table8 : classTable 8 2 19 cond8 = true := by decide +kernel

theorem table32 : classTable 32 2 33751571 cond32 = true := by decide +kernel

/-- If one of `m`, `n` is `2` modulo `4` and the other is odd, no rational `t` makes the three
diagonal cells of the G1 family squares. -/
theorem no_point_of_two_mod_four {m n : ℕ}
    (h : (m % 4 = 2 ∧ n % 2 = 1) ∨ (n % 4 = 2 ∧ m % 2 = 1)) (t : ℚ) :
    ¬DiagonalSquares m n t :=
  of_classTable Nat.prime_two (by norm_num) (by norm_num)
    (fun x hx => coprime_pow_of_not_dvd Nat.prime_two 3 x hx) (by decide +kernel) table8
    (by simp only [cond8, decide_eq_true_eq]; omega) t

/-- If `m` and `n` are odd, one of them `±1` and the other `±3` modulo `8`, no rational `t`
makes the three diagonal cells of the G1 family squares. -/
theorem no_point_of_odd {m n : ℕ}
    (h : ((m % 8 = 1 ∨ m % 8 = 7) ∧ (n % 8 = 3 ∨ n % 8 = 5)) ∨
      ((m % 8 = 3 ∨ m % 8 = 5) ∧ (n % 8 = 1 ∨ n % 8 = 7))) (t : ℚ) :
    ¬DiagonalSquares m n t :=
  of_classTable Nat.prime_two (by norm_num) (by norm_num)
    (fun x hx => coprime_pow_of_not_dvd Nat.prime_two 5 x hx) (by decide +kernel) table32
    (by simp only [cond32, decide_eq_true_eq]; omega) t

/-! ## The primes 7, 23 and 47 -/

/-- `a² ≡ q·b²` or `q·a² ≡ b²` modulo `M`, with `a`, `b` not both `0` modulo `M`. -/
def condRatio (M q a b : ℕ) : Bool :=
  decide ((a * a % M = q * (b * b) % M ∨ q * (a * a) % M = b * b % M) ∧
    ¬(a % M = 0 ∧ b % M = 0))

theorem table7 : classTable 7 7 23 (condRatio 7 4) = true := by decide +kernel

theorem table23 : classTable 23 23 340831 (condRatio 23 4) = true := by decide +kernel

theorem table47 : classTable 47 47 4626135339999 (condRatio 47 25) = true := by decide +kernel

private theorem of_ratio {m n M p S q : ℕ} (hp : p.Prime) (hM1 : 1 < M) (hpM : p ∣ M)
    (hcopM : ∀ x : ℕ, ¬p ∣ x → Nat.Coprime x M) (hS : S = sqMask M)
    (htab : classTable M p S (condRatio M q) = true)
    (h : m * m % M = q * (n * n) % M ∨ q * (m * m) % M = n * n % M)
    (h0 : ¬(m % M = 0 ∧ n % M = 0)) (t : ℚ) : ¬DiagonalSquares m n t := by
  have em : m % M * (m % M) ≡ m * m [MOD M] := (Nat.mod_modEq m M).mul (Nat.mod_modEq m M)
  have en : n % M * (n % M) ≡ n * n [MOD M] := (Nat.mod_modEq n M).mul (Nat.mod_modEq n M)
  have e1 : m % M * (m % M) % M = m * m % M := em
  have e2 : n % M * (n % M) % M = n * n % M := en
  have e3 : q * (m % M * (m % M)) % M = q * (m * m) % M := (Nat.ModEq.refl q).mul em
  have e4 : q * (n % M * (n % M)) % M = q * (n * n) % M := (Nat.ModEq.refl q).mul en
  refine of_classTable hp hM1 hpM hcopM hS htab ?_ t
  simp only [condRatio, decide_eq_true_eq, Nat.mod_mod]
  rw [e1, e2, e3, e4]
  exact ⟨h, h0⟩

/-- If `m² ≡ 4n²` or `4m² ≡ n²` modulo `7`, and `7` does not divide both, no rational `t`
makes the three diagonal cells of the G1 family squares. -/
theorem no_point_mod_seven {m n : ℕ}
    (h : m * m % 7 = 4 * (n * n) % 7 ∨ 4 * (m * m) % 7 = n * n % 7)
    (h0 : ¬(m % 7 = 0 ∧ n % 7 = 0)) (t : ℚ) : ¬DiagonalSquares m n t :=
  of_ratio (p := 7) (by norm_num) (by norm_num) (dvd_refl 7)
    (fun x hx => coprime_pow_of_not_dvd (by norm_num) 1 x hx) (by decide +kernel) table7 h h0 t

/-- The same statement modulo `23`. -/
theorem no_point_mod_twenty_three {m n : ℕ}
    (h : m * m % 23 = 4 * (n * n) % 23 ∨ 4 * (m * m) % 23 = n * n % 23)
    (h0 : ¬(m % 23 = 0 ∧ n % 23 = 0)) (t : ℚ) : ¬DiagonalSquares m n t :=
  of_ratio (p := 23) (by norm_num) (by norm_num) (dvd_refl 23)
    (fun x hx => coprime_pow_of_not_dvd (by norm_num) 1 x hx) (by decide +kernel) table23 h h0 t

/-- If `m² ≡ 25n²` or `25m² ≡ n²` modulo `47`, and `47` does not divide both, no rational `t`
makes the three diagonal cells of the G1 family squares. -/
theorem no_point_mod_forty_seven {m n : ℕ}
    (h : m * m % 47 = 25 * (n * n) % 47 ∨ 25 * (m * m) % 47 = n * n % 47)
    (h0 : ¬(m % 47 = 0 ∧ n % 47 = 0)) (t : ℚ) : ¬DiagonalSquares m n t :=
  of_ratio (p := 47) (by norm_num) (by norm_num) (dvd_refl 47)
    (fun x hx => coprime_pow_of_not_dvd (by norm_num) 1 x hx) (by decide +kernel) table47 h h0 t

end MagicSquareOfSquares.G1
