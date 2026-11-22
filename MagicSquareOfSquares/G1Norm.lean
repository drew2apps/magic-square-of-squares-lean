import MagicSquareOfSquares.G1Bridge
import MagicSquareOfSquares.G1Primes
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic

/-!
# Conditions on the number `m² + n²`

The ratio `m : n` of the G1 family enters the congruence obstructions mostly through the single
number `m² + n²`. This file records three facts.

* `G1.isSquare_norm_mod`: if the five forms are squares for some rational `t`, then `m² + n²` is a
  square modulo `7`, `17`, `31` and `47`. For these four primes *every* class `m : n` with
  `m² + n²` a non-residue is obstructed; the check is evaluated by the Lean kernel.
* `IsMSS.isSquare_ratio_norm`, `IsMSS.even_padicValNat_ratio_norm`, `IsMSS.edge_norm`: the same, and the statement about
  primes `≡ 5 (mod 8)` from `G1Primes.lean`, read as conditions on the edge roots of a magic
  square of squares: the number `(b + d)² + (h + f)²`, divided by the square of
  `gcd(b + d, h + f)`, is a square modulo `7`, `17`, `31`, `47`, and every prime `≡ 5 (mod 8)`
  divides `(b + d)² + (h + f)²` to an even power.
* `G1.exists_point_mod_of_not_isSquare_two`: if `2` is not a square modulo `p` then for every
  `m`, `n` the five forms are squares modulo `p` at `(u, v) = (1, 1)` or at `(u, v) = (0, 1)`.
  So a congruence modulo a prime `p ≡ ±3 (mod 8)` can never exclude a ratio. This explains why
  all 34 obstructing primes are `±1 (mod 8)`.

The values of `(b + d)² + (h + f)²` at `(1, 1)` and `(0, 1)` are `m² + n²` and `2(m² + n²)` up to
squares, so an obstruction at `p` needs both to be non-residues.
-/

namespace MagicSquareOfSquares.G1

/-- `normOk (p, S)`: `S` is the mask of squares modulo `p`, and every class `r` modulo `p` for
which `r² + 1` is not a square modulo `p` is obstructed. -/
def normOk (e : ℕ × ℕ) : Bool :=
  decide (e.2 = sqMask e.1) &&
    (List.range e.1).all fun r =>
      ((List.range e.1).any fun x => x * x % e.1 == (r * r + 1) % e.1) ||
        obstructsFull r 1 e.1 e.1 e.2

/-- The four primes at which every non-residue class is obstructed, with their masks. -/
def normPrimes : List (ℕ × ℕ) :=
  [(7, 23), (17, 107287), (31, 303908791), (47, 4626135339999)]

/-- Evaluated by the Lean kernel. -/
theorem normPrimes_valid : normPrimes.all normOk = true := by decide +kernel

theorem normPrimes_prime : ∀ e ∈ normPrimes, e.1.Prime := by
  intro e he
  simp only [normPrimes, List.mem_cons, List.not_mem_nil, or_false] at he
  rcases he with rfl | rfl | rfl | rfl <;> norm_num

/-- If every non-residue class modulo the prime `p` is obstructed, a solution of the five-form
system has `m² + n²` a square modulo `p`. -/
theorem isSquare_norm_of_normOk {p S : ℕ} (hp : p.Prime) (hok : normOk (p, S) = true)
    {m n : ℕ} {t : ℚ} (h : FullSquares m n t) : IsSquare (((m ^ 2 + n ^ 2 : ℕ) : ZMod p)) := by
  have : Fact p.Prime := ⟨hp⟩
  have : NeZero p := ⟨hp.ne_zero⟩
  obtain ⟨hS, hchk⟩ := Bool.and_eq_true_iff.mp hok
  have hS' : S = sqMask p := of_decide_eq_true hS
  by_cases hn : (n : ZMod p) = 0
  · refine ⟨(m : ZMod p), ?_⟩
    push_cast
    rw [hn]
    ring
  · obtain ⟨r, hr⟩ : ∃ r : ℕ, r = ((m : ZMod p) * (n : ZMod p)⁻¹).val := ⟨_, rfl⟩
    have hrlt : r < p := hr ▸ ZMod.val_lt _
    have hrcast : (r : ZMod p) = (m : ZMod p) * (n : ZMod p)⁻¹ := by
      rw [hr, ZMod.natCast_zmod_val]
    have hm : (m : ZMod p) = (r : ZMod p) * (n : ZMod p) := by
      rw [hrcast, mul_assoc, inv_mul_cancel₀ hn, mul_one]
    have hmn : m % p = r * n % p := by
      refine (ZMod.natCast_eq_natCast_iff' m (r * n) p).mp ?_
      rw [Nat.cast_mul]
      exact hm
    have hnp : n % p ≠ 0 := fun h0 =>
      hn ((CharP.cast_eq_zero_iff (ZMod p) p n).mpr (Nat.dvd_of_mod_eq_zero h0))
    have hall := List.all_eq_true.mp hchk r (List.mem_range.mpr hrlt)
    rcases Bool.or_eq_true_iff.mp hall with hx | hob
    · obtain ⟨x, -, hx⟩ := List.any_eq_true.mp hx
      have hx' : x * x % p = (r * r + 1) % p := by simpa using hx
      have hz : ((x * x : ℕ) : ZMod p) = ((r * r + 1 : ℕ) : ZMod p) :=
        (ZMod.natCast_eq_natCast_iff' _ _ p).mpr hx'
      push_cast at hz
      refine ⟨(n : ZMod p) * x, ?_⟩
      push_cast
      rw [hm]
      linear_combination (-(n : ZMod p) ^ 2) * hz
    · exact absurd h (no_solution_of_ratio hp hS' hob hmn hnp t)

/-- **If the five forms of the G1 family are squares for some rational `t`, then `m² + n²` is a
square modulo `7`, `17`, `31` and `47`.** -/
theorem isSquare_norm_mod {m n : ℕ} {t : ℚ} (h : FullSquares m n t) {p : ℕ}
    (hp : p = 7 ∨ p = 17 ∨ p = 31 ∨ p = 47) : IsSquare (((m ^ 2 + n ^ 2 : ℕ) : ZMod p)) := by
  have hv := List.all_eq_true.mp normPrimes_valid
  rcases hp with rfl | rfl | rfl | rfl
  · exact isSquare_norm_of_normOk (normPrimes_prime (7, 23) (by simp [normPrimes]))
      (hv (7, 23) (by simp [normPrimes])) h
  · exact isSquare_norm_of_normOk (normPrimes_prime (17, 107287) (by simp [normPrimes]))
      (hv (17, 107287) (by simp [normPrimes])) h
  · exact isSquare_norm_of_normOk (normPrimes_prime (31, 303908791) (by simp [normPrimes]))
      (hv (31, 303908791) (by simp [normPrimes])) h
  · exact isSquare_norm_of_normOk (normPrimes_prime (47, 4626135339999) (by simp [normPrimes]))
      (hv (47, 4626135339999) (by simp [normPrimes])) h

/-- **Why only primes `±1 (mod 8)` obstruct.** If `2` is not a square modulo the prime `p`, then
for all `m`, `n` the five forms are squares modulo `p` at `(u, v) = (1, 1)` or at
`(u, v) = (0, 1)`: the first works when `m² + n²` is a square, the second when it is not. -/
theorem exists_point_mod_of_not_isSquare_two {p : ℕ} [Fact p.Prime]
    (h2 : ¬IsSquare (2 : ZMod p)) (m n : ZMod p) :
    ∃ u v : ZMod p, (u ≠ 0 ∨ v ≠ 0) ∧
      IsSquare (m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2) ∧
      IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2)) ∧
      IsSquare (n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2) ∧
      IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) - 8 * m * n * u * v) ∧
      IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) + 8 * m * n * u * v) := by
  by_cases hN : IsSquare (m ^ 2 + n ^ 2)
  · obtain ⟨s, hs⟩ := hN
    exact ⟨1, 1, Or.inl one_ne_zero, ⟨s, by linear_combination hs⟩,
      ⟨2 * s, by linear_combination 4 * hs⟩, ⟨s, by linear_combination hs⟩,
      ⟨2 * (m - n), by ring⟩, ⟨2 * (m + n), by ring⟩⟩
  · have hχ2 : quadraticChar (ZMod p) 2 = -1 := quadraticChar_neg_one_iff_not_isSquare.mpr h2
    have hχN : quadraticChar (ZMod p) (m ^ 2 + n ^ 2) = -1 :=
      quadraticChar_neg_one_iff_not_isSquare.mpr hN
    have hχ : quadraticChar (ZMod p) (2 * (m ^ 2 + n ^ 2)) = 1 := by
      rw [map_mul, hχ2, hχN]
      norm_num
    have hne : 2 * (m ^ 2 + n ^ 2) ≠ 0 := by
      intro h0
      rw [h0, quadraticChar_zero] at hχ
      exact zero_ne_one hχ
    obtain ⟨s, hs⟩ := (quadraticChar_one_iff_isSquare hne).mp hχ
    exact ⟨0, 1, Or.inr one_ne_zero, ⟨m, by ring⟩, ⟨s, by linear_combination hs⟩, ⟨n, by ring⟩,
      ⟨s, by linear_combination hs⟩, ⟨s, by linear_combination hs⟩⟩

/-- The same for a prime `p ≡ 3` or `5 (mod 8)`. -/
theorem exists_point_mod_of_three_five_mod_eight {p : ℕ} [Fact p.Prime]
    (hp8 : p % 8 = 3 ∨ p % 8 = 5) (m n : ZMod p) :
    ∃ u v : ZMod p, (u ≠ 0 ∨ v ≠ 0) ∧
      IsSquare (m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2) ∧
      IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2)) ∧
      IsSquare (n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2) ∧
      IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) - 8 * m * n * u * v) ∧
      IsSquare (2 * (m ^ 2 + n ^ 2) * (u ^ 2 + v ^ 2) + 8 * m * n * u * v) := by
  refine exists_point_mod_of_not_isSquare_two ?_ m n
  rw [ZMod.exists_sq_eq_two_iff (by omega)]
  omega

end MagicSquareOfSquares.G1

namespace MagicSquareOfSquares

open G1

/-- **The congruence obstruction as a statement about all magic squares of squares.** If the
class `r` modulo the prime `p` is obstructed (as are the 484 classes of
`G1.fullTable`), then no magic square of squares has `(b + d) : (h + f) = m : n` with
`m ≡ r·n (mod p)`. There is no bound on `m` and `n`. -/
theorem IsMSS.ratio_not_in_class {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hn : 0 < n) (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) (hbd : b ≠ d)
    {p r S : ℕ} (hp : p.Prime) (hS : S = sqMask p) (hob : obstructsFull r 1 p p S = true)
    (hnp : n % p ≠ 0) : m % p ≠ r * n % p :=
  fun hmnp => no_solution_of_ratio hp hS hob hmnp hnp _ (hs.fullSquares hn hmn hne hbd)

/-- In a magic square of squares with `(b + d)·n = (h + f)·m`, the number `m² + n²` is a square
modulo `7`, `17`, `31` and `47`. -/
theorem IsMSS.isSquare_ratio_norm {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    {m n : ℕ} (hn : 0 < n) (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) (hbd : b ≠ d)
    {p : ℕ} (hp : p = 7 ∨ p = 17 ∨ p = 31 ∨ p = 47) :
    IsSquare (((m ^ 2 + n ^ 2 : ℕ) : ZMod p)) :=
  isSquare_norm_mod (hs.fullSquares hn hmn hne hbd) hp

/-- In a magic square of squares with `(b + d)·n = (h + f)·m` and `m`, `n` coprime, every prime
`p ≡ 5 (mod 8)` divides `m² + n²` to an even power. -/
theorem IsMSS.even_padicValNat_ratio_norm {a b c d e f g h i T : ℤ}
    (hs : IsMSS a b c d e f g h i T) {m n : ℕ} (hn : 0 < n) (hcop : Nat.Coprime m n)
    (hmn : (b + d) * n = (h + f) * m) (hne : h + f ≠ 0) (hbd : b ≠ d)
    {p : ℕ} (hp : p.Prime) (hp8 : p % 8 = 5) : Even (padicValNat p (m * m + n * n)) := by
  by_contra hodd
  exact no_point_of_five_mod_eight_of_coprime hp hp8 (Nat.not_even_iff_odd.mp hodd) hcop _
    (hs.fullSquares hn hmn hne hbd).diagonal

/-- The signs of the edge roots can be chosen so that `b + d ≥ 0` and `h + f ≥ 0`. -/
private theorem IsMSS.exists_abs {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hbd : b ≠ d) :
    ∃ b' d' f' h' : ℤ, IsMSS a b' c d' e f' g h' i T ∧ b' + d' = |b + d| ∧
      h' + f' = |h + f| ∧ b' ≠ d' := by
  have nb : IsMSS a (-b) c (-d) e f g h i T := by simpa [IsMSS] using hs
  have nh : IsMSS a b c d e (-f) g (-h) i T := by simpa [IsMSS] using hs
  have nbh : IsMSS a (-b) c (-d) e (-f) g (-h) i T := by simpa [IsMSS] using hs
  have hbd' : -b ≠ -d := fun h0 => hbd (neg_injective h0)
  rcases le_total 0 (b + d) with h1 | h1 <;> rcases le_total 0 (h + f) with h2 | h2
  · exact ⟨b, d, f, h, hs, (abs_of_nonneg h1).symm, (abs_of_nonneg h2).symm, hbd⟩
  · exact ⟨b, d, -f, -h, nh, (abs_of_nonneg h1).symm, by rw [abs_of_nonpos h2]; ring, hbd⟩
  · exact ⟨-b, -d, f, h, nb, by rw [abs_of_nonpos h1]; ring, (abs_of_nonneg h2).symm, hbd'⟩
  · exact ⟨-b, -d, -f, -h, nbh, by rw [abs_of_nonpos h1]; ring, by rw [abs_of_nonpos h2]; ring,
      hbd'⟩

/-- **In a magic square of squares with roots `a b c / d e f / g h i`, write
`b + d = ±g·m`, `h + f = ±g·n` with `g = gcd(b + d, h + f)`. Then `m² + n²` is a square modulo
`7`, `17`, `31` and `47`, and every prime `p ≡ 5 (mod 8)` divides it to an even power.** The
hypotheses are `b ≠ d` and `h + f ≠ 0`; no distinctness is assumed and the signs of the roots are
arbitrary. -/
theorem IsMSS.edge_norm {a b c d e f g h i T : ℤ} (hs : IsMSS a b c d e f g h i T)
    (hne : h + f ≠ 0) (hbd : b ≠ d) :
    ∃ m n : ℕ, Nat.Coprime m n ∧ (b + d).natAbs = m * Nat.gcd (b + d).natAbs (h + f).natAbs ∧
      (h + f).natAbs = n * Nat.gcd (b + d).natAbs (h + f).natAbs ∧
      (∀ p : ℕ, p = 7 ∨ p = 17 ∨ p = 31 ∨ p = 47 → IsSquare (((m ^ 2 + n ^ 2 : ℕ) : ZMod p))) ∧
      ∀ p : ℕ, p.Prime → p % 8 = 5 → Even (padicValNat p (m * m + n * n)) := by
  obtain ⟨b', d', f', h', hs', hX, hY, hbd'⟩ := hs.exists_abs hbd
  obtain ⟨m, n, hcop, hm, hn⟩ := Nat.exists_coprime (b + d).natAbs (h + f).natAbs
  have hY0 : 0 < (h + f).natAbs := Int.natAbs_pos.mpr hne
  have hg0 : 0 < Nat.gcd (b + d).natAbs (h + f).natAbs := Nat.gcd_pos_of_pos_right _ hY0
  have hn0 : 0 < n := by
    rcases Nat.eq_zero_or_pos n with h0 | h0
    · rw [h0, zero_mul] at hn
      omega
    · exact h0
  have hne' : h' + f' ≠ 0 := by
    rw [hY]
    exact abs_ne_zero.mpr hne
  have hmn : (b' + d') * (n : ℤ) = (h' + f') * (m : ℤ) := by
    rw [hX, hY, ← Int.natCast_natAbs, ← Int.natCast_natAbs]
    have e1 : (((b + d).natAbs : ℕ) : ℤ) =
        (m : ℤ) * (Nat.gcd (b + d).natAbs (h + f).natAbs : ℤ) := by exact_mod_cast hm
    have e2 : (((h + f).natAbs : ℕ) : ℤ) =
        (n : ℤ) * (Nat.gcd (b + d).natAbs (h + f).natAbs : ℤ) := by exact_mod_cast hn
    rw [e1, e2]
    ring
  exact ⟨m, n, hcop, hm, hn, fun p hp => hs'.isSquare_ratio_norm hn0 hmn hne' hbd' hp,
    fun p hp hp8 => hs'.even_padicValNat_ratio_norm hn0 hcop hmn hne' hbd' hp hp8⟩

end MagicSquareOfSquares
