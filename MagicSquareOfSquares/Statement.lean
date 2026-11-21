import MagicSquareOfSquares.Residues
import MagicSquareOfSquares.Progressions
import MagicSquareOfSquares.Primitive

/-!
# The open problem as stated in `formal-conjectures`

`OpenStatement` is the open problem in the form used by
`FormalConjectures/Wikipedia/MagicSquares.lean` in
[google-deepmind/formal-conjectures](https://github.com/google-deepmind/formal-conjectures):
a function `m : Fin 3 → Fin 3 → ℕ` of distinct positive squares with equal line sums.

This file restates the known necessary conditions for such an `m`, so that they apply
directly to any witness of that statement.

## Main results

* `IsMagicSquareOfSquares.sum_eq`: the magic sum is three times the centre.
* `IsMagicSquareOfSquares.odd`, `emod_24`, `sum_emod_72`: residues in the primitive case.
* `IsMagicSquareOfSquares.not_dvd_of_three_mod_eight`, `center_prime_factors`, `not_dvd_edge`,
  `dvd_edges_of_dvd_corner`, `dvd_center_of_dvd_corner`: prime divisors.
* `IsMagicSquareOfSquares.pos`: distinct entries are automatically positive.
* `openStatement_iff`: the positivity hypothesis in the open statement is redundant.
* `openStatement_iff_exists_isMSS`: the open statement in terms of nine integer roots.
* `openStatement_iff_three_progressions`: Robertson's reformulation of the open statement.
* `openStatement_iff_primitive`: if the open statement has a solution, it has a primitive one.
-/

namespace MagicSquareOfSquares

/-- The open problem: is there a 3×3 magic square of distinct positive squares?
This is the statement of `MagicSquares.exists_magic_square_squares` in `formal-conjectures`. -/
def OpenStatement : Prop :=
  ∃ m : Fin 3 → Fin 3 → ℕ, ∃ t : ℕ,
    m.Injective2 ∧
    (∀ i j, 0 < m i j ∧ IsSquare (m i j)) ∧
    (∀ i, ∑ j, m i j = t) ∧
    (∀ j, ∑ i, m i j = t) ∧
    m 0 0 + m 1 1 + m 2 2 = t ∧
    m 0 2 + m 1 1 + m 2 0 = t

/-- `m` is a 3×3 magic square of squares with sum `t`. Distinctness and positivity are not
assumed. -/
structure IsMagicSquareOfSquares (m : Fin 3 → Fin 3 → ℕ) (t : ℕ) : Prop where
  isSquare : ∀ i j, IsSquare (m i j)
  row : ∀ i, ∑ j, m i j = t
  col : ∀ j, ∑ i, m i j = t
  diag : m 0 0 + m 1 1 + m 2 2 = t
  anti : m 0 2 + m 1 1 + m 2 0 = t

/-- `m` is primitive if no prime divides all nine entries. -/
def IsPrimitive (m : Fin 3 → Fin 3 → ℕ) : Prop :=
  ∀ p : ℕ, p.Prime → ∃ i j, ¬p ∣ m i j

private theorem forall_of_nine {P : Fin 3 → Fin 3 → Prop}
    (h : P 0 0 ∧ P 0 1 ∧ P 0 2 ∧ P 1 0 ∧ P 1 1 ∧ P 1 2 ∧ P 2 0 ∧ P 2 1 ∧ P 2 2) :
    ∀ i j, P i j := by
  obtain ⟨h00, h01, h02, h10, h11, h12, h20, h21, h22⟩ := h
  intro i j
  fin_cases i <;> fin_cases j <;> assumption

private theorem nine_of_exists {P : Fin 3 → Fin 3 → Prop} (h : ∃ i j, P i j) :
    P 0 0 ∨ P 0 1 ∨ P 0 2 ∨ P 1 0 ∨ P 1 1 ∨ P 1 2 ∨ P 2 0 ∨ P 2 1 ∨ P 2 2 := by
  obtain ⟨i, j, hij⟩ := h
  fin_cases i <;> fin_cases j
  · exact Or.inl hij
  · exact Or.inr (Or.inl hij)
  · exact Or.inr (Or.inr (Or.inl hij))
  · exact Or.inr (Or.inr (Or.inr (Or.inl hij)))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hij))))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hij)))))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hij))))))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hij)))))))
  · exact Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr hij)))))))

private theorem cast_sum {x y z t : ℕ} (h : x * x + y * y + z * z = t) :
    (x : ℤ) ^ 2 + (y : ℤ) ^ 2 + (z : ℤ) ^ 2 = (t : ℤ) := by
  rw [← h]
  push_cast
  ring

private theorem dvd_iff {p x y : ℕ} (hp : p.Prime) (h : x = y ^ 2) :
    p ∣ x ↔ (p : ℤ) ∣ (y : ℤ) := by
  subst h
  rw [Int.natCast_dvd_natCast]
  exact ⟨fun h => hp.dvd_of_dvd_pow h, fun h => Dvd.dvd.pow h two_ne_zero⟩

private theorem natAbs_line {x y z e T : ℤ} (hT : T = 3 * e ^ 2)
    (h : x ^ 2 + y ^ 2 + z ^ 2 = T) :
    x.natAbs ^ 2 + y.natAbs ^ 2 + z.natAbs ^ 2 = 3 * e.natAbs ^ 2 := by
  zify
  simp only [sq_abs]
  linarith

namespace IsMagicSquareOfSquares

variable {m : Fin 3 → Fin 3 → ℕ} {t : ℕ}

/-- The entries have natural-number roots, and the roots satisfy `IsMSS`. -/
theorem exists_roots (H : IsMagicSquareOfSquares m t) :
    ∃ r : Fin 3 → Fin 3 → ℕ, (∀ i j, m i j = r i j ^ 2) ∧
      IsMSS (r 0 0) (r 0 1) (r 0 2) (r 1 0) (r 1 1) (r 1 2) (r 2 0) (r 2 1) (r 2 2) t := by
  choose r hr using H.isSquare
  refine ⟨r, fun i j => by rw [hr i j, sq], ?_⟩
  have row := H.row
  have col := H.col
  have dg := H.diag
  have an := H.anti
  simp only [Fin.sum_univ_three, hr] at row col dg an
  exact IsMagic.mk (cast_sum (row 0)) (cast_sum (row 1)) (cast_sum (row 2)) (cast_sum (col 0))
    (cast_sum (col 1)) (cast_sum (col 2)) (cast_sum dg) (cast_sum an)

private theorem primitive_roots {r : Fin 3 → Fin 3 → ℕ} (hr : ∀ i j, m i j = r i j ^ 2)
    (hprim : IsPrimitive m) :
    Primitive (r 0 0) (r 0 1) (r 0 2) (r 1 0) (r 1 1) (r 1 2) (r 2 0) (r 2 1) (r 2 2) := by
  intro p hp hall
  obtain ⟨i, j, hij⟩ := hprim p hp
  exact hij ((dvd_iff hp (hr i j)).mpr
    (forall_of_nine (P := fun i j => (p : ℤ) ∣ (r i j : ℤ)) hall i j))

/-- The magic sum is three times the centre entry. -/
theorem sum_eq (H : IsMagicSquareOfSquares m t) : t = 3 * m 1 1 := by
  have row := H.row
  have col := H.col
  simp only [Fin.sum_univ_three] at row col
  have r0 := row 0
  have r1 := row 1
  have r2 := row 2
  have c0 := col 0
  have c1 := col 1
  have c2 := col 2
  have dg := H.diag
  have an := H.anti
  omega

/-- In a primitive magic square of squares every entry is odd. -/
theorem odd (H : IsMagicSquareOfSquares m t) (hprim : IsPrimitive m) : ∀ i j, Odd (m i j) := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have key := forall_of_nine (P := fun i j => Odd (r i j : ℤ))
    (HR.odd (primitive_roots hr hprim))
  intro i j
  rw [hr i j]
  exact ((Int.odd_coe_nat _).mp (key i j)).pow

/-- In a primitive magic square of squares every entry is `≡ 1 (mod 24)`. -/
theorem emod_24 (H : IsMagicSquareOfSquares m t) (hprim : IsPrimitive m) :
    ∀ i j, m i j % 24 = 1 := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have key := forall_of_nine (P := fun i j => (r i j : ℤ) ^ 2 % 24 = 1)
    (HR.sq_emod_24 (primitive_roots hr hprim))
  intro i j
  have h : (r i j : ℤ) ^ 2 % 24 = 1 := key i j
  rw [hr i j]
  exact_mod_cast h

/-- In a primitive magic square of squares the magic sum is `≡ 3 (mod 72)`. -/
theorem sum_emod_72 (H : IsMagicSquareOfSquares m t) (hprim : IsPrimitive m) : t % 72 = 3 := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have h := HR.sum_emod_72 (primitive_roots hr hprim)
  omega

/-- In a primitive magic square of squares no prime `p ≡ 3 (mod 8)` divides any entry. -/
theorem not_dvd_of_three_mod_eight (H : IsMagicSquareOfSquares m t) (hprim : IsPrimitive m)
    {p : ℕ} (hp : p.Prime) (h8 : p % 8 = 3) : ∀ i j, ¬p ∣ m i j := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have : Fact p.Prime := ⟨hp⟩
  have key := forall_of_nine (P := fun i j => ¬(p : ℤ) ∣ (r i j : ℤ))
    (HR.not_dvd_of_three_mod_eight (primitive_roots hr hprim) h8)
  exact fun i j hd => key i j ((dvd_iff hp (hr i j)).mp hd)

/-- In a primitive magic square of squares every prime factor of the centre entry is
`≡ 1 (mod 4)`. -/
theorem center_prime_factors (H : IsMagicSquareOfSquares m t) (hprim : IsPrimitive m)
    {p : ℕ} (hp : p.Prime) (hd : p ∣ m 1 1) : p % 4 = 1 := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  exact HR.center_prime_factors (primitive_roots hr hprim) hp ((dvd_iff hp (hr 1 1)).mp hd)

/-- In a primitive magic square of squares no prime `p ≡ 3, 5 (mod 8)` divides an edge
entry. -/
theorem not_dvd_edge (H : IsMagicSquareOfSquares m t) (hprim : IsPrimitive m)
    {p : ℕ} (hp : p.Prime) (h8 : p % 8 = 3 ∨ p % 8 = 5) :
    ¬p ∣ m 0 1 ∧ ¬p ∣ m 1 0 ∧ ¬p ∣ m 1 2 ∧ ¬p ∣ m 2 1 := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have : Fact p.Prime := ⟨hp⟩
  obtain ⟨hb, hd, hf, hh⟩ := HR.not_dvd_edge (primitive_roots hr hprim) h8
  exact ⟨fun h => hb ((dvd_iff hp (hr 0 1)).mp h), fun h => hd ((dvd_iff hp (hr 1 0)).mp h),
    fun h => hf ((dvd_iff hp (hr 1 2)).mp h), fun h => hh ((dvd_iff hp (hr 2 1)).mp h)⟩

/-- A prime `p ≡ 3 (mod 4)` dividing a corner entry divides the two edge entries that do
not touch that corner. -/
theorem dvd_edges_of_dvd_corner (H : IsMagicSquareOfSquares m t) {p : ℕ} (hp : p.Prime)
    (h4 : p % 4 = 3) :
    (p ∣ m 0 0 → p ∣ m 1 2 ∧ p ∣ m 2 1) ∧ (p ∣ m 0 2 → p ∣ m 1 0 ∧ p ∣ m 2 1) ∧
    (p ∣ m 2 0 → p ∣ m 0 1 ∧ p ∣ m 1 2) ∧ (p ∣ m 2 2 → p ∣ m 0 1 ∧ p ∣ m 1 0) := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have : Fact p.Prime := ⟨hp⟩
  have D := fun i j => dvd_iff hp (hr i j)
  obtain ⟨k1, k2, k3, k4⟩ := HR.dvd_edges_of_dvd_corner h4
  exact ⟨fun h => ⟨(D 1 2).mpr (k1 ((D 0 0).mp h)).1, (D 2 1).mpr (k1 ((D 0 0).mp h)).2⟩,
    fun h => ⟨(D 1 0).mpr (k2 ((D 0 2).mp h)).1, (D 2 1).mpr (k2 ((D 0 2).mp h)).2⟩,
    fun h => ⟨(D 0 1).mpr (k3 ((D 2 0).mp h)).1, (D 1 2).mpr (k3 ((D 2 0).mp h)).2⟩,
    fun h => ⟨(D 0 1).mpr (k4 ((D 2 2).mp h)).1, (D 1 0).mpr (k4 ((D 2 2).mp h)).2⟩⟩

/-- A prime `p ≡ 3, 5 (mod 8)` dividing a corner entry divides the centre entry and the
opposite corner entry. -/
theorem dvd_center_of_dvd_corner (H : IsMagicSquareOfSquares m t) {p : ℕ} (hp : p.Prime)
    (h8 : p % 8 = 3 ∨ p % 8 = 5) :
    (p ∣ m 0 0 → p ∣ m 1 1 ∧ p ∣ m 2 2) ∧ (p ∣ m 0 2 → p ∣ m 1 1 ∧ p ∣ m 2 0) ∧
    (p ∣ m 2 0 → p ∣ m 1 1 ∧ p ∣ m 0 2) ∧ (p ∣ m 2 2 → p ∣ m 1 1 ∧ p ∣ m 0 0) := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have : Fact p.Prime := ⟨hp⟩
  have D := fun i j => dvd_iff hp (hr i j)
  obtain ⟨k1, k2, k3, k4⟩ := HR.dvd_center_of_dvd_corner h8
  exact ⟨fun h => ⟨(D 1 1).mpr (k1 ((D 0 0).mp h)).1, (D 2 2).mpr (k1 ((D 0 0).mp h)).2⟩,
    fun h => ⟨(D 1 1).mpr (k2 ((D 0 2).mp h)).1, (D 2 0).mpr (k2 ((D 0 2).mp h)).2⟩,
    fun h => ⟨(D 1 1).mpr (k3 ((D 2 0).mp h)).1, (D 0 2).mpr (k3 ((D 2 0).mp h)).2⟩,
    fun h => ⟨(D 1 1).mpr (k4 ((D 2 2).mp h)).1, (D 0 0).mpr (k4 ((D 2 2).mp h)).2⟩⟩

/-- If one entry is `0` then all nine entries are `0`. -/
theorem eq_zero_of_exists_eq_zero (H : IsMagicSquareOfSquares m t) (h0 : ∃ i j, m i j = 0) :
    ∀ i j, m i j = 0 := by
  obtain ⟨r, hr, HR⟩ := H.exists_roots
  have h0' : ∃ i j, (r i j : ℤ) = 0 := by
    obtain ⟨i, j, hij⟩ := h0
    rw [hr i j] at hij
    exact ⟨i, j, by exact_mod_cast (pow_eq_zero_iff two_ne_zero).mp hij⟩
  have key := forall_of_nine (P := fun i j => (r i j : ℤ) = 0)
    (HR.eq_zero_of_exists_eq_zero (nine_of_exists (P := fun i j => (r i j : ℤ) = 0) h0'))
  intro i j
  have h : (r i j : ℤ) = 0 := key i j
  have h' : r i j = 0 := by exact_mod_cast h
  rw [hr i j, h']
  rfl

/-- Distinct entries are automatically positive. -/
theorem pos (H : IsMagicSquareOfSquares m t) (hinj : m.Injective2) : ∀ i j, 0 < m i j := by
  intro i j
  by_contra hne
  have hz := H.eq_zero_of_exists_eq_zero ⟨i, j, by omega⟩
  have h01 : m 0 0 = m 0 1 := by rw [hz 0 0, hz 0 1]
  exact absurd (hinj h01).2 (by decide)

end IsMagicSquareOfSquares

/-- The positivity hypothesis in the open statement is redundant. -/
theorem openStatement_iff :
    OpenStatement ↔
      ∃ m : Fin 3 → Fin 3 → ℕ, ∃ t : ℕ, m.Injective2 ∧ IsMagicSquareOfSquares m t := by
  constructor
  · rintro ⟨m, t, hinj, hsq, hrow, hcol, hd, ha⟩
    exact ⟨m, t, hinj, ⟨fun i j => (hsq i j).2, hrow, hcol, hd, ha⟩⟩
  · rintro ⟨m, t, hinj, H⟩
    exact ⟨m, t, hinj, fun i j => ⟨H.pos hinj i j, H.isSquare i j⟩, H.row, H.col, H.diag,
      H.anti⟩

/-- The open statement in terms of nine integer roots with distinct squares. -/
theorem openStatement_iff_exists_isMSS :
    OpenStatement ↔ ∃ a b c d e f g h i T : ℤ, IsMSS a b c d e f g h i T ∧
      [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup := by
  rw [openStatement_iff]
  constructor
  · rintro ⟨m, t, hinj, H⟩
    obtain ⟨r, hr, HR⟩ := H.exists_roots
    refine ⟨_, _, _, _, _, _, _, _, _, _, HR, ?_⟩
    have hnd : ([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)] :
        List (Fin 3 × Fin 3)).Nodup := by decide
    have hinj' : Function.Injective (fun q : Fin 3 × Fin 3 => ((m q.1 q.2 : ℕ) : ℤ)) := by
      rintro ⟨u, v⟩ ⟨u', v'⟩ hq
      have hq' : ((m u v : ℕ) : ℤ) = (m u' v' : ℕ) := hq
      obtain ⟨h1, h2⟩ := hinj (Nat.cast_injective hq')
      rw [h1, h2]
    have hmap := hnd.map hinj'
    simpa [hr] using hmap
  · rintro ⟨a, b, c, d, e, f, g, h, i, T, H, hnd⟩
    have hT := H.sum_eq
    have HM : IsMagic (a ^ 2) (b ^ 2) (c ^ 2) (d ^ 2) (e ^ 2) (f ^ 2) (g ^ 2) (h ^ 2) (i ^ 2) T :=
      H
    let ρ : Fin 3 → Fin 3 → ℤ := ![![a, b, c], ![d, e, f], ![g, h, i]]
    refine ⟨fun u v => (ρ u v).natAbs ^ 2, 3 * e.natAbs ^ 2, ?_, ?_⟩
    · intro u₁ u₂ v₁ v₂ heq
      have hmem : ∀ q : Fin 3 × Fin 3, q ∈ ([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2),
          (2, 0), (2, 1), (2, 2)] : List (Fin 3 × Fin 3)) := by decide
      have hnd' : (([(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2)] :
          List (Fin 3 × Fin 3)).map (fun q => ρ q.1 q.2 ^ 2)).Nodup := hnd
      have hz : ρ u₁ v₁ ^ 2 = ρ u₂ v₂ ^ 2 := by
        have hc := congrArg (Nat.cast : ℕ → ℤ) heq
        simpa only [Nat.cast_pow, Int.natCast_natAbs, sq_abs] using hc
      exact Prod.mk.inj (List.inj_on_of_nodup_map hnd' (hmem (u₁, v₁)) (hmem (u₂, v₂)) hz)
    · refine ⟨fun u v => ⟨(ρ u v).natAbs, sq _⟩, ?_, ?_, natAbs_line hT HM.diag,
        natAbs_line hT HM.anti⟩
      · intro u
        rw [Fin.sum_univ_three]
        fin_cases u
        · exact natAbs_line hT HM.row₁
        · exact natAbs_line hT HM.row₂
        · exact natAbs_line hT HM.row₃
      · intro v
        rw [Fin.sum_univ_three]
        fin_cases v
        · exact natAbs_line hT HM.col₁
        · exact natAbs_line hT HM.col₂
        · exact natAbs_line hT HM.col₃

/-- Robertson's reformulation: the open statement holds if and only if there are three
3-term arithmetic progressions of squares with the same common difference whose middle
terms form an arithmetic progression, with all nine squares distinct. -/
theorem openStatement_iff_three_progressions :
    OpenStatement ↔ ∃ x y z n : ℤ, x + z = 2 * y ∧
      IsSquare (x - n) ∧ IsSquare x ∧ IsSquare (x + n) ∧
      IsSquare (y - n) ∧ IsSquare y ∧ IsSquare (y + n) ∧
      IsSquare (z - n) ∧ IsSquare z ∧ IsSquare (z + n) ∧
      [y + n, x - n, z, z - n, y, x + n, x, z + n, y - n].Nodup :=
  openStatement_iff_exists_isMSS.trans exists_isMSS_iff_three_progressions

/-- If the open problem has a solution then it has a primitive solution, so the necessary
conditions for primitive squares apply without loss of generality. -/
theorem openStatement_iff_primitive :
    OpenStatement ↔ ∃ a b c d e f g h i T : ℤ, IsMSS a b c d e f g h i T ∧
      [a ^ 2, b ^ 2, c ^ 2, d ^ 2, e ^ 2, f ^ 2, g ^ 2, h ^ 2, i ^ 2].Nodup ∧
      Primitive a b c d e f g h i := by
  rw [openStatement_iff_exists_isMSS]
  constructor
  · rintro ⟨a, b, c, d, e, f, g, h, i, T, H, hnd⟩
    exact H.exists_primitive hnd
  · rintro ⟨a, b, c, d, e, f, g, h, i, T, H, hnd, -⟩
    exact ⟨a, b, c, d, e, f, g, h, i, T, H, hnd⟩

end MagicSquareOfSquares
