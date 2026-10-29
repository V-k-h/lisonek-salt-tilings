import Mathlib

open Matrix

variable (p : ℕ) [Fact p.Prime]

def row12_orig (x : ZMod p) : Fin 4 → ZMod p :=
  ![(x^2 + x) * 2⁻¹, -(x^2 + 2*x), 0, 0]

def row34_orig (n y : ZMod p) : Fin 4 → ZMod p :=
  ![(y^2 - y) * (2 * n)⁻¹, (y - y^2) * n⁻¹, 1 - y, y]

def M_canon (n : ZMod p) : Matrix (Fin 4) (Fin 4) (ZMod p) :=
  !![4, 1, 0, n⁻¹; -2, -1, 0, 0; 0, 0, 1, 1; 0, 0, 0, 1]

/-- Lemma 1, block 12: the rows r_x^{12} go to the parabola (x², x, 0, 0). -/
theorem M_canon_row12 (n x : ZMod p) (h2 : (2 : ZMod p) ≠ 0) :
    (M_canon p n).mulVec (row12_orig p x) = ![x^2, x, 0, 0] := by
  ext i
  fin_cases i <;> simp [M_canon, row12_orig, mulVec, dotProduct, Fin.sum_univ_four] <;>
    field_simp <;> ring

/-- Lemma 1, block 34: the rows r_y^{34} go to the parabola (y²/n, 0, 1, y). -/
theorem M_canon_row34 (n y : ZMod p) (h2 : (2 : ZMod p) ≠ 0) (hn : n ≠ 0) :
    (M_canon p n).mulVec (row34_orig p n y) = ![y^2 * n⁻¹, 0, 1, y] := by
  ext i
  fin_cases i <;> simp [M_canon, row34_orig, mulVec, dotProduct, Fin.sum_univ_four] <;>
    field_simp <;> ring

/-- M_canon is invertible: det = -2. -/
theorem M_canon_det (n : ZMod p) : (M_canon p n).det = -2 := by
  simp [M_canon, det_succ_row_zero, Fin.sum_univ_succ, Fin.succAbove]
  ring

/-- Key step in Theorem 3: for a nonsquare n, w₄²/n − w₂² = 0 forces w₂ = w₄ = 0. -/
theorem quad_coeff_ne_zero (n w2 w4 : ZMod p) (hn : ¬ IsSquare n)
    (hw : (w2, w4) ≠ (0, 0)) : w4^2 * n⁻¹ - w2^2 ≠ 0 := by
  intro h
  have hn0 : n ≠ 0 := by rintro rfl; exact hn ⟨0, by simp⟩
  have h' : w4^2 = n * w2^2 := by field_simp at h; linear_combination h
  by_cases h2 : w2 = 0
  · subst h2
    have : w4 = 0 := by simpa using h'
    exact hw (by simp [this])
  · apply hn
    refine ⟨w4 / w2, ?_⟩
    field_simp
    linear_combination -h'
