/-
  Seed constants on Mathlib `ℝ` — residual / chaos / positivity lemmas.
-/
import Mathlib.Data.Real.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Exp
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Analysis.Real.Pi.Bounds
import Mathlib.Tactic
import FSOTGenetics.Seeds

namespace FSOTGenetics

open Real

/-- Golden ratio φ = (1+√5)/2. -/
noncomputable def phiR : ℝ := (1 + sqrt 5) / 2

/-- Euler's number. -/
noncomputable def eR : ℝ := exp 1

/-- π as Mathlib Real.pi. -/
noncomputable def piR : ℝ := Real.pi

theorem phiR_pos : 0 < phiR := by
  unfold phiR
  positivity

theorem phiR_gt_one : 1 < phiR := by
  unfold phiR
  have h5 : (1 : ℝ) < 5 := by norm_num
  have hsqrt : 1 < sqrt 5 := by
    have := Real.sqrt_lt_sqrt (by norm_num : (0:ℝ) ≤ 1) h5
    simpa using this
  -- 1 < (1+√5)/2 ↔ 2 < 1+√5 ↔ 1 < √5
  linarith

theorem eR_pos : 0 < eR := exp_pos 1

theorem eR_gt_one : 1 < eR := one_lt_exp_iff.mpr (by norm_num)

theorem piR_pos : 0 < piR := Real.pi_pos

/-- Residual factor (1 + |S| · P_new) ≥ 1 when P_new ≥ 0. -/
theorem residual_factor_ge_one (S Pnew : ℝ) (hP : 0 ≤ Pnew) :
    1 ≤ 1 + |S| * Pnew := by
  have : 0 ≤ |S| * Pnew := mul_nonneg (abs_nonneg _) hP
  linarith

/-- Chaos factor at D = 25 is identity. -/
theorem chaos_at_D25 (c : ℝ) : 1 + c * ((25 : ℝ) - 25) / 25 = 1 := by
  ring

/-- For D > 25 and chaos seed c < 0, chaos factor is strictly below 1. -/
theorem chaos_factor_lt_one_of_D_gt_25 {c D : ℝ}
    (hc : c < 0) (hD : 25 < D) :
    1 + c * (D - 25) / 25 < 1 := by
  have hden : (0 : ℝ) < 25 := by norm_num
  have hnum : 0 < D - 25 := sub_pos.mpr hD
  have hquot : 0 < (D - 25) / 25 := div_pos hnum hden
  have hprod : c * ((D - 25) / 25) < 0 := mul_neg_of_neg_of_pos hc hquot
  have : c * (D - 25) / 25 = c * ((D - 25) / 25) := by ring
  linarith

/-- (γ/e)·√2 > 0 (P_NEW shape). -/
theorem pnew_style_pos (γ e : ℝ) (hγ : 0 < γ) (he : 0 < e) :
    0 < (γ / e) * sqrt 2 := by
  positivity

/-- η_eff = 1/(π-1). Uses `pi > 3` so π − 1 > 0. -/
noncomputable def etaEff : ℝ := 1 / (piR - 1)

theorem etaEff_pos : 0 < etaEff := by
  unfold etaEff piR
  have h : (3 : ℝ) < Real.pi := Real.pi_gt_three
  have : (0 : ℝ) < Real.pi - 1 := by linarith
  positivity

/-- Long-range gate Nat fact (matches Python ⌈η·13⌉ = 7). -/
theorem long_range_gate_nat : longRangeGateBiochem = 7 := long_range_gate_biochem

end FSOTGenetics
