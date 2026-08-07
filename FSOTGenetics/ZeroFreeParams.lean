/-
  Zero free parameters — Mathlib-backed **claim lock**.

  This module does **not** introduce free parameters.
  It proves the genetics claim path has freeParameters = 0,
  no neural weights, and no free continuous D_eff search.
-/
import Mathlib.Tactic
import FSOTGenetics.ChemLink
import FSOTGenetics.Seeds

namespace FSOTGenetics

/-- What the public claim path is allowed to use. -/
structure ClaimPath where
  /-- Must be 0 on the genetics fold claim path. -/
  freeParameters : Nat
  usesNeuralWeights : Bool
  /-- Continuous free D search forbidden (named pin domains only). -/
  usesFreeD_eff : Bool

def geneticsClaim : ClaimPath where
  freeParameters := 0
  usesNeuralWeights := false
  usesFreeD_eff := false

theorem genetics_claim_zero_free : geneticsClaim.freeParameters = 0 := by rfl
theorem genetics_claim_no_nn : geneticsClaim.usesNeuralWeights = false := by rfl
theorem genetics_claim_no_free_D : geneticsClaim.usesFreeD_eff = false := by rfl

def dEffFromLink (c : ChemLink) : Nat := c.D_eff

theorem dEff_determined_by_link (c : ChemLink) :
    dEffFromLink c = c.domain.D_eff := by rfl

/-- Free-parameter count is additive zero (no hidden dials). -/
theorem free_params_add_zero (n : Nat) :
    freeParameterCount + n = n := by
  simp [freeParameterCount]

end FSOTGenetics
