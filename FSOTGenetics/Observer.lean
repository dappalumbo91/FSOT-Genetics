/-
  Observer policy by chemical connection (Mathlib).
-/
import Mathlib.Tactic
import FSOTGenetics.ChemLink

namespace FSOTGenetics

/-- Whether the fold treats this link as an *observation* of structure. -/
def ChemLink.observed : ChemLink → Bool
  | .backbone           => false
  | .disulfide          => true
  | .saltBridge         => true
  | .hydrophobicPack    => true
  | .hbondSecondary     => true
  | .molecularSidechain => true
  | .tertiaryBiochem    => true

theorem backbone_unobserved : ChemLink.backbone.observed = false := by rfl
theorem disulfide_observed : ChemLink.disulfide.observed = true := by rfl
theorem salt_observed : ChemLink.saltBridge.observed = true := by rfl
theorem pack_observed : ChemLink.hydrophobicPack.observed = true := by rfl
theorem tertiary_observed : ChemLink.tertiaryBiochem.observed = true := by rfl

def observerApplies (c : ChemLink) : Prop := c.observed = true

theorem observer_on_pack : observerApplies .hydrophobicPack := by rfl

theorem observer_off_backbone : ¬ observerApplies .backbone := by
  simp [observerApplies, ChemLink.observed]

/-- Only backbone is unobserved among the seven chem links. -/
theorem observed_iff_not_backbone (c : ChemLink) :
    c.observed = true ↔ c ≠ .backbone := by
  cases c <;> simp [ChemLink.observed]

end FSOTGenetics
