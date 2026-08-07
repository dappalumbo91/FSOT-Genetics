import Lake
open Lake DSL

/-!
# FSOT-Genetics — formal verification with Mathlib

Sibling of [FSOT-2.1-Lean](https://github.com/dappalumbo91/FSOT-2.1-Lean).

**Mathlib is required** (`mathlib4 @ v4.31.0`, Lean `v4.31.0`).

Gates: chem-link → D_eff, observer policy, zero free parameters, residual lemmas on ℝ.

Note: On Windows, linking a Mathlib-heavy `lean_exe` can hit the PE export limit
(~65k symbols). The formal gate is **`lake build`** (typecheck all modules), not a
linked CLI. CI runs `lake build` only.
-/

package «FSOT-Genetics» where
  -- package config

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.31.0"

@[default_target]
lean_lib FSOTGenetics where
  roots := #[`FSOTGenetics]
  globs := #[.submodules `FSOTGenetics]
