{- |
  Chemical connection → D_eff (must match Python full_scalar_law.chem_link_domain
  and Lean FSOTGenetics.ChemLink).

  If someone tries to invent a free continuous D, they have to leave this module —
  there is no @Double@ D_eff constructor, only named links.
-}
module FSOT.Genetics.ChemLink
  ( ChemLink (..)
  , DomainName (..)
  , domainOf
  , dEff
  , observed
  , longRangeGateBiochem
  ) where

-- | Connecting chemical systems on the fold path.
data ChemLink
  = Backbone
  | Disulfide
  | SaltBridge
  | HydrophobicPack
  | HBondSecondary
  | MolecularSidechain
  | TertiaryBiochem
  deriving (Eq, Ord, Show, Enum, Bounded)

-- | Pin-table domain names (35-domain FSOT table subset).
data DomainName
  = PhysicalChemistry
  | AtomicPhysics
  | Electromagnetism
  | CondensedMatter
  | Chemistry
  | MolecularChemistry
  | Biochemistry
  deriving (Eq, Ord, Show)

-- | Chem link → domain (interacting chemical systems).
domainOf :: ChemLink -> DomainName
domainOf Backbone           = PhysicalChemistry
domainOf Disulfide          = AtomicPhysics
domainOf SaltBridge         = Electromagnetism
domainOf HydrophobicPack    = CondensedMatter
domainOf HBondSecondary     = Chemistry
domainOf MolecularSidechain = MolecularChemistry
domainOf TertiaryBiochem    = Biochemistry

-- | Effective dimension from pin table (Int, not free Double).
dEff :: ChemLink -> Int
dEff c = case domainOf c of
  PhysicalChemistry  -> 8
  AtomicPhysics      -> 7
  Electromagnetism   -> 9
  CondensedMatter    -> 14
  Chemistry          -> 8
  MolecularChemistry -> 9
  Biochemistry       -> 13

-- | Observer flag by chemical system (backbone = geometry, not measurement).
observed :: ChemLink -> Bool
observed Backbone = False
observed _        = True

-- | F13 gate ⌈η·13⌉ = 7 for biochemistry long-range.
longRangeGateBiochem :: Int
longRangeGateBiochem = 7
