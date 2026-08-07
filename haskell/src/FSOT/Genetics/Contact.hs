{- |
  Long-range contact ranking types (next error-margin mode).

  Ranking is driven by evidence tags, not free trained weights.
  A @ContactScore@ is only built from closed seed-side contributions.
-}
module FSOT.Genetics.Contact
  ( ResidueIndex
  , ContactPair (..)
  , Evidence (..)
  , ContactScore (..)
  , scoreContact
  , rankContacts
  ) where

import Data.List (sortOn)
import Data.Ord (Down (..))

type ResidueIndex = Int

data ContactPair = ContactPair
  { cpI   :: !ResidueIndex
  , cpJ   :: !ResidueIndex
  , cpSep :: !Int
  } deriving (Eq, Show)

-- | Closed-form evidence channels (no free weights — combine by sum of seeds later in Python).
data Evidence
  = FromDistogramProximity Double  -- M_ij from F15
  | FromHydrophobicCore Double     -- SMILES KD product channel
  | FromSaltBridge Double          -- opposite charge
  | FromDisulfide Double           -- F18 gate
  | FromSecondaryRegister Double   -- F16/F17
  | FromPolarizability Double      -- SMILES §26
  deriving (Eq, Show)

-- | Total score = sum of evidence magnitudes (compile-time: only Evidence constructors).
newtype ContactScore = ContactScore { unContactScore :: Double }
  deriving (Eq, Ord, Show)

evidenceValue :: Evidence -> Double
evidenceValue (FromDistogramProximity x) = x
evidenceValue (FromHydrophobicCore x)    = x
evidenceValue (FromSaltBridge x)         = x
evidenceValue (FromDisulfide x)          = x
evidenceValue (FromSecondaryRegister x)  = x
evidenceValue (FromPolarizability x)     = x

scoreContact :: [Evidence] -> ContactScore
scoreContact = ContactScore . sum . map evidenceValue

-- | Rank pairs high-score first (for top-L selection).
rankContacts :: [(ContactPair, ContactScore)] -> [(ContactPair, ContactScore)]
rankContacts = sortOn (Down . unContactScore . snd)
