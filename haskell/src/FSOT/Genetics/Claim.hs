{- |
  Claim path lock — free parameters are a *forbidden count*, not a feature.

  Same spirit as Lean ZeroFreeParams: this module will not type-check a claim
  with freeParameters /= 0 if you use @validClaim@.
-}
module FSOT.Genetics.Claim
  ( ClaimPath (..)
  , geneticsClaim
  , validClaim
  ) where

data ClaimPath = ClaimPath
  { freeParameters    :: !Int
  , usesNeuralWeights :: !Bool
  , usesFreeDEff      :: !Bool  -- continuous free D search
  } deriving (Eq, Show)

geneticsClaim :: ClaimPath
geneticsClaim = ClaimPath
  { freeParameters    = 0
  , usesNeuralWeights = False
  , usesFreeDEff      = False
  }

-- | Only the zero-free, no-NN, no-free-D path is valid.
validClaim :: ClaimPath -> Bool
validClaim c =
  freeParameters c == 0
    && not (usesNeuralWeights c)
    && not (usesFreeDEff c)
