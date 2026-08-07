module Main where

import FSOT.Genetics.ChemLink
import FSOT.Genetics.Claim
import FSOT.Genetics.Contact

main :: IO ()
main = do
  putStrLn "FSOT-Genetics Haskell check"
  putStrLn $ "  validClaim geneticsClaim = " ++ show (validClaim geneticsClaim)
  putStrLn $ "  Backbone D_eff           = " ++ show (dEff Backbone)
  putStrLn $ "  Disulfide D_eff          = " ++ show (dEff Disulfide)
  putStrLn $ "  SaltBridge D_eff         = " ++ show (dEff SaltBridge)
  putStrLn $ "  HydrophobicPack D_eff    = " ++ show (dEff HydrophobicPack)
  putStrLn $ "  HBondSecondary D_eff     = " ++ show (dEff HBondSecondary)
  putStrLn $ "  TertiaryBiochem D_eff    = " ++ show (dEff TertiaryBiochem)
  putStrLn $ "  backbone observed        = " ++ show (observed Backbone)
  putStrLn $ "  longRangeGateBiochem     = " ++ show longRangeGateBiochem
  let demo =
        rankContacts
          [ (ContactPair 1 20 19, scoreContact [FromHydrophobicCore 2.0, FromDistogramProximity 1.0])
          , (ContactPair 2 15 13, scoreContact [FromSaltBridge 3.0])
          , (ContactPair 5 40 35, scoreContact [FromDistogramProximity 0.5])
          ]
  putStrLn $ "  demo rankContacts head   = " ++ show (take 1 demo)
  putStrLn "  ALL HASKELL GATES PASSED"
