/-
  Generated catalog of Genetics solves (Nat milliscale).
  Do not edit by hand — python verification/export_obligations.py
  Pin D1D38A. Law S = K(T1+T2+T3). 0 free parameters.
-/
import FSOTGenetics.Seeds

namespace FSOTGenetics

def freeParameters : Nat := 0
def phiMilli : Nat := 1618
def leftoverMilli : Nat := 382
def closeHomologMilli : Nat := 618
def backboneD : Nat := 8
def disulfideD : Nat := 7
def saltD : Nat := 9
def packD : Nat := 14
def hbondD : Nat := 8
def molecularD : Nat := 9
def tertiaryD : Nat := 13
def longRangeGate : Nat := 7
def chemLinkCard : Nat := 7
def productN : Nat := 10
def productSub2A : Nat := 10
def productMedianMilliA : Nat := 133
def alphafoldMedianMilliA : Nat := 471
def productLtAf : Nat := 133
def productLtAfRhs : Nat := 471
def productLtBulk : Nat := 133
def productLtBulkRhs : Nat := 13574
def homologMeasured : Nat := 26
def homologFolded : Nat := 26
def homologTrueMiss : Nat := 1
def homologClose : Nat := 10
def homologCloseRhs : Nat := 24
def analogJobs : Nat := 2
def maleOlfVsVnc : Nat := 6
def maleOlfVsVncRhs : Nat := 107402
def maleOlfVsJo : Nat := 6
def maleOlfVsJoRhs : Nat := 77729
def bancOlfVsVnc : Nat := 8
def bancOlfVsVncRhs : Nat := 101619
def bancOlfVsJo : Nat := 8
def bancOlfVsJoRhs : Nat := 66369
def larvaOlfVsMech : Nat := 4292
def larvaOlfVsMechRhs : Nat := 166016
def wormHermVsMaleSex : Nat := 745
def wormHermVsMaleSexRhs : Nat := 77442
def plantArabidopsisFolded : Nat := 6
def plantRiceMaizeFolded : Nat := 12
def plantRiceMaizeMiss : Nat := 0
def plantMinId : Nat := 618
def plantMinIdRhs : Nat := 741
def wormCells : Nat := 0
def wormCellsRhs : Nat := 453
def cionaCells : Nat := 0
def cionaCellsRhs : Nat := 205
def platynereisCells : Nat := 0
def platynereisCellsRhs : Nat := 1720
def maleNeurons : Nat := 100000
def maleNeuronsRhs : Nat := 165122
def bancNeurons : Nat := 100000
def bancNeuronsRhs : Nat := 175401
def cionaGabaFlag : Nat := 0
def platynereisGabaFlag : Nat := 0

theorem okFreeParametersZero : freeParameters = 0 := by decide

theorem okPhiMilli : phiMilli = 1618 := by decide

theorem okLeftoverFloorMilli : leftoverMilli = 382 := by decide

theorem okCloseHomologMilli : closeHomologMilli = 618 := by decide

theorem okChemBackboneD : backboneD = 8 := by decide

theorem okChemDisulfideD : disulfideD = 7 := by decide

theorem okChemSaltD : saltD = 9 := by decide

theorem okChemPackD : packD = 14 := by decide

theorem okChemHbondD : hbondD = 8 := by decide

theorem okChemMolecularD : molecularD = 9 := by decide

theorem okChemTertiaryD : tertiaryD = 13 := by decide

theorem okLongRangeGate : longRangeGate = 7 := by decide

theorem okChemLinkCard : chemLinkCard = 7 := by decide

theorem okProductN : productN = 10 := by decide

theorem okProductSub2A : productSub2A = 10 := by decide

theorem okProductMedianMilliA : productMedianMilliA = 133 := by decide

theorem okAlphafoldMedianMilliA : alphafoldMedianMilliA = 471 := by decide

theorem okProductLtAlphafold : productLtAf < productLtAfRhs := by decide

theorem okProductLtBulk : productLtBulk < productLtBulkRhs := by decide

theorem okHomologMeasured : homologMeasured = 26 := by decide

theorem okHomologFoldedEqMeasured : homologFolded = 26 := by decide

theorem okHomologTrueMiss : homologTrueMiss = 1 := by decide

theorem okHomologCloseCount : homologClose ≤ homologCloseRhs := by decide

theorem okAnalogJobs : analogJobs = 2 := by decide

theorem okMaleVncGtOlf : maleOlfVsVnc < maleOlfVsVncRhs := by decide

theorem okMaleJoGtOlf : maleOlfVsJo < maleOlfVsJoRhs := by decide

theorem okBancVncGtOlf : bancOlfVsVnc < bancOlfVsVncRhs := by decide

theorem okBancJoGtOlf : bancOlfVsJo < bancOlfVsJoRhs := by decide

theorem okLarvaMechGtOlf : larvaOlfVsMech < larvaOlfVsMechRhs := by decide

theorem okWormMaleSexGtHerm : wormHermVsMaleSex < wormHermVsMaleSexRhs := by decide

theorem okPlantArabidopsisFolded : plantArabidopsisFolded = 6 := by decide

theorem okPlantRiceMaizeFolded : plantRiceMaizeFolded = 12 := by decide

theorem okPlantRiceMaizeMiss : plantRiceMaizeMiss = 0 := by decide

theorem okPlantMinIdClose : plantMinId ≤ plantMinIdRhs := by decide

theorem okWormCells : wormCells < wormCellsRhs := by decide

theorem okCionaCells : cionaCells < cionaCellsRhs := by decide

theorem okPlatynereisCells : platynereisCells < platynereisCellsRhs := by decide

theorem okMaleNeurons : maleNeurons < maleNeuronsRhs := by decide

theorem okBancNeurons : bancNeurons < bancNeuronsRhs := by decide

theorem okCionaUnsignedGaba : cionaGabaFlag = 0 := by decide

theorem okPlatynereisUnsignedGaba : platynereisGabaFlag = 0 := by decide

end FSOTGenetics