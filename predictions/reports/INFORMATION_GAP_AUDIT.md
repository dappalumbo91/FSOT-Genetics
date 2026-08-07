# FSOT Structure Information-Gap Audit

## Conclusion

The current approximately 8 A RMSD wall is not an intrinsic limit of classical
MDS and is not evidence that FSOT must become machine learning. It is caused by
information missing before and during coordinate construction:

1. Predicted pair distances are inaccurate and mutually inconsistent.
2. F13 activates same-kind secondary regions without selecting a physical
   interaction graph.
3. The objective is reflection-invariant and cannot determine protein
   handedness.
4. Pair terms do not compete for saturable interactions or respond to the
   evolving three-dimensional environment.

## Measured Localization

The audit uses the frozen five-protein development set. It does not fit a
constant or alter production output.

| Measurement | Result |
| --- | ---: |
| Baseline complete-fold median RMSD | 8.590 A |
| F12c complete-fold median RMSD | 10.172 A |
| Direct MDS median RMSD from predicted distances | 9.180 A |
| Median predicted-distance MAE | 5.252 A |
| Median predicted/native distance Pearson | 0.416 |
| Median negative Gram-eigenvalue mass | 0.200 |
| F12c same-kind region edges that are native | 8 / 18 |
| Mirror stress difference | exactly zero on 5 / 5 proteins |

Classical MDS applied to the exact native distance matrix reconstructs four
proteins at 0.026-0.042 A RMSD. Engrailed is reconstructed as the opposite
enantiomer at 10.589 A under proper-rotation RMSD, but at 0.009 A after testing
the mirror. Therefore exact distances are sufficient up to reflection, while
the current FSOT distances are neither accurate nor Euclidean enough.

## Dataset Versus Machine Learning

Using PDB, PISCES, DSSP, Kaggle, or another empirical collection as a blind
test does not make the method machine learning. The boundary is whether labels
determine fitted weights, lookup values, thresholds, or model selection.

The Kaggle/PISCES validation measured H/E/C classification on 6,483 eligible
chains. F12c improved macro recall and beta recall, but those labels did not
train a neural network. The disclosed four-value development comparison did,
however, select one formula variant. It must therefore be treated as model
selection and followed by untouched holdouts, as was done here.

`free_parameters: 0` currently means no learned model weights. It does not mean
the runtime contains no engineering choices. Examples include top-L and L/2
contact caps, 24 refinement rounds, force weights, learning-rate clipping,
all-to-all same-kind region coupling, and fixed sparse-pair counts. These
choices must be derived from geometry, classified as numerical solver controls,
or removed before making a stronger zero-arbitrary-parameter claim.

## DNA Is Not the Primary Missing Input

DNA determines the translated amino-acid sequence, but synonymous codons
usually map to the same equilibrium protein sequence. The current fold runtime
does not consume codons; it consumes amino acids. That is the correct primary
level for predicting an isolated canonical protein fold.

Genome context matters when it changes the proteoform or folding process:

- splice isoform and start/stop choice;
- signal peptides and membrane targeting;
- post-translational modifications;
- cotranslational folding and codon-dependent translation rate;
- oligomer partners, ligands, metals, and cofactors;
- redox state and disulfide formation.

Those effects require more than raw DNA. They require organism, expression,
cellular-compartment, and experimental-condition metadata. Genome data should
therefore enter through an explicit context layer, not replace the protein
sequence used by the structure engine.

## Missing Physical State

### Oriented backbone geometry

A symmetric distance objective determines coordinates only up to reflection.
The runtime has no signed pseudo-dihedral, peptide-plane orientation, or
L-amino-acid chirality term. This is a proven non-identifiability, not a tuning
problem.

### Competitive and saturable interactions

Hydrogen bonds, salt bridges, disulfides, and packing contacts cannot all be
simultaneously assigned to every attractive pair. Current region and contact
terms are scored independently, so adding correct secondary regions creates
false global constraints. Interaction valence, hydrogen-bond directionality,
steric exclusion, and competition are absent.

### Environment-dependent energy

The pair matrix is built before coordinates exist. It therefore cannot use
solvent exposure, dielectric screening, burial, side-chain orientation,
neighbor density, or clashes in the emerging fold. Dipole attraction is always
attractive in F06 even though dipole interactions are orientation-dependent.

### Many-body topology

Sequence separation and residue-pair identity do not uniquely select a contact
map. Sheet order, parallel versus antiparallel orientation, helix-bundle
packing, loop closure, and domain arrangement are collective constraints. F17
scores register after two beta regions are assumed to interact; it does not
decide whether they should interact.

## Next Engineering Targets

1. Add an oriented internal-coordinate backbone state using virtual bond
   angles, signed pseudo-dihedrals, peptide handedness, and exact C-alpha bond
   geometry. The first gate is resolving the oracle-distance mirror without
   using native labels.
2. Replace all-to-all F13 coupling with a competitive region-interaction graph.
   Hydrogen-bond/register candidates must obey direction, saturation, and
   steric compatibility. The first gate is improved native-edge precision
   without losing edge recall.
3. Evaluate interactions in the evolving coordinates, including excluded
   volume, solvent exposure, electrostatic screening, and side-chain-facing
   orientation. Physical conditions are explicit inputs, not learned weights.
4. Produce a geometrically realizable distance matrix or optimize internal
   coordinates directly. The negative eigenvalue-mass gate must fall from the
   current median 0.200.
5. Freeze the resulting physics and solver before a new nonredundant PDB/CATH
   holdout. Report C-alpha RMSD together with TM-score, GDT, and lDDT so domain
   motion and local accuracy are not hidden by one global number.

The 1-2 A objective remains a legitimate research target, but crystallographic
resolution is not prediction RMSD and does not imply that a sequence method can
reach the same number. The next measurable milestone is not 2 A directly: it is
to remove the demonstrated reflection ambiguity, raise topology-edge precision,
and reduce distance inconsistency on untouched structures.

Machine-readable measurements are in `data/information_gap_audit.json`.
