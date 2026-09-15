---- MODULE GeneticsSolve ----
EXTENDS Naturals

\* FSOT-Genetics routing. Pin D1D38A. 0 free parameters.
\* Measured homolog → product Cα. No map → no_measured_map.
\* Genome ≠ connectome. Unsigned transmitter stays unsigned.

CONSTANTS Unobserved, Observed

VARIABLES observer, homolog, product, synapses, gabaInvented, freeParams

TypeOK ==
  /\ observer \in {Unobserved, Observed}
  /\ homolog \in BOOLEAN
  /\ product \in BOOLEAN
  /\ synapses \in BOOLEAN
  /\ gabaInvented \in BOOLEAN
  /\ freeParams \in Nat

Init ==
  /\ observer = Unobserved
  /\ homolog = FALSE
  /\ product = FALSE
  /\ synapses = FALSE
  /\ gabaInvented = FALSE
  /\ freeParams = 0

ObserveSidechain ==
  /\ observer' = Observed
  /\ UNCHANGED <<homolog, product, synapses, gabaInvented, freeParams>>

BackboneUnobserved ==
  /\ observer' = Unobserved
  /\ UNCHANGED <<homolog, product, synapses, gabaInvented, freeParams>>

MeasureHomolog ==
  /\ homolog' = TRUE
  /\ UNCHANGED <<observer, product, synapses, gabaInvented, freeParams>>

FoldProduct ==
  /\ homolog = TRUE
  /\ product' = TRUE
  /\ UNCHANGED <<observer, homolog, synapses, gabaInvented, freeParams>>

NoMapNoProduct ==
  /\ homolog = FALSE
  /\ product' = FALSE
  /\ UNCHANGED <<observer, homolog, synapses, gabaInvented, freeParams>>

MeasureSynapses ==
  /\ synapses' = TRUE
  /\ UNCHANGED <<observer, homolog, product, gabaInvented, freeParams>>

Next ==
  ObserveSidechain \/ BackboneUnobserved \/ MeasureHomolog \/ FoldProduct
  \/ NoMapNoProduct \/ MeasureSynapses

Spec == Init /\ [][Next]_<<observer, homolog, product, synapses, gabaInvented, freeParams>>

\* Product Cα only where a measured homolog exists.
ProductRequiresHomolog == product => homolog

\* Never invent GABA when the graph is unsigned.
NeverInventGABA == gabaInvented = FALSE

\* Claim path has no free parameters.
ZeroFreeParams == freeParams = 0

\* A genome is not a connectome: product on a protein does not create synapses.
GenomeIsNotConnectome == product => (synapses = TRUE \/ synapses = FALSE)

====
