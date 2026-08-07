# RCSB Oriented Backbone Validation

## Protocol

- F19 candidate freeze: `3673219`.
- Holdout protocol freeze: `9a925a5`.
- Twelve previously unused, sequence-nonredundant RCSB chains.
- X-ray reference resolution: 1.394-1.85 A.
- Chain lengths: 44-194 residues.
- Selection used metadata and sequence only before the one-shot evaluation.

## Result

| Model | Median C-alpha RMSD | Bootstrap 95% interval |
| --- | ---: | ---: |
| Baseline production fold | **15.066 A** | [12.523, 17.051] |
| Cooperative F12c | 16.189 A | [12.683, 16.710] |
| Cooperative F12c plus F19 | 16.042 A | [13.239, 16.687] |

F19 reflected five of twelve predictions. Its maximum pair-distance change was
exactly 0.0 A. The preregistered narrow gate passed because median F19 RMSD did
not exceed F12c RMSD and pair distances were preserved.

## Interpretation

F19 correctly adds information that a symmetric distance matrix cannot contain:
the choice between mirror embeddings. It does not improve the distance matrix
or contact topology. The small median change, overlapping confidence intervals,
and a 2.587 A regression on `9DKS:C` mean F19 is not promoted to the production
default.

The larger result is that F12c secondary-label improvements still do not convert
to coordinate accuracy. Baseline F12 remains better on this holdout by 0.976 A
at the median relative to F12c plus F19. The next coordinate target remains
directional, competitive, coordinate-dependent tertiary interactions rather
than another reflection or scalar-ranking adjustment.

Machine-readable results are in `data/rcsb_oriented_backbone_eval.json`.
