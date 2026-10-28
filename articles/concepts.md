# How the paradigms work (animated)

Five short animations, one per idea behind the integration paradigms.
Each is a small vector graphic, so it stays crisp at any size. The
colours follow the package palette: teal and amber for two sources, navy
for the result.

## Direct integration

The sources are concatenated into one wide matrix and clustered once.
Fast and transparent (`ADC`, `ADEC`).

![Two source blocks slide together into one, then are
clustered.](../anim/direct.svg)

## Similarity fusion

Each source gives an object-by-object similarity network. The networks
reinforce each other through cross-diffusion and converge to one fused
network (`SNF`).

![Two similarity networks exchange edges and fuse into
one.](../anim/fusion.svg)

## Voting consensus

Many partitions vote on whether each pair of objects belongs together.
The votes accumulate into a consensus that is then clustered (`M_ABC`,
`CEC`, `EvidenceAccumulation`).

![Three partitions of the same objects vote, and a consensus partition
emerges.](../anim/voting.svg)

## Hierarchy

Per-source hierarchies are consensed into a single dendrogram. The tree
is read from the leaves upward to the root
(`HierarchicalEnsembleClustering`, `EHC`).

![A dendrogram draws itself from the leaves up to the
root.](../anim/hierarchy.svg)

## Data nuggets

A large set of observations is compressed into a few weighted
representatives (nuggets). Feature weights and clustering then run on
the small set, which is why the scheme scales to large samples.

![Many scattered points collapse into a few larger weighted
nuggets.](../anim/nuggets.svg)

Each animation loops. For the full method details, see the [method
comparison](https://lucp12891.github.io/MosaiClusteR/articles/method-comparison.md)
and the guides.
