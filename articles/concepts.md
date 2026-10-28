# How the methods work (animated)

Short animations for the ideas behind MosaiClusteR. Each is a small
vector graphic that loops, so it stays crisp at any size. Teal and amber
mark two sources; blue and navy mark objects and results.

## Integration paradigms

The five ways MosaiClusteR combines sources.

### Direct integration

The sources are concatenated into one wide matrix and clustered once.
Fast and transparent (`ADC`, `ADEC`).

![Two source blocks slide together into one, then are
clustered.](../anim/direct.svg)

### Similarity fusion

Each source gives an object-by-object similarity network. The networks
reinforce each other through cross-diffusion and converge to one fused
network (`SNF`).

![Two similarity networks exchange edges and fuse into
one.](../anim/fusion.svg)

### Voting consensus

Many partitions vote on whether each pair of objects belongs together.
The votes accumulate into a consensus that is then clustered (`M_ABC`,
`CEC`, `EvidenceAccumulation`).

![Three partitions of the same objects vote, and a consensus partition
emerges.](../anim/voting.svg)

### Hierarchy

Per-source hierarchies are consensed into a single dendrogram, read from
the leaves up to the root (`HierarchicalEnsembleClustering`, `EHC`).

![A dendrogram draws itself from the leaves up to the
root.](../anim/hierarchy.svg)

## Inside the algorithms

A closer look at two mechanisms that recur across the methods.

### Feature bundling (the ABC family)

For each source, M-ABC repeatedly draws a weighted subsample of features
(a “bundle”), clusters on it, and records the result. Signal that hides
in a subset of features is found across many bundles.

![Feature columns are repeatedly resampled into
bundles.](../anim/bundles.svg)

### Consensus co-association

Every base clustering casts a vote for each pair of objects. The votes
build an object-by-object co-association matrix, and its block structure
is the consensus the final clustering recovers.

![A co-association matrix fills in, and three cluster blocks emerge on
the diagonal.](../anim/coassoc.svg)

## Choosing and checking the clusters

How many clusters are there? The average silhouette width rises to a
peak at the number of clusters the data supports.

![A silhouette-width curve draws itself and peaks at k =
3.](../anim/selectk.svg)

### Data nuggets

For very large samples, the objects are first compressed into a few
weighted representatives (nuggets). Feature weighting and clustering
then run on the small set.

![Many scattered points collapse into a few larger weighted
nuggets.](../anim/nuggets.svg)

For the full method details, see the [method
comparison](https://lucp12891.github.io/MosaiClusteR/articles/method-comparison.md)
and the
[guides](https://lucp12891.github.io/MosaiClusteR/articles/index.md).
