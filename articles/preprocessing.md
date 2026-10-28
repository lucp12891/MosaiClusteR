# Preprocessing and source alignment

Integration is only as good as the inputs it receives. This page covers
the two preprocessing decisions that matter most in MosaiClusteR: how to
make heterogeneous sources comparable, and how the sources are aligned
so that every method sees the same objects.

## The data model

Every source is a numeric matrix with **objects (samples) in rows and
features in columns**. All sources in a `List` describe the **same
objects**; they may differ in the number of features. This single
convention holds for every method, including the M-ABC family, and it is
the most common thing to get wrong. If a result looks off, check that a
source was not transposed.

``` r

library(MosaiClusteR)
data(mosaic_toy)
L <- mosaic_toy$List
vapply(L, dim, integer(2))   # 60 objects in rows; 120 and 150 features
```

## Normalisation

[`Normalization()`](https://lucp12891.github.io/MosaiClusteR/reference/Normalization.md)
rescales a source so that features sit on a comparable scale before a
distance is computed. Four strategies are available:

| `method` | What it does | Use when |
|----|----|----|
| `"Standardize"` | centre and scale each feature to unit variance | features are on different numeric scales |
| `"Range"` | rescale each feature to a fixed range | you want bounded values and to preserve zeros |
| `"Quantile"` | map features onto a common quantile distribution | distributions differ in shape, not only scale |
| `"Fisher-Yates"` | rank-based normal-scores transform | heavy tails or outliers dominate a raw scale |

``` r

Xn <- Normalization(L[[1]], method = "Quantile")
```

Whether to normalise depends on the source. Expression-style continuous
data usually benefits from standardisation or quantile normalisation.
Binary data (fingerprints, presence or absence) should not be
standardised; pick a distance made for binary data instead (see below).

## Distances

[`Distance()`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)
turns a source into an object-by-object dissimilarity. Choose the
measure to match the data:

| `distmeasure` | For | Notes |
|----|----|----|
| `"euclidean"` | continuous | true Euclidean distance via [`stats::dist`](https://rdrr.io/r/stats/dist.html); applies feature normalisation when `normalize =` is set |
| `"manhattan"` | continuous | sum of absolute differences |
| `"tanimoto"`, `"jaccard"` | binary | overlap of present features |
| `"hamming"` | binary or categorical | count of disagreements |
| `"gower"` | mixed-type | handles continuous and categorical columns together |
| `"cont tanimoto"`, `"MCA_coord"` | specialised | continuous Tanimoto; coordinates from multiple correspondence analysis |

``` r

D1 <- Distance(Xn, distmeasure = "euclidean")       # from already-normalised data
round(as.matrix(D1)[1:4, 1:4], 2)
```

Two equivalent paths give the same result for `"euclidean"`: normalise
first and pass the normalised matrix to
[`Distance()`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md),
or let
[`Distance()`](https://lucp12891.github.io/MosaiClusteR/reference/Distance.md)
normalise internally with its `normalize =` argument. Do one or the
other, not both.

## Preprocess once, or let the method do it

Integration methods accept per-source preprocessing arguments, so you
rarely have to preprocess by hand. A method call carries one entry per
source for `distmeasure`, `normalize` and `method`:

``` r

S  <- length(L)
DM <- rep("euclidean", S)       # one distance per source
NM <- rep(FALSE, S)             # normalise? one flag per source
ME <- rep(list(NULL), S)        # normalisation method per source (NULL = none)

SNF(L, type = "data", distmeasure = DM, normalize = NM, method = ME)
```

If you have already built dissimilarities yourself, pass them with
`type = "dist"` and the method skips its own distance step.

## How sources are aligned

Before combining sources, MosaiClusteR matches them object by object:

- **Row names present**: sources are aligned by object name. Rows in a
  different order are realigned to the first source, so order does not
  matter.
- **No row names**: sources are combined positionally (row 1 with row 1,
  and so on), so they must already be in the same order.
- **Mismatched sets, different row counts, or duplicated names**: these
  raise a clear error that names the problem, rather than failing later
  with a cryptic message.

The practical rule: give every source the same object row names. Then
you can reorder rows freely and trust the alignment.

## Samples missing from a source

If some samples are measured in only part of the sources, use
[`NEMO()`](https://lucp12891.github.io/MosaiClusteR/reference/NEMO.md).
It is the one method here that does not require every object in every
source: it averages neighbourhood affinities over the omics that
measured each pair, with no imputation. Every other method expects the
complete, aligned set.

## Mixing data types across sources

Because the methods accept precomputed dissimilarities, you can
integrate sources of different kinds. Compute a binary distance on a
fingerprint source and a continuous distance on an expression source,
then pass both as distance matrices:

``` r

Dfp <- Distance(fingerprints, distmeasure = "tanimoto")
Dex <- Distance(Normalization(expression, "Standardize"), distmeasure = "euclidean")
WeightedClust(list(Dfp, Dex), type = "dist", weight = seq(1, 0, -0.25))
```
