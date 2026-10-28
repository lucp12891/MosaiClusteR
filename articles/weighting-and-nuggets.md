# Feature weighting and data nuggets

Not every feature carries signal, and not every source deserves equal
weight. MosaiClusteR lets you weight both. This page covers feature
weighting inside the M-ABC family, the data-nugget scheme behind it, and
source weighting across the integration.

## Feature weighting in M-ABC

The M-ABC family repeats, per source: bootstrap the samples, draw a
**weighted subsample of features**, cluster, and record co-clustering.
The probability that a feature is drawn comes from a weighting scheme:

| `weighting` | Score behind the selection probability | Best when |
|----|----|----|
| `"var"` | feature variance (the classic choice) | signal lives in high-variance features |
| `"cv"` | coefficient of variation | you want a scale-free measure of spread |
| `"nugget"` | data-nugget between-group weighted variance | large N, heavy tails, outliers |
| `"equal"` | flat weights | an ablation, or no prior on features |

``` r

library(MosaiClusteR)
data(mosaic_toy)
L <- mosaic_toy$List

M_ABC(L, NC = 3, weighting = c("var", "var"))                 # variance weighting
M_ABC(L, NC = 3, weighting = c("nugget", "nugget"), seed = 1) # nugget weighting
```

Pass one weighting per source. Variance weighting is a sound default;
switch to `"nugget"` when a few extreme observations would otherwise
dominate the variance.

## Data nuggets

A data nugget summarises a group of observations by three numbers: a
**centre**, a **weight** (how many observations it stands for), and a
**scale** (its internal spread). Compressing the objects into nuggets
and asking which features separate the nugget centres gives a
feature-importance score that rewards real between-group signal and
discounts within-nugget noise. Because it works on a small set of
representatives, it scales to very large N.

``` r

dn <- create_data_nuggets(L[[1]], max_nuggets = 25)   # objects are in rows
dn                                                    # reports the compression
w  <- nugget_feature_weights(dn, type = "between")    # one weight per feature
head(sort(w, decreasing = TRUE))
```

[`nugget_feature_weights()`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_feature_weights.md)
offers different notions of importance through `type`: `"between"`
rewards separation of the nugget centres, while other options weigh
within-group spread. The `"between"` weights are the ones M-ABC uses
when `weighting = "nugget"`.

Two engines build the nuggets. The built-in native engine needs no extra
package. If the reference `datanugget` package is installed, pass
`engine = "datanugget"` to use it instead:

``` r

create_data_nuggets(L[[1]], max_nuggets = 25, engine = "datanugget")
```

You can also pass nugget options straight into M-ABC with `nugget_args`:

``` r

M_ABC(L, NC = 3, weighting = c("nugget", "nugget"),
      nugget_args = list(max_nuggets = 25, seed = 1))
```

## Clustering the nuggets directly

For big-data clustering, cluster the compressed representatives and map
the labels back to the objects.
[`nugget_cluster()`](https://lucp12891.github.io/MosaiClusteR/reference/nugget_cluster.md)
does the whole round trip;
[`Wkmeans()`](https://lucp12891.github.io/MosaiClusteR/reference/Wkmeans.md)
and
[`Whclust()`](https://lucp12891.github.io/MosaiClusteR/reference/Whclust.md)
are the weighted k-means and weighted Ward primitives underneath.

``` r

nugget_cluster(do.call(cbind, L), k = 3, max_nuggets = 40)$cluster
Wkmeans(dn, k = 3)
Whclust(dn, k = 3)
```

[`simulate_weighting_regime()`](https://lucp12891.github.io/MosaiClusteR/reference/simulate_weighting_regime.md)
generates data under a known weighting regime, which is useful for
checking how each scheme behaves before applying it to real data.

## Source weighting

Feature weighting decides which columns matter; source weighting decides
which sources matter. Two routes:

- **A weight grid.**
  [`WeightedClust()`](https://lucp12891.github.io/MosaiClusteR/reference/WeightedClust.md)
  forms a convex combination of the per-source distances over a grid of
  weights, so you can see how the solution changes as one source gains
  influence.
- **Data-driven weights.**
  [`DetermineWeight_SilClust()`](https://lucp12891.github.io/MosaiClusteR/reference/DetermineWeight_SilClust.md)
  and
  [`DetermineWeight_SimClust()`](https://lucp12891.github.io/MosaiClusteR/reference/DetermineWeight_SimClust.md)
  pick source weights from silhouette or similarity criteria, so the
  data chooses the balance.

``` r

WeightedClust(L, type = "data", distmeasure = rep("euclidean", length(L)),
              normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L)),
              weight = seq(1, 0, -0.1))
DetermineWeight_SilClust(L, type = "data")
```
