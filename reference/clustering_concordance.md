# Compare clustering solutions to each other (no ground truth needed)

Given several labellings of the same objects (e.g. different weightings
or methods), compute their pairwise agreement (ARI) and each object's
stability – how consistently it is co-clustered across the solutions.
Useful for real data where no ground truth is available.

## Usage

``` r
clustering_concordance(labels_list)
```

## Arguments

- labels_list:

  A named list of label vectors, all of the same length.

## Value

A list with `agreement` (pairwise ARI matrix), `stability` (mean
pairwise co-membership per object) and `changed` (objects whose label is
not constant across solutions).

## See also

[`cluster_agreement`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md),
[`compare_clusterings`](https://lucp12891.github.io/MosaiClusteR/reference/compare_clusterings.md)

## Examples

``` r
data(mosaic_toy)
v <- mosaic_labels(M_ABC(mosaic_toy$List, numsim = 40, NC = 3,
                           weighting = c("var", "var")), 3)
#> M_ABC: source 1 of 2 (weighting = var) ...
#> M_ABC: source 2 of 2 (weighting = var) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
n <- mosaic_labels(M_ABC(mosaic_toy$List, numsim = 40, NC = 3,
                           weighting = c("nugget", "nugget")), 3)
#> M_ABC: source 1 of 2 (weighting = nugget) ...
#> M_ABC: source 2 of 2 (weighting = nugget) ...
#> M_ABC: building consensus dissimilarity (C++ kernel) ...
clustering_concordance(list(var = v, nugget = n))$agreement
#>        var nugget
#> var      1      1
#> nugget   1      1
```
