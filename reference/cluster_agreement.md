# Agreement between two partitions

Computes the adjusted Rand index (ARI), normalised mutual information
(NMI) and pairwise Jaccard between two label vectors of equal length.

## Usage

``` r
cluster_agreement(a, b)
```

## Arguments

- a, b:

  Integer/character/factor label vectors of equal length, without
  missing values.

## Value

A named numeric vector `c(ARI, NMI, Jaccard)`.

## Details

ARI rule of thumb: `>0.90` excellent; `0.75-0.90` good; `0.50-0.75`
moderate; `<0.50` weak.

Edge cases: identical partitions (up to relabelling) return 1 for all
three measures, including the degenerate single-cluster and
all-singleton cases where the ARI and NMI formulae are otherwise
undefined (0/0). `NA` labels are not allowed and raise a clear error
rather than propagating.

## Examples

``` r
cluster_agreement(c(1, 1, 2, 2), c(2, 2, 1, 1))   # identical up to relabel
#>     ARI     NMI Jaccard 
#>       1       1       1 
cluster_agreement(rep(1, 5), rep(1, 5))           # identical, single cluster -> 1
#>     ARI     NMI Jaccard 
#>       1       1       1 
```
