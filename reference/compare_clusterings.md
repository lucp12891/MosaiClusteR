# Compare two clusterings derived from distance matrices

Clusters two dissimilarity matrices with the same linkage and number of
clusters and reports their agreement (ARI, NMI, pairwise Jaccard) and
the contingency table.

## Usage

``` r
compare_clusterings(
  D1,
  D2,
  NC = NULL,
  linkage = "ward.D2",
  labels = c("D1", "D2"),
  verbose = TRUE
)
```

## Arguments

- D1, D2:

  Square dissimilarity matrices over the same objects.

- NC:

  Number of clusters (default `floor(sqrt(n))`).

- linkage:

  Linkage method. Default `"ward.D2"`.

- labels:

  Length-2 labels for the two solutions.

- verbose:

  Logical; print a short report.

## Value

Invisibly, a list with the two label vectors, `ARI`, `NMI`,
`pair_Jaccard` and the `contingency` table.

## Examples

``` r
set.seed(1)
D1 <- as.matrix(dist(matrix(rnorm(60), 20)))
D2 <- as.matrix(dist(matrix(rnorm(60), 20)))
compare_clusterings(D1, D2, NC = 3, verbose = FALSE)$ARI
#>       ARI 
#> 0.2024809 
```
