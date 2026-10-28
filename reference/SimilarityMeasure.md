# Similarity of a set of clusterings to a reference

Given a colour/label matrix (clusterings in rows, objects in columns) or
a list of clustering outputs, the proportion of objects sharing the
label of the first (reference) row is computed for every row. When a
list is supplied it is first aligned with `ReorderToReference`;
supplying the label matrix directly avoids that dependency.

## Usage

``` r
SimilarityMeasure(
  List,
  nrclusters = NULL,
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL
)
```

## Arguments

- List:

  A label matrix (rows are clusterings, columns are objects) or a list
  of clustering outputs to be aligned via `ReorderToReference`.

- nrclusters:

  Number of clusters used when `List` is a list. Default `NULL`.

- fusionsLog, weightclust:

  Logical controls passed to `ReorderToReference` when `List` is a list.

- names:

  Optional method names passed to `ReorderToReference`.

## Value

A numeric vector of similarities in `[0, 1]`, one per row, where the
first entry is the reference compared to itself (always `1`).

## References

Fodeh, S. J., Brandt, C., Luong, T. B., Haddad, A., Schultz, M., Murphy,
T. and Krauthammer, M. (2013). Complementary ensemble clustering of
biomedical data. Journal of Biomedical Informatics, 46(3), 436-443.

## Examples

``` r
if (FALSE) { # \dontrun{
M <- rbind(c(1, 1, 2, 2, 3, 3),
           c(1, 1, 2, 3, 3, 3),
           c(2, 2, 1, 1, 3, 3))
SimilarityMeasure(M)
} # }
```
