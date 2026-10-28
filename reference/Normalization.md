# Normalisation of features

When data of different scales are combined it is recommended to
normalise so that the structures are comparable.

## Usage

``` r
Normalization(
  Data,
  method = c("Quantile", "Fisher-Yates", "Standardize", "Range", "Q", "q", "F", "f", "S",
    "s", "R", "r")
)
```

## Arguments

- Data:

  A data matrix; rows are assumed to correspond to the objects.

- method:

  One of `"Quantile"`, `"Fisher-Yates"`, `"Standardize"`, `"Range"` (or
  any unambiguous first letter).

## Value

The normalised matrix.

## Details

`"Quantile"` performs quantile normalisation (common in omics);
`"Fisher-Yates"` uses normal scores based only on the number of rows;
`"Standardize"` centres and scales each column (base
[`scale`](https://rdrr.io/r/base/scale.html)); `"Range"` maps the matrix
to `[0, 1]`.

## Examples

``` r
x <- matrix(rnorm(100), 10, 10)
Norm_x <- Normalization(x, method = "R")
```
