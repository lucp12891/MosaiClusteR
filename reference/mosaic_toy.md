# Toy multi-omics example data

A small two-source example produced by
[`mosaic_sim`](https://lucp12891.github.io/MosaiClusteR/reference/mosaic_sim.md)
with three planted groups, bundled for examples and tests.

## Usage

``` r
data(mosaic_toy)
```

## Format

A list with two elements:

- List:

  A list of two numeric matrices with the 60 shared objects in rows and
  features in columns (120 and 150 features).

- truth:

  Integer vector of the 60 ground-truth group labels.

## Examples

``` r
data(mosaic_toy)
str(mosaic_toy, max.level = 2)
#> List of 2
#>  $ List :List of 2
#>   ..$ source1: num [1:60, 1:120] 2.47 2.92 1.97 3.35 2.69 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   ..$ source2: num [1:60, 1:150] 2.358 2.819 0.784 2.313 0.958 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>  $ truth: Named int [1:60] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "names")= chr [1:60] "obj01" "obj02" "obj03" "obj04" ...
```
