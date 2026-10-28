# Intersection over resulting gene sets of `PathwaysIter` function for a selection of objects

Internal function of `Geneset.intersect`.

## Usage

``` r
Geneset.intersectSelection(
  list.output,
  sign = 0.05,
  names = NULL,
  seperatetables = FALSE,
  separatepvals = FALSE
)
```

## Arguments

- list.output:

  The output of the `PathwaysIter` function.

- sign:

  The significance level to be handled for cutting of the pathways.
  Default is 0.05.

- names:

  Optional. Names of the methods. Default is NULL.

- seperatetables:

  Logical. If TRUE, a separate element is created per cluster containing
  the pathways for each iteration. Default is FALSE.

- separatepvals:

  Logical. If TRUE, the p-values of the each iteration of each pathway
  in the intersection is given. If FALSE, only the mean p-value is
  provided. Default is FALSE.

## Examples

``` r
if (FALSE) { # \dontrun{
Geneset.intersectSelection(MCF7_Paths_FandT,sign=0.05)
} # }
```
