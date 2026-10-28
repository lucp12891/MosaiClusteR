# Find the genes shared across methods

The `FindGenes` function collects the differentially expressed genes per
cluster per method (the output of `DiffGenes`) and determines which
genes are found across the methods.

## Usage

``` r
FindGenes(dataLimma, names = NULL)
```

## Arguments

- dataLimma:

  The output of the `DiffGenes` function.

- names:

  Optional. Names of the methods. Default is NULL.

## Value

The result is a list with for each cluster a list of the found genes and
the methods that found them.

## Examples

``` r
if (FALSE) { # \dontrun{
MCF7_SharedGenes=FindGenes(dataLimma=MCF7_FT_DE,names=c("FP","TP"))
} # }
```
