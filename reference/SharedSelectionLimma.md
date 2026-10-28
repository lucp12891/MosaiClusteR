# Intersection of genes over multiple methods for a selection of objects.

Internal function of `SharedGenesPathsFeat`.

## Usage

``` r
SharedSelectionLimma(DataLimma = NULL, names = NULL)
```

## Arguments

- DataLimma:

  Optional. The output of a `DiffGenes` function. Default is NULL.

- names:

  Optional. Names of the methods or "Selection" if it only considers a
  selection of objects. Default is NULL.

## Value

A list with a Table and a Which element for the shared genes.

## Examples

``` r
if (FALSE) { # \dontrun{
SharedSelectionLimma(DataLimma=MCF7_FT_DE)
} # }
```
