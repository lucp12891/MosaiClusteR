# Intersection of genes and pathways over multiple methods for a selection of objects.

Internal function of `SharedGenesPathsFeat`.

## Usage

``` r
SharedSelection(
  DataLimma = NULL,
  DataMLP = NULL,
  DataFeat = NULL,
  names = NULL
)
```

## Arguments

- DataLimma:

  Optional. The output of a `DiffGenes` function. Default is NULL.

- DataMLP:

  Optional. The output of `Geneset.intersect` function. Default is NULL.

- DataFeat:

  Optional. The output of `CharacteristicFeatures` function. Default is
  NULL.

- names:

  Optional. Names of the methods or "Selection" if it only considers a
  selection of objects. Default is NULL.

## Value

A list with a Table and a Which element for the selection.

## Examples

``` r
if (FALSE) { # \dontrun{
SharedSelection(DataLimma=MCF7_FT_DE,names="Selection")
} # }
```
