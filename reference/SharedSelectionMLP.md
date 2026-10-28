# Intersection of pathways over multiple methods for a selection of objects.

Internal function of `SharedGenesPathsFeat`.

## Usage

``` r
SharedSelectionMLP(DataMLP = NULL, names = NULL)
```

## Arguments

- DataMLP:

  Optional. The output of `Geneset.intersect` function. Default is NULL.

- names:

  Optional. Names of the methods or "Selection" if it only considers a
  selection of objects. Default is NULL.

## Value

A list with a Table and a Which element for the shared pathways.

## Examples

``` r
if (FALSE) { # \dontrun{
SharedSelectionMLP(DataMLP=MCF7_Paths_intersection)
} # }
```
