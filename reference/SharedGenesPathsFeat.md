# Shared genes, pathways and features across methods

The `SharedGenesPathsFeat` function takes the intersection of the
differentially expressed genes, the pathways and the characteristic
features over the methods per cluster.

## Usage

``` r
SharedGenesPathsFeat(
  DataLimma = NULL,
  DataMLP = NULL,
  DataFeat = NULL,
  names = NULL,
  Selection = FALSE
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

  Optional. Names of the methods. Default is NULL.

- Selection:

  Logical. Whether a selection of objects was considered. Default is
  FALSE.

## Value

A list with a Table summarizing the shared elements and a Which element
detailing the shared genes, pathways and features per cluster.

## Examples

``` r
if (FALSE) { # \dontrun{
SharedGenesPathsFeat(DataLimma=MCF7_FT_DE,DataMLP=MCF7_Paths_intersection,names=c("FP","TP"))
} # }
```
