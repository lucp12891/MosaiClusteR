# Differential expression for a selection of objects

Internal function of `DiffGenes`.

## Usage

``` r
DiffGenesSelection(
  List,
  Selection,
  geneExpr = NULL,
  nrclusters = NULL,
  method = "limma",
  sign = 0.05,
  topG = NULL,
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL
)
```

## Arguments

- List:

  A list of the clustering outputs to be compared. The first element of
  the list will be used as the reference in `ReorderToReference`.

- Selection:

  If differential gene expression should be investigated for a specific
  selection of objects, this selection can be provided here. Selection
  can be of the type "character" (names of the objects) or "numeric"
  (the number of specific cluster). Default is NULL.

- geneExpr:

  The gene expression matrix or ExpressionSet of the objects. The rows
  should correspond with the genes.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in. The number
  of clusters should not be specified if the interest lies only in a
  specific selection of objects which is known by name. Otherwise, it is
  required. Default is NULL.

- method:

  The method to applied to look for DE genes. For now, only the limma
  method is available. Default is "limma".

- sign:

  The significance level to be handled. Default is 0.05.

- topG:

  Overrules sign. The number of top genes to be shown. Default is NULL.

- fusionsLog:

  Logical. To be handed to `ReorderToReference`: indicator for the
  fusion of clusters. Default is TRUE

- weightclust:

  Logical. To be handed to `ReorderToReference`: to be used for the
  outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only the
  result of the Clust element is considered. Default is TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

## Examples

``` r
if (FALSE) { # \dontrun{
MCF7_FT_DE = DiffGenes(List=L,Selection=c("Cpd1","Cpd2"),geneExpr=geneMat,method="limma")
} # }
```
