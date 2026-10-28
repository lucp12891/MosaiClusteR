# Pathway analysis with intersection over iterations

The `PathwayAnalysis` function performs pathway analysis multiple times
via `PathwaysIter` and takes the intersection over the iterations via
`Geneset.intersect`.

## Usage

``` r
PathwayAnalysis(
  List,
  Selection = NULL,
  geneExpr = NULL,
  nrclusters = NULL,
  method = c("limma", "MLP"),
  geneInfo = NULL,
  geneSetSource = "GOBP",
  topP = NULL,
  topG = NULL,
  GENESET = NULL,
  sign = 0.05,
  niter = 10,
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL,
  seperatetables = FALSE,
  separatepvals = FALSE
)
```

## Arguments

- List:

  A list of clustering outputs or output of the`DiffGenes` function.

- Selection:

  If pathway analysis should be conducted for a specific selection of
  objects, this selection can be provided here. Default is NULL.

- geneExpr:

  The gene expression matrix or ExpressionSet of the objects.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in. Default is
  NULL.

- method:

  The methods for gene and pathway analysis. Default is
  c("limma","MLP").

- geneInfo:

  A data frame with at least the columns ENTREZID and SYMBOL. Default is
  NULL.

- geneSetSource:

  The source for the getGeneSets function, defaults to "GOBP".

- topP:

  Overrules sign. The number of pathways to display for each cluster.
  Default is NULL.

- topG:

  Overrules sign. The number of top genes to be returned. Default is
  NULL.

- GENESET:

  Optional. Can provide own candidate gene sets. Default is NULL.

- sign:

  The significance level to be handled. Default is 0.05.

- niter:

  The number of times to perform pathway analysis. Default is 10.

- fusionsLog:

  Logical. To be handed to `ReorderToReference`. Default is TRUE

- weightclust:

  Logical. To be handed to `ReorderToReference`. Default is TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

- seperatetables:

  Logical. Default is FALSE.

- separatepvals:

  Logical. Default is FALSE.

## Value

The intersection over the iterations of the pathway analysis.

## Details

This function relies on the suggested Bioconductor packages 'MLP',
'biomaRt' and 'org.Hs.eg.db'.

## Examples

``` r
if (FALSE) { # \dontrun{
MCF7_PathsFandT=PathwayAnalysis(List=L, geneExpr = geneMat, nrclusters = 7,
method = c("limma","MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP",niter=2)
} # }
```
