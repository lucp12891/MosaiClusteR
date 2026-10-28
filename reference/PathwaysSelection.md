# Pathway analysis for a selection of objects

Internal function of `Pathways`.

## Usage

``` r
PathwaysSelection(
  List = NULL,
  Selection,
  geneExpr = NULL,
  nrclusters = NULL,
  method = c("limma", "MLP"),
  geneInfo = NULL,
  geneSetSource = "GOBP",
  topP = NULL,
  topG = NULL,
  GENESET = NULL,
  sign = 0.05,
  fusionsLog = TRUE,
  weightclust = TRUE,
  names = NULL
)
```

## Arguments

- List:

  A list of clustering outputs or output of the`DiffGenes` function. The
  first element of the list will be used as the reference in
  `ReorderToReference`. The output of `ChooseFeatures` is also accepted.

- Selection:

  If pathway analysis should be conducted for a specific selection of
  objects, this selection can be provided here. Selection can be of the
  type "character" (names of the objects) or "numeric" (the number of
  specific cluster). Default is NULL.

- geneExpr:

  The gene expression matrix or ExpressionSet of the objects. The rows
  should correspond with the genes.

- nrclusters:

  Optional. The number of clusters to cut the dendrogram in. The number
  of clusters should not be specified if the interest lies only in a
  specific selection of objects which is known by name. Otherwise, it is
  required. Default is NULL.

- method:

  The method to applied to look for differentially expressed genes and
  related pathways. For now, only the limma method is available for gene
  analysis and the MLP method for pathway analysis. Default is
  c("limma","MLP").

- geneInfo:

  A data frame with at least the columns ENTREZID and SYMBOL. This is
  necessary to connect the symbolic names of the genes with their
  EntrezID in the correct order. The order of the gene is here not in
  the order of the rownames of the gene expression matrix but in the
  order of their significance. Default is NULL.

- geneSetSource:

  The source for the getGeneSets function, defaults to "GOBP".

- topP:

  Overrules sign. The number of pathways to display for each cluster. If
  not specified, only the significant genes are shown. Default is NULL.

- topG:

  Overrules sign. The number of top genes to be returned in the result.
  If not specified, only the significant genes are shown. Defaults is
  NULL.

- GENESET:

  Optional. Can provide own candidate gene sets. Default is NULL.

- sign:

  The significance level to be handled. Default is 0.05.

- fusionsLog:

  Logical. To be handed to `ReorderToReference`: indicator for the
  fusion of clusters. Default is TRUE

- weightclust:

  Logical. To be handed to `ReorderToReference`: to be used for the
  outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only the
  result of the Clust element is considered. Default is TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

## Details

This function relies on the suggested Bioconductor packages 'MLP',
'biomaRt' and 'org.Hs.eg.db'.

## Examples

``` r
if (FALSE) { # \dontrun{
PathwaysSelection(List=L,Selection=c("Cpd1","Cpd2"),geneExpr=geneMat,geneInfo=GeneInfo)
} # }
```
