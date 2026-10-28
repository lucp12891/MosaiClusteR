# Iterations of the pathway analysis

The MLP method to perform pathway analysis is based on resampling of the
data. Therefore it is recommended to perform the pathway analysis
multiple times to observe how much the results are influenced by a
different resample. The function `PathwaysIter` performs the pathway
analysis as described in `Pathways` a specified number of times. The
input can be one data set or a list as in `Pathway.2` and `Pathways`.

## Usage

``` r
PathwaysIter(
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

  The gene expression matrix of the objects. The rows should correspond
  with the genes.

- nrclusters:

  The number of clusters to cut the dendrogram in. Default is NULL.

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

  The source for the getGeneSets function ("GOBP", "GOMF","GOCC", "KEGG"
  or "REACTOME"). Default is "GOBP".

- topP:

  Overrules sign. The number of pathways to display for each cluster. If
  not specified, only the significant genes are shown. Default is NULL.

- topG:

  Overrules sign. The number of top genes to be returned in the result.
  If not specified, only the significant genes are shown. Default is
  NULL.

- GENESET:

  Optional. Can provide own candidate gene sets. Default is NULL.

- sign:

  The significance level to be handled. Default is 0.05.

- niter:

  The number of times to perform pathway analysis. Default is 10.

- fusionsLog:

  Logical. To be handed to `ReorderToReference`: indicator for the
  fusion of clusters. Default is TRUE

- weightclust:

  Logical. To be handed to `ReorderToReference`: to be used for the
  outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only the
  result of the Clust element is considered. Default is TRUE.

- names:

  Optional. Names of the methods. Default is NULL.

## Value

A list with an element per iteration, each being the output of
`Pathways`.

## Details

This function relies on the suggested Bioconductor packages 'MLP',
'biomaRt' and 'org.Hs.eg.db'.

## Examples

``` r

if (FALSE) { # \dontrun{
data(fingerprintMat)
data(targetMat)
data(geneMat)
data(GeneInfo)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr = geneMat, nrclusters = 7, method =
c("limma", "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE,
weightclust = TRUE, names =names)
} # }
```
