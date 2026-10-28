# Pathway analysis for multiple clustering results

A pathway analysis per cluster per method is conducted.

## Usage

``` r
Pathways(
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

## Value

The returned value is a list with an element per cluster per method.

## Details

After finding differently expressed genes, it can be investigated
whether pathways are related to those genes. This can be done with the
help of the function `Pathways` which makes use of the `MLP` function of
the MLP package. Given the output of a method, the cutree function is
performed which results into a specific number of clusters. For each
cluster, the limma method is performed comparing this cluster to the
other clusters. This to obtain the necessary p-values of the genes.
These are used as the input for the `MLP` function to find interesting
pathways. By default the candidate gene sets are determined by the
`AnnotateEntrezIDtoGO` function. The default source will be GOBP, but
this can be altered. Further, it is also possible to provide own
candidate gene sets in the form of a list of pathway categories in which
each component contains a vector of Entrez Gene identifiers related to
that particular pathway. The default values for the minimum and maximum
number of genes in a gene set for it to be considered were used. For MLP
this is respectively 5 and 100. If a list of outputs of several methods
is provided as data input, the cluster numbers are rearranged according
to a reference method. The first method is taken as the reference and
ReorderToReference is applied to get the correct ordering. When the
clusters haven been re-appointed, the pathway analysis as described
above is performed for each cluster of each method.

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

MCF7_PathsFandT=Pathways(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma",
"MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
topG = NULL, GENESET = NULL, sign = 0.05,fusionsLog = TRUE, weightclust = TRUE,
 names =names)
 } # }
```
