# Intersection over resulting gene sets of `PathwaysIter` function

The function `Geneset.intersect` collects the results of the
`PathwaysIter` function per method for each cluster and takes the
intersection over the iterations per cluster per method. This is to see
if over the different resamplings of the data, similar pathways were
discovered.

## Usage

``` r
Geneset.intersect(
  PathwaysOutput,
  Selection = FALSE,
  sign = 0.05,
  names = NULL,
  seperatetables = FALSE,
  separatepvals = FALSE
)
```

## Arguments

- PathwaysOutput:

  The output of the `PathwaysIter` function.

- Selection:

  Logical. Indicates whether or not the output of the pathways function
  were concentrated on a specific selection of objects. If this was the
  case, Selection should be put to TRUE. Otherwise, it should be put to
  FALSE. Default is TRUE.

- sign:

  The significance level to be handled for cutting of the pathways.
  Default is 0.05.

- names:

  Optional. Names of the methods. Default is NULL.

- seperatetables:

  Logical. If TRUE, a separate element is created per cluster containing
  the pathways for each iteration. Default is FALSE.

- separatepvals:

  Logical. If TRUE, the p-values of the each iteration of each pathway
  in the intersection is given. If FALSE, only the mean p-value is
  provided. Default is FALSE.

## Value

The output is a list with an element per method. For each method, it is
portrayed per cluster which pathways belong to the intersection over all
iterations and their corresponding mean p-values.

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

MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr = geneMat, nrclusters = 7, method =
c("limma", "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE,
weightclust = TRUE, names =names)

MCF7_Paths_intersection=Geneset.intersect(PathwaysOutput=MCF7_Paths_FandT,
sign=0.05,names=c("FP","TP"),seperatetables=FALSE,separatepvals=FALSE)

str(MCF7_Paths_intersection)
} # }
```
