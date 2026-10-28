# Prepare data for pathway analysis

The `PreparePathway` function prepares the input for the pathway
analysis by computing differentially expressed genes (via limma) for a
selection of objects if these were not provided.

## Usage

``` r
PreparePathway(Object, geneExpr, topG, sign)
```

## Arguments

- Object:

  A clustering output or DiffGenes output element.

- geneExpr:

  The gene expression matrix or ExpressionSet of the objects.

- topG:

  Overrules sign. The number of top genes. Default is NULL.

- sign:

  The significance level to be handled. Default is 0.05.

## Value

A list with the p-values of the genes, the objects and the genes.

## Details

This function may rely on the suggested Bioconductor package 'limma'.

## Examples

``` r
if (FALSE) { # \dontrun{
PreparePathway(Object=MCF7_F,geneExpr=geneMat,topG=NULL,sign=0.05)
} # }
```
