# MosaiClusteR: biological interpretation (genes and pathways)
# https://lucp12891.github.io/MosaiClusteR/
#
# TEMPLATE. This layer needs an expression matrix and gene annotation, which the
# bundled toy data does not provide, so fill in your own data where marked. It
# also needs Bioconductor packages, declared under Suggests.
#
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("limma", "MLP", "biomaRt", "org.Hs.eg.db"))

library(MosaiClusteR)

## ---- your data -------------------------------------------------------------
# L        : list of source matrices, objects (samples) in rows
# geneExpr : a genes-by-samples expression matrix aligned to the objects in L
# GeneInfo : gene annotation for the pathway step (see ?PathwayAnalysis)
#
# L       <- my_sources
# geneExpr <- my_expression
# GeneInfo <- my_gene_info

## Cluster, then confirm the clusters are stable -----------------------------
fit <- SNF(L, type = "data", distmeasure = rep("euclidean", length(L)),
           normalize = rep(FALSE, length(L)), method = rep(list(NULL), length(L)))
lab <- mosaic_labels(fit, 3)

## Differential genes per cluster -------------------------------------------
dg <- DiffGenes(List = L, geneExpr = geneExpr, nrclusters = 3, method = "limma")

## Pathway enrichment per cluster -------------------------------------------
pa <- PathwayAnalysis(List = L, geneExpr = geneExpr, nrclusters = 3,
                      method = c("limma", "MLP"), geneSetSource = "GOBP")

plot_pathways(pa, style = "circular")       # also "dotplot" and "lollipop"

## What the solutions agree on biologically ---------------------------------
shared <- Geneset.intersect(PathwaysOutput = pa)
