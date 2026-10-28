# MosaiClusteR: 60-second quick start
# https://lucp12891.github.io/MosaiClusteR/
#
# install.packages("MosaiClusteR")   # or remotes::install_github("lucp12891/MosaiClusteR")

library(MosaiClusteR)

data(mosaic_toy)                 # two sources, 60 shared objects, 3 planted groups
L     <- mosaic_toy$List         # a list of matrices, objects in rows
truth <- mosaic_toy$truth

## Integrate: every method shares one interface, so swapping it changes the paradigm
fit <- SNF(L, type = "data", distmeasure = c("euclidean", "euclidean"),
           normalize = c(FALSE, FALSE), method = list(NULL, NULL))

labels <- mosaic_labels(fit, k = 3)      # pull a flat partition from any fitted result
cluster_agreement(labels, truth)         # ARI, NMI, Jaccard, purity vs ground truth

## Same score from a different paradigm, e.g. a neighbourhood method
cluster_agreement(NEMO(L, k = 3)$cluster, truth)

## Compare two integrations directly
compare_clusterings(fit$DistM,
                    ADC(L, distmeasure = "euclidean", normalize = FALSE)$DistM, NC = 3)
