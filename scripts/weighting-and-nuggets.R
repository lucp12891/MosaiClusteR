# MosaiClusteR: feature weighting and data nuggets
# https://lucp12891.github.io/MosaiClusteR/
#
# Shows the feature-weighting schemes in the M-ABC family, the data-nugget
# compression behind "nugget" weighting, and clustering the nuggets directly.

library(MosaiClusteR)

data(mosaic_toy)
L     <- mosaic_toy$List
truth <- mosaic_toy$truth

## Feature weighting inside M-ABC -------------------------------------------
# One weighting per source: "var", "cv", "nugget" or "equal".
abc_var    <- M_ABC(L, NC = 3, weighting = c("var", "var"),       seed = 1)
abc_nugget <- M_ABC(L, NC = 3, weighting = c("nugget", "nugget"), seed = 1)

c(var    = cluster_agreement(mosaic_labels(abc_var, 3),    truth)[["ARI"]],
  nugget = cluster_agreement(mosaic_labels(abc_nugget, 3), truth)[["ARI"]])

## Data nuggets: compression and feature importance -------------------------
# Objects are in rows, so pass the source directly (no transpose).
dn <- create_data_nuggets(L[[1]], max_nuggets = 25)
dn                                                   # reports the compression
w  <- nugget_feature_weights(dn, type = "between")   # one weight per feature
head(sort(w, decreasing = TRUE))

## Cluster the compressed representatives directly (big-data clustering) -----
nugget_cluster(do.call(cbind, L), k = 3, max_nuggets = 40)$cluster
Wkmeans(dn, k = 3)
Whclust(dn, k = 3)

## Source weighting ----------------------------------------------------------
S <- length(L)
WeightedClust(L, type = "data", distmeasure = rep("euclidean", S),
              normalize = rep(FALSE, S), method = rep(list(NULL), S),
              weight = seq(1, 0, -0.1))
DetermineWeight_SilClust(L, type = "data")
