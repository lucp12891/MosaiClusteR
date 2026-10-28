# Column names referenced inside ggplot2 aes() mappings (non-standard evaluation).
# Declared here so R CMD check does not report them as undefined global variables.
utils::globalVariables(c(
  "cluster", "term", "neglogp", "genes", "ratio",   # pathway/target enrichment tables
  "value", "rid", "ang", "hj", "start", "end",      # plot_pathways(style = "circular") helpers
  "status", "y", "k", "xintercept"                  # cluster_biology_panel / k-sweep / selection
))
