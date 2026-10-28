## Submission — MosaiClusteR 0.1.4

### R CMD check results

0 errors | 0 warnings | 0 notes

* Checked locally with `R CMD check --as-cran` on Windows 11, R 4.6.0: Status OK.

### Changes in this version

* Fixed a bug in `CharacteristicFeatures()`: when called with `binData` (as in
  the documented example) an inner loop index overwrote the cluster index, so
  only the first cluster's characteristic features were returned, and a
  `datanames` indexing error left the per-source result list unnamed. Features
  are now computed and named for every cluster.
* Added `enriched_features()`, a helper that turns `CharacteristicFeatures()`
  output into per-cluster tables ready for `plot_pathways()`.
* Added `plot_pathways(style = "circular")` (radial enrichment plot) with new
  `cluster_cols` and `wrap` arguments.
* Added `cexB`/`cexR`/`cexN` label-size arguments to `FeaturesOfCluster()` and
  `BinFeaturesPlot_SingleData()`; `Cyclogram()` now uses a larger default
  `cex.labels`.

### Maintainer note

* The maintainer email address has been updated to b.osangir@gmail.com. The
  confirmation request for this submission will be sent to that address.

### Reverse dependencies

There are no reverse dependencies on CRAN at this time.
