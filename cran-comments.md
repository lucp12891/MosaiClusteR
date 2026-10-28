## Submission: MosaiClusteR 1.0.0

### R CMD check results

0 errors | 0 warnings | 0 notes

Checked locally with `R CMD check --as-cran` on Windows 11, R 4.6.0: Status OK.

### Changes in this version

This is the first stable (1.0.0) release. It is a correctness release: there are
no new exported functions and no change to the package interface or to the
declared dependencies.

* `Cluster()` now forwards its normalisation method to `Distance()`. Previously
  the `method` argument was accepted but had no effect.
* `cluster_agreement()` returns 1 for identical partitions (including the
  single-cluster and all-singleton cases) and raises a clear error when labels
  contain `NA`.
* `ADC()`, `ADEC()`, `M_ABC()` and `CEC()` validate and align their multi-source
  inputs by object identifier, with clear errors for mismatched, reordered or
  duplicated objects, and accept sources without row names (combined
  positionally).
* `mosaic_labels()` validates the number of clusters `k`.
* `M_ABC()` and `ADEC()` gained a `seed` argument for reproducible results.
* Vignette and documentation examples were corrected to match the above.

### Reverse dependencies

There are no reverse dependencies on CRAN at this time.
