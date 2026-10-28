# FAQ and troubleshooting

Short answers to the questions that come up most often. If something
here does not cover your case, please open an issue at
<https://github.com/lucp12891/MosaiClusteR/issues>.

## Installing and loading

**Do I need a compiler?** The M-ABC family ships C++ kernels, so
building from source needs a C/C++ toolchain:
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows, the
Xcode command-line tools on macOS. The CRAN binaries already include the
compiled code, so `install.packages("MosaiClusteR")` needs no toolchain.
Every method also has a pure-R path, so the package still loads if a
kernel is unavailable.

**A method says a package is “required here”.** Heavier and specialised
methods declare their dependencies under *Suggests* and load them only
when called. The message names the exact install command, for example
`install.packages("SNFtool")` or `BiocManager::install("MLP")`. Install
the named package and run the call again.

**`could not find function "M_ABCpp"`.** The function was renamed to
[`M_ABC()`](https://lucp12891.github.io/MosaiClusteR/reference/M_ABC.md)
in 0.1.4. Update `M_ABCpp(...)` to `M_ABC(...)`; the arguments and
behaviour are unchanged.

## Data shape and alignment

**Which way round do my matrices go?** Objects (samples) in **rows**,
features in **columns**, for every source and every method, including
the M-ABC family. This is the single most common mistake. If a result
looks wrong or an example compresses the wrong dimension, check that you
did not transpose a source.

``` r

# right: 60 objects in rows, 120 features in columns
dim(mosaic_toy$List[[1]])   #> 60 120
```

**“the sources do not describe the same objects” (or a subscript error
on older versions).** Integration needs every source to cover the same
objects. MosaiClusteR aligns the sources by row name and reports clearly
when they do not match: different object sets, different row counts, or
duplicated identifiers. Give each source row names, and make sure the
sets agree. Reordered rows are fine; they are realigned for you.

**My samples are not all measured in every source.** Use
[`NEMO()`](https://lucp12891.github.io/MosaiClusteR/reference/NEMO.md).
It is the one method here that handles samples missing a modality, by
averaging affinities over the omics that measured each pair, with no
imputation. The other methods need the same objects across all sources.

## Reproducibility

**I get a different answer each run.** Several methods are stochastic:
the M-ABC family,
[`ADEC()`](https://lucp12891.github.io/MosaiClusteR/reference/ADEC.md),
[`intNMF()`](https://lucp12891.github.io/MosaiClusteR/reference/intNMF.md),
and any step that uses k-means initialisation
([`spectral_clustering()`](https://lucp12891.github.io/MosaiClusteR/reference/spectral_clustering.md),
[`NEMO()`](https://lucp12891.github.io/MosaiClusteR/reference/NEMO.md),
[`HBGF()`](https://lucp12891.github.io/MosaiClusteR/reference/HBGF.md)).
Fix the result with the method’s `seed =` argument where it has one, or
with [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call.

``` r

set.seed(1)
M_ABC(mosaic_toy$List, NC = 3, seed = 1)   # seed = also fixes it directly
```

## Choosing the number of clusters

**[`SelectnrClusters()`](https://lucp12891.github.io/MosaiClusteR/reference/SelectnrClusters.md)
returns a table, not a number.** It gives one recommendation per
criterion (a one-row data frame). Reduce it to a single integer, for
example by taking the most frequent recommendation, before using it as a
cluster count.

``` r

sel  <- SelectnrClusters(DistM, maxK = 8)
recs <- as.integer(unlist(sel$Optimal_Nr_of_CLusters))
k    <- as.integer(names(sort(table(recs), decreasing = TRUE))[1])
```

**How do I read an agreement score?** For
[`cluster_agreement()`](https://lucp12891.github.io/MosaiClusteR/reference/cluster_agreement.md)
ARI: above 0.90 excellent, 0.75 to 0.90 good, 0.50 to 0.75 moderate,
below 0.50 weak. ARI, NMI and Jaccard all return 1 for identical
partitions (up to relabelling), including the single-cluster and
all-singleton cases.

## Methods that need external tools

**METIS or MATLAB.** A few methods use an external graph-partitioning
executable:
[`EnsembleClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/EnsembleClustering.md),
the METIS path of
[`EHC()`](https://lucp12891.github.io/MosaiClusteR/reference/EHC.md),
and `LinkBasedClustering(linkBasedMethod = "asrs")` (MATLAB). They raise
an informative error when the tool is absent. Pure-R alternatives cover
the same need:
[`HBGF()`](https://lucp12891.github.io/MosaiClusteR/reference/HBGF.md),
[`EvidenceAccumulation()`](https://lucp12891.github.io/MosaiClusteR/reference/EvidenceAccumulation.md),
and the CTS or SRS paths of
[`LinkBasedClustering()`](https://lucp12891.github.io/MosaiClusteR/reference/LinkBasedClustering.md).

## Omics interpretation

**Pathway and gene-set functions need Bioconductor.**
[`DiffGenes()`](https://lucp12891.github.io/MosaiClusteR/reference/DiffGenes.md),
[`PathwayAnalysis()`](https://lucp12891.github.io/MosaiClusteR/reference/PathwayAnalysis.md),
[`Pathways()`](https://lucp12891.github.io/MosaiClusteR/reference/Pathways.md)
and the related functions use Bioconductor packages (`limma`, `MLP`,
`biomaRt`, `org.Hs.eg.db`). Install them with
`BiocManager::install(...)` as the error message names.

## Performance

**M-ABC is slow on my data.** Runtime grows with the number of
simulations (`numsim`), the number of sources, and the sample size.
Start with a modest `numsim`, and for large N or heavy-tailed data use
`weighting = "nugget"`, which derives feature weights from a compact set
of representatives rather than from every observation.
