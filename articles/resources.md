# Resources and attribution

## Where to find MosaiClusteR

- **CRAN** (stable release):
  <https://CRAN.R-project.org/package=MosaiClusteR>
- **Source and development** (GitHub):
  <https://github.com/lucp12891/MosaiClusteR>
- **Report a bug or request a feature**:
  <https://github.com/lucp12891/MosaiClusteR/issues>
- **Documentation site**: <https://lucp12891.github.io/MosaiClusteR/>

## Worked example scripts

Runnable `.R` versions of the worked examples, ready to download and
adapt:

- [quickstart.R](https://lucp12891.github.io/MosaiClusteR/scripts/quickstart.R):
  the 60-second quick start.
- [end-to-end.R](https://lucp12891.github.io/MosaiClusteR/scripts/end-to-end.R):
  a complete analysis, from preprocessing to characterisation.
- [compare-methods.R](https://lucp12891.github.io/MosaiClusteR/scripts/compare-methods.R):
  run one method per paradigm and compare them by ARI.

Each script runs on the bundled `mosaic_toy` data; swap in your own list
of matrices (objects in rows) to use it on real data.

## How to cite

Please cite the package when it contributes to published work.

``` r

citation("MosaiClusteR")
```

``` bibtex
@Manual{MosaiClusteR2026,
  title  = {MosaiClusteR: An Umbrella Framework for Multi-Source and
            Multi-Omics Clustering},
  author = {Osang'ir, Bernard Isekah and Van Moerbeke, Marijke and
            Claesen, J\"urgen and Gupta, Surya and Cabrera, Javier and
            Amaratunga, Dhammika and Shkedy, Ziv},
  year   = {2026},
  note   = {R package version 1.0.0},
  url    = {https://github.com/lucp12891/MosaiClusteR}
}
```

If you use a specific method, please also cite the original work behind
it, listed in the references below.

## Packages MosaiClusteR builds on

MosaiClusteR implements many methods natively, and interoperates with or
wraps a number of other packages. These are declared under *Suggests*
and loaded only when the relevant method is called. Credit for the
wrapped functionality belongs to the authors of these packages.

| Package | Role in MosaiClusteR | Source |
|----|----|----|
| **SNFtool** | Similarity Network Fusion (`SNF`) | <https://CRAN.R-project.org/package=SNFtool> |
| **datanugget** | Reference data-nugget engine (`create_data_nuggets(engine = "datanugget")`) | <https://CRAN.R-project.org/package=datanugget> |
| **WCluster** | Weighted clustering primitives used with data nuggets | <https://CRAN.R-project.org/package=WCluster> |
| **LUCIDus** | Outcome-guided latent clustering (`LUCID`) | <https://CRAN.R-project.org/package=LUCIDus> |
| **FactoMineR**, **lsa**, **analogue** | Extra distances (MCA, cosine, chi-squared) | CRAN |
| **circlize**, **plotrix**, **gplots**, **dendextend** | Dendrograms, cyclograms, heatmaps | CRAN |
| **uwot**, **Rtsne** | UMAP and t-SNE embeddings for visualisation | CRAN |
| **limma**, **MLP**, **biomaRt**, **org.Hs.eg.db** | Differential expression and pathway analysis | <https://bioconductor.org> |

The core of the package (the integration paradigms, the M-ABC family,
the native data-nugget engine, evaluation and the
`List -> {DistM, Clust}` interface) is pure R and compiled C++, with no
dependency beyond base R and a small set of imports.

## Method references

- Amaratunga, Cabrera and Shkedy (2008). *Exploration and Analysis of
  DNA Microarray and Other High-Dimensional Data.* Wiley. (Aggregating
  Bundles of Clusters, the basis of the M-ABC family.)
- Cherasia, Cabrera, Fernholz and Fernholz (2023). *Data Nuggets in
  Supervised Learning.* In *Robust and Multivariate Statistical
  Methods*, Springer, pages 429 to 449.
- Wang et al. (2014). *Similarity Network Fusion for aggregating data
  types on a genomic scale.* Nature Methods 11, 333 to 337.
- Rappoport and Shamir (2019). *NEMO: cancer subtyping by integration of
  partial multi-omic data.* Bioinformatics.
- Perualila-Tan et al. (2016). *Weighted similarity-based clustering of
  chemical structures and bioactivity data in early drug discovery.*
- Fodeh, Punch and Tan (2013). *On unifying multi-source information for
  biomedical text and data mining.*
- Strehl and Ghosh (2002). *Cluster ensembles: a knowledge reuse
  framework for combining multiple partitions.* (CSPA, HGPA, MCLA.)
- Fern and Brodley (2004). *Solving cluster ensemble problems by
  bipartite graph partitioning.* (HBGF.)
- Gionis, Mannila and Tsaparas (2007). *Clustering aggregation.*
- Nguyen and Caruana (2007). *Consensus clusterings.* (IVC, IPVC, IPC.)
- Fred and Jain (2005). *Combining multiple clusterings using evidence
  accumulation.*
- Iam-On et al. (2010). *Link-based cluster ensembles.* (CTS, SRS,
  ASRS.)
- Saeed et al. (2012, 2014). *Cumulative voting consensus.* (CVAA,
  W-CVAA.)
- Zheng, Chen and Zheng (2014). *A novel approach for hierarchical
  ensemble clustering.* IEEE TPAMI.
- Hossain et al. (2012). *Ensemble for hierarchical clustering.*
- Chalise and Fridley (2017). *Integrative clustering of multi-level
  omic data based on non-negative matrix factorisation.* (intNMF.)
- Zhao, Jia, Goodrich and Conti (2024). *LUCIDus: an R package for
  latent unknown clustering with integrated data.* The R Journal 16(2).
- Zhang, Zhou, Xu and Liu (2022). *Integrative clustering methods for
  multi-omics data.* WIREs Computational Statistics 14, e1553. (Review
  covering SNF, NEMO, the iCluster family, ab-SNF and CIMLR.)

## Licence

MosaiClusteR is released under the **GPL-3** licence.
