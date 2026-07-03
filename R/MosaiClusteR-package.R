#' MosaiClusteR: An Umbrella Framework for Multi-Source and Multi-Omics Clustering
#'
#' \strong{MosaiClusteR} ("MoSaIC" = Multi-Omics Source-Agnostic Integration
#' Clustering in R) unifies a large family of multi-source clustering
#' methodologies behind a single, consistent list-of-matrices interface, and
#' wraps them in a complete analysis framework: preprocessing, single-source
#' baselines, integrative clustering across five paradigms, and downstream
#' comparison/evaluation.
#'
#' @section The MoSaIC framework (how to use the package):
#' \enumerate{
#'   \item \strong{Preprocess} each source with \code{\link{Normalization}} and
#'         \code{\link{Distance}}.
#'   \item \strong{Baseline} each source on its own (any hierarchical clusterer)
#'         to motivate integration.
#'   \item \strong{Integrate} with one or more methods, e.g. \code{\link{M_ABCpp}},
#'         \code{\link{M_ABCdist}}, \code{\link{WeightedClust}}, \code{\link{SNF}},
#'         \code{\link{CEC}}, \code{\link{HierarchicalEnsembleClustering}}.
#'   \item \strong{Evaluate} and compare solutions with
#'         \code{\link{compare_clusterings}} and \code{\link{cluster_agreement}}.
#' }
#' Every clustering method returns a list containing at least \code{DistM} (a
#' dissimilarity matrix) and \code{Clust} (a hierarchical-clustering object),
#' so results are directly comparable and composable.
#'
#' @section Feature weighting with data nuggets:
#' The M-ABC family weights features for subsampling. Beyond the classic
#' variance (\code{"var"}) and coefficient-of-variation (\code{"cv"}) schemes,
#' MosaiClusteR adds \code{weighting = "nugget"}, a robust, Big-Data-friendly
#' alternative based on \code{\link{create_data_nuggets}}. See
#' \code{\link{nugget_feature_weights}} and \code{\link{ABCpp.SingleInMultiple}}.
#'
#' @keywords internal
#' @aliases MosaiClusteR-package
#' @name MosaiClusteR
#' @import stats
#' @import cluster
#' @importFrom matrixStats rowVars rowSds rowMeans2 colVars
#' @importFrom Rdpack reprompt
#' @importFrom Rcpp evalCpp
#' @useDynLib MosaiClusteR, .registration = TRUE
"_PACKAGE"
