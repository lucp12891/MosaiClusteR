# =============================================================================
# ABCpp.R
# Single-source ABC inner loop -- the canonical, C++-accelerated ABC engine.
# Main function: ABCpp.SingleInMultiple
#
# This retains the ORIGINAL implementation (Osang'ir, based on Amaratunga,
# Cabrera & Kovtun 2008): the per-iteration weighted feature subsample is drawn
# in pure R with sample.int(prob = zf); clustering is done with
# fastcluster::hclust; the C++ acceleration lives in the consensus step
# (count_coclustering_cpp, called by f.clustABC.MultiSource in M_ABCpp.R).
#
# Corrections carried over from the original engine:
#   B1  c solved adaptively via uniroot (.solve_weight_c), not hardcoded nf/183
#   B2  feature subsampling (Step 3) always runs; bag only controls Step 2
#   B3  NC defaults to floor(sqrt(N))
#   B4  g = floor(sqrt(G)) as an explicit integer
#   B5  as.dist() applied before hclust
#
# NEW in 0.1.1:
#   * weighting = "nugget" -- data-nugget feature weighting (see the weighting
#     block below); the classic "var"/"cv" schemes are unchanged.
#   * nfeat -- optional override of the feature-subsample size g.
# =============================================================================

#' Single-source Aggregating Bundles of Clusters (ABC), C++-accelerated engine
#'
#' The canonical single-source ABC engine of
#' \insertCite{Amaratunga2008}{MosaiClusteR}. Over \code{numsim} iterations it
#' (optionally) bootstraps the objects, draws a weighted subsample of features,
#' clusters the sub-matrix with \code{fastcluster::hclust}, and records each
#' selected object's base-cluster label. \code{\link{M_ABCpp}} turns the stacked
#' label matrix into a consensus dissimilarity using the C++ co-clustering
#' kernel.
#'
#' @section Feature weighting (incl. data nuggets):
#' Features are subsampled with probability proportional to a rank-transformed
#' importance score:
#' \itemize{
#'   \item \code{"var"} -- feature variance (the classic ABC choice);
#'   \item \code{"cv"}  -- coefficient of variation;
#'   \item \code{"nugget"} -- \emph{NEW}: the data-nugget between-group weighted
#'     variance (\code{\link{create_data_nuggets}}), a robust, Big-Data-friendly
#'     alternative to the raw variance;
#'   \item anything else -- equal weights.
#' }
#'
#' @param data Numeric matrix with \strong{objects (samples) in rows} and
#'   features in columns (the MosaiClusteR convention; transposed internally to
#'   the features-by-objects layout of the original engine).
#' @param distmeasure Distance measure passed to \code{\link{Distance}}.
#' @param weighting \code{"var"}, \code{"cv"}, \code{"nugget"} or equal.
#' @param normalize Optional normalisation (see \code{\link{Normalization}}).
#' @param gr Unused; kept for interface compatibility.
#' @param bag Logical; bootstrap the objects each iteration (Step 2).
#' @param numsim Number of iterations \eqn{R}.
#' @param numvar Kept for backward compatibility; superseded by \code{nfeat}.
#' @param linkage Agglomeration method for \code{fastcluster::hclust}.
#' @param NC Base-cluster count (\code{NULL} = \eqn{\lfloor\sqrt{N}\rfloor}).
#' @param mds Unused placeholder.
#' @param nfeat Optional feature-subsample size. \code{NULL} (default) uses the
#'   classic \eqn{g=\lfloor\sqrt{G}\rfloor}; a value \eqn{\ge 1} sets an absolute
#'   count (raise it for low-dimensional sources where most features matter).
#' @param nugget_type,nugget_args Passed to the data-nugget weighting when
#'   \code{weighting = "nugget"}.
#' @return A data frame, \code{numsim} rows by one column per object; entries are
#'   base-cluster labels (\code{0} = object not selected that iteration).
#' @references \insertAllCited{}
#' @seealso \code{\link{M_ABCpp}}, \code{\link{create_data_nuggets}}
#' @examples
#' data(mosaic_toy)
#' lab <- ABCpp.SingleInMultiple(mosaic_toy$List[[1]], numsim = 40, NC = 3,
#'                               weighting = "nugget")
#' dim(lab)
#' @export
ABCpp.SingleInMultiple <- function(data,
                                   distmeasure = "euclidean",
                                   weighting   = "var",
                                   normalize   = NULL,
                                   gr          = c(),
                                   bag         = TRUE,
                                   numsim      = 1000,
                                   numvar      = 100,
                                   linkage     = "ward.D2",
                                   NC          = NULL,
                                   mds         = FALSE,
                                   nfeat       = NULL,
                                   nugget_type = "between",
                                   nugget_args = list()) {

  # ---- input orientation -----------------------------------------------------
  # MosaiClusteR feeds objects-in-rows; the original engine works on a
  # features-by-objects matrix, so we transpose once here and keep the original
  # internal logic unchanged below.
  data <- as.matrix(data)
  samp_names <- rownames(data)
  X  <- t(data)        # X: features (rows) x objects (columns)
  ns <- ncol(X)        # N objects
  nf <- nrow(X)        # G features

  # ---- Step 4: base-cluster count  n = floor(sqrt(N))                 [B3] ----
  if (is.null(NC)) NC <- max(2L, floor(sqrt(ns)))
  NC <- max(2L, min(as.integer(NC), ns - 1L))

  # ---- Step 3: feature count  g = floor(sqrt(G))                      [B4] ----
  # (nfeat overrides g when supplied; useful for low-dimensional sources.)
  g <- if (is.null(nfeat)) max(2L, floor(sqrt(nf))) else max(2L, min(nf, as.integer(nfeat)))

  # ---- Step 1: feature weights (computed once, invariant over iterations)[B1] ----
  zf <- rep(1.0, nf)
  if (weighting == "var") {
    # classic ABC: rank features by variance; top features get more mass
    vars <- matrixStats::rowVars(X, na.rm = TRUE)
    cc   <- .solve_weight_c(nf)
    zf   <- 1.0 / (rank(-vars, ties.method = "first") + cc)
  } else if (weighting == "cv") {
    # coefficient of variation
    mn   <- matrixStats::rowMeans2(X, na.rm = TRUE)
    mn[abs(mn) < .Machine$double.eps] <- .Machine$double.eps
    cvs  <- matrixStats::rowSds(X, na.rm = TRUE) / abs(mn)
    cc   <- .solve_weight_c(nf)
    zf   <- 1.0 / (rank(-cvs, ties.method = "first") + cc)
  } else if (weighting == "nugget") {
    # NEW (0.1.1): data-nugget feature weighting. Compress the SAMPLES into
    # weighted data nuggets, then score each feature by how strongly it separates
    # the nugget centres (size-weighted between-nugget variance). Robust to
    # outliers and scales to large N. Drops straight into the same rank scheme.
    dn    <- do.call(create_data_nuggets, c(list(x = t(X)), nugget_args))  # t(X): objects x features
    score <- nugget_feature_weights(dn, type = nugget_type)
    cc    <- .solve_weight_c(nf)
    zf    <- 1.0 / (rank(-score, ties.method = "first") + cc)
  }
  zf <- zf / sum(zf)   # normalise to a probability vector over features

  ok_lnk <- c("single", "complete", "average", "mcquitty",
              "ward", "ward.D", "ward.D2", "centroid", "median")
  if (!(linkage %in% ok_lnk))
    stop("'linkage' must be one of: ", paste(ok_lnk, collapse = ", "))

  # label matrix: rows = iterations, cols = objects; 0 = object not selected
  res <- matrix(0L, nrow = numsim, ncol = ns,
                dimnames = list(paste0("iter_", seq_len(numsim)), samp_names))

  for (i in seq_len(numsim)) {

    # ---- Step 2: bootstrap sample selection                           [B2] ----
    if (bag) {
      ids <- sort(unique(sample.int(ns, size = ns, replace = TRUE)))
    } else {
      ids <- seq_len(ns)
    }
    dat_s <- X[, ids, drop = FALSE]

    # ---- Step 3: weighted feature subsample (pure R, prob = zf) -------------
    idf <- sort(sample.int(nf, size = g, replace = FALSE, prob = zf))
    dat <- dat_s[idf, , drop = FALSE]

    # ---- Step 4: cluster the g x N* sub-matrix                        [B5] ----
    DistM <- as.matrix(Distance(t(dat), distmeasure = distmeasure, normalize = normalize))
    # Robustness: the Gower-based distance is undefined (NaN) when a sampled
    # feature is constant across the selected objects (common with sparse data).
    # Fall back to a plain Euclidean distance for that iteration so the run never
    # aborts; the fallback is only reached when the primary distance fails.
    if (anyNA(DistM) || any(!is.finite(DistM))) DistM <- as.matrix(stats::dist(t(dat)))
    HClust <- fastcluster::hclust(stats::as.dist(DistM), method = linkage)
    labels <- as.integer(stats::cutree(HClust, k = NC))

    # ---- Step 5: record the labels of the selected objects ------------------
    res[i, ids] <- labels
  }

  as.data.frame(res, stringsAsFactors = FALSE)
}
