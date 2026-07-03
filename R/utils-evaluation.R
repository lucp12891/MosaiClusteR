# =============================================================================
# utils-evaluation.R
# Partition extraction, agreement metrics and method comparison.
# =============================================================================

#' Extract a flat partition from a MosaiClusteR result
#'
#' Cuts the hierarchical clustering carried in a MosaiClusteR result object into
#' \code{k} groups. Works with any method returning a \code{Clust} element (an
#' \code{agnes}/\code{hclust} object) and/or a \code{DistM}.
#'
#' @param object A MosaiClusteR result (e.g. from \code{\link{M_ABCpp}}), an
#'   \code{agnes}/\code{hclust} object, or a distance matrix.
#' @param k Number of clusters.
#' @param linkage Linkage used when \code{object} is a bare distance matrix.
#' @return An integer vector of cluster labels, named by object.
#' @examples
#' set.seed(1)
#' s <- paste0("o", 1:20)
#' L <- list(matrix(rnorm(20 * 40), 20, dimnames = list(s, NULL)),
#'           matrix(rnorm(20 * 50), 20, dimnames = list(s, NULL)))
#' fit <- M_ABCpp(L, numsim = 30)
#' mosaic_labels(fit, k = 3)
#' @export
mosaic_labels <- function(object, k, linkage = "ward") {
  cl <- NULL
  if (is.list(object) && !is.null(object$Clust)) {
    cl <- object$Clust
    if (is.list(cl) && !is.null(cl$Clust)) cl <- cl$Clust  # WeightedClust nesting
  } else if (inherits(object, c("agnes", "hclust", "twins"))) {
    cl <- object
  } else if (is.matrix(object)) {
    cl <- cluster::agnes(as.matrix(object), diss = TRUE, method = linkage)
  } else {
    stop("Cannot extract a clustering from 'object'.")
  }
  stats::cutree(stats::as.hclust(cl), k = k)
}

#' Agreement between two partitions
#'
#' Computes the adjusted Rand index (ARI), normalised mutual information (NMI)
#' and pairwise Jaccard between two label vectors of equal length.
#'
#' @param a,b Integer/character/factor label vectors of equal length.
#' @return A named numeric vector \code{c(ARI, NMI, Jaccard)}.
#' @details ARI rule of thumb: \code{>0.90} excellent; \code{0.75-0.90} good;
#'   \code{0.50-0.75} moderate; \code{<0.50} weak.
#' @examples
#' cluster_agreement(c(1, 1, 2, 2), c(2, 2, 1, 1))   # identical up to relabel
#' @export
cluster_agreement <- function(a, b) {
  a <- as.integer(as.factor(a)); b <- as.integer(as.factor(b))
  if (length(a) != length(b)) stop("'a' and 'b' must have equal length.")
  n  <- length(a)
  ct <- table(a, b)
  rs <- rowSums(ct); cs <- colSums(ct)

  sum_a <- sum(choose(rs, 2)); sum_b <- sum(choose(cs, 2))
  sum_ab <- sum(choose(ct, 2)); sum_n <- choose(n, 2)
  expected <- sum_a * sum_b / sum_n
  max_val  <- 0.5 * (sum_a + sum_b)
  ari <- if (max_val == expected) 0 else (sum_ab - expected) / (max_val - expected)

  p <- ct / n; p1 <- rowSums(p); p2 <- colSums(p)
  H1 <- -sum(p1[p1 > 0] * log(p1[p1 > 0]))
  H2 <- -sum(p2[p2 > 0] * log(p2[p2 > 0]))
  pm <- outer(p1, p2); mask <- p > 0 & pm > 0
  MI <- sum(p[mask] * log(p[mask] / pm[mask]))
  nmi <- if (H1 + H2 == 0) 0 else 2 * MI / (H1 + H2)

  same1 <- outer(a, a, "=="); same2 <- outer(b, b, "==")
  lt <- lower.tri(same1)
  A  <- sum(same1[lt] & same2[lt]); B <- sum(same1[lt] & !same2[lt])
  C  <- sum(!same1[lt] & same2[lt])
  jac <- if (A + B + C == 0) 1 else A / (A + B + C)

  c(ARI = ari, NMI = nmi, Jaccard = jac)
}

#' Compare two clusterings derived from distance matrices
#'
#' Clusters two dissimilarity matrices with the same linkage and number of
#' clusters and reports their agreement (ARI, NMI, pairwise Jaccard) and the
#' contingency table.
#'
#' @param D1,D2 Square dissimilarity matrices over the same objects.
#' @param NC Number of clusters (default \code{floor(sqrt(n))}).
#' @param linkage Linkage method. Default \code{"ward.D2"}.
#' @param labels Length-2 labels for the two solutions.
#' @param verbose Logical; print a short report.
#' @return Invisibly, a list with the two label vectors, \code{ARI},
#'   \code{NMI}, \code{pair_Jaccard} and the \code{contingency} table.
#' @examples
#' set.seed(1)
#' D1 <- as.matrix(dist(matrix(rnorm(60), 20)))
#' D2 <- as.matrix(dist(matrix(rnorm(60), 20)))
#' compare_clusterings(D1, D2, NC = 3, verbose = FALSE)$ARI
#' @export
compare_clusterings <- function(D1, D2, NC = NULL, linkage = "ward.D2",
                                labels = c("D1", "D2"), verbose = TRUE) {
  stopifnot(identical(dim(D1), dim(D2)))
  ns <- nrow(D1)
  if (is.null(NC)) NC <- max(2L, floor(sqrt(ns)))
  h1 <- fastcluster::hclust(stats::as.dist(D1), method = linkage)
  h2 <- fastcluster::hclust(stats::as.dist(D2), method = linkage)
  c1 <- as.integer(stats::cutree(h1, k = NC))
  c2 <- as.integer(stats::cutree(h2, k = NC))
  names(c1) <- names(c2) <- rownames(D1)

  ag <- cluster_agreement(c1, c2)
  ct <- table(c1, c2)
  dimnames(ct) <- list(sprintf("%s_c%d", labels[1], seq_len(nrow(ct))),
                       sprintf("%s_c%d", labels[2], seq_len(ncol(ct))))
  if (verbose) {
    cat(sprintf("NC = %d, linkage = %s, n = %d\n", NC, linkage, ns))
    cat(sprintf("Adjusted Rand Index : %.4f\n", ag["ARI"]))
    cat(sprintf("Normalized Mut. Info: %.4f\n", ag["NMI"]))
    cat(sprintf("Pair-Jaccard (same) : %.4f\n", ag["Jaccard"]))
    print(ct)
  }
  invisible(list(clusters_1 = c1, clusters_2 = c2, ARI = ag["ARI"],
                 NMI = ag["NMI"], pair_Jaccard = ag["Jaccard"], contingency = ct))
}
