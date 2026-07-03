# =============================================================================
# distance.R  - distance/dissimilarity computation dispatcher
# =============================================================================

#' Compute a distance matrix
#'
#' Dispatcher computing an object-by-object distance matrix under a range of
#' measures suitable for continuous, binary and mixed-type data.
#'
#' @param Data A data matrix; rows are the objects.
#' @param distmeasure One of \code{"tanimoto"}, \code{"jaccard"},
#'   \code{"euclidean"} (Gower), \code{"hamming"}, \code{"cont tanimoto"},
#'   \code{"MCA_coord"}, \code{"gower"}, \code{"chi.squared"}, \code{"cosine"}.
#' @param normalize Optional normalisation method (see \code{\link{Normalization}}),
#'   applied for the \code{"euclidean"} / \code{"cont tanimoto"} measures.
#' @param method Ignored; retained for backward compatibility.
#' @return A symmetric numeric distance matrix.
#' @details \code{"cosine"}, \code{"chi.squared"} and \code{"MCA_coord"} require
#'   the suggested packages \pkg{lsa}, \pkg{analogue} and \pkg{FactoMineR}
#'   respectively.
#' @examples
#' m <- matrix(rbinom(40, 1, 0.5), nrow = 8)
#' D <- Distance(m, distmeasure = "tanimoto")
#' dim(D)
#' @export
Distance <- function(Data,
                     distmeasure = c("tanimoto", "jaccard", "euclidean",
                                     "hamming", "cont tanimoto", "MCA_coord",
                                     "gower", "chi.squared", "cosine"),
                     normalize = NULL, method = NULL) {
  distmeasure <- match.arg(distmeasure)

  if (distmeasure == "euclidean" && !is.null(normalize))
    Data <- Normalization(Data, method = normalize)

  tanimoto <- function(m) {
    m <- as.matrix(m); storage.mode(m) <- "numeric"; n <- nrow(m)
    if (!anyNA(m)) {
      N.C <- m %*% t(m); N.A <- m %*% (1 - t(m)); N.B <- (1 - m) %*% t(m)
      denom <- N.A + N.B + N.C
      S <- N.C / denom; S[denom == 0] <- 0
    } else {
      S <- matrix(NA_real_, n, n)
      for (i in seq_len(n)) for (j in i:n) {
        xi <- m[i, ]; xj <- m[j, ]; ok <- !is.na(xi) & !is.na(xj)
        if (!any(ok)) { S[i, j] <- S[j, i] <- 0; next }
        xi <- xi[ok]; xj <- xj[ok]
        c <- sum(xi == 1 & xj == 1); u <- sum(xi == 1 | xj == 1)
        S[i, j] <- S[j, i] <- if (u == 0) 0 else c / u
      }
    }
    D <- 1 - S; diag(D) <- 0; (D + t(D)) / 2
  }

  if (distmeasure == "jaccard") {
    dist <- as.matrix(ade4::dist.binary(Data, method = 1))
  } else if (distmeasure == "tanimoto") {
    dist <- as.matrix(tanimoto(Data)); rownames(dist) <- rownames(Data)
  } else if (distmeasure %in% c("euclidean", "gower")) {
    dist <- as.matrix(FD::gowdis(Data))
  } else if (distmeasure == "hamming") {
    dist <- as.matrix(e1071::hamming.distance(Data))
  } else if (distmeasure == "MCA_coord") {
    if (!requireNamespace("FactoMineR", quietly = TRUE))
      stop("distmeasure='MCA_coord' requires the 'FactoMineR' package.")
    Data <- as.data.frame(Data); Data <- Data[, sapply(Data, nlevels) > 1]
    MCA_result <- FactoMineR::MCA(Data, graph = FALSE)
    dist <- Distance(MCA_result$ind$coord, distmeasure = "euclidean", normalize = NULL)
  } else if (distmeasure == "chi.squared") {
    if (!requireNamespace("analogue", quietly = TRUE))
      stop("distmeasure='chi.squared' requires the 'analogue' package.")
    Temp <- ade4::acm.disjonctif(Data)
    dist <- analogue::distance(Temp, method = "chi.distance")
  } else if (distmeasure == "cosine") {
    if (!requireNamespace("lsa", quietly = TRUE))
      stop("distmeasure='cosine' requires the 'lsa' package.")
    Temp <- lsa::cosine(t(as.matrix(Data)))
    dist <- 1 - (1 - (2 * acos(Temp)) / pi)
  } else if (distmeasure == "cont tanimoto") {
    if (!is.null(normalize)) Data <- Normalization(Data, method = normalize)
    m <- as.matrix(Data); X_A_X_B <- m %*% t(m); temp <- diag(X_A_X_B)
    X_A <- matrix(rep(temp, nrow(Data)), nrow(Data), nrow(Data), byrow = FALSE)
    X_B <- t(X_A); Denom <- X_A + X_B - X_A_X_B
    dist <- 1 - (X_A_X_B / Denom)
  } else {
    stop("Incorrect choice of distmeasure.")
  }
  colnames(dist) <- rownames(dist)
  dist
}

#' @rdname Distance
#' @export
Distance_v2 <- function(Data, distmeasure = "tanimoto", normalize = NULL)
  Distance(Data, distmeasure = distmeasure, normalize = normalize)
