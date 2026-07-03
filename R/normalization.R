# =============================================================================
# normalization.R  - feature normalisation utilities
# =============================================================================

#' Normalisation of features
#'
#' When data of different scales are combined it is recommended to normalise so
#' that the structures are comparable.
#'
#' @param Data A data matrix; rows are assumed to correspond to the objects.
#' @param method One of \code{"Quantile"}, \code{"Fisher-Yates"},
#'   \code{"Standardize"}, \code{"Range"} (or any unambiguous first letter).
#' @details \code{"Quantile"} performs quantile normalisation (common in omics);
#'   \code{"Fisher-Yates"} uses normal scores based only on the number of rows;
#'   \code{"Standardize"} centres and scales each column (base
#'   \code{\link[base]{scale}}); \code{"Range"} maps the matrix to \code{[0, 1]}.
#' @return The normalised matrix.
#' @examples
#' x <- matrix(rnorm(100), 10, 10)
#' Norm_x <- Normalization(x, method = "R")
#' @export
Normalization <- function(Data,
                          method = c("Quantile", "Fisher-Yates", "Standardize",
                                     "Range", "Q", "q", "F", "f", "S", "s",
                                     "R", "r")) {
  method <- substring(method[1], 1, 1)   # accept the first letter of each method
  method <- match.arg(method)

  if (method %in% c("S", "s")) {
    return(scale(Data, center = TRUE, scale = TRUE))
  }
  if (method %in% c("R", "r")) {
    minc <- min(as.vector(Data)); maxc <- max(as.vector(Data))
    return((Data - minc) / (maxc - minc))
  }
  if (method %in% c("Q", "q")) {
    DataColSorted <- apply(Data, 2, sort)
    NValues <- apply(DataColSorted, 1, stats::median, na.rm = TRUE)
  } else {  # Fisher-Yates
    NValues <- stats::qnorm(seq_len(nrow(Data)) / (nrow(Data) + 1))
  }
  DataRanked <- c(apply(Data, 2, rank))
  array(stats::approx(seq_len(nrow(Data)), NValues, DataRanked)$y,
        dim(Data), dimnames(Data))
}
