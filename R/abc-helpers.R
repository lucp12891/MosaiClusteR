# =============================================================================
# abc-helpers.R  - internal helpers shared by the ABC / M-ABC family.
# =============================================================================

# Solve for c so the top `top_frac` of features carry `target` of the combined
# selection probability under w = 1 / (rank + c). (Amaratunga et al. 2008.)
.solve_weight_c <- function(G, top_frac = 0.01, target = 0.20) {
  n_top <- max(1L, floor(G * top_frac))
  if (n_top >= G) return(1e-6)
  obj <- function(cc) {
    w <- 1.0 / (seq_len(G) + cc)
    sum(w[seq_len(n_top)]) / sum(w) - target
  }
  lo <- 1e-6; hi <- G * 1e4
  if (obj(lo) <= 0) return(lo)
  if (obj(hi) >= 0) return(hi)
  stats::uniroot(obj, lower = lo, upper = hi, tol = 1e-10)$root
}

# Feature-subsample size per iteration. nfeat = NULL -> floor(sqrt(G)) (the ABC
# default, best for high-dimensional sparse-signal data); a value in (0,1) is a
# fraction of G; a value >= 1 is an absolute count. Raising it helps when G is
# small and most features are informative (feature subsampling would otherwise
# discard signal).
.abc_gsize <- function(nf, nfeat = NULL) {
  if (is.null(nfeat))            return(max(2L, floor(sqrt(nf))))
  if (nfeat > 0 && nfeat < 1)    return(max(2L, min(nf, as.integer(round(nfeat * nf)))))
  max(2L, min(nf, as.integer(nfeat)))
}

# Per-feature selection-probability vector for the ABC family. X is
# features-by-objects (internal layout). Implements the data-nugget weighting
# option alongside the classic variance / coefficient-of-variation schemes.
.abc_zf <- function(X, weighting, nugget_type = "between", nugget_args = list()) {
  nf <- nrow(X)
  score <- switch(as.character(weighting),
    var = matrixStats::rowVars(X, na.rm = TRUE),
    cv  = {
      mn <- matrixStats::rowMeans2(X, na.rm = TRUE)
      mn[abs(mn) < .Machine$double.eps] <- .Machine$double.eps
      matrixStats::rowSds(X, na.rm = TRUE) / abs(mn)
    },
    nugget = {
      dn <- do.call(create_data_nuggets, c(list(x = t(X)), nugget_args))
      nugget_feature_weights(dn, type = nugget_type)
    },
    NULL)
  if (is.null(score)) return(rep(1.0 / nf, nf))   # equal weights
  score[!is.finite(score)] <- 0
  cc <- .solve_weight_c(nf)
  zf <- 1.0 / (rank(-score, ties.method = "first") + cc)
  zf / sum(zf)
}
