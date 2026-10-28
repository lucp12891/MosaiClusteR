# =============================================================================
# simulate-weighting.R
# Simulation designed to evaluate FEATURE-WEIGHTING schemes fairly, in
# particular the coefficient-of-variation (CV) weighting. Two regimes:
#
#   "variance" - the discriminative features carry high ABSOLUTE variance and
#                the decoys are low-variance. Ranking by variance and by CV
#                coincide, so variance / CV / nugget weighting all succeed.
#
#   "cv"       - the discriminative features have a LOW mean with moderate
#                spread (high CV) while the decoys have a HIGH mean and high
#                absolute variance but are NOT discriminative (variance grows
#                with the mean, as in count / intensity omics data). Variance
#                weighting is then dominated by the high-abundance decoys and
#                misses the signal, whereas CV weighting surfaces it.
#
# This lets the package demonstrate *when* CV weighting is advantageous, rather
# than only the variance-dominated regime in which CV appears inferior.
# =============================================================================

#' Simulate data for evaluating feature-weighting schemes (variance vs CV)
#'
#' Generates a multi-source data set with a known cluster structure in which the
#' discriminative signal is carried either by high-absolute-variance features
#' (\code{regime = "variance"}) or by high-\emph{coefficient-of-variation}
#' features buried under high-abundance, high-variance decoys
#' (\code{regime = "cv"}). It is intended for benchmarking the \code{"var"},
#' \code{"cv"} and \code{"nugget"} weightings of the ABC family (and any other
#' feature-weighting method) on the same footing.
#'
#' @param regime \code{"cv"} (default) or \code{"variance"}; see Details.
#' @param n_samples Number of objects (rows), split evenly across clusters.
#' @param n_clusters Number of clusters.
#' @param n_sources Number of data sources (independent draws of the same design).
#' @param n_signal,n_decoy,n_noise Numbers of discriminative, decoy
#'   (non-discriminative but high-variance in the CV regime) and pure-noise
#'   features per source.
#' @param effect Between-cluster mean shift on the signal features (defaults:
#'   1.4 for the variance regime, 0.9 for the cv regime).
#' @param seed Optional RNG seed.
#'
#' @return A list with \code{List} (the sources, objects in rows), \code{truth}
#'   (cluster labels), \code{signal_features} (indices of the discriminative
#'   features), \code{regime} and \code{params}.
#'
#' @seealso \code{\link{M_ABC}} (the \code{weighting} argument),
#'   \code{\link{mosaic_sim}}
#' @examples
#' cv  <- simulate_weighting_regime("cv", seed = 1)
#' # variance weighting is misled by high-abundance decoys, CV recovers the truth:
#' fit_var <- M_ABC(cv$List, numsim = 80, NC = 3, weighting = c("var", "var"))
#' fit_cv  <- M_ABC(cv$List, numsim = 80, NC = 3, weighting = c("cv",  "cv"))
#' c(var = cluster_agreement(mosaic_labels(fit_var, 3), cv$truth)[["ARI"]],
#'   cv  = cluster_agreement(mosaic_labels(fit_cv,  3), cv$truth)[["ARI"]])
#' @export
simulate_weighting_regime <- function(regime      = c("cv", "variance"),
                                      n_samples   = 120,
                                      n_clusters  = 3,
                                      n_sources   = 2,
                                      n_signal    = 40,
                                      n_decoy     = 500,
                                      n_noise     = 60,
                                      effect      = NULL,
                                      seed        = NULL) {
  regime <- match.arg(regime)
  if (!is.null(seed)) set.seed(seed)
  z  <- rep(seq_len(n_clusters), length.out = n_samples)   # balanced clusters
  dir <- as.integer(factor(z)) - (n_clusters + 1) / 2      # per-cluster direction
  if (is.null(effect)) effect <- if (regime == "variance") 1.4 else 0.9

  one_source <- function() {
    if (regime == "variance") {
      # signal: cluster mean-shift on unit-scale features -> high absolute variance
      Xs <- vapply(seq_len(n_signal), function(j)
        stats::rnorm(n_samples, 5 + effect * dir, 1.2), numeric(n_samples))
      Xd <- matrix(stats::rnorm(n_samples * n_decoy, 5, 0.3), n_samples)  # low-variance decoys
      Xn <- matrix(stats::rnorm(n_samples * n_noise, 5, 1.0), n_samples)
    } else {
      # signal: LOW mean, moderate spread, cluster mean-shift -> high CV, modest abs. variance
      Xs <- vapply(seq_len(n_signal), function(j)
        abs(stats::rnorm(n_samples, 3 + effect * dir, 0.7)), numeric(n_samples))
      # decoys: HIGH mean, variance proportional to mean (count/intensity-like),
      # NOT cluster-discriminative -> high absolute variance, LOW CV
      Xd <- vapply(seq_len(n_decoy), function(j) {
        mu <- stats::runif(1, 25, 60); stats::rnorm(n_samples, mu, sqrt(mu * 2)) },
        numeric(n_samples))
      Xn <- matrix(stats::rnorm(n_samples * n_noise, 8, 1.0), n_samples)
    }
    X <- cbind(Xs, Xd, Xn)
    rownames(X) <- sprintf("obj%03d", seq_len(n_samples))
    colnames(X) <- sprintf("f%04d", seq_len(ncol(X)))
    X
  }

  List <- lapply(seq_len(n_sources), function(i) one_source())
  names(List) <- sprintf("source%d", seq_len(n_sources))
  list(List = List,
       truth = stats::setNames(z, sprintf("obj%03d", seq_len(n_samples))),
       signal_features = seq_len(n_signal),
       regime = regime,
       params = list(n_samples = n_samples, n_clusters = n_clusters,
                     n_sources = n_sources, n_signal = n_signal,
                     n_decoy = n_decoy, n_noise = n_noise, effect = effect))
}
