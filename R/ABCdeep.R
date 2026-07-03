# =============================================================================
# ABCdeep.R
# ABC with a deep-learning (autoencoder) Step 4 -- the ABCdeep engine.
# Main function: ABCdeep.SingleInMultiple
#
# This follows the ABC algorithm of Amaratunga, Cabrera & Kovtun (2008) with a
# deep Step 4:
#   Step 1  feature weights zf (var / cv / nugget)              -- as in ABCpp
#   Step 2  bootstrap the objects                               -- as in ABCpp
#   Step 3  DEEP feature focus: rank features by zf and keep the top fraction
#           (topk_frac); the autoencoder learns on those informative features
#           (this replaces ABCpp's per-iteration g=floor(sqrt(G)) subsample,
#           which does not suit a trained encoder -- see below)
#   Step 4  DEEP: a neural autoencoder embeds the objects into a latent space,
#           then the latent object-vectors are clustered into NC groups
#           (replaces ABCpp's distance -> hclust on the raw sub-matrix)
#   Step 5  record the base-cluster labels                      -- as in ABCpp
# The stacked label matrix is fused into a consensus dissimilarity by
# f.clustABC.MultiSource (the same C++-accelerated kernel used by M_ABCpp).
#
# DIGEST OF THE ORIGINAL ABCdeep + WHAT CHANGED HERE
# --------------------------------------------------
# The original ABCdeep (keras/tensorflow) trained a fresh autoencoder on EVERY
# iteration (slow) and silently fell back to PCA when TF was unavailable. This
# version removes both, and -- crucially -- gets clustering quality on par with
# (or better than) M_ABCpp:
#   * WHY NOT PER-ITERATION FEATURE SUBSAMPLING? A subsample of only
#     g = floor(sqrt(G)) features pushed through an encoder trained on all G
#     features is a train/inference mismatch and gives a degraded embedding on
#     dense-signal data. Instead we focus the encoder ONCE on the top-`topk_frac`
#     most-informative features (by the same var/cv/nugget weighting): this drops
#     the noise features that dilute a full-feature embedding on sparse-signal
#     data, yet keeps enough features to embed dense-signal data well. Empirically
#     this matches M_ABCpp on sparse signals and beats it on dense signals.
#   * SUPER FAST -- the autoencoder is trained ONCE per source (on the top-K
#     features); each of the numsim iterations does only a forward pass over the
#     bootstrap. The network is a self-contained, vectorised R autoencoder
#     (He-init, ReLU, Adam, MSE) -- no Python, no TensorFlow, no torch.
#   * NO FALLBACK -- exactly one code path. No PCA, no random labels, no
#     availability check; a degenerate input errors honestly.
# =============================================================================

## ---- self-contained autoencoder (vectorised, He-init, ReLU, Adam, MSE) ------

.deep_sigmoid <- function(z) 1 / (1 + exp(-z))
.deep_relu    <- function(z) { z[z < 0] <- 0; z }

# min-max scale each column to [0, 1] (autoencoder targets a sigmoid output).
.deep_minmax <- function(M) {
  cr   <- matrixStats::colRanges(M, na.rm = TRUE)
  span <- cr[, 2] - cr[, 1]; span[span == 0 | !is.finite(span)] <- 1
  X <- sweep(sweep(M, 2, cr[, 1], "-"), 2, span, "/")
  X[!is.finite(X)] <- 0
  X
}

# one Adam update over a named list of parameters (matrices and/or vectors).
.deep_adam <- function(par, grad, state, lr, t,
                       b1 = 0.9, b2 = 0.999, eps = 1e-8) {
  for (nm in names(par)) {
    state$m[[nm]] <- b1 * state$m[[nm]] + (1 - b1) * grad[[nm]]
    state$v[[nm]] <- b2 * state$v[[nm]] + (1 - b2) * grad[[nm]]^2
    mhat <- state$m[[nm]] / (1 - b1^t)
    vhat <- state$v[[nm]] / (1 - b2^t)
    par[[nm]] <- par[[nm]] - lr * mhat / (sqrt(vhat) + eps)
  }
  list(par = par, state = state)
}

# Train a 4-layer autoencoder  X (n x m, in [0,1])
#   encoder:  m --W1--> hidden --ReLU--> --W2--> latent (linear)
#   decoder:  latent --W3--> hidden --ReLU--> --W4--> m (sigmoid)
# Full-batch Adam on the MSE reconstruction loss. Returns the encoder weights.
.deep_ae_fit <- function(X, hidden, latent, epochs, lr) {
  n <- nrow(X); m <- ncol(X)
  h <- as.integer(hidden); d <- as.integer(latent)
  he  <- function(nr, nc, fan) matrix(stats::rnorm(nr * nc, 0, sqrt(2 / fan)), nr, nc)
  par <- list(W1 = he(m, h, m), b1 = numeric(h),
              W2 = he(h, d, h), b2 = numeric(d),
              W3 = he(d, h, d), b3 = numeric(h),
              W4 = he(h, m, h), b4 = numeric(m))
  zero  <- lapply(par, function(p) p * 0)
  state <- list(m = zero, v = zero)
  denom <- n * m

  for (t in seq_len(epochs)) {
    # ---- forward ----
    z1  <- sweep(X   %*% par$W1, 2, par$b1, "+"); a1  <- .deep_relu(z1)
    enc <- sweep(a1  %*% par$W2, 2, par$b2, "+")                 # linear latent
    z3  <- sweep(enc %*% par$W3, 2, par$b3, "+"); a3  <- .deep_relu(z3)
    z4  <- sweep(a3  %*% par$W4, 2, par$b4, "+"); out <- .deep_sigmoid(z4)
    # ---- backward (MSE, mean over all n*m elements) ----
    dz4 <- (2 * (out - X) / denom) * out * (1 - out)
    gW4 <- crossprod(a3, dz4);              gb4 <- colSums(dz4)
    dz3 <- (dz4 %*% t(par$W4)) * (z3 > 0)
    gW3 <- crossprod(enc, dz3);             gb3 <- colSums(dz3)
    denc <- dz3 %*% t(par$W3)                                    # latent is linear
    gW2 <- crossprod(a1, denc);             gb2 <- colSums(denc)
    dz1 <- (denc %*% t(par$W2)) * (z1 > 0)
    gW1 <- crossprod(X, dz1);               gb1 <- colSums(dz1)
    grad <- list(W1 = gW1, b1 = gb1, W2 = gW2, b2 = gb2,
                 W3 = gW3, b3 = gb3, W4 = gW4, b4 = gb4)
    upd   <- .deep_adam(par, grad, state, lr, t)
    par   <- upd$par; state <- upd$state
  }
  list(W1 = par$W1, b1 = par$b1, W2 = par$W2, b2 = par$b2)       # encoder only
}

# Encode a batch of objects through the (already trained) encoder. Fast: two
# matrix products. `A` is objects x (top-K features), matching the AE input.
.deep_encode <- function(A, P) {
  z1 <- sweep(A %*% P$W1, 2, P$b1, "+")
  sweep(.deep_relu(z1) %*% P$W2, 2, P$b2, "+")
}

# Cluster the latent object-vectors into NC groups (Step 4b).
.deep_cluster <- function(Z, NC, method, linkage) {
  if (method == "kmeans")
    return(as.integer(stats::kmeans(Z, centers = NC, nstart = 5L, iter.max = 50L)$cluster))
  hc <- fastcluster::hclust(stats::dist(Z), method = linkage)
  as.integer(stats::cutree(hc, k = NC))
}

#' Single-source ABC with a deep-learning (autoencoder) Step 4 (ABCdeep)
#'
#' A deep-learning variant of \code{\link{ABCpp.SingleInMultiple}} that keeps the
#' ABC algorithm of \insertCite{Amaratunga2008}{MosaiClusteR} intact and replaces
#' \strong{only Step 4}: instead of clustering the raw feature sub-matrix, a
#' neural autoencoder embeds the objects into a latent space and the latent
#' vectors are clustered into \code{NC} base clusters. The stacked label matrix
#' is fused across sources by \code{\link{f.clustABC.MultiSource}}, exactly as in
#' \code{\link{M_ABCpp}}.
#'
#' @section Speed (trained once, forward pass per iteration):
#' The autoencoder is trained a single time on the source's top-\code{topk_frac}
#' most-informative features; each of the \code{numsim} iterations then performs
#' only a forward pass over the bootstrapped objects. The network is a
#' self-contained, vectorised R autoencoder (He initialisation, ReLU, Adam, MSE)
#' -- there is no Python, TensorFlow or torch dependency, and \strong{no PCA or
#' other fallback}: a single code path.
#'
#' @param data Numeric matrix, \strong{objects (samples) in rows}, features in
#'   columns (the MosaiClusteR convention).
#' @param weighting Feature weighting used to rank features for the top-K focus:
#'   \code{"var"}, \code{"cv"}, \code{"nugget"} or equal (as in
#'   \code{\link{ABCpp.SingleInMultiple}}).
#' @param normalize Optional normalisation applied before the internal min-max
#'   scaling (see \code{\link{Normalization}}); \code{NULL}/\code{FALSE} = none.
#' @param gr Unused; kept for interface compatibility.
#' @param bag Logical; bootstrap the objects each iteration (Step 2).
#' @param numsim Number of iterations \eqn{R}.
#' @param numvar Kept for backward compatibility; unused.
#' @param NC Base-cluster count (\code{NULL} = \eqn{\lfloor\sqrt{N}\rfloor}).
#' @param topk_frac Fraction of the highest-weighted features the autoencoder is
#'   trained/embedded on (default 0.2). Focusing on the informative features
#'   drops the noise that would dilute a full-feature embedding on sparse signals,
#'   while keeping enough features for dense signals; the count is floored at
#'   \code{min(G, 20)} and capped at \code{min(G, 4000)}.
#' @param latent_dim Autoencoder bottleneck size: an integer, or \code{"auto"}
#'   (default) = \code{max(NC + 2, 8)}, capped at the hidden width and \eqn{K-1}.
#' @param hidden_units Hidden-layer width of the encoder/decoder.
#' @param ae_epochs Training epochs for the one-time autoencoder fit.
#' @param ae_lr Adam learning rate.
#' @param cluster_in_latent Latent clustering (Step 4b): \code{"ward"} (default,
#'   hierarchical) or \code{"kmeans"}.
#' @param linkage Agglomeration method when \code{cluster_in_latent = "ward"}.
#' @param nugget_type,nugget_args Passed to the nugget weighting when
#'   \code{weighting = "nugget"}.
#' @param min_samples Emit a warning when \code{nrow(data)} is below this (the
#'   autoencoder is underdetermined at very small \eqn{N}); default 40.
#' @param seed Optional RNG seed (controls the autoencoder init and the
#'   bootstrap).
#' @return A data frame, \code{numsim} rows by one column per object; entries are
#'   base-cluster labels (\code{0} = object not selected that iteration) -- the
#'   same contract as \code{\link{ABCpp.SingleInMultiple}}.
#' @references \insertAllCited{}
#' @seealso \code{\link{M_ABCdeep}}, \code{\link{ABCpp.SingleInMultiple}},
#'   \code{\link{f.clustABC.MultiSource}}
#' @examples
#' data(mosaic_toy)
#' lab <- ABCdeep.SingleInMultiple(mosaic_toy$List[[1]], numsim = 40, NC = 3,
#'                                 ae_epochs = 40, seed = 1)
#' dim(lab)
#' @export
ABCdeep.SingleInMultiple <- function(data,
                                     weighting         = "var",
                                     normalize         = NULL,
                                     gr                = c(),
                                     bag               = TRUE,
                                     numsim            = 1000,
                                     numvar            = 100,
                                     NC                = NULL,
                                     topk_frac         = 0.2,
                                     latent_dim        = "auto",
                                     hidden_units      = 32L,
                                     ae_epochs         = 80L,
                                     ae_lr             = 1e-3,
                                     cluster_in_latent = "ward",
                                     linkage           = "ward.D2",
                                     nugget_type       = "between",
                                     nugget_args       = list(),
                                     min_samples       = 40L,
                                     seed              = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (!cluster_in_latent %in% c("ward", "kmeans"))
    stop("'cluster_in_latent' must be 'ward' or 'kmeans'.")

  data <- as.matrix(data)
  samp_names <- rownames(data)
  n <- nrow(data)   # N objects
  m <- ncol(data)   # G features

  # Deep learning needs samples: an autoencoder trained on very few objects is
  # underdetermined and the latent embedding is unreliable. Warn (do not fail)
  # so tiny-n callers know to prefer ABCpp; the demo skips ABCdeep below n.
  if (n < min_samples)
    warning("ABCdeep: n = ", n, " objects is small for an autoencoder (< ",
            min_samples, "); the deep embedding may be unreliable -- consider ",
            "ABCpp.SingleInMultiple / M_ABCpp for datasets this small.",
            call. = FALSE)

  # optional normalisation, then min-max scaling for the sigmoid autoencoder
  if (!is.null(normalize) && !isFALSE(normalize))
    data <- as.matrix(Normalization(data, method = normalize))
  Xs <- .deep_minmax(data)

  if (is.null(NC)) NC <- max(2L, floor(sqrt(n)))
  NC <- max(2L, min(as.integer(NC), n - 1L))

  # ---- Step 1: feature weights (var / cv / nugget), reusing the ABC helper ----
  # Rank features on the RAW (normalised) data, exactly as ABCpp does -- NOT on
  # the min-max-scaled matrix, which flattens per-feature variance and would make
  # the "var" ranking pick the wrong (noise) features. .abc_zf wants
  # features-by-objects, so transpose the objects-in-rows input.
  zf <- .abc_zf(t(data), weighting, nugget_type, nugget_args)

  # ---- Step 3 (deep focus): keep the top-K most-informative features ---------
  # A per-iteration sqrt(G) subsample mismatches an encoder trained on all G
  # features; instead we focus the encoder on the top-K informative features so
  # it learns signal, not noise. K is floored at 20 and capped at 4000.
  K    <- min(m, 4000L, max(as.integer(round(topk_frac * m)), 20L))
  keep <- order(zf, decreasing = TRUE)[seq_len(K)]
  Xk   <- Xs[, keep, drop = FALSE]

  # ---- autoencoder bottleneck size ------------------------------------------
  if (identical(latent_dim, "auto")) latent_dim <- max(NC + 2L, 8L)
  latent_dim <- max(2L, min(as.integer(latent_dim), K - 1L, as.integer(hidden_units)))

  # ---- train the autoencoder ONCE on the top-K features ---------------------
  message("ABCdeep: training autoencoder (G=", m, ", top-K=", K, ", N=", n,
          ", hidden=", hidden_units, ", latent=", latent_dim,
          ", epochs=", ae_epochs, ") ...")
  P <- .deep_ae_fit(Xk, hidden_units, latent_dim, as.integer(ae_epochs), ae_lr)

  # ---- iterate: bootstrap objects -> encode -> cluster ----------------------
  res <- matrix(0L, nrow = numsim, ncol = n,
                dimnames = list(paste0("iter_", seq_len(numsim)), samp_names))

  for (i in seq_len(numsim)) {
    ids <- if (bag) sort(unique(sample.int(n, size = n, replace = TRUE)))
           else      seq_len(n)
    if (length(ids) < max(4L, NC)) next            # not enough objects to cut NC

    Z <- .deep_encode(Xk[ids, , drop = FALSE], P)                    # Step 4a
    res[i, ids] <- .deep_cluster(Z, NC, cluster_in_latent, linkage)  # Step 4b
  }

  as.data.frame(res, stringsAsFactors = FALSE)
}
