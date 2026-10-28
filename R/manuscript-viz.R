# =============================================================================
# manuscript-viz.R
# Publication helpers: justify the number of clusters (silhouette / elbow),
# visualise clusters in 2-D (PCA / UMAP / t-SNE), and compare clustering
# solutions to each other when no ground truth exists.
# =============================================================================

# distinct, colour-blind-friendly cluster palette (ColorPalette is a sequential
# ramp; for categorical cluster colours we want maximally distinct hues).
.mc_pal <- function(n) {
  base <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F",
            "#EDC948", "#B07AA1", "#9C755F", "#BAB0AC", "#D37295")
  if (n <= length(base)) base[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

# Shared publication theme: white background, sharp bold axes, grey facet strips
# and (by default) no grid lines. `grid` = "none"/"x"/"y"/"both" adds faint
# guide lines only where they help read values.
.mc_theme <- function(base_size = 13, grid = c("none", "x", "y", "both")) {
  grid <- match.arg(grid)
  ink <- "#111827"
  th <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background    = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border       = ggplot2::element_blank(),
      axis.line          = ggplot2::element_line(colour = ink, linewidth = 0.7),
      axis.ticks         = ggplot2::element_line(colour = ink, linewidth = 0.7),
      axis.ticks.length  = grid::unit(3, "pt"),
      axis.title         = ggplot2::element_text(face = "bold", colour = ink),
      axis.text          = ggplot2::element_text(face = "bold", colour = ink),
      plot.title         = ggplot2::element_text(face = "bold", colour = ink, hjust = 0.5),
      strip.background   = ggplot2::element_rect(fill = "#E5E7EB", colour = NA),
      strip.text         = ggplot2::element_text(face = "bold", colour = ink, size = base_size),
      legend.title       = ggplot2::element_text(face = "bold", colour = ink),
      legend.key         = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major   = ggplot2::element_blank())
  if (grid %in% c("x", "both"))
    th <- th + ggplot2::theme(panel.grid.major.x =
      ggplot2::element_line(colour = "grey88", linewidth = 0.4, linetype = "22"))
  if (grid %in% c("y", "both"))
    th <- th + ggplot2::theme(panel.grid.major.y =
      ggplot2::element_line(colour = "grey88", linewidth = 0.4, linetype = "22"))
  th
}

# coerce a fit / dist / list to an object-by-object distance matrix
.as_distm <- function(x, distmeasure = "euclidean") {
  if (is.list(x) && !is.null(x$DistM)) return(as.matrix(x$DistM))
  if (is.list(x) && is.null(x$DistM)) {                       # a List of sources
    Z <- do.call(cbind, lapply(x, function(m) { z <- scale(as.matrix(m)); z[is.na(z)] <- 0; z }))
    return(as.matrix(stats::dist(Z)))
  }
  as.matrix(x)
}

#' Justify the number of clusters (silhouette + elbow / scree)
#'
#' For a fused dissimilarity (or a fitted result, or a \code{List} of sources)
#' compute, over a range of \eqn{k}, the average silhouette width and the total
#' within-cluster dispersion (an elbow / scree curve), and plot them so the
#' chosen \eqn{k} can be justified in the manuscript.
#'
#' @param x A fitted result (with \code{$DistM}), a dissimilarity matrix, or a
#'   \code{List} of data sources (objects in rows).
#' @param k_range Integer vector of candidate cluster counts.
#' @param linkage Agglomeration method used to cut the dissimilarity.
#' @param select Which criterion to show \emph{and} use for \code{best_k}:
#'   \code{"both"} (default) draws both panels and picks \eqn{k} from the elbow;
#'   \code{"elbow"} draws only the within-cluster dispersion (scree) curve and
#'   uses its knee; \code{"silhouette"} draws only the silhouette curve and uses
#'   its maximum. Passing a single criterion therefore shows a single panel.
#' @param plot Logical; draw the figure.
#' @param point_size Marker size for the curve points.
#' @param base_size Base font size for the \pkg{ggplot2} theme.
#' @return Invisibly a data frame with \code{k}, \code{silhouette} and
#'   \code{within_dispersion}, plus attributes \code{"best_k"}, \code{"elbow_k"},
#'   \code{"silhouette_k"} and \code{"plot"} (the \pkg{ggplot} object).
#' @seealso \code{\link{embedding_plot}}, \code{\link{SelectnrClusters}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#' cluster_selection_plot(fit, k_range = 2:6, select = "elbow")   # one panel
#' @export
cluster_selection_plot <- function(x, k_range = 2:8, linkage = "ward.D2",
                                   select = c("both", "elbow", "silhouette"),
                                   plot = TRUE, point_size = 2.6, base_size = 13) {
  select <- match.arg(select)
  D  <- .as_distm(x); d <- stats::as.dist(D)
  hc <- fastcluster::hclust(d, method = linkage)
  res <- do.call(rbind, lapply(k_range, function(k) {
    cl  <- stats::cutree(hc, k)
    sil <- mean(cluster::silhouette(cl, d)[, "sil_width"])
    wss <- sum(vapply(split(seq_along(cl), cl), function(idx)
      if (length(idx) < 2) 0 else sum(as.matrix(D)[idx, idx][lower.tri(diag(length(idx)))]),
      numeric(1)))
    data.frame(k = k, silhouette = sil, within_dispersion = wss)
  }))
  # ---- elbow (kneedle): the k where the scree curve starts to flatten. Scale
  # k and dispersion to [0,1] and take the point furthest below the straight
  # line joining the first and last point (largest drop in marginal gain).
  kk <- res$k; w <- res$within_dispersion
  kn <- (kk - min(kk)) / (max(kk) - min(kk))
  wn <- (w  - min(w))  / (max(w)  - min(w) + 1e-12)
  chord   <- 1 - kn
  elbow_k <- res$k[which.max(chord - wn)]
  sil_k   <- res$k[which.max(res$silhouette)]
  best    <- if (select == "silhouette") sil_k else elbow_k
  gp <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("ggplot2", quietly = TRUE))
      stop("cluster_selection_plot() needs the 'ggplot2' package.")
    EL <- "Within-cluster dispersion (elbow)"; SI <- "Average silhouette width"
    long <- rbind(data.frame(k = res$k, value = res$within_dispersion, metric = EL),
                  data.frame(k = res$k, value = res$silhouette,        metric = SI))
    show <- switch(select, elbow = EL, silhouette = SI, both = c(EL, SI))
    long <- long[long$metric %in% show, ]
    long$metric <- factor(long$metric, levels = c(EL, SI))
    vl <- data.frame(metric = factor(c(EL, SI), levels = c(EL, SI)),
                     xintercept = c(elbow_k, sil_k))
    vl <- vl[as.character(vl$metric) %in% show, ]
    gp <- ggplot2::ggplot(long, ggplot2::aes(x = k, y = value)) +
      ggplot2::geom_vline(data = vl, ggplot2::aes(xintercept = xintercept),
                          linetype = "dashed", colour = "#E15759", linewidth = 0.8) +
      ggplot2::geom_line(colour = "#4E79A7", linewidth = 0.9) +
      ggplot2::geom_point(colour = "#4E79A7", size = point_size) +
      ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = if (length(show) == 1) 1 else 2) +
      ggplot2::scale_x_continuous(breaks = res$k) +
      ggplot2::labs(x = "Number of clusters (k)", y = NULL) +
      .mc_theme(base_size, grid = "none")
    print(gp)
  }
  attr(res, "best_k")       <- best
  attr(res, "elbow_k")      <- elbow_k
  attr(res, "silhouette_k") <- sil_k
  attr(res, "plot")         <- gp
  invisible(res)
}

#' Show the clustering at every k (see how it deteriorates)
#'
#' Rather than only reporting the best \eqn{k}, project the objects to a
#' \emph{fixed} 2-D embedding and recolour them by the \eqn{k}-cut for each
#' candidate \eqn{k} in a grid, annotating the average silhouette per panel. This
#' makes it visible how the structure splits and degrades as \eqn{k} grows.
#'
#' The panels are drawn with \pkg{ggplot2} for a clean, high-resolution result
#' whose quality is preserved when the figure is rescaled -- save the returned
#' object to a vector device (PDF/SVG) or a high-DPI raster for publication, e.g.
#' \code{ggplot2::ggsave("ksweep.pdf", cluster_k_sweep(fit), width = 9, height = 6)}.
#'
#' @param x A fitted result (\code{$DistM}), a dissimilarity, or a \code{List}.
#' @param k_range Candidate cluster counts to display.
#' @param linkage Agglomeration method used to cut the fused dissimilarity.
#' @param embed Embedding for the (fixed) 2-D coordinates: \code{"pca"} (default),
#'   \code{"umap"} or \code{"tsne"}.
#' @param ncol Columns in the panel grid.
#' @param point_size Diameter of the plotted points (default \code{2.2}); raise it
#'   for larger, clearer markers.
#' @param alpha Point opacity (default \code{0.9}).
#' @param palette Discrete \code{\link[grDevices]{hcl.colors}} palette used to
#'   colour the clusters (default \code{"Dark 3"}).
#' @param base_size Base font size passed to the \pkg{ggplot2} theme.
#' @return Invisibly the \pkg{ggplot} object (also printed, so it draws on an open
#'   device); the average silhouette per \eqn{k} is attached as the attribute
#'   \code{"silhouette"}.
#' @seealso \code{\link{cluster_selection_plot}}, \code{\link{embedding_plot}}
#' @examples
#' \donttest{
#' data(mosaic_toy)
#' fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#' cluster_k_sweep(fit, k_range = 2:6, point_size = 3)
#' }
#' @export
cluster_k_sweep <- function(x, k_range = 2:6, linkage = "ward.D2",
                            embed = c("pca", "umap", "tsne"), ncol = 3,
                            point_size = 2.2, alpha = 0.9, palette = "Dark 3",
                            base_size = 13) {
  embed <- match.arg(embed)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("cluster_k_sweep() needs the 'ggplot2' package.")
  D <- .as_distm(x); d <- stats::as.dist(D); n <- nrow(D)
  hc <- fastcluster::hclust(d, method = linkage)
  # fixed 2-D embedding + axis labels (with % variance for the PCoA case)
  xlab <- paste0(toupper(embed), " 1"); ylab <- paste0(toupper(embed), " 2")
  Y <- switch(embed,
    pca  = { mds <- stats::cmdscale(d, k = 2, eig = TRUE)
             ev <- mds$eig; pv <- 100 * pmax(ev, 0) / sum(pmax(ev, 0))
             xlab <- sprintf("Coordinate 1 (%.1f%%)", pv[1])
             ylab <- sprintf("Coordinate 2 (%.1f%%)", pv[2]); mds$points },
    umap = { if (!requireNamespace("uwot", quietly = TRUE)) stop("embed='umap' needs 'uwot'.")
             uwot::umap(d, n_neighbors = min(15, n - 1), n_components = 2) },
    tsne = { if (!requireNamespace("Rtsne", quietly = TRUE)) stop("embed='tsne' needs 'Rtsne'.")
             Rtsne::Rtsne(d, is_distance = TRUE, dims = 2,
                          perplexity = max(5, min(30, floor((n - 1) / 3))))$Y })
  # tidy: one row per (object, k); facet label carries the silhouette
  sils <- stats::setNames(numeric(length(k_range)), k_range)
  dat  <- do.call(rbind, lapply(k_range, function(k) {
    cl  <- stats::cutree(hc, k)
    sil <- mean(cluster::silhouette(cl, d)[, "sil_width"])
    sils[as.character(k)] <<- sil
    data.frame(x = Y[, 1], y = Y[, 2], cluster = factor(cl),
               panel = sprintf("k = %d   (silhouette %.2f)", k, sil),
               stringsAsFactors = FALSE)
  }))
  dat$panel <- factor(dat$panel, levels = unique(dat$panel))
  cols <- grDevices::hcl.colors(max(k_range), palette)
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y, colour = cluster)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::facet_wrap(~ panel, ncol = ncol) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::labs(x = xlab, y = ylab) +
    .mc_theme(base_size, grid = "none") +
    ggplot2::theme(panel.spacing = grid::unit(1.0, "lines"), aspect.ratio = 1)
  print(p)
  attr(p, "silhouette") <- sils
  invisible(p)
}

#' Two-dimensional cluster visualisation (PCA / UMAP / t-SNE)
#'
#' Project the objects to two dimensions and colour them by cluster. \code{"pca"}
#' uses classical MDS of the fused dissimilarity (always available); \code{"umap"}
#' needs \pkg{uwot} and \code{"tsne"} needs \pkg{Rtsne}.
#'
#' @param x A fitted result (\code{$DistM}), a dissimilarity, or a \code{List}.
#' @param labels Cluster label of each object.
#' @param method \code{"pca"} (default), \code{"umap"} or \code{"tsne"}.
#' @param col Optional cluster colours.
#' @param main Plot title.
#' @return Invisibly the 2-column coordinate matrix.
#' @seealso \code{\link{cluster_selection_plot}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#' embedding_plot(fit, mosaic_labels(fit, 3), method = "pca")
#' @export
embedding_plot <- function(x, labels, method = c("pca", "umap", "tsne"),
                           col = NULL, main = NULL) {
  method <- match.arg(method)
  D <- .as_distm(x); n <- nrow(D)
  Y <- switch(method,
    pca  = stats::cmdscale(stats::as.dist(D), k = 2),
    umap = { if (!requireNamespace("uwot", quietly = TRUE))
               stop("method = 'umap' needs the 'uwot' package.")
             uwot::umap(stats::as.dist(D), n_neighbors = min(15, n - 1), n_components = 2) },
    tsne = { if (!requireNamespace("Rtsne", quietly = TRUE))
               stop("method = 'tsne' needs the 'Rtsne' package.")
             Rtsne::Rtsne(as.dist(D), is_distance = TRUE, dims = 2,
                          perplexity = max(5, min(30, floor((n - 1) / 3))))$Y })
  cl <- factor(labels)
  if (is.null(col)) col <- .mc_pal(nlevels(cl))
  if (is.null(main)) main <- paste0(toupper(method), " coloured by cluster")
  plot(Y, col = col[as.integer(cl)], pch = 19, xlab = paste(toupper(method), "1"),
       ylab = paste(toupper(method), "2"), main = main)
  graphics::legend("topright", levels(cl), col = col[seq_len(nlevels(cl))],
                   pch = 19, bty = "n", title = "cluster", cex = 0.8)
  invisible(Y)
}

#' Visually compare clustering solutions with ComparePlot (no ground truth)
#'
#' Aligns several clustering solutions of the same objects and draws them with
#' \code{\link{ComparePlot}} so their agreement can be read off directly -- which
#' objects the methods place together and where they disagree -- without any
#' ground truth. Each solution may be a fitted result (with \code{$Clust} /
#' \code{$DistM}) or a plain label vector; label vectors are turned into a
#' co-membership dissimilarity so every solution is comparable.
#'
#' @param solutions A named list of fitted results and/or label vectors, all over
#'   the same objects in the same order.
#' @param nrclusters Number of clusters to cut each solution into.
#' @param cols Optional cluster colours; defaults to a distinct palette.
#' @param cex.labels,margins Passed to \code{\link{ComparePlot}}.
#' @param ... Further arguments to \code{\link{ComparePlot}}.
#' @return Invisibly \code{NULL}; called for the plot. Needs \pkg{circlize} and
#'   \pkg{plotrix}.
#' @seealso \code{\link{clustering_concordance}}, \code{\link{ComparePlot}}
#' @examples
#' data(mosaic_toy)
#' sols <- list(
#'   var    = M_ABC(mosaic_toy$List, numsim = 40, NC = 3, weighting = c("var", "var")),
#'   nugget = M_ABC(mosaic_toy$List, numsim = 40, NC = 3, weighting = c("nugget", "nugget")),
#'   NEMO   = NEMO(mosaic_toy$List, k = 3)$cluster)
#' \dontrun{compare_methods_plot(sols, nrclusters = 3)}
#' @export
compare_methods_plot <- function(solutions, nrclusters, cols = NULL,
                                 cex.labels = 0.5, margins = c(5, 11, 3, 2), ...) {
  if (!requireNamespace("circlize", quietly = TRUE) ||
      !requireNamespace("plotrix", quietly = TRUE))
    stop("compare_methods_plot() needs the 'circlize' and 'plotrix' packages.")
  # object names: from the first fit's DistM, else first labelled vector, else index
  first <- solutions[[1]]
  samp  <- if (is.list(first) && !is.null(first$DistM)) rownames(as.matrix(first$DistM))
           else names(first)
  n <- if (is.list(first) && !is.null(first$DistM)) nrow(as.matrix(first$DistM)) else length(first)
  if (is.null(samp)) samp <- sprintf("obj%03d", seq_len(n))

  # build a hierarchy + co-membership distance from *any* solution shape:
  # a fitted result with a usable $Clust, or a plain label vector, and fall back
  # to the object's labels (via mosaic_labels) if $Clust cannot be coerced.
  from_labels <- function(cl) {
    cl <- as.integer(factor(cl)); if (length(cl) != length(samp)) cl <- rep(1L, length(samp))
    D  <- outer(cl, cl, function(a, b) as.numeric(a != b)); dimnames(D) <- list(samp, samp)
    hc <- stats::as.hclust(cluster::agnes(stats::as.dist(D), method = "average"))
    hc$labels <- samp; list(hc = hc, D = D)
  }
  wrap <- function(x) {
    hc <- NULL; D <- NULL
    if (is.list(x) && !is.null(x$Clust)) {
      hc <- tryCatch(stats::as.hclust(x$Clust), error = function(e) NULL)
      if (!is.null(hc)) {
        D <- tryCatch(as.matrix(x$DistM), error = function(e) NULL)
        if (is.null(D) || length(D) == 0L || nrow(D) != length(samp)) {
          cl <- stats::cutree(hc, min(nrclusters, length(hc$order)))
          D  <- outer(cl, cl, function(a, b) as.numeric(a != b))
        }
        if (is.null(rownames(D))) dimnames(D) <- list(samp, samp)
        if (is.null(hc$labels)) hc$labels <- samp
      }
    }
    if (is.null(hc)) {                       # label vector, or a $Clust we could not coerce
      cl <- if (is.list(x)) tryCatch(mosaic_labels(x, nrclusters), error = function(e) NULL) else x
      fl <- from_labels(if (is.null(cl)) rep(1L, length(samp)) else cl)
      hc <- fl$hc; D <- fl$D
    }
    hc$order.lab <- hc$labels[hc$order]
    el <- list(Clust = hc, DistM = D); attr(el, "method") <- "Single Clustering"; el
  }
  res <- lapply(solutions, wrap)
  # ComparePlot may need a few extra colours when a method has a cluster that
  # does not match the reference during alignment (otherwise those cells render
  # white). Give the palette a small buffer.
  if (is.null(cols)) cols <- .mc_pal(nrclusters + 3L)
  ComparePlot(res, nrclusters = nrclusters, cols = cols, names = names(solutions),
              circle = FALSE, plottype = "sweave", cex.labels = cex.labels,
              margins = margins, ...)
  invisible(NULL)
}

#' Compare clustering solutions to each other (no ground truth needed)
#'
#' Given several labellings of the same objects (e.g. different weightings or
#' methods), compute their pairwise agreement (ARI) and each object's stability
#' -- how consistently it is co-clustered across the solutions. Useful for real
#' data where no ground truth is available.
#'
#' @param labels_list A named list of label vectors, all of the same length.
#' @return A list with \code{agreement} (pairwise ARI matrix), \code{stability}
#'   (mean pairwise co-membership per object) and \code{changed} (objects whose
#'   label is not constant across solutions).
#' @seealso \code{\link{cluster_agreement}}, \code{\link{compare_clusterings}}
#' @examples
#' data(mosaic_toy)
#' v <- mosaic_labels(M_ABC(mosaic_toy$List, numsim = 40, NC = 3,
#'                            weighting = c("var", "var")), 3)
#' n <- mosaic_labels(M_ABC(mosaic_toy$List, numsim = 40, NC = 3,
#'                            weighting = c("nugget", "nugget")), 3)
#' clustering_concordance(list(var = v, nugget = n))$agreement
#' @export
clustering_concordance <- function(labels_list) {
  stopifnot(length(labels_list) >= 2)
  nm <- names(labels_list); if (is.null(nm)) nm <- paste0("sol", seq_along(labels_list))
  m  <- length(labels_list); n <- length(labels_list[[1]])
  A  <- matrix(1, m, m, dimnames = list(nm, nm))
  for (i in seq_len(m)) for (j in seq_len(m)) if (i < j) {
    a <- cluster_agreement(labels_list[[i]], labels_list[[j]])[["ARI"]]
    A[i, j] <- A[j, i] <- a
  }
  # per-object stability: average pairwise co-membership agreement across solutions
  co <- matrix(0, n, n)
  for (l in labels_list) co <- co + outer(l, l, function(a, b) as.numeric(a == b))
  co <- co / m
  stability <- (rowSums(co) - 1) / (n - 1)
  changed <- which(apply(vapply(labels_list, function(l)
    match(l, unique(l)), integer(n)), 1, function(r) length(unique(r)) > 1))
  list(agreement = A, stability = stability, changed = changed)
}
