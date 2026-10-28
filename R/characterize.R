# =============================================================================
# characterize.R
# Post-clustering biological / clinical characterisation. Given cluster labels
# and a table of external variables (measured on the same objects), summarise
# how the clusters differ -- with effect sizes and gradient detection, not just
# p-values -- so that unsupervised clusters can be interpreted biologically.
# =============================================================================

# effect sizes ---------------------------------------------------------------
.kw_eta2 <- function(x, g) {                 # rank eta^2 for Kruskal-Wallis
  ok <- !is.na(x) & !is.na(g); x <- x[ok]; g <- factor(g[ok])
  n <- length(x); k <- nlevels(g)
  if (n <= k || k < 2) return(NA_real_)
  H <- suppressWarnings(stats::kruskal.test(x, g)$statistic)
  max(0, (H - k + 1) / (n - k))
}
.cramers_v <- function(tab) {
  n <- sum(tab); if (n == 0) return(NA_real_)
  chi <- suppressWarnings(stats::chisq.test(tab)$statistic)
  sqrt(as.numeric(chi) / (n * (min(dim(tab)) - 1)))
}

#' Characterise clusters with external (biological / clinical) variables
#'
#' For a clustering, summarise how each external variable differs across the
#' clusters. Continuous variables get \eqn{n}, mean, SD, median, IQR and range
#' \emph{per cluster}, an overall Kruskal-Wallis test, a rank \eqn{\eta^2} effect
#' size, and \strong{gradient detection} (are the cluster means monotonically
#' ordered?). Categorical variables get per-cluster counts and column
#' percentages, a chi-squared / Fisher test and Cramer's V. The emphasis is on
#' effect sizes and descriptive summaries rather than p-values alone.
#'
#' @param labels Cluster label of each object (vector; \code{NA} allowed).
#' @param metadata A data frame of external variables, one row per object, in the
#'   same order as \code{labels}.
#' @param continuous,categorical Optional character vectors naming the columns to
#'   treat as continuous / categorical. If \code{NULL}, numeric columns are taken
#'   as continuous and the rest as categorical.
#'
#' @return A list with \code{continuous} (per-cluster summary, long format),
#'   \code{continuous_tests} (one row per variable: \code{KW_p}, \code{eta2},
#'   \code{gradient} -- the clusters ordered by mean -- and a monotone-trend
#'   \code{trend_rho}/\code{trend_p}), \code{categorical} (counts and column
#'   percentages) and \code{categorical_tests} (\code{p}, \code{cramers_v}).
#'
#' @seealso \code{\link{plot_cluster_profile}}, \code{\link{cluster_biology_panel}}
#' @examples
#' data(mosaic_toy)
#' fit <- M_ABC(mosaic_toy$List, numsim = 40, NC = 3)
#' cl  <- mosaic_labels(fit, 3)
#' meta <- data.frame(score = rnorm(length(cl)), grp = sample(c("A", "B"), length(cl), TRUE))
#' ch  <- characterize_clusters(cl, meta)
#' ch$continuous_tests
#' @export
characterize_clusters <- function(labels, metadata,
                                  continuous = NULL, categorical = NULL) {
  metadata <- as.data.frame(metadata)
  if (length(labels) != nrow(metadata))
    stop("length(labels) must equal nrow(metadata).")
  cl <- factor(labels)
  keep <- !is.na(cl); cl <- droplevels(cl[keep]); metadata <- metadata[keep, , drop = FALSE]
  levs <- levels(cl)

  if (is.null(continuous) && is.null(categorical)) {
    num <- vapply(metadata, is.numeric, logical(1))
    continuous  <- names(metadata)[num]
    categorical <- names(metadata)[!num]
  }
  continuous  <- intersect(continuous,  names(metadata))
  categorical <- intersect(categorical, names(metadata))

  ## ---- continuous ----------------------------------------------------------
  cont_rows <- list(); cont_tests <- list()
  for (v in continuous) {
    x <- suppressWarnings(as.numeric(metadata[[v]]))
    per <- lapply(levs, function(l) {
      xi <- x[cl == l & !is.na(x)]
      data.frame(variable = v, cluster = l, n = length(xi),
                 mean = mean(xi), sd = stats::sd(xi),
                 median = stats::median(xi),
                 q25 = stats::quantile(xi, .25, names = FALSE),
                 q75 = stats::quantile(xi, .75, names = FALSE),
                 min = suppressWarnings(min(xi)), max = suppressWarnings(max(xi)),
                 stringsAsFactors = FALSE)
    })
    cont_rows[[v]] <- do.call(rbind, per)
    # gradient: order clusters by their mean, test a monotone trend across it
    cmean <- vapply(levs, function(l) mean(x[cl == l], na.rm = TRUE), numeric(1))
    ord   <- order(cmean)
    grad  <- paste(sprintf("%s(%.3g)", levs[ord], cmean[ord]), collapse = " < ")
    rankpos <- stats::setNames(order(ord), levs)          # position of each cluster in the mean-ordering
    obj_rank <- rankpos[as.character(cl)]
    tr <- suppressWarnings(stats::cor.test(x, obj_rank, method = "spearman"))
    kw <- suppressWarnings(stats::kruskal.test(x, cl))
    cont_tests[[v]] <- data.frame(variable = v, KW_p = kw$p.value,
                                  eta2 = .kw_eta2(x, cl), gradient = grad,
                                  trend_rho = unname(tr$estimate), trend_p = tr$p.value,
                                  stringsAsFactors = FALSE)
  }

  ## ---- categorical ---------------------------------------------------------
  cat_rows <- list(); cat_tests <- list()
  for (v in categorical) {
    f <- factor(metadata[[v]])
    tab <- table(f, cl)
    pct <- sweep(tab, 2, colSums(tab), "/") * 100
    df  <- as.data.frame.matrix(tab); df$level <- rownames(df); df$variable <- v
    cat_rows[[v]] <- df
    p <- tryCatch(suppressWarnings(
      if (any(tab < 5)) stats::fisher.test(tab, simulate.p.value = TRUE)$p.value
      else stats::chisq.test(tab)$p.value), error = function(e) NA_real_)
    cat_tests[[v]] <- data.frame(variable = v, p = p, cramers_v = .cramers_v(tab),
                                 stringsAsFactors = FALSE)
  }

  list(continuous       = if (length(cont_rows)) do.call(rbind, cont_rows) else NULL,
       continuous_tests = if (length(cont_tests)) do.call(rbind, cont_tests) else NULL,
       categorical      = cat_rows,
       categorical_tests = if (length(cat_tests)) do.call(rbind, cat_tests) else NULL,
       n_clusters = nlevels(cl))
}

#' Plot one external variable by cluster
#'
#' Boxplot (continuous) or grouped bar plot of column proportions (categorical)
#' of an external variable across clusters -- the per-cluster distribution and
#' gradient the manuscript should show.
#'
#' @param labels Cluster labels.
#' @param x The variable (numeric -> boxplot; factor/character -> barplot).
#' @param name Axis / title label for the variable.
#' @param col Optional vector of cluster colours.
#' @return Invisibly \code{NULL}; called for the plot.
#' @seealso \code{\link{characterize_clusters}}, \code{\link{cluster_biology_panel}}
#' @export
plot_cluster_profile <- function(labels, x, name = deparse(substitute(x)), col = NULL) {
  df <- data.frame(x = x); names(df) <- name
  cluster_biology_panel(labels, df, ncol = 1)
}

#' Panel of cluster characterisation plots
#'
#' Lay out \code{\link{plot_cluster_profile}} for every external variable in a
#' single figure -- boxplots for continuous variables, bar plots for categorical
#' ones -- ready for the manuscript.
#'
#' @param labels Cluster labels.
#' @param metadata Data frame of external variables (rows = objects).
#' @param vars Optional subset of columns to show (default: all).
#' @param ncol Number of columns in the panel grid.
#' @param base_size Base font size passed to the plot theme.
#' @return Invisibly \code{NULL}.
#' @seealso \code{\link{characterize_clusters}}
#' @export
cluster_biology_panel <- function(labels, metadata, vars = NULL, ncol = 3,
                                  base_size = 13) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("cluster_biology_panel() needs the 'ggplot2' package.")
  metadata <- as.data.frame(metadata)
  if (is.null(vars)) vars <- names(metadata)
  vars <- intersect(vars, names(metadata))
  cl <- factor(labels)
  is_num <- vapply(metadata[vars], is.numeric, logical(1))
  catv <- vars[!is_num]; conv <- vars[is_num]

  # categorical variables -> column-percentage grouped bars
  catdf <- if (length(catv)) do.call(rbind, lapply(catv, function(v) {
    x <- factor(metadata[[v]]); ok <- !is.na(x) & !is.na(cl)
    pct <- sweep(table(x[ok], cl[ok]), 2, table(cl[ok])[colnames(table(x[ok], cl[ok]))], "/") * 100
    d <- as.data.frame(as.table(pct)); names(d) <- c("status", "cluster", "value")
    d$variable <- v; d
  })) else NULL
  # continuous variables -> boxplots
  condf <- if (length(conv)) {
    d <- do.call(rbind, lapply(conv, function(v)
      data.frame(cluster = cl, value = metadata[[v]], variable = v)))
    d[!is.na(d$value) & !is.na(d$cluster), ]
  } else NULL

  lev <- c(catv, conv)
  if (!is.null(catdf)) catdf$variable <- factor(catdf$variable, levels = lev)
  if (!is.null(condf)) condf$variable <- factor(condf$variable, levels = lev)
  ns   <- if (!is.null(catdf)) nlevels(factor(catdf$status)) else 0L
  spal <- if (ns <= 2) c("#4E79A7", "#F28E2B")[seq_len(max(ns, 1L))] else .mc_pal(ns)

  p <- ggplot2::ggplot()
  if (!is.null(catdf))
    p <- p + ggplot2::geom_col(data = catdf,
             ggplot2::aes(x = cluster, y = value, fill = status),
             position = ggplot2::position_dodge2(preserve = "single"), width = 0.72)
  if (!is.null(condf))
    p <- p + ggplot2::geom_boxplot(data = condf,
             ggplot2::aes(x = cluster, y = value), fill = "#9CA3AF",
             colour = "#374151", width = 0.6, outlier.size = 0.5)
  p <- p +
    ggplot2::facet_wrap(~ variable, ncol = ncol, scales = "free_y") +
    ggplot2::scale_fill_manual(values = spal, name = "Status") +
    ggplot2::labs(x = "Cluster", y = if (length(conv)) NULL else "Percentage (%)") +
    .mc_theme(base_size, grid = "none") +
    ggplot2::theme(panel.spacing = grid::unit(0.9, "lines"))
  print(p)
  invisible(p)
}
