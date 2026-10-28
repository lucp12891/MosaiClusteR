#' Modern Gene Ontology / pathway-enrichment plots
#'
#' Publication-quality replacements for the classic green GO DAG and the
#' horizontal bar plot. Given the output of [PathwayAnalysis()] (or a single
#' MLP result table) it draws either a \emph{dot plot} that compares the enriched
#' terms across clusters (the enrichment-map / \code{clusterProfiler::dotplot}
#' idiom) or a \emph{lollipop} plot for a single cluster. Dot area encodes the
#' number of genes in the term and colour encodes \eqn{-\log_{10}} p on a
#' perceptually-uniform viridis scale.
#'
#' @param pathways Output of [PathwayAnalysis()] (a list of sources, each a list
#'   of per-cluster results carrying a `$Pathways` table), a single source's
#'   per-cluster list, or a single MLP/`data.frame` result table.
#' @param style `"dotplot"` to compare clusters side by side, `"lollipop"`
#'   for one cluster, or `"circular"` for a radial (sunburst-style) bar plot
#'   that groups the terms by cluster around a central hub; bar length encodes
#'   \eqn{-\log_{10}} p, colour encodes the cluster, and each term is a
#'   single-line radial label that stays legible as the number of terms grows.
#' @param top Number of top terms (by p-value) to show per cluster.
#' @param cluster For `"lollipop"` (or to subset a dot plot): which cluster to
#'   show, by name (e.g. `"Cluster 1"`) or index. `NULL` shows all.
#' @param source Which source (index or name) of a full [PathwayAnalysis()]
#'   object to plot. Default `1`.
#' @param title Optional plot title.
#' @param palette viridis option passed to [ggplot2::scale_colour_viridis_c()]
#'   (`"D"` viridis, `"C"` plasma, `"B"` inferno, `"G"` mako, ...).
#' @param max_chars Term descriptions longer than this are wrapped onto two
#'   lines for legibility (dot plot / lollipop).
#' @param wrap For `style = "circular"`, the character width at which long term
#'   labels are wrapped onto successive lines, so long names stay compact and
#'   legible around the ring. Smaller values wrap more aggressively.
#' @param cluster_cols For `style = "circular"`, an optional vector of fill
#'   colours for the clusters (cluster bars and legend). Supply the same palette
#'   used elsewhere - e.g. the `cols` passed to [Cyclogram()] - so the clusters
#'   share one colour scheme across figures. Positional by cluster order, or a
#'   vector named by cluster. Default `NULL` uses the package palette.
#' @return A \pkg{ggplot} object (invisibly); it is also printed so the function
#'   works directly on an open graphics device.
#' @examples
#' \dontrun{
#' P <- PathwayAnalysis(List = L, geneExpr = ge, nrclusters = 3,
#'                      method = c("limma", "MLP"), geneInfo = geneInfo,
#'                      geneSetSource = "GOBP", sign = 0.05, niter = 2)
#' plot_pathways(P, style = "dotplot", top = 8)          # compare clusters
#' plot_pathways(P, style = "lollipop", cluster = 1)     # one cluster
#' }
#' @seealso [PathwayAnalysis()], [PlotPathways()]
#' @export
plot_pathways <- function(pathways, style = c("dotplot", "lollipop", "circular"), top = 8,
                          cluster = NULL, source = 1, title = NULL,
                          palette = "D", max_chars = 42, wrap = 20,
                          cluster_cols = NULL) {
  style <- match.arg(style)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_pathways() needs the 'ggplot2' package.")

  ## ---- 1. coerce any accepted input into a named list of $Pathways tables ----
  as_table <- function(x) {
    if (is.data.frame(x)) return(x)
    if (is.list(x) && !is.null(x$Pathways)) return(as.data.frame(x$Pathways))
    as.data.frame(x)
  }
  is_cluster_list <- function(x)
    is.list(x) && length(x) && all(vapply(x, function(e)
      is.data.frame(e) || (is.list(e) && !is.null(e$Pathways)), logical(1)))

  if (is.data.frame(pathways)) {                      # a single result table
    tabs <- list("Cluster 1" = as.data.frame(pathways))
  } else if (is_cluster_list(pathways)) {             # one source: clusters -> results
    tabs <- lapply(pathways, as_table)
  } else if (is.list(pathways)) {                     # full PathwayAnalysis object
    src <- if (is.character(source)) source else min(source, length(pathways))
    inner <- pathways[[src]]
    if (!is_cluster_list(inner)) stop("Could not locate per-cluster pathway tables in 'pathways'.")
    tabs <- lapply(inner, as_table)
  } else stop("Unrecognised 'pathways' input.")
  if (is.null(names(tabs))) names(tabs) <- paste("Cluster", seq_along(tabs))

  if (!is.null(cluster)) {
    sel <- if (is.numeric(cluster)) cluster else which(names(tabs) %in% cluster)
    tabs <- tabs[sel]
  }

  ## ---- 2. tidy the tables: term, p, gene count, cluster ---------------------
  pick <- function(d, patterns) {
    for (p in patterns) { j <- grep(p, colnames(d), ignore.case = TRUE); if (length(j)) return(colnames(d)[j[1]]) }
    NULL
  }
  wrapfun <- function(s) vapply(s, function(t) {
    t <- as.character(t)
    if (nchar(t) <= max_chars) return(t)
    words <- strsplit(t, " ")[[1]]; line <- ""; out <- character(0)
    for (w in words) { if (nchar(paste(line, w)) > max_chars) { out <- c(out, line); line <- w }
                       else line <- if (nzchar(line)) paste(line, w) else w }
    paste(c(out, line), collapse = "\n")
  }, character(1))

  rows <- list()
  for (cl in names(tabs)) {
    d <- tabs[[cl]]; if (!nrow(d)) next
    pcol <- pick(d, c("pvalue", "p\\.val", "^p$", "mean.*p"))
    dcol <- pick(d, c("geneSetDescription", "description", "term", "Term"))
    tcol <- pick(d, c("testedGeneSetSize", "tested", "Count", "count"))
    ncol_ <- pick(d, c("totalGeneSetSize", "total", "Size"))
    term  <- if (!is.null(dcol)) as.character(d[[dcol]]) else rownames(d)
    p     <- if (!is.null(pcol)) suppressWarnings(as.numeric(d[[pcol]])) else NA_real_
    size  <- if (!is.null(tcol)) suppressWarnings(as.numeric(d[[tcol]])) else 1
    tot   <- if (!is.null(ncol_)) suppressWarnings(as.numeric(d[[ncol_]])) else NA_real_
    ord   <- order(p, na.last = TRUE)
    keep  <- utils::head(ord, top)
    rows[[cl]] <- data.frame(
      cluster = cl, term = term[keep],
      neglogp = -log10(pmax(p[keep], .Machine$double.xmin)),
      genes = size[keep], ratio = size[keep] / tot[keep],
      stringsAsFactors = FALSE)
  }
  dat <- do.call(rbind, rows)
  if (is.null(dat) || !nrow(dat)) stop("No enriched terms to plot.")
  dat$term <- wrapfun(dat$term)
  # order terms by their strongest signal (most significant at the top)
  ord_term <- stats::aggregate(neglogp ~ term, dat, max)
  dat$term <- factor(dat$term, levels = ord_term$term[order(ord_term$neglogp)])
  dat$cluster <- factor(dat$cluster, levels = names(tabs))

  gg  <- ggplot2::ggplot
  aes <- ggplot2::aes
  vir <- ggplot2::scale_colour_viridis_c
  plab <- expression(-log[10] * "(" * italic(p) * "-value)")   # clear axis / legend label

  ## ---- 3. draw -------------------------------------------------------------
  if (style == "dotplot") {
    p <- gg(dat, aes(x = cluster, y = term)) +
      ggplot2::geom_point(aes(size = genes, colour = neglogp)) +
      vir(option = palette, name = plab, direction = -1,
          guide = ggplot2::guide_colourbar(barheight = 7, barwidth = 0.9)) +
      ggplot2::scale_size_continuous(range = c(2.5, 10), name = "Gene count") +
      ggplot2::labs(x = NULL, y = NULL, title = title) +
      .mc_theme(12, grid = "both") +
      ggplot2::theme(axis.line = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_text(face = "plain", size = 10),
                     legend.position = "right")
  } else if (style == "circular") {
    ## radial (sunburst-style) bar plot: one spoke per term, grouped by cluster
    D <- dat
    D$value <- D$neglogp
    clusters <- levels(dat$cluster); clusters <- clusters[clusters %in% D$cluster]
    D <- D[order(D$cluster, -D$value), ]
    nbar <- nrow(D)
    gap  <- max(1L, round(nbar / 20))                  # inter-cluster spacer, scales with size
    # wrap long names onto several lines; wrap less as spokes get thin, so the
    # (tangentially stacked) multi-line labels do not crowd their neighbours
    wrap_eff <- max(wrap, round(wrap * nbar / 18))
    D$term <- vapply(D$term, function(t) {
      if (is.na(t)) return(NA_character_)
      paste(strwrap(gsub("\n", " ", t), width = wrap_eff), collapse = "\n")
    }, character(1))
    DD <- NULL; id <- 0L
    for (cl in clusters) {
      dd <- D[D$cluster == cl, , drop = FALSE]; if (!nrow(dd)) next
      dd$rid <- id + seq_len(nrow(dd)); id <- id + nrow(dd)
      sp <- dd[rep(1L, gap), , drop = FALSE]
      sp$term <- NA; sp$value <- NA; sp$genes <- NA; sp$rid <- id + seq_len(gap); id <- id + gap
      DD <- rbind(DD, dd, sp)
    }
    n    <- nrow(DD)
    ymax <- max(DD$value, na.rm = TRUE)
    hole <- ymax * 0.55; pad <- ymax * 0.06            # centre hole + label offset
    txt  <- max(2.8, 4.0 - 0.03 * nbar)                # shrink text as spokes grow
    lines   <- unlist(strsplit(DD$term[!is.na(DD$term)], "\n"))
    maxline <- max(nchar(lines), 1L)                   # longest single label line
    outer   <- ymax * (1.35 + 0.02 * maxline)          # outer radius tracks label length
    DD$ang <- 90 - 360 * (DD$rid - 0.5) / n
    DD$hj  <- ifelse(DD$ang < -90, 1, 0)
    DD$ang <- ifelse(DD$ang < -90, DD$ang + 180, DD$ang)
    DD$cluster <- factor(DD$cluster, levels = clusters)
    arc <- do.call(rbind, lapply(clusters, function(cl) {
      ids <- DD$rid[DD$cluster == cl & !is.na(DD$term)]
      data.frame(cluster = cl, start = min(ids) - 0.45, end = max(ids) + 0.45)
    }))
    arc$cluster <- factor(arc$cluster, levels = clusters)
    cols <- if (is.null(cluster_cols)) .mc_pal(length(clusters))
            else if (!is.null(names(cluster_cols))) unname(cluster_cols[clusters])
            else rep_len(cluster_cols, length(clusters))
    names(cols) <- clusters
    yarc <- -hole * 0.30
    p <- gg(DD[!is.na(DD$value), , drop = FALSE], aes(x = rid, y = value, fill = cluster)) +
      ggplot2::geom_col(width = 1, na.rm = TRUE, colour = "white", linewidth = 0.9) +
      ggplot2::geom_segment(data = arc,
                            aes(x = start, xend = end, y = yarc, yend = yarc, colour = cluster),
                            inherit.aes = FALSE, linewidth = 3, lineend = "round", show.legend = FALSE) +
      ggplot2::geom_text(data = DD[!is.na(DD$term), , drop = FALSE],
                         aes(x = rid, y = value + pad, label = term, angle = ang, hjust = hj),
                         inherit.aes = FALSE, size = txt, lineheight = 0.82) +
      ggplot2::scale_fill_manual(values = cols, name = NULL) +
      ggplot2::scale_colour_manual(values = cols, guide = "none") +
      ggplot2::ylim(-hole, outer) +
      ggplot2::coord_polar(start = 0) +
      ggplot2::labs(title = title) +
      ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 5))) +
      ggplot2::theme_void() +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 13),
                     legend.position = "bottom",
                     legend.text = ggplot2::element_text(size = 11, face = "bold"),
                     plot.margin = ggplot2::margin(0, 2, 0, 2))
  } else {
    rng <- max(dat$neglogp, na.rm = TRUE)
    p <- gg(dat, aes(x = neglogp, y = term)) +
      ggplot2::geom_segment(aes(x = 0, xend = neglogp, yend = term),
                            colour = "grey78", linewidth = 0.7) +
      ggplot2::geom_point(aes(size = genes, colour = neglogp)) +
      ggplot2::geom_text(aes(label = sprintf("%.1f", neglogp)),
                         hjust = -0.45, size = 3.4, colour = "#111827") +
      vir(option = palette, name = plab, direction = -1) +
      ggplot2::scale_size_continuous(range = c(3, 10), name = "Gene count") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.12)),
                                  limits = c(0, rng * 1.12)) +
      ggplot2::labs(x = plab, y = NULL, title = title) +
      .mc_theme(12, grid = "x") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(face = "plain", size = 10),
                     legend.position = "right")
    if (length(unique(dat$cluster)) > 1)
      p <- p + ggplot2::facet_wrap(~ cluster, scales = "free_y")
  }
  print(p)
  invisible(p)
}
