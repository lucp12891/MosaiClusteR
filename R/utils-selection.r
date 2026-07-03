# =============================================================================
# utils-selection.r  - cluster-number selection and cluster-set similarity
#
# Utility functions ported from the legacy MoSaIC code base. Inputs follow the
# package convention: a `List` of data / distance / clustering objects, with
# objects in the rows of each data matrix.
# =============================================================================

#' Select the number of clusters via silhouette widths
#'
#' For each element of \code{List}, partitions around medoids
#' (\code{cluster::pam}) are computed over a range of cluster numbers and the
#' average silhouette width is recorded. The number of clusters maximising the
#' average silhouette width is reported per source and for the across-source
#' average. The analytical result (the silhouette table and the chosen number
#' of clusters) is returned and is testable without producing any plot.
#'
#' @param List A list of data matrices (objects in rows), distance matrices, or
#'   \code{cluster::pam} outputs, depending on \code{type}.
#' @param type One of \code{"data"} (compute a distance first via
#'   \code{\link{Distance}}), \code{"dist"} (elements are dissimilarity
#'   matrices) or \code{"pam"} (elements are \code{pam} objects).
#' @param distmeasure Distance measure(s) used when \code{type = "data"}.
#' @param normalize,method Normalisation controls passed to \code{\link{Distance}}.
#' @param nrclusters Sequence of cluster numbers to evaluate.
#' @param names Optional names for the sources; used in output labels.
#' @param StopRange Logical; if \code{FALSE} and dissimilarities fall outside
#'   \code{[0, 1]} they are range-normalised for comparability.
#' @param plottype One of \code{"none"} (default; no plot, safe for headless
#'   runs), \code{"new"}, \code{"pdf"} or \code{"sweave"}. Any value other than
#'   \code{"none"} attempts to draw the silhouette curve; drawing is guarded so
#'   that a missing graphics device does not abort the computation.
#' @param location File location (without extension) used when
#'   \code{plottype = "pdf"}. Default \code{NULL}.
#' @return A list with two elements: \code{Silhoutte_Widths} (a matrix of
#'   average silhouette widths per source and their mean) and
#'   \code{Optimal_Nr_of_CLusters} (a one-row data frame with the silhouette
#'   optimal number of clusters per source and overall).
#' @references
#' Kaufman, L. and Rousseeuw, P. J. (1990). Finding Groups in Data: An
#' Introduction to Cluster Analysis. Wiley, New York.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' sel <- SelectnrClusters(mosaic_toy$List, type = "data",
#'                         distmeasure = c("euclidean", "euclidean"),
#'                         nrclusters = seq(2, 5), plottype = "none")
#' sel$Optimal_Nr_of_CLusters
#' }
#' @export
SelectnrClusters <- function(List, type = c("data", "dist", "pam"),
                             distmeasure = c("tanimoto", "tanimoto"),
                             normalize = c(FALSE, FALSE), method = c(NULL, NULL),
                             nrclusters = seq(5, 25, 1), names = NULL,
                             StopRange = FALSE, plottype = "none", location = NULL) {

	type = match.arg(type)
	avsilwidth <- matrix(0, ncol = length(List), nrow = length(nrclusters))
	pamfunction <- function(DistM, nrclusters) {
		asw = sapply(nrclusters, function(x) cluster::pam(DistM, x)$silinfo$avg.width)
		return(asw)
	}

	CheckDist <- function(Dist, StopRange) {
		if (StopRange == FALSE & !(0 <= min(Dist) & max(Dist) <= 1)) {
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist = Normalization(Dist, method = "Range")
		}
		else {
			Dist = Dist
		}
	}


	if (type == "data") {
		OrderNames = rownames(List[[1]])
		for (i in 1:length(List)) {
			List[[i]] = List[[i]][OrderNames, ]
		}
		Dist = lapply(seq(length(List)), function(i) Distance(List[[i]], distmeasure[i], normalize[i], method[i]))
		Dist = lapply(seq(length(Dist)), function(i) CheckDist(Dist[[i]], StopRange))

		avsilwidth = sapply(Dist, function(x) pamfunction(x, nrclusters = nrclusters))
		rownames(avsilwidth) = nrclusters
	}
	else if (type == "dist") {
		OrderNames = rownames(List[[1]])
		for (i in 1:length(List)) {
			List[[i]] = List[[i]][OrderNames, OrderNames]
		}
		Dist = List
		Dist = lapply(seq(length(Dist)), function(i) CheckDist(Dist[[i]], StopRange))

		avsilwidth = sapply(Dist, function(x) pamfunction(x, nrclusters = nrclusters))
		rownames(avsilwidth) = nrclusters
	}
	else {
		avsilwidth = sapply(List, function(x) return(x$silinfo$avg.width))
	}

	plottypein <- function(plottype, location) {
		if (plottype == "pdf" & !(is.null(location))) {
			grDevices::pdf(paste(location, ".pdf", sep = ""))
		}
		if (plottype == "new") {
			try(grDevices::dev.new(), silent = TRUE)
		}
		if (plottype == "sweave") {

		}
	}
	plottypeout <- function(plottype) {
		if (plottype == "pdf") {
			grDevices::dev.off()
		}
	}


	if (is.null(names)) {
		names1 = c()
		names2 = c()
		for (i in 1:length(List)) {
			names1 = c(names1, paste("Silhouette widths for Data", i, sep = " "))
			names2 = c(names2, paste("Nr Clusters for Data", i, sep = ' '))
		}
		names1 = c(names1, "Average Silhoutte Widths")
		names2 = c(names2, "Optimal nr of clusters")

	}
	else {
		names1 = c()
		names2 = c()
		for (i in 1:length(List)) {
			names1 = c(names1, paste("Silhouette widths for", names[i], sep = " "))
			names2 = c(names2, paste("Nr Clusters for", names[i], sep = ' '))
		}
		names1 = c(names1, "Average Silhoutte Widths")
		names2 = c(names2, "Optimal nr of clusters")
	}


	rownames(avsilwidth) = nrclusters


	avsil = apply(avsilwidth, 1, mean)
	avsilwidth = cbind(avsilwidth, avsil)
	colnames(avsilwidth) = names1


	plotsil <- function(sils, plottype, location, name) {
		k.best = as.numeric(names(sils)[which.max(sils)])
		cat("silhouette-optimal number of clusters:", k.best, "\n")
		plottypein(plottype, location)
		try({
			graphics::plot(nrclusters, sils, type = "h", main = name,
					xlab = "k  (# clusters)", ylab = "average silhouette width")
			graphics::axis(1, k.best, paste("best", k.best, sep = "\n"), col = "red", col.axis = "red")
		}, silent = TRUE)
		plottypeout(plottype)
	}

	# Plotting is skippable: only attempt to draw when explicitly requested.
	# The analytical result below is always computed and returned.
	if (!is.null(plottype) && plottype != "none") {
		try(
			sapply(c(1:ncol(avsilwidth)), function(x) plotsil(avsilwidth[, x], plottype, location, names1[x])),
			silent = TRUE
		)
	}

	Output = list()
	Output[[1]] = avsilwidth
	nrclusters = apply(avsilwidth, 2, function(x) return(as.numeric(names(x)[which.max(x)])))
	nrclusters = as.data.frame(t(nrclusters))
	colnames(nrclusters) = names2
	rownames(nrclusters) = "NrClusters"


	Output[[2]] = nrclusters

	names(Output) = c("Silhoutte_Widths", "Optimal_Nr_of_CLusters")
	return(Output)
}


#' Similarity of a set of clusterings to a reference
#'
#' Given a colour/label matrix (clusterings in rows, objects in columns) or a
#' list of clustering outputs, the proportion of objects sharing the label of
#' the first (reference) row is computed for every row. When a list is supplied
#' it is first aligned with \code{ReorderToReference}; supplying the label
#' matrix directly avoids that dependency.
#'
#' @param List A label matrix (rows are clusterings, columns are objects) or a
#'   list of clustering outputs to be aligned via \code{ReorderToReference}.
#' @param nrclusters Number of clusters used when \code{List} is a list. Default
#'   \code{NULL}.
#' @param fusionsLog,weightclust Logical controls passed to
#'   \code{ReorderToReference} when \code{List} is a list.
#' @param names Optional method names passed to \code{ReorderToReference}.
#' @return A numeric vector of similarities in \code{[0, 1]}, one per row, where
#'   the first entry is the reference compared to itself (always \code{1}).
#' @references
#' Fodeh, S. J., Brandt, C., Luong, T. B., Haddad, A., Schultz, M., Murphy, T.
#' and Krauthammer, M. (2013). Complementary ensemble clustering of biomedical
#' data. Journal of Biomedical Informatics, 46(3), 436-443.
#' @examples
#' \dontrun{
#' M <- rbind(c(1, 1, 2, 2, 3, 3),
#'            c(1, 1, 2, 3, 3, 3),
#'            c(2, 2, 1, 1, 3, 3))
#' SimilarityMeasure(M)
#' }
#' @export
SimilarityMeasure <- function(List, nrclusters = NULL, fusionsLog = TRUE,
                              weightclust = TRUE, names = NULL) {

	if (!methods::is(List, "list")) {
		MatrixColors = List
	}
	else {
		MatrixColors = ReorderToReference(List, nrclusters, fusionsLog, weightclust, names)
	}

	#Compare every row to the first row
	Similarity = c()
	for (i in 1:dim(MatrixColors)[1]) {
		Shared = 0
		for (j in 1:dim(MatrixColors)[2])
			if (MatrixColors[i, j] == MatrixColors[1, j]) {
				Shared = Shared + 1
			}
		Similarity = c(Similarity, Shared / ncol(MatrixColors))


	}

	return(Similarity)
}
