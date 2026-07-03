
#' @title Link-based cluster ensembles
#'
#' @description Link-based cluster ensembles combine multiple base clusterings
#' into a single, refined similarity matrix by exploiting the link structure
#' between clusters of the ensemble. Three refinement schemes are available:
#' Connected-Triple-Based Similarity ("cts"), SimRank-Based Similarity ("srs")
#' and Approximate SimRank-Based Similarity ("asrs"). The refined object-by-object
#' similarity matrix is turned into a dissimilarity matrix on which a final
#' agglomerative hierarchical clustering is performed.
#' @export
#' @param List A list of data matrices, distance matrices or clustering results.
#'   It is assumed the rows are corresponding with the objects.
#' @param type Indicates whether the provided matrices in "List" are either data
#'   matrices, distance matrices or clustering results obtained from the data.
#'   Should be one of "data", "dist" or "clust".
#' @param distmeasure A vector of the distance measures to be used on each data
#'   matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming".
#'   Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices
#'   or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if
#'   different distance types are used. More details on normalization in
#'   \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile",
#'   "Fisher-Yates", "standardize","Range" or any of the first letters of these
#'   names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data
#'   set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the
#'   agnes function. Defaults to 0.625 and is only used if the linkage is set to
#'   "flexible".
#' @param nrclusters A vector with the number of clusters to cut each base
#'   dendrogram into. If NULL and gap=FALSE, the function stops.
#' @param gap Logical. Whether to use the gap statistic to determine the number
#'   of clusters per data set. Defaults to FALSE.
#' @param maxK The maximal number of clusters to consider when gap=TRUE.
#' @param linkBasedMethod The link-based refinement scheme. One of "cts", "srs"
#'   or "asrs". Defaults to "cts".
#' @param decayfactor The decay (confidence) factor used by the refinement scheme.
#'   Defaults to 0.8.
#' @param niter The number of iterations for the SimRank-based scheme ("srs").
#'   Defaults to 5.
#' @param linkBasedLinkage The linkage used in the final agglomerative clustering
#'   of the refined similarity matrix. Defaults to "ward".
#' @param waitingtime Retained for backward compatibility with the external
#'   executable / MATLAB workflow; ignored by the pure-R schemes.
#' @param file_number Retained for backward compatibility with the external
#'   executable / MATLAB workflow; ignored by the pure-R schemes.
#' @param executable Logical. If TRUE, the refined similarity matrix is computed
#'   by an external executable / MATLAB. The "cts" and "srs" schemes are pure R
#'   and work without an executable. The "asrs" scheme is only available through
#'   the external executable; if executable=FALSE it stops with an informative
#'   message.
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting (refined) distance matrix}
#' \item{Clust}{The resulting hierarchical clustering}
#' The value has class 'LinkBased'.
#' @details Each base clustering of the ensemble is represented as a hard
#' partition. The clusters of all partitions form the vertices of a weighted
#' graph; the link weights between clusters are derived from shared objects. The
#' refined object-by-object similarity is obtained by propagating these
#' cluster-to-cluster links ("cts", "srs", "asrs"). See Iam-On et al. (2010).
#' @references
#' Iam-On, N., Boongoen, T., Garrett, S. and Price, C. (2010). Link-based cluster
#' ensembles for the combination of multiple clusterings. IEEE Transactions on
#' Knowledge and Data Engineering, 23(12), 1839-1853.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res <- LinkBasedClustering(List = L, type = "data",
#'   distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
#'   linkBasedMethod = "cts", decayfactor = 0.8, niter = 5,
#'   linkBasedLinkage = "ward", nrclusters = c(3, 3))
#' }
LinkBasedClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,linkBasedMethod=c("cts","srs","asrs"),decayfactor=0.8,niter=5,linkBasedLinkage="ward",waitingtime=300,file_number=00,executable=FALSE){

	linkBasedMethod <- match.arg(linkBasedMethod)

	# Guard: the "asrs" scheme is only available through the external
	# executable / MATLAB. Refuse early if it is requested without one.
	if(linkBasedMethod=="asrs" && !executable){
		stop("The 'asrs' (Approximate SimRank-Based Similarity) scheme requires an external executable or MATLAB. Set executable=TRUE and make the executable available, or use linkBasedMethod='cts' or 'srs' (pure R).")
	}

	# Step 1: Generate several clustering results on the objects
	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))

		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)

		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

		Clusters=lapply(seq(length(Clusterings)),function(i) stats::cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
	}

	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))

		Dist=List

		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters of put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=clusters
			}
		}

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

		Clusters=lapply(seq(length(Clusterings)),function(i) stats::cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))


	}
	else if(type=="clust"){

		Clusterings=List

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}

		Clusters=lapply(seq(length(Clusterings)),function(i) stats::cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))

	}
	#### Step 1 Complete

	# Step 2: Get the ensemble in matrix format: one row per clustering result,
	# one column per object.
	if(type=="data" | type=="dist"){
		nc=nrow(List[[1]])
	}
	else if(type=="clust"){
		nc=nrow(List[[1]]$DistM)
	}
	matlabdata=t(matrix(unlist(Clusters),ncol=nc,byrow=TRUE))

	if(executable){
		# Legacy path: hand the ensemble to an external executable / MATLAB,
		# which writes the refined similarity matrix to S_<file_number>.csv.
		utils::write.table(matlabdata,file=paste("matlabdata_",file_number,".csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)

		if(linkBasedMethod=="cts"){
			system(paste("./LinkBasedClusteringcts ",file_number,sep=""),intern=TRUE)
		}
		if(linkBasedMethod=="srs"){
			system(paste("./LinkBasedClusteringsrs ",file_number,sep=""),intern=TRUE)
		}
		if(linkBasedMethod=="asrs"){
			system(paste("./LinkBasedClusteringasrs ",file_number,sep=""),intern=TRUE)
		}

		Continue=FALSE
		time=0
		while(Continue==FALSE){
			Sys.sleep(15)
			time=time+15
			Continue=file.exists(paste("S_",file_number,".csv",sep=""))
			if(time>waitingtime & Continue==FALSE){
				stop(paste("Waited",waitingtime, "seconds for completion of the ensemble clustering procedure. Increase waiting time to continue.",sep=" "))
			}
		}

		SimM=as.matrix(utils::read.table(paste("S_",file_number,".csv",sep=""),sep=","))

		file.remove(paste("matlabdata_",file_number,".csv",sep=""))
		file.remove(paste("S_",file_number,".csv",sep=""))
	}
	else{
		# Pure-R path: compute the refined object-by-object similarity matrix
		# directly. matlabdata has one column per base clustering (objects in rows).
		SimM=.linkBasedSimilarity(matlabdata,linkBasedMethod=linkBasedMethod,decayfactor=decayfactor,niter=niter)
	}

	SimM=as.matrix(SimM)
	DistM=as.matrix(1-SimM)
	rownames(DistM)=rownames(Clusterings[[1]]$DistM)
	colnames(DistM)=rownames(Clusterings[[1]]$DistM)

	LinkBasedCluster=cluster::agnes(DistM,diss=TRUE,method=linkBasedLinkage,par.method=alpha)

	Out=list("DistM"=DistM,"Clust"=LinkBasedCluster)
	attr(Out,"method")="LinkBased"

	return(Out)
}


# -----------------------------------------------------------------------------
# Internal helpers: pure-R link-based similarity schemes (Iam-On et al., 2010)
# -----------------------------------------------------------------------------

# Compute the refined object-by-object similarity matrix from an ensemble.
# Ensemble: a matrix with one row per object and one column per base clustering;
# entries are (arbitrary) cluster labels.
.linkBasedSimilarity <- function(Ensemble, linkBasedMethod = c("cts", "srs", "asrs"),
                                 decayfactor = 0.8, niter = 5) {
	linkBasedMethod <- match.arg(linkBasedMethod)
	Ensemble <- as.matrix(Ensemble)
	N <- nrow(Ensemble)   # number of objects
	M <- ncol(Ensemble)   # number of base clusterings

	# Build a binary object-by-cluster membership matrix B (N x sum_of_clusters)
	# and record, for each global cluster, which base clustering it belongs to.
	blocks <- vector("list", M)
	partOf <- integer(0)
	for (m in seq_len(M)) {
		labs <- Ensemble[, m]
		ul <- sort(unique(labs))
		Bm <- matrix(0, N, length(ul))
		for (j in seq_along(ul)) Bm[, j] <- as.numeric(labs == ul[j])
		blocks[[m]] <- Bm
		partOf <- c(partOf, rep(m, length(ul)))
	}
	B <- do.call(cbind, blocks)   # N x C
	C <- ncol(B)

	# Cluster-to-cluster shared-object overlap (used by CTS).
	# Compute cluster-to-cluster similarity S_cluster (C x C).
	if (linkBasedMethod == "cts") {
		S_cluster <- .ctsClusterSim(B, partOf, decayfactor)
	} else { # srs (and asrs falls back here only when forced, but is guarded out)
		S_cluster <- .srsClusterSim(B, decayfactor, niter)
	}

	# Refined object similarity: average cluster similarity of the clusters that
	# each pair of objects belongs to, across the M base clusterings.
	# membership index of object i in base clustering m -> column of B.
	# Build, per object, the list of its C-space cluster column for each m.
	colOf <- matrix(0L, N, M)
	offset <- 0L
	for (m in seq_len(M)) {
		Bm <- blocks[[m]]
		# column (within block) that each object belongs to
		within <- max.col(Bm, ties.method = "first")
		colOf[, m] <- offset + within
		offset <- offset + ncol(Bm)
	}

	Sim <- matrix(0, N, N)
	for (m in seq_len(M)) {
		idx <- colOf[, m]
		Sim <- Sim + S_cluster[idx, idx, drop = FALSE]
	}
	Sim <- Sim / M
	diag(Sim) <- 1
	rownames(Sim) <- rownames(Ensemble)
	colnames(Sim) <- rownames(Ensemble)
	Sim
}

# Connected-Triple-Based Similarity between clusters.
# B: N x C binary membership matrix; partOf: base-clustering index per cluster.
.ctsClusterSim <- function(B, partOf, decayfactor) {
	C <- ncol(B)
	# weighted graph of clusters: edge weight = number of shared objects (only
	# meaningful across different base clusterings).
	W <- crossprod(B)           # C x C, W[p,q] = #objects shared by clusters p,q
	diag(W) <- 0
	# zero out edges within the same base clustering (clusters of one partition
	# never share objects anyway, but keep it explicit).
	same <- outer(partOf, partOf, FUN = "==")
	W[same] <- 0

	# Connected-triple count between two clusters via common neighbours.
	wMax <- max(W)
	if (wMax == 0) {
		S <- diag(1, C)
		return(S)
	}
	# Number of connected triples shared: for clusters i,j the sum over common
	# neighbour k of min(W[i,k], W[k,j]) (simple weighted triple count).
	S <- diag(1, C)
	for (i in seq_len(C)) {
		for (j in seq_len(C)) {
			if (i < j) {
				wc <- sum(pmin(W[i, ], W[, j]))
				S[i, j] <- wc
				S[j, i] <- wc
			}
		}
	}
	maxS <- max(S[upper.tri(S)])
	if (maxS > 0) S <- S / maxS * decayfactor
	diag(S) <- 1
	S
}

# SimRank-Based Similarity between clusters (iterative).
# B: N x C binary membership matrix.
.srsClusterSim <- function(B, decayfactor, niter) {
	C <- ncol(B)
	# bipartite weighted adjacency cluster<->object: B is N x C, so the
	# cluster-to-object weights are columns of B.
	# SimRank on the cluster side, where the "neighbours" of a cluster are its
	# member objects, and the neighbours of an object are its clusters.
	# Normalised neighbour weights.
	clusterDeg <- colSums(B)               # number of objects per cluster
	objectDeg  <- rowSums(B)               # number of clusters per object
	clusterDeg[clusterDeg == 0] <- 1
	objectDeg[objectDeg == 0]   <- 1

	# Transition: from cluster to object (C x N), row-normalised
	P_co <- t(B) / clusterDeg              # C x N
	# Transition: from object to cluster (N x C), row-normalised
	P_oc <- B / objectDeg                  # N x C

	Sc <- diag(1, C)                       # cluster-cluster similarity
	So <- diag(1, nrow(B))                 # object-object similarity
	for (it in seq_len(niter)) {
		Sc_new <- decayfactor * (P_co %*% So %*% t(P_co))
		So_new <- decayfactor * (P_oc %*% Sc %*% t(P_oc))
		diag(Sc_new) <- 1
		diag(So_new) <- 1
		Sc <- Sc_new
		So <- So_new
	}
	maxS <- max(Sc[upper.tri(Sc)])
	if (is.finite(maxS) && maxS > 0) {
		offdiag <- Sc
		diag(offdiag) <- 0
		Sc <- offdiag / maxS * decayfactor
	} else {
		Sc <- matrix(0, C, C)
	}
	diag(Sc) <- 1
	Sc
}
