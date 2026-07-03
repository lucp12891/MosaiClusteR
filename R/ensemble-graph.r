
#' @title Ensemble clustering via graph partitioning (CSPA / HGPA / MCLA)
#'
#' @description Cluster-based ensemble clustering that combines several
#' single-source partitions into a single consensus partition through the
#' graph-partitioning algorithms of Strehl and Ghosh: the Cluster-based
#' Similarity Partitioning Algorithm (CSPA), the HyperGraph Partitioning
#' Algorithm (HGPA) and the Meta-CLustering Algorithm (MCLA), as well as a
#' "Best" selection strategy.
#'
#' \strong{External executable required.} This method does not run in pure R.
#' It writes the stacked label matrix to disk and delegates the actual graph
#' partitioning to an \strong{external graph-partitioning executable}
#' (e.g. a compiled CSPA/HGPA/MCLA binary backed by METIS / \code{gpmetis} /
#' \code{shmetis}) or to a \strong{MATLAB} session. These binaries are
#' \strong{not bundled} with the package; the user must supply and configure
#' them. The path to the executable is selected through the \code{executable}
#' argument. When \code{executable=FALSE} (the default), or when the configured
#' executable cannot be located, the function stops with an informative error
#' rather than attempting to launch an unavailable program.
#'
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clust" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clust".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param ensembleMethod A character string indicating the consensus function to use. One of "CSPA", "HGPA", "MCLA" or "Best". Defaults to "CSPA".
#' @param waitingtime The maximum number of seconds to wait for the external program to write its result before stopping. Defaults to 300.
#' @param file_number An identifier appended to the temporary files exchanged with the external program. Defaults to 00.
#' @param executable Either FALSE (the default) or the path/flag identifying the external graph-partitioning executable (or MATLAB) to use. When FALSE the function stops, since no consensus can be computed in pure R. Precompiled binaries are not bundled with the package.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting consensus clustering}
#' The value has class 'Ensemble'.
#' @references Strehl A. and Ghosh J. (2002). Cluster ensembles - a knowledge reuse framework for combining multiple partitions. Journal of Machine Learning Research, 3, 583-617.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#'
#' # Requires an external graph-partitioning executable / MATLAB to be
#' # configured through the 'executable' argument.
#' res <- EnsembleClustering(List = L, type = "data",
#'   distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
#'   method = c(NULL, NULL), clust = "agnes",
#'   linkage = c("flexible", "flexible"), nrclusters = c(7, 7),
#'   gap = FALSE, maxK = 15, ensembleMethod = "CSPA",
#'   executable = "./EnsembleClusteringC")
#' }
EnsembleClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,ensembleMethod=c("CSPA","HGPA","MCLA","Best"),waitingtime=300,file_number=00,executable=FALSE){

	# This method relies on an external graph-partitioning executable (METIS /
	# gpmetis / shmetis based) or a MATLAB session to compute the consensus
	# partition. Such binaries are not bundled with the package. Refuse to run
	# (with a clear message) when no usable executable has been configured.
	if(isFALSE(executable) || is.null(executable)){
		stop(paste0("EnsembleClustering requires an external graph-partitioning executable ",
			"(e.g. a compiled CSPA/HGPA/MCLA binary backed by METIS / gpmetis / shmetis) ",
			"or MATLAB to compute the consensus partition. These binaries are not bundled ",
			"with the package. Configure one through the 'executable' argument; ",
			"it is currently FALSE."))
	}
	if(is.character(executable) && !nzchar(Sys.which(executable)) && !file.exists(executable)){
		stop(paste0("EnsembleClustering could not locate the configured graph-partitioning ",
			"executable '", executable, "'. Provide a valid path to a compiled CSPA/HGPA/MCLA ",
			"binary (METIS / gpmetis / shmetis) or to MATLAB. Precompiled binaries are not ",
			"bundled with the package."))
	}

	ensembleMethod <- match.arg(ensembleMethod)

	# Step 1: Generate several clustering results on the objects
	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
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

		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
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
				nrclusters=ceiling(mean(clusters))
			}
		}

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))


	}
	else if(type=="clust"){

		Clusterings=List

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}

		Clusters=lapply(seq(length(Clusterings)),function(i) cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))

	}
	#### Step 1 Complete

	# Step 2: Get data in format to be given to matlab =>  matrix with a row per clustering result
	if(type=="data" | type=="dist"){
		nc=nrow(List[[1]])
	}
	else if(type=="clust"){
		nc=nrow(List[[1]]$DistM)
	}

	matlabdata=matrix(unlist(Clusters),ncol=nc,byrow=TRUE)
	utils::write.table(matlabdata,file=paste("matlabdata_",file_number,".csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)

	if(ensembleMethod=="CSPA"){
		system(paste(executable," ",file_number,sep=""),intern=TRUE)
	}
	if(ensembleMethod=="HGPA"){
		system(paste(executable," ",file_number,sep=""),intern=TRUE)
	}
	if(ensembleMethod=="MCLA"){
		system(paste(executable," ",file_number,sep=""),intern=TRUE)
	}
	if(ensembleMethod=="Best"){
		system(paste(executable," -nodisplay -r \"run('EnsembleClusteringB ",file_number,"'); exit\"",sep=""),intern=TRUE)
	}

	Continue=FALSE
	time=0
	while(Continue==FALSE){
		Sys.sleep(15)
		time=time+15
		Continue=file.exists(paste("ClusterEnsemble_",file_number,".csv",sep=""))
		if(time>waitingtime & Continue==FALSE){
			stop(paste("Waited",waitingtime, "seconds for completion of the ensemble clustering procedure. Increase waiting time to continue.",sep=" "))
		}
	}

	ClusterEnsemble <- utils::read.table(paste("ClusterEnsemble_",file_number,".csv",sep=""),sep=",")
	ClusterEnsemble=as.vector(as.matrix(ClusterEnsemble))
	names(ClusterEnsemble)=rownames(List[[1]])


	DistM=NULL


	Clust=list()
	Clust$order=sort(ClusterEnsemble,index=TRUE)$ix
	Clust$order.lab=names(ClusterEnsemble[sort(ClusterEnsemble,index=TRUE)$ix])
	Clust$Clusters=ClusterEnsemble

	Out=list("DistM"=DistM,"Clust"=Clust)
	attr(Out,"method")="Ensemble"

	return(Out)
}


#' @title Ensemble hierarchical clustering via graph partitioning (METIS / MST)
#'
#' @description Ensemble Hierarchical Clustering (EHC) builds a cluster
#' association graph from the agglomerative clustering results of each data
#' source. Strengths of association between pairs of objects are accumulated
#' across the dendrograms (weighted by cluster depth and intra-cluster
#' proximity), yielding a weighted adjacency matrix that is then partitioned.
#'
#' Two partitioning strategies are available through \code{graphPartitioning}:
#' \describe{
#'   \item{"METIS"}{Delegates the partitioning of the association graph to an
#'     \strong{external graph-partitioning executable} (METIS / \code{gpmetis} /
#'     \code{shmetis}, or the corresponding MATLAB \code{hmetis} routine). These
#'     binaries are \strong{not bundled} with the package and must be configured
#'     through the \code{executable} argument. When \code{executable=FALSE}
#'     (the default), or when the configured executable cannot be located, the
#'     function stops with an informative error.}
#'   \item{"MST"}{A pure-R minimum spanning tree partitioning that needs no
#'     external executable. It does require the \pkg{igraph} package; if
#'     \pkg{igraph} is not installed the function stops with an informative
#'     error.}
#' }
#'
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clust" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clust".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param graphPartitioning A character string indicating the graph partitioning strategy. One of "METIS" (external executable) or "MST" (pure R, requires igraph). Defaults to "METIS".
#' @param optimalk An estimate of the final optimal number of clusters. Default is 7.
#' @param waitingtime The maximum number of seconds to wait for the external program to write its result before stopping. Defaults to 300.
#' @param file_number An identifier appended to the temporary files exchanged with the external program. Defaults to 00.
#' @param executable Either FALSE (the default) or the path/flag identifying the external graph-partitioning executable (METIS / gpmetis / shmetis) or MATLAB to use for the "METIS" strategy. When FALSE the "METIS" path stops, since no partition can be computed in pure R. Precompiled binaries are not bundled with the package.
#' @return The returned value is a list of two elements:
#' \item{DistM}{The weighted cluster association matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references Mirzaei A. and Rahmati M. (2010). A novel hierarchical-clustering-combination scheme based on fuzzy-similarity relations. IEEE Transactions on Fuzzy Systems, 18(1), 27-39.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#'
#' # METIS path requires an external graph-partitioning executable / MATLAB
#' # configured through the 'executable' argument.
#' res <- EHC(List = L, type = "data",
#'   distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
#'   method = c(NULL, NULL), clust = "agnes",
#'   linkage = c("flexible", "flexible"), alpha = 0.625, gap = FALSE,
#'   maxK = 15, graphPartitioning = "METIS", optimalk = 7,
#'   executable = "./MetisAlgorithm")
#'
#' # MST path is pure R but needs the igraph package.
#' res2 <- EHC(List = L, type = "data", graphPartitioning = "MST")
#' }
EHC<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,gap = FALSE, maxK = 15,graphPartitioning=c("METIS","MST"),optimalk=7,waitingtime=300,file_number=00,executable=FALSE){

	graphPartitioning <- match.arg(graphPartitioning)

	# Guard the executable-dependent METIS path: refuse to run with a clear
	# message when no usable external graph-partitioning executable / MATLAB is
	# configured. Such binaries are not bundled with the package.
	if(graphPartitioning=="METIS"){
		if(isFALSE(executable) || is.null(executable)){
			stop(paste0("EHC with graphPartitioning='METIS' requires an external graph-partitioning ",
				"executable (METIS / gpmetis / shmetis) or MATLAB to partition the association ",
				"graph. These binaries are not bundled with the package. Configure one through ",
				"the 'executable' argument (it is currently FALSE), or use graphPartitioning='MST' ",
				"for the pure-R minimum spanning tree alternative."))
		}
		if(is.character(executable) && !nzchar(Sys.which(executable)) && !file.exists(executable)){
			stop(paste0("EHC could not locate the configured graph-partitioning executable '",
				executable, "'. Provide a valid path to a METIS / gpmetis / shmetis binary or to ",
				"MATLAB. Precompiled binaries are not bundled with the package."))
		}
	}

	# The cluster association graph is built and (for MST) partitioned with
	# igraph, which is an optional dependency. Fail clearly if it is missing.
	if(!requireNamespace("igraph", quietly=TRUE)){
		stop(paste0("EHC requires the 'igraph' package to build (and, for graphPartitioning='MST', ",
			"partition) the cluster association graph. Please install it with ",
			"install.packages('igraph')."))
	}

	## Step 1: perfom aggl clustering on each data set

	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

	}

	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

	}
	else if(type=="clust"){

		Clusterings=List

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}


	}

	## Step 2: compute strengths of association
	theta_ab=matrix(0,nrow=nrow(Clusterings[[1]]$DistM),ncol=nrow(Clusterings[[1]]$DistM))
	rownames(theta_ab)=rownames(Clusterings[[1]]$DistM)
	colnames(theta_ab)=rownames(Clusterings[[1]]$DistM)

	for(i in 1:length(Clusterings)){
		Dend=Clusterings[[i]]$Clust
		nclusters=nrow(Dend$merge)
		MergingList=list()

		for(a in 1:nclusters){
			NewMerge=Dend$merge[a,]

			if(sign(NewMerge[1])==-1 & sign(NewMerge[2])==-1){
				MergingList[[length(MergingList)+1]]=abs(NewMerge)
			}

			else if(sign(NewMerge[1])==1 & sign(NewMerge[2])==-1){
				MergingList[[length(MergingList)+1]]=c(MergingList[[NewMerge[1]]],abs(NewMerge[2]))
			}

			else if(sign(NewMerge[1])==-1 & sign(NewMerge[2])==1){
				MergingList[[length(MergingList)+1]]=c(abs(NewMerge[1]),MergingList[[NewMerge[2]]])
			}

			else if(sign(NewMerge[1])==1 & sign(NewMerge[2])==1){
				MergingList[[length(MergingList)+1]]=c(MergingList[[NewMerge[1]]],MergingList[[NewMerge[2]]])
			}
		}

		#depths
		depths=c()
		for(b in 1:length(MergingList)){
			SubList=MergingList[[b]]
			depth=0
			for(c in 1:length(MergingList)){
				if(all(SubList%in%MergingList[[c]])){
					depth=depth+1
				}
			}
			depths=c(depths,depth)
		}

		depths=depths-1 #root cluster has depth 0

		#max_depth
		max_depth=max(depths)

		#intra_cluster proximity values: based on cd ; compare with heights
		values=sort(Dend$height)

		#values=as.matrix(cophenetic(as.hclust(Clusterings[[i]]$Clust)))

		CheckDist<-function(Dist,StopRange){
			if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
				#message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				Dist=Normalization(Dist,method="Range")
			}
			else{
				Dist=Dist
			}
		}

		values=CheckDist(values,StopRange=FALSE)  #use 1-height as the intra-cluster proximity values after normalizing height between 0 and 1
		proximityvalues=1-values  #one for each cluster

		for(j in 1:length(MergingList)){
			Pairs=utils::combn(MergingList[[j]],2)
			for(k in 1:ncol(Pairs)){
				ab=Pairs[,k]
				theta_ab[ab[1],ab[2]]=theta_ab[ab[1],ab[2]]+depths[j]*proximityvalues[j]/max_depth
				theta_ab[ab[2],ab[1]]=theta_ab[ab[1],ab[2]]
			}
		}
	}



	## Step 3: generate cluster association graph
	## theta_ab can be seen as an adjacency matrix of the graph: the higher the value, the closer the objects
	net=igraph::graph_from_adjacency_matrix(theta_ab,mode="undirected",weighted=TRUE,diag=FALSE)
	igraph::plot.igraph(net,vertex.label=igraph::V(net)$name,layout=igraph::layout_with_fr, edge.color="black",edge.width=igraph::E(net)$weight)


	if(graphPartitioning=="METIS"){

		matlabdata=theta_ab
		utils::write.table(matlabdata,file=paste("matlabdata_",file_number,".csv",sep=""),sep=",",col.names=FALSE,row.names=FALSE)

		system(paste(executable," ",optimalk," ",file_number,sep=""),intern=TRUE)

		Continue=FALSE
		time=0
		while(Continue==FALSE){
			Sys.sleep(15)
			time=time+15
			Continue=file.exists(paste("Partition_",file_number,".csv",sep=""))
			if(time>waitingtime & Continue==FALSE){
				stop(paste("Waited",waitingtime, "seconds for completion of the ensemble clustering procedure. Increase waiting time to continue.",sep=" "))
			}
		}


		Partition <- utils::read.table(paste("Partition_",file_number,".csv",sep=""),sep=",")
		Partition=as.vector(as.matrix(Partition))
		names(Partition)=rownames(theta_ab)

		order=sort(Partition,index=TRUE)$ix
		order.lab=names(Partition[sort(Partition,index=TRUE)$ix])
		Clusters=Partition

		file.remove(paste("Partition_",file_number,".csv",sep=""))
		file.remove(paste("matlabdata_",file_number,".csv",sep=""))

	}

	else if(graphPartitioning=="MST"){

		#are assigning to the first encounter to break ties and get clusters: if changed  to join all if one is in common, all end up in 1 cluster

		Graph=igraph::graph_from_adjacency_matrix(adjmatrix=theta_ab, mode=c( "undirected"), weighted=TRUE, diag=TRUE,add.colnames=NULL)
		MST_Graph=igraph::as_data_frame(igraph::mst(Graph))

		Partition=list()
		Partition[[1]]=c(MST_Graph[1,1],MST_Graph[1,2])
		Placed=rep(FALSE,nrow(theta_ab))
		Placed[Partition[[1]]]=TRUE
		for(j in 2:nrow(MST_Graph)){
			Edge=as.numeric(MST_Graph[j,c(1:2)])
			k=1

			while(any(Placed[Edge]==FALSE)){

				if(k>length(Partition)){
					Partition[[length(Partition)+1]]=Edge
					Placed[Edge]=TRUE
				}

				else if(Edge[1]%in%Partition[[k]] | Edge[2]%in%Partition[[k]]){
					Partition[[k]]=c(unlist(Partition[[k]]),Edge)
					Partition[[k]]=unique(Partition[[k]])
					Placed[Edge]=TRUE
				}
				k=k+1
			}


		}

		order=unlist(Partition)
		order.lab=rownames(theta_ab)[order]
		t1<-sapply(1:length(Partition),function(i) rep(i,length(Partition[[i]])))
		t2<-sapply(1:length(Partition),function(i) rownames(theta_ab)[Partition[[i]]])
		t3<-cbind(unlist(t1),unlist(t2))
		clus=as.numeric(t3[,1])
		names(clus)=t3[,2]
		Clusters=clus[rownames(theta_ab)]


	}

	Out=list(DistM=theta_ab,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)

}
