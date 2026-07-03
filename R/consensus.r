#' @title Voting-based consensus clustering
#'
#' @description Consensus Clustering performs voting-based consensus on a set of
#' single-source partitions. Each data matrix is clustered separately and the
#' resulting label matrix is reconciled into a single consensus partition by an
#' iterative voting scheme. Three voting methods are available: Iterative Voting
#' Consensus (IVC), Iterative Probabilistic Voting Consensus (IPVC) and Iterative
#' Pairwise Consensus (IPC).
#'
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type Indicates whether the provided matrices in "List" are either data matrices ("data"),
#' distance matrices ("dist") or clustering results obtained from the data ("clust").
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param votingMethod The voting method to be performed: one of "IVC", "IPVC" or "IPC".
#' @param optimalk The number of clusters for the final consensus partition. Defaults to 7.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references Nguyen N. and Caruana R. (2007). Consensus Clusterings. Seventh IEEE International Conference on Data Mining (ICDM 2007), 607-612.
#' @export
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' fit <- ConsensusClustering(List = mosaic_toy$List, type = "data",
#'   distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
#'   method = c(NULL, NULL), clust = "agnes", linkage = c("ward", "ward"),
#'   nrclusters = c(3, 3), gap = FALSE, votingMethod = "IVC", optimalk = 3)
#' }
ConsensusClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,votingMethod=c("IVC","IPVC","IPC"),optimalk=7){

	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,StopRange=TRUE))

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

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,StopRange=TRUE))

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


	Y=Reduce(cbind,Clusters)
	rownames(Y)=rownames(List[[1]])
	colnames(Y)=seq(1,ncol(Y))


	if(votingMethod=="IVC"){
		target=sample(optimalk,nrow(Y),replace=TRUE)
		update=rep(0,nrow(Y))
		Continue=TRUE
		while(Continue){

			#representation of each cluster
			P=list()
			y_P=list()
			for(i in unique(target)){
				P_i=which(target==i)
				y_Pi=apply(Y,2,function(k) as.numeric(names(which.max(table(k[P_i])))))
				P[[i]]=P_i
				y_P[[i]]=y_Pi

			}

			#update y
			for(y in 1:length(target)){
				distances=c()
				for(a in 1:length(y_P)){
					distances=c(distances,sum(Y[y,] != y_P[[a]])	)
				}

				update[y]=which.min(distances)

			}
			if(all(update==target)){
				Continue=FALSE
				#print("here")
			}
			else{
				target=update
			}

		}

	}
	else if(votingMethod=="IPVC"){
		target=sample(optimalk,nrow(Y),replace=TRUE)
		update=rep(0,nrow(Y))
		Continue=TRUE

		while(Continue){

			#representation of each cluster
			P=list()
			n_P=list()
			for(i in unique(target)){
				P_i=which(target==i)
				P[[i]]=P_i
				n_P[[i]]=length(P_i)

			}

			#update y
			for(y in 1:length(target)){
				distances=c()

				for(a in 1:length(P)){
					dist=0
					for(b in 1:ncol(Y)){
						dist=dist+sum(Y[y,b] != Y[P[[a]],b])/n_P[[a]]
					}
					distances=c(distances,dist)
				}

				update[y]=which.min(distances)

			}
			if(all(update==target)){
				Continue=FALSE
				#print("here")
			}
			else{
				target=update
			}

		}



	}
	else if(votingMethod=="IPC"){


		#similarity matrix
		H=list()
		for(i in 1:length(List)){
			Bin=matrix(0,nrow=max(Clusters[[i]]),ncol=length(Clusters[[i]]))
			for(j in 1:nrow(Bin)){
				Bin[j,which(Clusters[[i]]==j)]=1
			}
			H[[i]]=as.data.frame(Bin)
		}
		T=data.table::rbindlist(H)
		S=t(as.matrix(T))%*%(as.matrix(T))
		S=S/length(List)

		if(type!="clust"){
			rownames(S)=rownames(List[[1]])
			colnames(S)=rownames(List[[1]])
		}
		else{
			rownames(S)=rownames(List[[1]]$DistM)
			colnames(S)=rownames(List[[1]]$DistM)
		}

		target=sample(optimalk,nrow(Y),replace=TRUE)
		update=rep(0,nrow(Y))
		Continue=TRUE

		while(Continue){

			#representation of each cluster
			P=list()
			n_P=list()
			for(i in unique(target)){
				P_i=which(target==i)
				P[[i]]=P_i
				n_P[[i]]=length(P_i)

			}

			#update x
			for(x in 1:length(target)){
				similarities=c()

				for(a in 1:length(P)){
					sim=sum(S[x,P[[a]]])/n_P[[a]]
					similarities=c(similarities,sim)
				}

				update[x]=which.max(similarities)
			}


			if(all(update==target)){
				Continue=FALSE
				#print("here")
			}
			else{
				target=update
			}

		}

	}



	Clusters=target

	if(type!="clust"){
		names(Clusters)=rownames(List[[1]])
	}
	else{
		names(Clusters)=rownames(List[[1]]$DistM)
	}

	order=sort(Clusters)
	order.lab=names(order)
	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)

}


#' @title Weighting on Membership clustering
#'
#' @description Weighting on Membership (WonM) is a hierarchy-based consensus
#' technique. Each data matrix is hierarchically clustered and the resulting
#' dendrogram is cut into a range of K values. For every cut a co-membership
#' matrix is built, recording which objects fall in the same cluster. These
#' matrices are averaged within and across data sets to produce an overall
#' consensus matrix, which is transformed into a dissimilarity and clustered a
#' final time.
#'
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type Indicates whether the provided matrices in "List" are either data matrices ("data"),
#' distance matrices ("dist") or clustering results obtained from the data ("clusters").
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param nrclusters A list with, for each data matrix, the sequence of numbers of clusters to cut the dendrogram in. Defaults to list(seq(5,25), seq(5,25)).
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting consensus dissimilarity matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'WonM'.
#' @references Saeed F., Salim N. and Abdo A. (2012). Voting-based consensus clustering for combining multiple clusterings of chemical structures. Journal of Cheminformatics, 4, 37.
#' @export
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' fit <- WonM(List = mosaic_toy$List, type = "data",
#'   distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
#'   method = c(NULL, NULL), nrclusters = list(seq(3, 6), seq(3, 6)),
#'   clust = "agnes", linkage = c("ward", "ward"))
#' }
WonM=function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),nrclusters=list(seq(5,25),seq(5,25)),clust="agnes",linkage=c("flexible","flexible"),alpha=0.625){

	type<-match.arg(type)


	#Step 1: Distance Matrices
	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,,drop=FALSE]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))

		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
	}


	#Step 2: perform hierarchical clustering on both distance matrices

	HClustering=lapply(seq(length(List)),function(i) cluster::agnes(Dist[[i]],diss=TRUE,method=linkage[i],par.method=alpha))


	#Step 3: cut the dendrograms into a range of K values

	#Give 0 to pair belonging together, give 1 to a pair not belonging together : ==> Distances created otherwise similarities.
	ClusterMembers<-function(HClust,nrclusters){
		Temp=lapply(seq(length(nrclusters)),function(i) stats::cutree(HClust,nrclusters[i]))
		CM=lapply(seq(length(nrclusters)),function(i) matrix(0,dim(List[[1]])[1],dim(List[[1]])[1]))

		clusters<-function(temp,cm){
			for(l in 1:length(temp)){
				label=temp[l]
				sameclust=which(temp==label)
				cm[l,sameclust]=1
			}
			return(cm)
		}

		CM2=lapply(seq(length(nrclusters)),function(i) clusters(temp=Temp[[i]],cm=CM[[i]]))
		Consensus2=Reduce("+",CM2)
		Consensus2=Consensus2/length(nrclusters)
		return(Consensus2)


	}

	Consensus=lapply(seq(length(List)), function(i) ClusterMembers(HClustering[[i]],nrclusters[[i]]))

	OverallConsensus=Reduce("+",Consensus)
	OverallConsensus=OverallConsensus/length(List)
	OverallConsensus=as.matrix(OverallConsensus)
	rownames(OverallConsensus)=rownames(Dist[[1]])
	colnames(OverallConsensus)=rownames(Dist[[1]])
	OverallConsensusD=1-OverallConsensus
	OverallClusteringR=cluster::agnes(OverallConsensusD,diss=TRUE,method="ward")

	out=list(DistM=OverallConsensusD,Clust=OverallClusteringR)
	attr(out,'method')<-'WonM'
	return(out)
}
