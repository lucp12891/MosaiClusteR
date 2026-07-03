
#' @title Hybrid Bipartite Graph Formulation
#'
#' @description The \code{HBGF} function implements the Hybrid Bipartite Graph
#' Formulation, a graph-based consensus clustering method. A bipartite graph is
#' constructed from the cluster ensemble (objects on one side, clusters on the
#' other) and partitioned with a spectral method based on the singular value
#' decomposition followed by k-means.
#'
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type Indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clust" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clust".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param graphPartitioning The graph partitioning method. Currently "Spec" (spectral) is supported.
#' @param optimalk An estimate of the final optimal number of clusters. Default is 7.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references Fern, X. Z. and Brodley, C. E. (2004). Solving cluster ensemble
#' problems by bipartite graph partitioning. Proceedings of the Twenty-First
#' International Conference on Machine Learning (ICML).
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' HBGF_toy <- HBGF(List = L, type = "data",
#'                  distmeasure = c("euclidean", "euclidean"),
#'                  normalize = c(FALSE, FALSE), method = c(NULL, NULL),
#'                  clust = "agnes", linkage = c("ward", "ward"),
#'                  nrclusters = c(7, 7), graphPartitioning = "Spec", optimalk = 7)
#' }
HBGF<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,graphPartitioning="Spec",optimalk=7){
	# Step 1: Generate several clustering results on the objects

	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
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

		Clusters=lapply(seq(length(Clusterings)),function(i) stats::cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
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

	# Step 2: construct a graph G from the cluster ensemble: based on matlab code from Brodley and Fern (http://web.engr.oregonstate.edu/~xfern/)

	ClusterEnsembles=as.data.frame(Clusters)
	rownames(ClusterEnsembles)=rownames(List[[1]])
	for(i in 1:ncol(ClusterEnsembles)){
		colnames(ClusterEnsembles)[i]=paste("Solution_",i,sep="")
	}


	A=lapply(seq(length(Clusters)),function(i) stats::model.matrix(~as.factor(Clusters[[i]])-1))
	A=Reduce(cbind,A)
	rownames(A)=rownames(List[[1]])
	for(c in 1:ncol(A)){
		colnames(A)[c]=paste("Cluster ",c,sep="")
	}



	#Need W or just the connectivity matrix A? Ferns Algorthym uses A
	if(graphPartitioning=="Spec"){
		D<-diag(sqrt(colSums(A)))
		L<-A%*%solve(D)
		SingularValues=svd(L,nu=optimalk,nv=optimalk)
		U=SingularValues$u
		V=SingularValues$v

		U=U/matrix(rep(sqrt(rowSums(U^2)),optimalk),nrow=nrow(A),ncol=optimalk,byrow=FALSE)
		V=V/matrix(rep(sqrt(rowSums(V^2)),optimalk),nrow=ncol(A),ncol=optimalk,byrow=FALSE)

		permutated=gtools::permute(1:nrow(A))
		centers=U[permutated[1],,drop=FALSE]

		c=rep(0,nrow(A))
		c[permutated[1]]=2*optimalk

		#finsih this last for loop

		for(j in 2:optimalk){

			c=c+abs(U%*%t(centers[j-1,,drop=FALSE]))
			m = which.min(c)
			centers = rbind(centers,U[m,])
			c[m] = 2*optimalk

		}

		#k-means clustering

		clusterid=try(stats::kmeans(x=rbind(U,V),centers=centers,iter.max=200),silent=TRUE)
		if(inherits(clusterid,"try-error")){
			clusterid=try(stats::kmeans(x=rbind(U,V),centers=length(centers),iter.max=200),silent=TRUE)
		}
		Clusters=clusterid$cluster[1:nrow(A)]
		names(Clusters)=rownames(A)

		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}

		order.lab=as.character(order)

	}

	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)

}


#' @title Clustering aggregation
#'
#' @description The \code{ClusteringAggregation} includes the ensemble clustering methods Balls, Agglomerative (Aggl.) and Furthest which are graph-based consensus methods.
#'
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clust" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clust".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7) for two data sets.
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Default is FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param agglMethod The method to be performed: "Balls","Aggl","Furthest" or "LocalSearch".
#' @param improve Logical. If TRUE, a local search is performed to improve the obtained results. Default is TRUE.
#' @param distThresh_B  A distance threshold for the Balls algoritme. Default is 0.5.
#' @param distThresh_A A distance threshold for the Aggl. algoritme. Default is 0.8.
#' @details Gionis, Mannila and Tsaparas (2007) propose heuristic algorithms in order to find a solution for the consensus problem. In a first step, a
#' weighted graph is built from the objects with weights between two vertices determined by the fraction of clusterings that place the two vertices
#' in different clusters. In a second step, an algorithm searches for the partition that minimizes the total number of disagreements with the given
#' partitions. The Balls algorithm is an iterative process which finds a cluster for the consensus partition in each iteration. For each object,
#' all objects at a distance of at most 0.5 are collected and the average distance of this set to the object of interest is calculated. If the
#' average distance is less or equal to a parameter the objects form a cluster; otherwise the object forms a singleton. The Agglomerative
#' (Aggl.) algorithm starts by considering every object as a singleton cluster. Next, the two closest clusters are merged if the average distance
#' between the clusters is less than 0.5. If there are no two clusters with an average distance smaller than 0.5, the algorithm stops and returns
#' the created clusters as a solution. The Furthest algorithm starts by placing all objects into a single cluster. In each iteration, the pair of
#' objects that are the furthest apart are considered as new cluster centers. The remaining objects are appointed to the center that increases the
#' cost of the partition the least and the new cost is computed. The cost is the sum of the all distances between the obtained partition and the
#' partitions in the ensemble. The iteration continues until the cost of the new partition is higher than the previous partition.
#' @return The returned value is a list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble'.
#' @references Gionis, A., Mannila, H. and Tsaparas, P. (2007). Clustering
#' aggregation. ACM Transactions on Knowledge Discovery from Data, 1(1), 4-es.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' Aggl_toy <- ClusteringAggregation(List = L, type = "data",
#'   distmeasure = c("euclidean", "euclidean"), normalize = c(FALSE, FALSE),
#'   method = c(NULL, NULL), clust = "agnes", linkage = c("ward", "ward"),
#'   alpha = 0.625, nrclusters = c(7, 7), agglMethod = "Aggl",
#'   improve = TRUE, distThresh_B = 0.5, distThresh_A = 0.8)
#' }
ClusteringAggregation<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,agglMethod=c("Balls","Aggl","Furthest","LocalSearch"),improve=TRUE,distThresh_B=0.5,distThresh_A=0.8){
	# Step 1: Generate several clustering results on the objects

	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
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

		Clusters=lapply(seq(length(Clusterings)),function(i) stats::cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
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

	# Step 2: Set up the distance matrix X
	H=list()
	for(i in 1:length(List)){
		Bin=matrix(0,nrow=max(Clusters[[i]]),ncol=length(Clusters[[i]]))
		for(j in 1:nrow(Bin)){
			Bin[j,which(Clusters[[i]]==j)]=1
		}
		H[[i]]=as.data.frame(Bin)
	}
	T=do.call(rbind,H)
	S=t(as.matrix(T))%*%(as.matrix(T))
	S=S/length(List)
	X=1-S
	if(type!="clust"){
		rownames(X)=rownames(List[[1]])
		colnames(X)=rownames(List[[1]])
	}
	else{
		rownames(X)=rownames(List[[1]]$DistM)
		colnames(X)=rownames(List[[1]]$DistM)
	}

	# Step 3: The specified algorithm on the distance matrix X

	if(agglMethod=="Balls"){
		AggClusters=list()
		SortedVertices=sort(rowSums(X))
		while(length(SortedVertices)>0){
			u=SortedVertices[1]
			B=which(SortedVertices[-c(1)]<=u+0.5)
			if(length(B)>0){
				AvDistBtoU=mean(X[which(rownames(X)==names(u)),which(colnames(X)%in%names(B))])

				if(AvDistBtoU<=distThresh_B){
					AggClusters[[length(AggClusters)+1]]=c(names(u),names(B))
					SortedVertices=SortedVertices[-which(names(SortedVertices)%in%c(names(u),names(B)))]
				}
				else{
					AggClusters[[length(AggClusters)+1]]=names(u)
					SortedVertices=SortedVertices[-which(names(SortedVertices)%in%c(names(u)))]
				}
			}
			else{
				AggClusters[[length(AggClusters)+1]]=names(u)
				SortedVertices=SortedVertices[-which(names(SortedVertices)%in%c(names(u)))]
			}
		}

		nrclusters=length(AggClusters)
		order.lab=unlist(AggClusters)

		if(type!="clust"){
			order=match(order.lab,rownames(List[[1]]))
		}
		else{
			order=match(order.lab,rownames(List[[1]]$DistM))
		}

		clusters=sapply(1:nrclusters, function(i) rep(i,length(AggClusters[[i]])))
		Clusters=unlist(clusters)
		names(Clusters)=order.lab

		if(type!="clust"){
			Clusters=Clusters[rownames(List[[1]])]
		}
		else{
			Clusters=Clusters[rownames(List[[1]]$DistM)]
		}

		if(improve==TRUE){agglMethod="LocalSearch"}
	}

	if(agglMethod=="Aggl"){
		AvDist<-function(DistMat, pairs,temp_clusters){
			Cluster1=names(temp_clusters)[which(temp_clusters==pairs[1])]
			Cluster2=names(temp_clusters)[which(temp_clusters==pairs[2])]
			AvDist1to2=mean(DistMat[which(rownames(DistMat)%in%Cluster2),which(colnames(DistMat)%in%Cluster1)])

			return(AvDist1to2)

		}

		temp_clusters=seq(1:length(rownames(X)))
		names(temp_clusters)=rownames(X)

		Continue=TRUE
		while(Continue==TRUE){

			Pairs=utils::combn(unique(temp_clusters),2)
			AvDistances=apply(Pairs,2,AvDist,DistMat=X,temp_clusters=temp_clusters)

			if(min(AvDistances)<distThresh_A){
				ChosenPair=Pairs[,which.min(AvDistances)]
				temp_clusters[which(temp_clusters==ChosenPair[2])]=ChosenPair[1]
			}
			else{
				Continue=FALSE
			}

			if(length(unique(temp_clusters))==1){
				Continue=FALSE # all have been put into the same cluster
			}

		}

		Clusters=as.numeric(as.factor(temp_clusters))
		names(Clusters)=names(temp_clusters)

		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}

		order.lab=as.character(order)
		if(improve==TRUE){agglMethod="LocalSearch"}
	}

	if(agglMethod=="Furthest"){

		clusters=rep(0,length(rownames(X)))
		names(clusters)=rownames(X)

		centers=c(rownames(X)[which(X==max(X),arr.ind = TRUE)[1,1]],colnames(X)[which(X==max(X),arr.ind = TRUE)[1,2]])
		for(i in 1:length(centers)){
			clusters[which(names(clusters)%in%centers[i])]=i
		}

		for(j in names(clusters)[-c(which(names(clusters)%in%centers))]){
			DistancesToCenters=X[which(rownames(X)==j),which(colnames(X)%in%centers)]
			AssignedCenter=names(DistancesToCenters)[which.min(DistancesToCenters)]
			clusters[which(names(clusters)==j)]=clusters[which(names(clusters)==AssignedCenter)]

		}

		Part_1=0
		Part_2=0

		Pairs=utils::combn(length(rownames(X)),2)
		for(k in 1:ncol(Pairs)){
			Pair=rownames(X)[Pairs[,k]]

			if(clusters[which(names(clusters)%in%Pair[1])]==clusters[which(names(clusters)%in%Pair[2])]){
				Part_1=Part_1+X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])]
			}
			else{
				Part_2=Part_2+(1-X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])])
			}

		}

		Cost_clusters=Part_1+Part_2
		Solution=clusters

		Continue=TRUE

		while(Continue==TRUE){

			AvDistToCenters=c()
			for(i in names(clusters)[-c(which(names(clusters)%in%centers))]){
				DistancesToCenters=X[which(rownames(X)==i),which(colnames(X)%in%centers)]
				AvDistToCenters=c(AvDistToCenters,mean(DistancesToCenters))
			}
			names(AvDistToCenters)=names(clusters)[-c(which(names(clusters)%in%centers))]

			new_center=names(AvDistToCenters)[which.max(AvDistToCenters)]

			centers<-c(centers,new_center)

			for(i in 1:length(centers)){
				clusters[which(names(clusters)%in%centers[i])]=i
			}

			for(j in names(clusters)[-c(which(names(clusters)%in%centers))]){
				DistancesToCenters=X[which(rownames(X)==j),which(colnames(X)%in%centers)]
				AssignedCenter=names(DistancesToCenters)[which.min(DistancesToCenters)]
				clusters[which(names(clusters)==j)]=clusters[which(names(clusters)==AssignedCenter)]
			}

			Part_1=0
			Part_2=0

			Pairs=utils::combn(length(rownames(X)),2)
			for(k in 1:ncol(Pairs)){
				Pair=rownames(X)[Pairs[,k]]

				if(clusters[which(names(clusters)%in%Pair[1])]==clusters[which(names(clusters)%in%Pair[2])]){
					Part_1=Part_1+X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])]
				}
				else{
					Part_2=Part_2+(1-X[which(rownames(X)%in%Pair[1]),which(colnames(X)%in%Pair[2])])
				}

			}

			Cost_clusters_new=Part_1+Part_2


			if(Cost_clusters_new<Cost_clusters){
				Cost_clusters=Cost_clusters_new
				Solution=clusters
			}
			else{
				Continue=FALSE
				clusters=Solution
			}


		}

		Clusters=clusters
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}

		order.lab=as.character(order)

		if(improve==TRUE){
			agglMethod="LocalSearch"
		}
	}

	if(agglMethod=="LocalSearch"){
		UpdateClustering=Clusters

		Stay=rep(0,length(UpdateClustering))
		names(Stay)=names(UpdateClustering)
		Moved=c()

		while(length(which(Stay=="Stay"))!=length(UpdateClustering)){  #continue untill all nodes want to stay where they are

			for(i in names(UpdateClustering)){   #Pick a node

				PresentCluster=UpdateClustering[which(names(UpdateClustering)==i)]  # Present cluster of the node
				ClusterSizes=table(UpdateClustering)  # Current cluster sizes

				Cost_Cj=c()   #Compute the cost of moving node i to any cluster
				for(j in sort(unique(UpdateClustering))){
					Cj=names(UpdateClustering)[which(UpdateClustering==j)]
					M_v_Cj=sum(X[which(rownames(X)==i),which(colnames(X)%in%Cj)])
					Cost_Cj=c(Cost_Cj,M_v_Cj)
				}
				names(Cost_Cj)=sort(unique(UpdateClustering))

				Dist_k=c()  #Compute the cost of moving node i to any cluster or a singleton
				for(k in sort(unique(UpdateClustering))){
					ClusterToMoveTo=k
					Cost_Singl=sum(ClusterSizes[which(names(ClusterSizes)!=ClusterToMoveTo)]-Cost_Cj[which(names(Cost_Cj)!=ClusterToMoveTo)])
					Dist_k=c(Dist_k,Cost_Cj[which(names(ClusterSizes)==ClusterToMoveTo)]+Cost_Singl)
				}
				PotentialClusters=which(Dist_k==min(Dist_k)) # To which clusters does the move hae minimal cost


				if(PresentCluster%in%PotentialClusters){
					Stay[which(names(UpdateClustering)==i)]="Stay"	#If present cluster is among the potential ones: stay
				}
				else{
					UpdateClustering[which(names(UpdateClustering)==i)]=PotentialClusters[1]  #else move and if mover before stay is changed to 0 again
					if(Stay[which(names(UpdateClustering)==i)]=="Stay"){
						Stay[which(names(UpdateClustering)==i)]=0
					}
					Moved=c(Moved,i)
					Moved=unique(Moved)
				}
			}

		}

		Clusters=UpdateClustering

		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}

		order.lab=as.character(order)

	}

	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)
}
