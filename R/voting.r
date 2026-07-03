
#' @title Cumulative Voting Aggregation (CVAA / W-CVAA)
#'
#' @description \code{CVAA} performs cumulative voting consensus clustering. A
#' reference partition is used to align the individual source partitions, which
#' are then aggregated into a consensus stochastic membership matrix. Two
#' variants are available: cumulative voting (\code{"CVAA"}) and the
#' information-weighted variant (\code{"W-CVAA"}).
#'
#' @export
#' @param Reference A reference partition (a clustering result, optionally an
#' object produced by a weighted/CEC method). Required.
#' @param nrclustersR The number of clusters of the reference partition. Default is 7.
#' @param List A list of data matrices. It is assumed the rows correspond with the objects.
#' @param typeL Indicates whether the provided matrices in "List" are either data
#' matrices ("data"), distance matrices ("dist") or clustering results ("clust").
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE).
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL).
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible","flexible").
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625.
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7).
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param votingMethod The method to be performed: "CVAA" or "W-CVAA".
#' @param optimalk An estimate of the final optimal number of clusters. Defaults to nrclustersR.
#' @return A list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering (order, order.lab and Clusters)}
#' The value has class 'Ensemble'.
#' @details Cumulative voting aggregation aligns each source membership matrix to
#' a running consensus by least-squares voting, then averages cumulatively. The
#' weighted variant orders the partitions by average information content and
#' weights their contribution accordingly. When \code{optimalk} is smaller than
#' \code{nrclustersR} the JS-ALink agglomeration is used to merge reference
#' clusters.
#' @references Ayad H.G. and Kamel M.S. (2008). Cumulative voting consensus
#' method for partitions with variable number of clusters. IEEE Transactions on
#' Pattern Analysis and Machine Intelligence, 30(1), 160-173.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' Ref <- Cluster(L[[1]], type = "data", distmeasure = "euclidean")
#' attr(Ref, "method") <- "Single"
#' res <- CVAA(Reference = Ref, nrclustersR = 5, List = L, typeL = "data",
#'   distmeasure = c("euclidean","euclidean"), votingMethod = "CVAA", optimalk = 5)
#' }
CVAA<-function(Reference=NULL,nrclustersR=7,List,typeL=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,votingMethod=c("CVAA","W-CVAA"),optimalk=nrclustersR){

	#needs a reference partition: which one to choose? Weighted? Or one of the single source clusterings? Let user decide but suggest Weighted as default
	#Reference can be a "method": List needs to be data or dist
	#if Reference is a clust, list can still be anything

	if(typeL=="data"){

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

		U_all=lapply(seq(length(Clusters)),function(i) stats::model.matrix(~as.factor(Clusters[[i]])-1))
		for(i in 1:length(U_all)){
			colnames(U_all[[i]])=seq(1,ncol(U_all[[i]]))
			rownames(U_all[[i]])=rownames(List[[i]])
		}

	}

	else if(typeL=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],typeL,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap,maxK,StopRange=TRUE))

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
		U_all=lapply(seq(length(Clusters)),function(i) stats::model.matrix(~as.factor(Clusters[[i]])-1))
		for(i in 1:length(U_all)){
			colnames(U_all[[i]])=seq(1,ncol(U_all[[i]]))
			rownames(U_all[[i]])=rownames(List[[i]])
		}


	}
	else if(typeL=="clust"){

		Clusterings=List

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}

		if(is.null(nrclusters)){
			stop("Please specify a number of clusters")
		}

		Clusters=lapply(seq(length(Clusterings)),function(i) stats::cutree(Clusterings[[i]]$Clust,k=nrclusters[i]))
		U_all=lapply(seq(length(Clusters)),function(i) stats::model.matrix(~as.factor(Clusters[[i]])-1))
		for(i in 1:length(U_all)){
			colnames(U_all[[i]])=seq(1,ncol(U_all[[i]]))
			rownames(U_all[[i]])=rownames(List[[i]]$DistM)

		}

	}

	if(!(is.null(Reference))){

		ListNew=list()
		element=0

		if(attributes(Reference)$method != "CEC" & attributes(Reference)$method != "Weighted" & attributes(Reference)$method!= "WeightedSim"){
			ResultsClust=list()
			ResultsClust[[1]]=list()
			ResultsClust[[1]][[1]]=Reference
			names(ResultsClust[[1]])[1]="Clust"
			element=element+1
			ListNew[[element]]=ResultsClust[[1]]
			#attr(ListNew[element],"method")="Weights"
		}
		else if(attributes(Reference)$method=="CEC" | attributes(Reference)$method=="Weighted" | attributes(Reference)$method == "WeightedSim"){
			ResultsClust=list()
			ResultsClust[[1]]=list()
			if(attributes(Reference)$method != "WeightedSim"){
				ResultsClust[[1]][[1]]=Reference$Clust
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1
				ListNew[[element]]=ResultsClust[[1]]
				attr(ListNew[[element]]$Clust,"method")="Weights"

			}
			else{
				for (j in 1:length(Reference$Results)){
					ResultsClust[[j]]=list()
					ResultsClust[[j]][[1]]=Reference$Results[[j]]
					names(ResultsClust[[j]])[1]="Clust"
					attr(ResultsClust[[1]],"method")="Weights"
					element=element+1
					ListNew[[element]]=ResultsClust[[j]]
					ListNew=unlist(ListNew,recursive=FALSE)
				}
			}
		}

		Reference=ListNew

		if(is.null(nrclustersR)){
			stop("Please specify a number of clusters for the reference partition")
		}

		CutTree<-function(Data,nrclusters){
			if(attributes(Data$Clust)$method == "Ensemble"){
				Clusters=Data$Clust$Clust$Clusters
				names(Clusters)=NULL
			}
			else{
				Clusters=stats::cutree(Data$Clust$Clust,k=nrclusters)
			}
			return(Clusters)
		}

		Ref_Clusters=CutTree(Reference[[1]],nrclusters=nrclustersR)

		U_0=stats::model.matrix(~as.factor(Ref_Clusters)-1)
		colnames(U_0)=seq(1,ncol(U_0))  #assign reference to U_0
		rownames(U_0)=rownames(Reference[[1]]$Clust$DistM)
		U_Ref=U_0

	}
	else{
		stop("No Reference specified")
	}

	if(votingMethod=="CVAA"){

		if(is.null(Reference)){
			random_part<-sample(1:length(U_all),size=1)
			U_0=U_all[[random_part]]
			U_all=U_all[-random_part]
		}

		for(i in 1:length(U_all)){ #for i in 2 to b to:

			W_i=solve((t(U_all[[i]])%*%U_all[[i]]))%*%t(U_all[[i]])%*%U_0

			V_i=U_all[[i]]%*%W_i

			U_0=((i-1)/i)*U_0+(1/i)*V_i
		}
	}

	else if(votingMethod=="W-CVAA"){

		#Preparation for the Weights for the partitions
		H_c=c()
		for(i in 1:length(U_all)){
			P_c=(colSums(U_all[[i]]))/nrow(U_all[[i]])
			H_c=c(H_c,-sum(P_c*log(P_c)))
		}

		#message("Reordering the clusterings in decreasing order of average amount of information")
		Order=sort(H_c,decreasing=TRUE,index.return=TRUE)$ix
		U_all=U_all[Order]

		if(is.null(Reference)){
			#message("Performing the Ada-cVote algorithm by Ayad and Kamel")
			U_0=U_all[Order[1]]
			U_all=U_all[Order[-c(1)]]
		}

		for(i in 1:length(U_all)){ #for i in 2 to b to:

			# The weight for partition i
			T_i=(H_c[i])/sum(H_c)


			W_i=solve((t(U_all[[i]])%*%U_all[[i]]))%*%t(U_all[[i]])%*%U_0

			V_i=U_all[[i]]%*%W_i

			U_0=((i-1)/i)*U_0+(T_i/i)*V_i
		}


	}

	U_consensus=U_0

	if(!is.null(optimalk) & optimalk==nrclustersR){

		Clusters=apply(U_consensus,1,function(i) which.max(i))
		names(Clusters)=rownames(U_consensus)
		clusters=unique(Clusters)
		order=c()
		for(j in clusters){
			order=c(order,which(Clusters==j))
		}

		order.lab=as.character(order)
	}

	else if(!is.null(optimalk)){

		p_hat_joint_c_x=(1/nrow(U_consensus))*U_consensus
		p_hat_x=1/nrow(U_consensus)

		#The JS-Alink Algorithm (applied to get the best partition out of the consensus partion U_consensus).
		#This does not do anything if the optimal number of clusters is equal to the number of clusters of the reference. It will when cluster of the reference are to be merged (less clusters than the reference)

		Pairs=utils::combn(ncol(U_consensus),2)

		JS=c()
		JS_ALink<-function(Pair,U_Ref,p_hat_joint_c_x){
			p_hat_cl=colSums(U_Ref)[Pair[1]]/nrow(U_Ref)    #the number of objects assigned to cluster l in the reference divided by the total
			p_hat_cm=colSums(U_Ref)[Pair[2]]/nrow(U_Ref)  	#the number of objects assigned to cluster m in the reference divided by the total

			p_hat_x_cl=p_hat_joint_c_x[,Pair[1]]/p_hat_cl
			p_hat_x_cm=p_hat_joint_c_x[,Pair[2]]/p_hat_cm

			alpha_1=p_hat_cl/(p_hat_cl+p_hat_cm)
			alpha_2=p_hat_cm/(p_hat_cl+p_hat_cm)


			temp1=alpha_1*p_hat_x_cl+alpha_2*p_hat_x_cm
			temp2=temp1*log(temp1)
			temp2[which(is.na(temp2))]=0
			part1=-sum(temp2)

			temp3=p_hat_x_cl*log(p_hat_x_cl)
			temp3[which(is.na(temp3))]=0
			part2=alpha_1*(-sum(temp3))

			temp4=p_hat_x_cm*log(p_hat_x_cm)
			temp4[which(is.na(temp4))]=0
			part3=alpha_2*(-sum(temp4))


			JS_cl_cm=part1-part2-part3
			JS=c(JS,JS_cl_cm)

		}
		JS=apply(Pairs,2,function(i) JS_ALink(Pair=i,U_Ref,p_hat_joint_c_x))

		Clust=cluster::agnes(JS,diss=TRUE,method="average")

		if(is.null(optimalk)){
			k_part=stats::cutree(stats::as.hclust(Clust),k=optimalk)

			priors=c()
			jointdistr=c()
			for(a in 1:optimalk){
				S_g=which(k_part==a)
				prior=sum(colSums(U_Ref)[S_g]/nrow(U_Ref))
				priors=c(priors,prior)

				joint_temp=rowSums(p_hat_joint_c_x[,S_g,drop=FALSE]/(colSums(U_Ref)[S_g,drop=FALSE]/nrow(U_Ref)))*colSums(U_Ref)[S_g,drop=FALSE]/nrow(U_Ref)
				jointdistr=cbind(jointdistr,joint_temp)

			}
			colnames(jointdistr)=c(1:ncol(jointdistr))

			U_hat=jointdistr/p_hat_x
			Clusters=apply(U_consensus,1,function(i) which.max(i))
			names(Clusters)=rownames(U_consensus)

			clusters=unique(Clusters)
			order=c()
			for(j in clusters){
				order=c(order,which(Clusters==j))
			}

			order.lab=as.character(order)

		}
	}
	else if(is.null(optimalk)){
		#message("No optimal number of clusters was specified. The consenus matrix will be returned but an optimal partition was not extracted by the JS-ALink algorithm")
		Out=U_consensus
		return(Out)

	}

	Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)

}


#' @title Evidence accumulation clustering (co-association)
#'
#' @description \code{EvidenceAccumulation} builds a co-association
#' (similarity) matrix from the individual source partitions and extracts a
#' consensus partition by graph partitioning. Three graph-partitioning options
#' are available: a minimum-spanning-tree procedure (\code{"MST"}, requires the
#' \pkg{igraph} package), a thresholded single-link procedure (\code{"SL"}), and
#' an agnes single-link procedure (\code{"SL_agnes"}).
#'
#' @export
#' @param List A list of data matrices. It is assumed the rows correspond with the objects.
#' @param type Indicates whether the provided matrices in "List" are either data
#' matrices ("data"), distance matrices ("dist") or clustering results ("clust").
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE).
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL).
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible","flexible").
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625.
#' @param nrclusters The number of clusters to divide each individual dendrogram in. Default is c(7,7).
#' @param gap Logical. Whether the optimal number of clusters should be determined with the gap statistic. Defaults to FALSE.
#' @param maxK The maximal number of clusters to investigate in the gap statistic. Default is 15.
#' @param graphPartitioning The graph partitioning method to be performed: "MST", "SL" or "SL_agnes".
#' @param t A threshold used by the "MST" and "SL" partitioning procedures. Defaults to NULL.
#' @return A list of two elements:
#' \item{DistM}{A NULL object}
#' \item{Clust}{The resulting clustering}
#' The value has class 'Ensemble' (or 'Single Clustering' for "SL_agnes").
#' @details The co-association matrix \code{S} records the fraction of source
#' partitions in which two objects are co-clustered. \code{"MST"} prunes a
#' minimum spanning tree of the similarity graph (optionally at threshold
#' \code{t}); \code{"SL"} performs a thresholded single-link grouping; and
#' \code{"SL_agnes"} runs agnes single linkage on the dissimilarity \code{1-S}.
#' The \code{"MST"} path requires the suggested package \pkg{igraph}.
#' @references Fred A.L.N. and Jain A.K. (2005). Combining multiple clusterings
#' using evidence accumulation. IEEE Transactions on Pattern Analysis and
#' Machine Intelligence, 27(6), 835-850.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res <- EvidenceAccumulation(List = L, type = "data",
#'   distmeasure = c("euclidean","euclidean"), nrclusters = c(5,5),
#'   graphPartitioning = "SL_agnes")
#' }
EvidenceAccumulation<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625,nrclusters=c(7,7),gap = FALSE, maxK = 15,graphPartitioning=c("MST","SL","SL_agnes"),t=NULL){

	#needs a reference partition: which one to choose? Weighted? Or one of the single source clusterings? Let user decide but suggest Weighted as default
	#Reference can be a "method": List needs to be data or dist
	#if Reference is a clust, list can still be anything

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

	#similarity matrix
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

	if(type!="clust"){
		rownames(S)=rownames(List[[1]])
		colnames(S)=rownames(List[[1]])
	}
	else{
		rownames(S)=rownames(List[[1]]$DistM)
		colnames(S)=rownames(List[[1]]$DistM)
	}


	if(graphPartitioning=="MST"){

		#are assigning to the first encounter to break ties and get clusters: if changed  to join all if one is in common, all end up in 1 cluster

		if (!requireNamespace("igraph", quietly = TRUE)) {
			stop("graphPartitioning = \"MST\" requires the suggested package 'igraph'. ",
				 "Install it or use graphPartitioning = \"SL\" / \"SL_agnes\".")
		}

		# Strip dimnames so igraph vertices are integer-indexed (1..n); modern
		# igraph::as_data_frame() otherwise returns edge endpoints as vertex NAMES
		# (character), which breaks the integer-based Placed[]/Partition[] logic.
		Su=S; dimnames(Su)=NULL
		Graph=igraph::graph_from_adjacency_matrix(adjmatrix=Su, mode="undirected", weighted=TRUE, diag=TRUE, add.colnames=NULL)
		MST_Graph=igraph::as_data_frame(igraph::mst(Graph))
		MST_Graph[,1]=as.integer(MST_Graph[,1]); MST_Graph[,2]=as.integer(MST_Graph[,2])

		if(!is.null(t)){
			MST_Graph=MST_Graph[-c(which(as.numeric(MST_Graph[,3])<t)),]
		}

		Partition=list()
		Partition[[1]]=c(MST_Graph[1,1],MST_Graph[1,2])
		Placed=rep(FALSE,nrow(S))
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
		order.lab=rownames(S)[order]
		t1<-sapply(1:length(Partition),function(i) rep(i,length(Partition[[i]])))
		t2<-sapply(1:length(Partition),function(i) rownames(S)[Partition[[i]]])
		t3<-cbind(unlist(t1),unlist(t2))
		clus=as.numeric(t3[,1])
		names(clus)=t3[,2]
		Clusters=clus[rownames(S)]
		if(any(is.na(Clusters))){
			Clusters[which(is.na(Clusters))]=c(1:length(which(is.na(Clusters))))
			names(Clusters)=rownames(S)
		}


	}

	else if(graphPartitioning=="SL"){
		if(is.null(t)){
			stop("Specify a treshold t for the Sl algorithm")
		}
		#better to let run for a number of thresholds t (default is 0.5)
		#Not discerning enough if only 2 datasets with 1 clustering each: need more to get more clusters

		Clusters=list()
		for(a in 1:nrow(S)){
			if(a==1){
				Clusters[[length(Clusters)+1]]=a
			}
			else{
				for(b in 1:(a-1)){

					if(S[a,b]>t){

						if(a%in%unlist(Clusters) & b%in%unlist(Clusters)){

							Clus1=which(sapply(1:length(Clusters),function(i) a%in%Clusters[[i]]))
							Clus2=which(sapply(1:length(Clusters),function(i) b%in%Clusters[[i]]))
							if(Clus1 != Clus2){
								Clusters[[Clus1]]=unique(c(unlist(Clusters[[Clus1]]),unlist(Clusters[[Clus2]])))
								Clusters[[Clus2]]=NULL
							}

						}
						else if(a%in%unlist(Clusters) & !(b%in%unlist(Clusters))){
							Clus1=which(sapply(1:length(Clusters),function(i) a%in%Clusters[[i]]))
							Clusters[[Clus1]]=unique(c(unlist(Clusters[[Clus1]]),b))

						}
						else if(!(a%in%unlist(Clusters)) & b%in%unlist(Clusters)){
							Clus1=which(sapply(1:length(Clusters),function(i) b%in%Clusters[[i]]))
							Clusters[[Clus1]]=unique(c(unlist(Clusters[[Clus1]]),a))
						}
						else{
							Clusters[[length(Clusters)+1]]=unique(c(a,b))
						}

					}

				}
			}
		}

		Singletons=which(!seq(1,nrow(S))%in%unlist(Clusters))
		if(length(Singletons)>0){

			for(l in Singletons){
				Clusters[[length(Clusters)+1]]=l
			}

		}


		order=unlist(Clusters)
		order.lab=rownames(S)[order]
		t1<-sapply(1:length(Clusters),function(i) rep(i,length(Clusters[[i]])))
		t2<-sapply(1:length(Clusters),function(i) rownames(S)[Clusters[[i]]])
		t3<-cbind(unlist(t1),unlist(t2))
		clus=as.numeric(t3[,1])
		names(clus)=t3[,2]
		Clusters=clus[rownames(S)]

		#Or perform single linkage on the similarity matrix S
		#Clusters=agnes(1-S,diss=TRUE,method="single")

	}

	else if(graphPartitioning=="SL_agnes"){
		Clust=cluster::agnes(1-S,diss=TRUE,method="single")
		order=Clust$order
		order.lab=rownames(S)[order]
		if(is.null(nrclusters)){
			Clusters=stats::cutree(stats::as.hclust(Clust),k=2)
		}
		else{
			Clusters=stats::cutree(stats::as.hclust(Clust),k=nrclusters[1])
		}
		names(Clusters)=rownames(S)
		Out=list(DistM=NULL,Clust=list(order=order,order.lab=order.lab,Clusters=Clusters))
		attr(Out,"method")="Single Clustering"
		return(Out)

	}

	Out=list(DistM=NULL,Clust=list(order=order,Clusters=Clusters))
	attr(Out,"method")="Ensemble"
	return(Out)

}
