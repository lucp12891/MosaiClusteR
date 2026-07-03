
#' @title Aggregated data clustering
#'
#' @description Aggregated Data Clustering (ADC) is a direct clustering multi-source technique. The data matrices are merged into a single fused data matrix on which a dissimilarity matrix is computed. Hierarchical clustering is then performed on this dissimilarity matrix.
#'
#' @export
#' @param List A list of data matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param distmeasure Choice of metric for the dissimilarity matrix (character). Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize Logical. Indicates whether to normalize the distance matrices or not, defaults to FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list with the following three elements.
#' \item{AllData}{Fused data matrix of the data matrices}
#' \item{DistM}{The resulting distance matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'ADC'. The Clust element will be of interest for further applications.
#' @references Fodeh, S. J., Punch, W. & Tan, P.-N. (2013). On unifying multi-source information from probabilistically generated functional data. Proteins 81:1298-1310.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res_ADC <- ADC(List=L, distmeasure="euclidean", normalize=FALSE, method=NULL,
#'                clust="agnes", linkage="flexible", alpha=0.625)
#' }
ADC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="flexible",alpha=0.625){

	#Checking required data types and methods:
	if(!is.list(List)){
		stop("Data must be of type list")
	}


	#Fuse variables into 1 data Matrix

	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}

	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}

	#Compute Distance Matrix on AllData

	AllDataDist=Distance(AllData,distmeasure,normalize,method)

	#Perform hierarchical clustering with ward link on distance matrix

	HClust = cluster::agnes(AllDataDist,diss=TRUE,method=linkage,par.method=alpha)


	out=list(AllData=AllData,DistM=AllDataDist,Clust=HClust)
	attr(out,'method')<-'ADC'
	return(out)

}


#' @title Aggregated data ensemble clustering
#'
#' @description Aggregated Data Ensemble Clustering (ADEC) is a direct clustering multi-source technique. ADEC is an iterative procedure which starts with the merging of the data sets. In each iteration, a random sample of the features is selected and/or a resulting dendrogram is divided into k clusters for a range of values of k.
#'
#' @export
#' @param List A list of data matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param distmeasure Choice of metric for the dissimilarity matrix (character). Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param t The number of iterations. Defaults to 10.
#' @param r The number of features to take for the random sample. If NULL (default), all features are considered.
#' @param nrclusters A sequence of numbers of clusters to cut the dendrogram in. If NULL (default), the function stops.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list with the following three elements.
#' \item{AllData}{Fused data matrix of the data matrices}
#' \item{DistM}{The resulting co-association matrix}
#' \item{Clust}{The resulting clustering}
#' The value has class 'ADEC'. The Clust element will be of interest for further applications.
#' @details If r is specified and nrclusters is a fixed number, only a random sampling of the features will be performed for the t iterations (ADECa). If r is NULL and the nrclusters is a sequence, the clustering is performed on all features and the dendrogam is divided into clusters for the values of nrclusters (ADECb). If both r is specified and nrclusters is a sequence, the combination is performed (ADECc).
#' After every iteration, either be random sampling, multiple divisions of the dendrogram or both, an incidence matrix is set up. All incidence matrices are summed and represent the distance matrix on which a final clustering is performed.
#' @references Fodeh, S. J., Punch, W. & Tan, P.-N. (2013). On unifying multi-source information from probabilistically generated functional data. Proteins 81:1298-1310.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res_ADEC <- ADEC(List=L, distmeasure="euclidean", normalize=FALSE, method=NULL,
#'                  t=5, r=50, nrclusters=5, clust="agnes", linkage="flexible", alpha=0.625)
#' }
ADEC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,t=10,r=NULL,nrclusters=NULL,clust="agnes",linkage="flexible",alpha=0.625){

	if(!is.list(List)){
		stop("Data must be of type lists")
	}

	if(is.null(nrclusters)){
		stop("Give a number of cluters to cut the dendrogram into.")
	}

	if(!is.null(r)&length(nrclusters)==1){
		message("Performing a random sampling of the features with a fixed number of clusters.")
	}else if(is.null(r)&length(nrclusters)>1){
		message("Dividing the dendrogram in k clusters for a range of values of k.")
		t=1
	}else if(!is.null(r)&length(nrclusters)>1){
		message("Performing a random sampling of the features and dividing the dendrogram in k clusters for a range of values of k.")
	}
	else{
		stop("Specify r and/or nrclusters in order to perform an ADEC method.")
	}

	#Fuse A1 and A2 into 1 Data Matrix

	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}

	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}


	#take random sample of features

	nc=ncol(AllData)
	evenn=function(x){if(x%%2!=0) x=x-1 else x}
	nc=evenn(nc)


	#Put up Incidence matrix
	Incidence=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
	rownames(Incidence)=rownames(AllData)
	colnames(Incidence)=rownames(AllData)


	#Repeat for t iterations


	for(g in 1:t){
		#message(g)

		#if r is not fixed: changes per iteration. Need 1 value for r.
		if(is.null(r)){
			r=ncol(AllData)
		}

		#take random sample:
		ZeroPresent=TRUE
		while(ZeroPresent){
			temp1=sample(ncol(AllData),r,replace=FALSE)

			A_prime=AllData[,temp1]

			if(all(rowSums(A_prime)!=0)){
				ZeroPresent=FALSE
			}

		}
		#Step 2: apply hierarchical clustering on A_prime  + cut tree into nrclusters

		DistM=Distance(A_prime,distmeasure,normalize,method)

		HClust_A_prime=cluster::agnes(DistM,diss=TRUE,method=linkage,par.method=alpha)


		for(k in 1:length(nrclusters)){
			#message(k)
			Temp=stats::cutree(HClust_A_prime,nrclusters[k])
			MembersofClust=matrix(1,dim(List[[1]])[1],dim(List[[1]])[1])

			for(l in 1:length(Temp)){
				label=Temp[l]
				sameclust=which(Temp==label)
				MembersofClust[l,sameclust]=0
			}
			Incidence=Incidence+MembersofClust
		}

	}

	Clust=cluster::agnes(Incidence,diss=TRUE,method="ward",par.method=alpha)

	out=list(AllData=AllData,DistM=Incidence,Clust=Clust)
	attr(out,'method')<-'ADEC'
	return(out)

}
