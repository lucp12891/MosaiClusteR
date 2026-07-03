
#' @title Similarity network fusion
#'
#' @description Similarity Network Fusion (SNF) is a similarity-based multi-source clustering technique. SNF consists of two steps. In the initial step a similarity network is set up for each data matrix. The network is the visualization of the similarity matrix as a weighted graph with the objects as vertices and the pairwise similarities as weights on the edges. In the network-fusion step, each network is iteratively updated with information of the other network which results in more alike networks every time. This eventually converges to a single network.
#' @export
#' @param List A list of data matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have so.
#' If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
#' This is recommended if different types of data are used such that these are comparable. Default is FALSE.
#' @param NN The number of neighbours to be used in the procedure. Defaults to 20.
#' @param mu The parameter epsilon. The value is recommended to be between 0.3 and 0.8. Defaults to 0.5.
#' @param T The number of iterations.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for the final clustering. Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @return The returned value is a list with the following three elements.
#' \item{FusedM }{The fused similarity matrix}
#' \item{DistM }{The distance matrix computed by subtracting FusedM from one}
#' \item{Clust}{The resulting clustering}
#' The value has class 'SNF'.
#' @details If r is specified and nrclusters is a fixed number, only a random sampling of the features will be performed for the t iterations (ADECa). If r is NULL and the nrclusters is a sequence, the clustering is performedon all features and the dendrogam is divided into clusters for the values of nrclusters (ADECb). If both r is specified and nrclusters is a sequence, the combination is performed (ADECc).
#' After every iteration, either be random sampling, multiple divisions of the dendrogram or both, an incidence matrix is set up. All incidence matrices are summed and represent the distance matrix on which a final clustering is performed.
#' @references
#' \insertRef{Wang2014a}{MosaiClusteR}
SNF=function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,NN=20,mu=0.5,T=20,clust="agnes",linkage="ward",alpha=0.625){
  
  #Checking required data types and methods:
  if (!requireNamespace("SNFtool", quietly = TRUE))
    stop("SNF() requires the suggested package 'SNFtool'. ",
         "Install it with install.packages('SNFtool').")
  if(!is.list(List)){
    stop("Data must be of type list")
  }
  
  if(mu<0.3 | mu >0.8){
    message("Warning: mu is recommended to be between 0.3 and 0.8 for the SNF method. Default is 0.5.")
  }
  
  
  CheckDist<-function(Dist,StopRange){
    if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
      #message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
      Dist=Normalization(Dist,method="Range")
    }
    else{
      Dist=Dist
    }
  }
  
  
  #STEP 1: Distance Matrices
  if(type=="data"){
    OrderNames=rownames(List[[1]])
    for(i in 1:length(List)){
      List[[i]]=List[[i]][OrderNames,,drop=FALSE]
    }
    DistM=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
    DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
  }
  else if(type=="dist"){
    OrderNames=rownames(List[[1]])
    for(i in 1:length(List)){
      List[[i]]=List[[i]][OrderNames,OrderNames]
    }
    DistM=List
    DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
  }
  else{
    DistM=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
    DistM=lapply(seq(length(DistM)),function(i) CheckDist(DistM[[i]],StopRange))
    
    OrderNames=rownames(DistM[[1]])
    for(i in 1:length(DistM)){
      DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
    }
  }
  
  
  #STEP 2: Affinity Matrices
  
  AffM=lapply(seq(length(List)), function(i) SNFtool::affinityMatrix(DistM[[i]], NN, mu))
  
  #STEP 3: Fuse Networks Into 1 Single Network
  
  SNF_FusedM=SNFtool::SNF(AffM, NN, T)
  rownames(SNF_FusedM)=rownames(List[[1]])
  colnames(SNF_FusedM)=rownames(List[[1]])
  Dist=1-SNF_FusedM
  
  #STEP 4: Perform Hierarchical Clustering with WARD Link
  
  HClust = cluster::agnes(Dist,diss=TRUE,method=linkage,par.method=alpha)
  
  
  #Output= list with the fused matrix and the performed clustering
  out=list(SNF_FusedM=SNF_FusedM,DistM=Dist,Clust=HClust)
  attr(out,'method')<-'SNF'
  return(out)
}
