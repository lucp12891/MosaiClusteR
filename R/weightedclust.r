#' @title Weighted clustering
#'
#' @description Weighted Clustering (Weighted) is a similarity-based multi-source clustering technique. Weighted clustering is performed with the function \code{WeightedClust}. Given a
#' list of the data matrices, a dissimilarity matrix is computed of each with
#' the provided distance measures. These matrices are then combined resulting
#' in a weighted dissimilarity matrix. Hierarchical clustering is performed
#' on this weighted combination with the agnes function and the ward link
#'
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one should be standardized to have values between zero and one.
#' If FALSE the range normalization is performed. See \code{Normalization}. If TRUE, the distance matrices are not changed.
#' This is recommended if different types of data are used such that these are comparable. Default is FALSE.
#' @param weight Optional. A list of different weight combinations for the data sets in List.
#' If NULL, the weights are determined to be equal for each data set.
#' It is further possible to fix weights for some data matrices and to
#'  let it vary randomly for the remaining data sets. Defaults to seq(1,0,-0.1).  An example is provided in the details.
#' @param weightclust A weight for which the result will be put aside of the other results. This was done for comparative reason and easy access. Defaults to 0.5 (two data sets)
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for the final clustering. Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list of four elements:
#' \item{DistM}{A list with the distance matrix for each data structure}
#' \item{WeightedDist }{A list with the weighted distance matrices}
#' \item{Results}{The hierarchical clustering result for each element in WeightedDist}
#' \item{Clust}{The result for the weight specified in Clustweight}
#' The value has class 'Weighted'.
#' @details The weight combinations should be provided as elements in a list. For three data
#' matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
#' a fixed weight for some data sets and let it vary randomly for others, the element "x"
#' indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7
#' is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
#' the other two data sets. This implies that every combination of the sequence from 0 to 0.3
#' with steps of 0.1 will be reported and clustering will be performed for each.
#' @references
#' \insertRef{PerualilaTan2016}{MosaiClusteR}
WeightedClust <- function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),StopRange=FALSE,weight=seq(1,0,-0.1),weightclust=0.5,clust="agnes",linkage="ward",alpha=0.625){ # weight = weight to data1
  
  
  #Step 1: compute distance matrices:
  type<-match.arg(type)
  
  CheckDist<-function(Dist,StopRange){
    if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
      #message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
      Dist=Normalization(Dist,method="Range")
    }
    else{
      Dist=Dist
    }
  }
  
  
  if(type=="data"){
    OrderNames=rownames(List[[1]])
    for(i in 1:length(List)){
      List[[i]]=List[[i]][OrderNames,,drop=FALSE]
    }
    Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
    Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
  }
  else if(type=="dist"){
    OrderNames=rownames(List[[1]])
    for(i in 1:length(List)){
      List[[i]]=List[[i]][OrderNames,OrderNames]
    }
    Dist=List
    Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
  }
  else{
    Dist=lapply(seq(length(List)),function(i) return(List[[i]]$Dist))
    Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
    OrderNames=rownames(DistM[[1]])
    for(i in 1:length(DistM)){
      DistM[[i]]=DistM[[i]][OrderNames,OrderNames]
    }
  }
  
  #Step 2: Weighted linear combination of the distance matrices:
  if(is.null(weight)){
    equalweights=1/length(List)
    weight=list(rep(equalweights,length(List)))
    
  }
  else if(is.list(weight) & length(weight[[1]])!=length(List)){
    stop("Give a weight for each data matrix or specify a sequence of weights")
  }

  if(!is.list(weight)){
    condition<-function(l){
      l=as.numeric(l)
      if( sum(l)==1 ){  #working with characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
        #return(row.match(l,t1))
        return(l)
      }
      else(return(0))
    }
    
    if(all(seq(1,0,-0.1)!=weight)){
      for(i in 1:length(weight)){
        rest=1-weight[i]
        if(!(rest%in%weight)){
          weight=c(weight,rest)
        }
      }
    }
    
    
    
    
    t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
    t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
    t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
    t4=t2[which(t3!=0)]
    weight=lapply(seq(length(t4)),function(i) rev(t4[[i]]))
    
  }
  if(is.list(weight) & "x" %in% weight[[1]]){ #x indicates a free weight
    newweight=list()
    for(k in 1:length(weight)){
      w=weight[[k]]
      weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
      givenweights=as.numeric(w[weightsfordata])
      
      stilltodistribute=1-sum(givenweights)
      
      newweights=seq(stilltodistribute,0,-0.1)
      
      t1=gtools::permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
      Input1=as.list(seq_len(nrow(t1)))
      Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
      return(Input1[[i]])})
      t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
      t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
      weightsforotherdata=t2[which(t3!=0)]
      
      new=list()
      for(i in 1:length(weightsforotherdata)){
        w1=weightsforotherdata[[i]]
        new[[i]]=rep(0,length(List))
        new[[i]][weightsfordata]=givenweights
        new[[i]][which(new[[i]]==0)]=w1
      }
      
      newweight[k]=new
    }
    
    weight=newweight
  }
  weightedcomb<-function(w,Dist){
    temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
    temp=Reduce("+",temp)
    return(temp)
  }
  
  DistClust=NULL
  Clust=NULL
  
  DistM=lapply(seq(length(weight)),function(i) weightedcomb(weight[[i]],Dist=Dist))
  namesweights=c()
  WeightedClust=lapply(seq(length(weight)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method=linkage,par.method=alpha))
  for(i in 1:length(WeightedClust)){
    namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
    if(length(weight[[i]])==length(weightclust) && all(weight[[i]]==weightclust)){
      Clust=WeightedClust[[i]]
      DistClust=DistM[[i]]
    }
  }

  if(is.null(DistClust)){
    # Focus weight for the returned 'Clust'. A scalar 'weightclust' only makes
    # sense for two sources; for a different number (or a non-matching scalar)
    # fall back to equal weights so weightedcomb() gets a per-source vector
    # (a scalar would recycle to NA for sources 2..K and break agnes).
    wf=if(length(weightclust)==length(Dist)) weightclust else rep(1/length(Dist),length(Dist))
    DistClust=weightedcomb(wf,Dist=Dist)
    Temp=cluster::agnes(DistClust,diss=TRUE,method=linkage,par.method=alpha)
    Clust=Temp
  }
  
  Results=lapply(seq(1,length(WeightedClust)),function(i) return(c("DistM"=DistM[i],"Clust"=WeightedClust[i])))
  names(Results)=namesweights
  
  # return list with objects
  out=list(Dist=Dist,Results=Results,Clust=list("DistM"=DistClust,"Clust"=Clust))
  attr(out,'method')<-'Weighted'
  return(out)
  
}
