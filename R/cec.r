
#' @title Complementary ensemble clustering
#'
#' @description Complementary Ensemble Clustering (CEC) Complementary Ensemble Clustering (CEC, \cite{Fodeh2013}) shows similarities with ADEC.
#' However, instead of merging the data matrices, ensemble clustering is performedon each data matrix separately. The resulting incidence
#' matrices for each data sets are combined in a weighted linear equation. The weighted incidence matrix is the input for the final clustering
#' algorithm. Similarly as ADEC, there are versions depending of the specification of the number of features to sample and the number of clusters.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param t The number of iterations. Defaults to 10.
#' @param r A vector with the number of features to take for the random sample for each element in List. If NULL (default), all features are considered.
#' @param nrclusters A list with a sequence of numbers of clusters to cut the dendrogram in for each element in List. If NULL (default), the function stops.
#' @param weight The weights for the weighted linear combination.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible"
#' @param weightclust A weight for which the result will be put aside of the other results. This was done for comparative reason and easy access.
#' @return The returned value is a list of four elements:
#' \item{DistM}{The resulting incidence matrix}
#' \item{Results}{The hierarchical clustering result for each element in WeightedDist}
#' \item{Clust}{The result for the weight specified in Clustweight}
#' The value has class 'CEC'.
#' @details If r is specified and nrclusters is a fixed number, only a random sampling of the features will be performed for the t iterations (CECa). If r is NULL and the nrclusters is a sequence, the clustering is performedon all features and the dendrogam is divided into clusters for the values of nrclusters (CECb). If both r is specified and nrclusters is a sequence, the combination is performed (CECc).
#' After every iteration, either be random sampling, multiple divisions of the dendrogram or both, an incidence matrix is set up. All incidence matrices are summed and represent the distance matrix on which a final clustering is performed.
#' @references
#' \insertRef{Fodeh2013}{MosaiClusteR}
CEC<-function(List,distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),t=10,r=NULL,nrclusters=NULL,weight=NULL,clust="agnes",linkage=c("flexible","flexible"),alpha=0.625,weightclust=0.5){
  
  if(!is.list(List)){
    stop("Data must be of type list")
  }
  
  if(is.null(nrclusters)){
    stop("Give a number of clusters to cut the dendrogram into for each data modality.")
  }
  
  if(!is.null(r)&length(nrclusters[[1]])==1){
    message("Performing a random sampling of the features with a fixed number of clusters.")
  }else if(is.null(r)&length(nrclusters[[1]])>1){
    message("Dividing the dendrogram in k clusters for a range of values of k.")
    t=1
  }else if(!is.null(r)&length(nrclusters[[1]])>1){
    message("Performing a random sampling of the features and dividing the dendrogram in k clusters for a range of values of k.")
  }
  else{
    stop("Specify r and/or nrclusters in order to perform an ADEC method.")
  }
  
  #Put all data in the same order
  OrderNames=rownames(List[[1]])
  for(i in 1:length(List)){
    List[[i]]=List[[i]][OrderNames,]
  }
  
  #Put up Incidence matrix for each data modality
  nc=c()
  Incidence=list()
  for (i in 1:length(List)){
    Incidence[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
    rownames(Incidence[[i]])=rownames(List[[i]])
    colnames(Incidence[[i]])=rownames(List[[i]])
    nc=c(nc,ncol(List[[i]]))
  }
  evenn=function(x){if(x%%2!=0)x=x-1 else x}
  nc=lapply(nc,FUN=evenn)
  nc=unlist(nc)
  
  #Repeat for t iterations
  for(g in 1:t){
    #message(g)
    if(is.null(r)){
      r=unlist(sapply(List,ncol))
    }
    
    #take random sample:
    A_prime=list()
    for(i in 1:length(r)){
      A=List[[i]]
      temp=sample(ncol(A),r[i],replace=FALSE)
      A_prime[[i]]=A[,temp]
      
      Ok=FALSE
      while(Ok==FALSE){
        if(is.numeric(A_prime[[i]])){
          if(any(rowSums(A_prime[[i]])==0)){
            temp=sample(ncol(A),r[i],replace=FALSE)
            A_prime[[i]]=A[,temp]
          }
          else{
            Ok=TRUE
          }
        }
        else{
          Ok=TRUE
        }
      }
    }
    
    #protect against zero rows:
    
    
    #Step 2: apply hierarchical clustering on each + cut tree into nrclusters
    
    DistM=lapply(seq(length(A_prime)),function(i) Distance(A_prime[[i]],distmeasure=distmeasure[i],normalize[i],method[i]))
    
    
    HClust_A_prime=lapply(seq(length(DistM)),function(i) cluster::agnes(DistM[[i]],diss=TRUE,method=linkage[i],par.method=alpha))
    
    
    for(o in seq(length(HClust_A_prime))){
      for(k in 1:length(nrclusters[[o]])){
        Temp=stats::cutree(HClust_A_prime[[o]],nrclusters[[o]][k])
        MembersofClust=matrix(1,dim(List[[o]])[1],dim(List[[o]])[1])
        
        for(l in 1:length(Temp)){
          label=Temp[l]
          sameclust=which(Temp==label)
          MembersofClust[l,sameclust]=0
        }
        Incidence[[o]]=Incidence[[o]]+MembersofClust
      }
    }
    
    
    
    
  }
  
  if(is.null(weight)){
    equalweights=1/length(List)
    weight=list(rep(equalweights,length(List)))
  }
  else if(is.list(weight) & length(weight[[1]])!=length(List)){
    stop("Give a weight for each data matrix or specify a sequence of weights")
  }
  else{
    #message('The weights are considered to be a sequence, each situation is investigated')
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
    for(i in 1:length(weight)){
      w=weight[[i]]
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
      
      weight=new
    }
  }
  
  weightedcomb<-function(w,Dist){
    temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
    temp=Reduce("+",temp)
    return(temp)
  }
  
  IncidenceComb=lapply(weight,weightedcomb,Incidence)
  namesweights=c()
  CEC=list()
  Clust=NULL; DistClust=NULL
  for (i in 1:length(IncidenceComb)){
    CEC[[i]]=cluster::agnes(IncidenceComb[[i]],diss=TRUE,method="ward")
    namesweights=c(namesweights,paste("Weight",weight[i],sep=" "))
    if(length(weight[[i]])==length(weightclust) && all(weight[[i]]==weightclust)){
      Clust=CEC[i]
      DistClust=IncidenceComb[i]
    }
  }
  # Fallback when the scalar 'weightclust' matches no generated combination
  # (e.g. the equal-weight focus does not sum to 1 for >2 sources): use the
  # equal-weights combination as the focus result.
  if(is.null(DistClust)){
    wf=rep(1/length(Incidence),length(Incidence))
    Dfocus=weightedcomb(wf,Incidence)
    Clust=list(cluster::agnes(Dfocus,diss=TRUE,method="ward"))
    DistClust=list(Dfocus)
  }

  Results=lapply(seq(1,length(CEC)),function(i) return(c("DistM"=IncidenceComb[i],"Clust"=CEC[i])))
  names(Results)=namesweights
  
  out=list(Incidence=Incidence,Results=Results,Clust=c("DistM"=DistClust,"Clust"=Clust))
  attr(out,'method')<-'CEC'
  return(out)
}
