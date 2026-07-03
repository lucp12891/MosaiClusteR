
#' @title Hierarchical ensemble clustering
#'
#' @description \insertCite{Zheng2014}{MosaiClusteR} proposed the Hierarchical Ensemble Clustering (HEC) algorithm. For each dendrogram, the cophenetic
#' distances between the object are calculated. The distances are aggregated across the data sets and an ultra-metric which is the closest to the
#' distance matrix is determined. A final hierarchical clustering is based on the ultra-metric values.
#' @export
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or"clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, default is FALSE. This is recommended if different distance types are used. More details on normalization in \code{Normalization}
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity (character) for each data set. Defaults to c("flexible", "flexible") for two data sets.
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625 and is only used if the linkage is set to "flexible".
#' @return The returned value is a list of two elements:
#' \item{DistM}{The resulting distance matrix}
#' \item{Clust}{The resulting hierarchical structure}
#' The value has class 'HEC'.
#' @references \insertRef{Zheng2014}{MosaiClusteR}
HierarchicalEnsembleClustering<-function(List,type=c("data","dist","clust"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),clust = "agnes", linkage = c("flexible","flexible"),alpha=0.625){
  
  ## Step 1: perfom aggl clustering on each data set
  
  if(type=="data"){
    
    if(is.null(rownames(List[[1]]))){
      for(i in 1:length(List)){
        rownames(List[[i]])=seq(1:nrow(List[[i]]))
      }
    }
    else{
      OrderNames=rownames(List[[1]])
      
      for(i in 1:length(List)){
        List[[i]]=List[[i]][OrderNames,,drop=FALSE]
      }
    }
    Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type="data",distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap=FALSE,maxK=5,StopRange=TRUE))
    #Clusterings=lapply(seq(length(List)),function(i) agnes(List[[i]],method=linkage[i],par.method=0.625))
    
    for(i in 1:length(Clusterings)){
      names(Clusterings)[i]=paste("Clust",i,sep=' ')
    }
    
  }
  
  else if(type=="dist"){
    OrderNames=rownames(List[[1]])
    for(i in 1:length(List)){
      List[[i]]=List[[i]][OrderNames,OrderNames]
    }
    
    Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize=FALSE,method=NULL,clust,linkage[i],alpha,gap=FALSE,maxK=5,StopRange=TRUE))
    
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
  
  ## Step 2: cophenetic distance of each object in each dendrogram
  
  CD=list()
  for(j in 1:length(Clusterings)){
    CD[[j]]=as.matrix(cophenetic(as.hclust(Clusterings[[j]]$Clust)))
  }
  
  ## Step 3: Aggregate the distance of the dendrograms and compute the matrix D
  ## ? what is meant by aggregate?
  ## Assumed: aggragete ==  added
  
  D_H=Reduce("+",CD)
  
  D=1/2*D_H
  
  ## Step 4: Find an ultra-metric distance T closest to D
  ## Apply the modified Floyd-Warshall algorithm on D
  
  Floyd_Warshall<-function(G){
    M=G
    N=nrow(G)
    for(k in c(1:N)){
      for(i in c(1:N)){
        for(j in c(1:N)){
          M[i,j]=min(M[i,j],max(M[i,k],M[k,j]))
          
        }
      }
    }
    return(M)
  }
  
  T=Floyd_Warshall(D)  # Distance matrix
  
  ## Step 5: Construct the final hierarchical clustering based on T
  ## the alpha-cut method form Meyer et al 2004
  
  ## ? defuzzify
  
  #	## Should
  #	test=cut(T, level = 1, type = c("alpha"), strict = FALSE)
  #
  #	#Testing
  #	A=matrix(c(0,1,14,18,18,18,1,0,14,18,18,18,14,14,0,18,18,18,18,18,18,0,15,15,18,18,18,15,0,2,18,18,18,15,2,0),ncol=6,nrow=6)
  #
  
  
  #try alpha cut for each value: denote cluster
  #build tree from this?
  
  alpha=sort(unique(as.vector(T)),decreasing = TRUE)
  C=list()
  #	for(i in alpha[-length(alpha)][c(1:11)]){
  #		C[[length(C)+1]]=archi(as.dist(T),i)[[1]]
  #
  #		if(length(C)>1){
  #
  #			if(max(C[[length(C)]])>max(C[[length(C)-1]])+1){
  #				Changed=which(C[[length(C)]]+ C[[length(C)-1]]<2*C[[length(C)]]-1)
  #				Groups=unique(C[[length(C)]][Changed])
  #
  #				if(length(Groups)==1){
  #					N=C[[length(C)]][-c(Changed)]
  #					O=C[[length(C)-1]][-c(Changed)]
  #					Joinedtemp=which(!as.vector(c(table(N)))==as.vector(c(table(O),0)))[c(2)]
  #					Joined=as.numeric(names(table(O))[Joinedtemp])
  #					New2=which(C[[length(C)]]==Joined[1])
  #					Groups=c(Groups,unique(C[[length(C)]][New2]))
  #					Changed=c(Changed,New2)
  #				}
  #
  #				Group1=Changed[which(C[[length(C)]][Changed]==Groups[1])]
  #				Group2=Changed[which(C[[length(C)]][Changed]==Groups[2])]
  #
  #				Clus=C[[length(C)-1]]
  #				Clusnext=C[[length(C)]]
  #				Clus1=Clusnext
  #				ChangeAfter=unique(Clus1[Group2])
  #				Clus1[Group2]=Clus[Group2]
  #				Clus1[which(Clus1>ChangeAfter)]=Clus1[which(Clus1>ChangeAfter)]-1
  #
  #				Clus2=Clusnext
  #
  #				C=C[-c(length(C))]
  #				C[[length(C)+1]]=Clus1
  #				C[[length(C)+1]]=Clus2
  #
  #			}
  #		}
  #	}
  
  for(k in c(1:nrow(List[[1]]))){
    Y=cutree(as.hclust(agnes(T,diss=TRUE)),k)
    names(Y)=NULL
    C[[length(C)+1]]=Y
  }
  #C[[length(C)+1]]=seq(1,nrow(T))
  
  Record=rep(0,length(C))
  names(Record)=seq(1,length(C))
  Merge=matrix(0,ncol=2,nrow=1)
  MergingList=list()
  
  for(i in length(C):1){
    if(i==length(C)){
      #All Leafs - Bottom of the tree
      Leafs=C[[length(C)]]
      Record=Leafs
    }
    
    else{
      
      NewMerge=C[[i]]
      if(i==(length(C)-1)){
        NewOnes=which(duplicated(NewMerge))
        for(n in 1:length(NewOnes)){
          MergedTogether=which(NewMerge==NewMerge[NewOnes][n])
          if(n==1){
            Merge[1,]=MergedTogether*(-1)
            
          }
          else{
            Merge=rbind(Merge,MergedTogether*(-1))
            
          }
          MergingList[[length(MergingList)+1]]=MergedTogether
          Record=NewMerge
        }
      }
      
      else{
        
        Changed=which(NewMerge!=Record)
        if(length(Changed)!=0){
          
          PositionChanged=Position(function(x) identical(x,Changed), MergingList, nomatch = 0)
          if(PositionChanged==0){
            
            formedafter=FALSE
            NewOnes=unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))[which(!unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))%in%unique(c(which(duplicated(Record)),which(duplicated(Record,fromLast=TRUE)))))]
            if(any(NewMerge[Changed]%in%NewMerge[NewOnes])){
              Temp=Changed[which(NewMerge[Changed]%in%NewMerge[NewOnes])]
              if(length(intersect(Temp,NewOnes))>0){
                NewOnes=Temp
                formedafter=TRUE
              }
            }
          }
          else{
            formedafter=FALSE
            NewOnes=which(NewMerge==unique(NewMerge[Changed]))[!which(NewMerge==unique(NewMerge[Changed]))%in%Changed]
          }
          if(length(NewOnes)==0){
            Joined=which(!as.vector(c(table(NewMerge),0))==as.vector(table(Record)))[c(1,2)]
            NewOnes=which(Record==Joined[1])
            
          }
          for(k in 1:length(unique(NewMerge[NewOnes]))){
            MergedTogether=which(NewMerge==unique(NewMerge[NewOnes])[k])
            
            Position1=Position(function(x) identical(x,(as.integer(MergedTogether[-which(MergedTogether%in%NewOnes)]))), MergingList, nomatch = 0)
            if(Position1==0){
              Merge=rbind(Merge,MergedTogether*(-1))  #they are both leaves
              
            }
            else{
              #First group was formed before
              if(formedafter==TRUE){
                TrulyNew=unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))[which(!unique(c(which(duplicated(NewMerge)),which(duplicated(NewMerge,fromLast=TRUE))))%in%unique(c(which(duplicated(Record)),which(duplicated(Record,fromLast=TRUE)))))]
                Group=NewOnes[-c(which(NewOnes%in%TrulyNew))]
                Set=list(TrulyNew,Group)
                for(s in Set){
                  if(length(s)>0){
                    Position2=Position(function(x) identical(x,(as.integer(MergedTogether[which(MergedTogether%in%s)]))), MergingList, nomatch = 0)
                    
                    if(Position2==0){  #first one was a group, newones is a leaf
                      Merge=rbind(Merge,c(Position1,(intersect(s,MergedTogether)*(-1))))
                      
                    }
                    else{  #both were grouped before
                      Merge=rbind(Merge,c(Position1,Position2))
                      
                    }
                  }
                }
              }
              else{
                Position2=Position(function(x) identical(x,(as.integer(MergedTogether[which(MergedTogether%in%NewOnes)]))), MergingList, nomatch = 0)
                
                if(Position2==0){  #first one was a group, newones is a leaf
                  Merge=rbind(Merge,c(Position1,(intersect(NewOnes,MergedTogether)*(-1))))
                  
                }
                else{  #both were grouped before
                  Merge=rbind(Merge,c(Position1,Position2))
                  
                }
              }
            }
            
            MergingList[[length(MergingList)+1]]=MergedTogether
          }
          Record=NewMerge
        }
      }
    }
  }
  
  
  Start=which(Merge==-1,arr.ind=TRUE)[1]
  Checked=rep(0,(nrow(Merge)))
  out=c()
  Row=Start
  while(!all(Checked==1)){
    Joined=Merge[Row,]
    #print(Row)
    if(sign(Joined[1])==-1 & sign(Joined[2])==-1){
      
      Checked[Row]=1
      
      if(Row==Start){
        
        out=abs(Joined)
        RNext=which(Merge==Row,arr.ind=TRUE)[1]
        
      }
      else{
        
        if(any(Checked>1)){
          out=c(out,abs(Joined))
          RNext=which.max(Checked)
        }
        
        else if(any(Checked==0)){
          out=c(out,abs(Joined))
          RNext=which(Merge==Row,arr.ind=TRUE)[1]
        }
      }
      
    }
    
    else if(sign(Joined[1])==1 & sign(Joined[2])==1){
      
      if(all(Checked[Joined]==1)){
        Checked[Row]=1
        
        if(any(Checked>1)){
          RNext=which.max(Checked)
        }
        else if(!(all(Checked==1))){
          RNext=which(Merge==Row,arr.ind=TRUE)[1]
        }
        
      }
      
      else{
        Checked[Row]=max(Checked)+1
        RNext=Joined[which(!Joined%in%which(Checked==1))]
        
        if(length(RNext)==2){
          Checked[RNext[2]]=max(Checked)+1
          RNext=RNext[1]
        }
      }
      
    }
    
    else if(sign(Joined[1])==1 & sign(Joined[2])==-1){
      
      Checked[Row]=1
      
      if(Checked[Joined[1]]==0){
        out=c(out,abs(Joined[2]))
        RNext=Joined[1]
      }
      else{
        if(!(abs(Joined)[2]%in%out)){
          out=c(out,abs(Joined[2]))
        }
        RNext=which(Merge==Row,arr.ind=TRUE)[1]
        
      }
    }
    
    else if(sign(Joined[1])==-1 & sign(Joined[2])==1){
      Checked[Row]=1
      
      if(Checked[Joined[2]]==0){
        out=c(out,abs(Joined[1])) #possibly to be shifted
        RNext=Joined[2]
        
      }
      else{
        if(!(abs(Joined)[1]%in%out)){
          out=c(out,abs(Joined[1]))
        }
        #					print("301")
        RNext=which(Merge==Row,arr.ind=TRUE)[1]
        
      }
      
    }
    
    Row=RNext
  }
  
  
  #	OrderObjects<-function(Row=Start,M=Merge,out=c(),Checked=rep(0,(nrow(T)-1))){
  #
  #		if(!all(Checked==1)){
  ##			print(out)
  #			print(Row)
  ##			print(Checked)
  ##			print("212")
  #
  #			Joined=M[Row,]
  #
  #			if(sign(Joined[1])==-1 & sign(Joined[2])==-1){
  #
  #				Checked[Row]=1
  #
  #				if(Row==Start){
  #
  #					out=abs(Joined)
  #					RNext=which(Merge==Row,arr.ind=TRUE)[1]
  ##					print("218")
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #
  #				}
  #				else{
  ##					print("222")
  #
  #
  #					if(any(Checked>1)){
  #						out=c(out,abs(Joined))
  #						#RNext=which(Merge==which.max(Checked),arr.ind=TRUE)[1]
  #						RNext=which.max(Checked)
  ##						print("227")
  #						#Checked[which.max(Checked)]=1
  #						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #					}
  #
  #					else if(any(Checked==0)){
  #						out=c(out,abs(Joined))
  #						RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #					}
  ##					return(out)
  #				}
  #
  #				if((all(Checked==1))){
  #					print(out)
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #				return(out)
  #			}
  #
  #
  #			else if(sign(Joined[1])==1 & sign(Joined[2])==1){
  #
  #				if(all(Checked[Joined]==1)){
  #					Checked[Row]=1
  ##					print("239")
  #					print(Checked)
  #
  #					if(any(Checked>1)){
  #						RNext=which.max(Checked)
  ##						print("242")
  #						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #					}
  #					else if(!(all(Checked==1))){
  #						RNext=which(Merge==Row,arr.ind=TRUE)[1]
  ##						print("247")
  #						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #					}
  #
  #				}
  #
  #				else{
  #					Checked[Row]=max(Checked)+1
  #					RNext=Joined[which(!Joined%in%which(Checked==1))]
  ##					print("259")
  #
  #					if(length(RNext)==2){
  #						Checked[RNext[2]]=max(Checked)+1
  #						RNext=RNext[1]
  ##						print("263")
  #						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #					}
  #					else{
  ##						print("269")
  #						out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #					}
  #				}
  #
  #				if((all(Checked==1))){
  #					print(out)
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #
  #				return(out)
  #			}
  #
  #			else if(sign(Joined[1])==1 & sign(Joined[2])==-1){
  #
  ##				print("276")
  #				Checked[Row]=1
  #
  #
  #				if(Checked[Joined[1]]==0){
  #					out=c(out,abs(Joined[2]))
  ##					print("280")
  #					RNext=Joined[1]
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #				else{
  ##					print("288")
  #					if(!(abs(Joined)[2]%in%out)){
  #						out=c(out,abs(Joined[2]))
  #					}
  #					RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #				if((all(Checked==1))){
  ##					print(out)
  ##					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #			}
  #				return(out)
  #			}
  #
  #			else if(sign(Joined[1])==-1 & sign(Joined[2])==1){
  #				Checked[Row]=1
  ##				print("293")
  #
  #				if(Checked[Joined[2]]==0){
  #					out=c(out,abs(Joined[1])) #possibly to be shifted
  ##					print("296")
  #					RNext=Joined[2]
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #				else{
  #					if(!(abs(Joined)[1]%in%out)){
  #						out=c(out,abs(Joined[1]))
  #					}
  ##					print("301")
  #					RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #
  #				if((all(Checked==1))){
  #					out=OrderObjects(Row=RNext,M=M,out=out,Checked=Checked)
  #				}
  #				return(out)
  #			}
  #
  #		}
  #
  #		else{
  ##			print("here")
  ##			print(out)
  ##			print(Checked)
  #			return(out)
  #		}
  #
  #
  #	}
  #	Order=OrderObjects(Row=Start,M=Merge,out=c(),Checked=rep(0,(nrow(T)-1)))
  #
  
  Order=out
  Heights=c()
  for(j in 1:nrow(Merge)){
    joined=Merge[j,]
    
    if(sign(joined[1])==-1){
      a=abs(joined[1])
    }
    else{
      a=MergingList[[joined[1]]]
      a=abs(a)
    }
    
    if(sign(joined[2])==-1){
      b=abs(joined[2])
    }
    else{
      b=MergingList[[joined[2]]]
      b=abs(b)
    }
    
    Heights=c(Heights,unique(as.vector(T[a,b]))) #should be put in the order of the names
    
  }
  
  
  #	HeightsObjects<-function(Row=Start,M=Merge,H=Heights,heights=c(),Checked=rep(0,(nrow(T)-1))){
  #
  #		if(!all(Checked==1)){
  #			print(heights)
  #			print(Row)
  #			print(Checked)
  #			print("212")
  #
  #			Joined=M[Row,]
  #
  #			if(sign(Joined[1])==-1 & sign(Joined[2])==-1){
  #
  #				Checked[Row]=1
  #
  #				if(Row==Start){
  #
  #					heights=c(heights,H[Row])
  #					RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #					print("218")
  #					heights=HeightsObjects(Row=RNext,M=M,H=H,heights=heights,Checked=Checked)
  #
  #				}
  #				else{
  #					print("222")
  #					heights=c(heights,H[Row])
  #
  #					if(any(Checked>1)){
  #						#RNext=which(Merge==which.max(Checked),arr.ind=TRUE)[1]
  #						RNext=which.max(Checked)
  #						print("227")
  #						#Checked[which.max(Checked)]=1
  #						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #					}
  ##					return(out)
  #				}
  #
  #				if((all(Checked==1))){
  #					print(heights)
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #				return(heights)
  #			}
  #
  #
  #			else if(sign(Joined[1])==1 & sign(Joined[2])==1){
  #
  #				if(all(Checked[Joined]==1)){
  #					Checked[Row]=1
  #					heights=c(heights,H[Row])
  #					print("239")
  #
  #					if(any(Checked>1)){
  #						RNext=which.max(Checked)
  #						print("242")
  #						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #					}
  #					else if(!(all(Checked==1))){
  #						RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #						print("247")
  #						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #					}
  #
  #				}
  #
  #				else{
  #					Checked[Row]=max(Checked)+1
  #					#Checked[Row]=1
  #					#heights=c(heights,H[Row])
  #					RNext=Joined[which(!Joined%in%which(Checked==1))]
  #					print("259")
  #
  #					if(length(RNext)==2){
  #						Checked[RNext[2]]=max(Checked)+1
  #						RNext=RNext[1]
  #						print("263")
  #						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #					}
  #					else{
  #						print("269")
  #						heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #					}
  #				}
  #
  #				if((all(Checked==1))){
  #					print(heights)
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #
  #				return(heights)
  #			}
  #
  #			else if(sign(Joined[1])==1 & sign(Joined[2])==-1){
  #
  #				print("276")
  #				Checked[Row]=1
  #				heights=c(heights,H[Row]) #possibly to be shifted
  #
  #				if(Checked[Joined[1]]==0){
  #					print("280")
  #					RNext=Joined[1]
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #				else{
  #					print("288")
  #					RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #				if((all(Checked==1))){
  #					print(heights)
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #				return(heights)
  #			}
  #
  #			else if(sign(Joined[1])==-1 & sign(Joined[2])==1){
  #				Checked[Row]=1
  #				print("293")
  #				heights=c(heights,H[Row]) #possibly to be shifted
  #				if(Checked[Joined[2]]==0){
  #					print("296")
  #					RNext=Joined[2]
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #				else{
  #					print("301")
  #					RNext=which(Merge==Row,arr.ind=TRUE)[1]
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #
  #				if((all(Checked==1))){
  #					print(heights)
  #					heights=HeightsObjects(Row=RNext,H=H,heights=heights,Checked=Checked)
  #				}
  #				return(heights)
  #			}
  #
  #		}
  #
  #		else{
  #			print("here")
  #			print(heights)
  #			return(heights)
  #		}
  #
  #
  #	}
  #	HeightsOrder<-HeightsObjects(Row=Start,M=Merge,H=Heights,heights=c(),Checked=rep(0,(nrow(T)-1)))
  
  Labels=rownames(T)
  
  out=list(height=Heights,merge=Merge,order=Order,labels=rownames(T))
  class(out)="hclust"
  plot(out)
  
  order.lab=rownames(T)[Order]
  Out=list(DistM=T,Clust=list(height=Heights,merge=Merge,order=Order,labels=rownames(T),order.lab=order.lab))
  attr(Out$Clust,"class")="hclust"
  attr(Out,"method")="HEC"
  return(Out)
  
}

