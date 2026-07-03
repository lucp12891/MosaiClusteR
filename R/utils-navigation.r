# =============================================================================
# utils-navigation.r - navigation / comparison utilities for clustering results
#
# Ported (verbatim) from the legacy MoSaiClusteR project. These functions
# operate on a `List` of clustering-result objects. Each element is a list with
# `$Clust` (an agnes/hclust object) and/or `$DistM` (a dissimilarity matrix), as
# produced by methods such as Cluster(), WeightedClust() and CEC(). Objects are
# assumed to be in the rows (the package convention).
# =============================================================================

#' @title Relabel and reorder clusterings to a reference
#'
#' @description \code{ReorderToReference} relabels and reorders the clusters of
#' a series of clustering results so that the cluster numbers of every method
#' are matched, via a Gale-Shapley stable-matching scheme, to those of the first
#' (reference) method. The result is a matrix in which the columns represent the
#' objects (ordered as for the reference method) and the rows represent the
#' methods. Each cell contains the number of the cluster the object belongs to
#' for that method, relabelled to the reference. It is a possibility that two or
#' more clusters are fused together compared to the reference; in that case the
#' function alerts the user and asks to set \code{fusionsLog} to \code{TRUE}.
#'
#' @param List A list of clustering outputs to be compared. The first element
#' of the list will be used as the reference.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param fusionsLog Logical. Indicator for the fusion of clusters. Default is FALSE.
#' @param weightclust Logical. To be used for the outputs of CEC, WeightedClust
#' or WeightedSimClust. If TRUE, only the result of the Clust element is
#' considered. Default is FALSE.
#' @param names Optional. A character vector with the names of the methods.
#' @return A matrix of which the cells indicate to what cluster the objects
#' belong to according to the rearranged methods.
#' @note The \code{ReorderToReference} function was optimized for the situations
#' presented by the data sets at hand. It is noted that the function might fail
#' in a particular situation which results in an infinite loop.
#' @references Van Moerbeke, M. (2017). MoSaiClusteR: integrative clustering of
#' multi-source data. Hasselt University.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res <- list(Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
#'             Cluster(L[[2]], type = "data", distmeasure = "euclidean"))
#' M <- ReorderToReference(res, nrclusters = 5, fusionsLog = FALSE,
#'                         weightclust = FALSE, names = c("S1", "S2"))
#' }
#' @export
ReorderToReference<-function(List,nrclusters=NULL,fusionsLog=FALSE,weightclust=FALSE,names=NULL){

	matequal <- function(x, y)
		is.matrix(x) && is.matrix(y) && all(dim(x) == dim(y)) && all(x == y)

	ListNew=list()
	element=0
	for(i in 1:length(List)){
		if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
			ResultsClust=list()
			ResultsClust[[1]]=list()
			ResultsClust[[1]][[1]]=List[[i]]
			names(ResultsClust[[1]])[1]="Clust"
			element=element+1
			ListNew[[element]]=ResultsClust[[1]]
			#attr(ListNew[element],"method")="Weights"
		}
		else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
			ResultsClust=list()
			if(weightclust==TRUE){
				ResultsClust[[1]]=list()
				if(attributes(List[[i]])$method != "WeightedSim"){

					ResultsClust[[1]][[1]]=List[[i]]$Clust
					names(ResultsClust[[1]])[1]="Clust"
					attr(ResultsClust[[1]]$Clust,"method")="Weights"
					element=element+1
					ListNew[[element]]=ResultsClust[[1]]

				}
				else{
					ResultsClust[[1]]=list()
					ResultsClust[[1]][[1]]=List[[i]]
					names(ResultsClust[[1]])[1]="Clust"
					attr(ResultsClust[[1]]$Clust,"method")="Weights"
					element=element+1
					ListNew[[element]]=ResultsClust[[1]]
				}
			}
			else{
				for (j in 1:length(List[[i]]$Results)){
					ResultsClust[[j]]=list()
					ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
					names(ResultsClust[[j]])[1]="Clust"
					attr(ResultsClust[[j]]$Clust,"method")="Weights"
					element=element+1
					ListNew[[element]]=ResultsClust[[j]]
					#attr(ListNew[[element]],"method")="Weights"
				}
			}
		}
	}

	if(is.null(names)){
		names=seq(1,length(ListNew),1)
		for(i in 1:length(ListNew)){
			names[i]=paste("Method",i,sep=" ")
		}
	}
	names(ListNew)=names
	List=ListNew

	Clusters=list()
	CutTree<-function(i,Data,nrclusters){
		if(attributes(Data$Clust)$method == "Ensemble"){
			Clusters=Data$Clust$Clust$Clusters
			names(Clusters)=NULL
		}
		else{
			Clusters=stats::cutree(Data$Clust$Clust,k=nrclusters)
		}
		return(Clusters)
	}
	Clusters=lapply(seq(1,length(List)),function(i) CutTree(i,Data=ListNew[[i]],nrclusters=nrclusters))


	xaxis=List[[1]]$Clust$Clust$order #order of the objects as for method 1.
	xaxis.names=List[[1]]$Clust$Clust$order.lab #might be that names of methods are not in the same order...

	ordercolors=Clusters[[1]][xaxis]
	order=seq(1,nrclusters)

	for (k in 1:length(unique(Clusters[[1]][xaxis]))){
		select=which(Clusters[[1]][xaxis]==unique(Clusters[[1]][xaxis])[k])
		ordercolors[select]=order[k]
	}

	cols=unique(ordercolors) #order of the colors as depicted by method 1

	Ordered=list()

	autograph=list()
	for(i in cols){
		autograph[[i]]=xaxis[which(ordercolors==i)]
	}

	#for(j in 1:length(List)){
	#		temp=Clusters[[j]][xaxis]  #put clusternumbers of the other method into the same order as those of method (1)
	#	clusternumbers=temp		   #problem:cutree is based on the ordering of the names as they are in the rownames not in the order of joined objects
	#	for(k in 1:length(cols)){
	#		change=which(temp==unique(temp)[k])
	#		clusternumbers[change]=cols[which(cols==unique(temp)[k])]
	#	}
	#	Ordered[[j]]=clusternumbers
	#}

	for (j in 1:length(List)){
		message(j)
		#ordercolorsj=Clusters[[j]][xaxis]
		if(attributes(List[[j]]$Clust)$method=="Ensemble"){
			DistM=matrix(0,ncol=length(List[[j]]$Clust$Clust$order),nrow=length(List[[j]]$Clust$Clust$order))
			if(is.null(names(List[[j]]$Clust$Clust$Clusters))){
				names(List[[j]]$Clust$Clust$Clusters)=paste("Comp", c(1:nrow(DistM)),sep=" ")
			}
			colnames(DistM)=names(List[[j]]$Clust$Clust$Clusters)
			rownames(DistM)=names(List[[j]]$Clust$Clust$Clusters)
			List[[j]]$Clust$DistM=DistM
		}
		ordercolorsj=Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]
		order=seq(1,nrclusters)

		#for (k in 1:length(unique(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]))){
		for (k in unique(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))])){
			#select=which(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]==unique(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))])[k])
           	if(k<=nrclusters){
				select=which(Clusters[[j]][match(xaxis.names,rownames(List[[j]]$Clust$DistM))]==k)
				ordercolorsj[select]=order[k]
			}
		}


		temp2=ordercolorsj
		#temp3=xaxis
		temp3=match(xaxis.names,rownames(List[[j]]$Clust$DistM))
		fan=list()
		for(i in cols){
			fan[[i]]=match(xaxis.names,rownames(List[[j]]$Clust$DistM))[which(temp2==i)]
		}

		favors=matrix(0,length(autograph),length(fan))
		rownames(favors)=seq(1,length(autograph))
		colnames(favors)=seq(1,length(fan))

		for(a in 1:length(autograph)){
			for (b in 1:length(fan)){
				favorab=length(which(rownames(List[[j]]$Clust$DistM)[fan[[b]]] %in% rownames(List[[1]]$Clust$DistM)[autograph[[a]]]))/length(autograph[[a]])
				favors[a,b]=favorab
			}
		}

		#See function woman and men CB (put back what has value replaced)

		tempfavors=favors

		matched=c(rep("Free",nrclusters))
		proposed=c(rep("No",nrclusters))
		Switches=c(rep("Open",nrclusters))

		proposals=matrix(0,length(autograph),length(fan))

		#First match does "fans" that only have 1 element in their column: only one choice
		for(a in 1:dim(tempfavors)[1]){
			for (b in 1:dim(tempfavors)[2]){
				if(favors[a,b]==1){
					matched[a]=b
					proposed[b]="Yes"
					proposals[a,b]=1
					col=a

					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[b]]])
					temp3[change]=col

					tempfavors[,b]=0
					tempfavors[a,]=0

					Switches[a]="Closed"
				}

			}
		}


		#OneLeftC=FALSE
		#OneLeftR=FALSE
		for(b in 1:dim(tempfavors)[2]){
			if(length(which(tempfavors[,b]!=0))==1){
				match=which(tempfavors[,b]!=0)
				test=which(tempfavors[match,]==max(tempfavors[match,]))[1]
				if(length(which(tempfavors[,test]!=0))!=1 | b %in% which(tempfavors[match,]==max(tempfavors[match,])) ){
					matched[match]=b
					proposed[b]="Yes"
					proposals[match,b]=1
					col=match

					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[b]]])
					temp3[change]=col

					tempfavors[,b]=0
					tempfavors[match,]=0

					Switches[match]="Closed"
				}
			}
			#Unneccesary?
			#if(length(which(tempfavors[,b]==1))>1){
			#	matches=which(tempfavors[,b]==1)
			#	matched[matches]=b
			#	proposed[b]="Yes"
			#	proposals[matches,b]=1
			#	col=matches[1]
			#
			#	change=which(xaxis %in% fan[[b]])
			#	temp3[change]=col
			#
			#	tempfavors[,b]=0
			#	tempfavors[matches,]=0
			#
			#	Switches[matches]="Closed"
			#
			#	OneLeftC=TRUE
			#}
		}

		for(a in 1:dim(tempfavors)[1]){
			if(length(which(tempfavors[a,]!=0))==1){
				propose=which(tempfavors[a,]!=0)
				test=which(tempfavors[,propose]==max(tempfavors[,propose]))[1]
				if(length(which(tempfavors[test,]!=0))!=1 | a %in% which(tempfavors[,propose]==max(tempfavors[,propose]))){

					matched[a]=propose
					proposed[propose]="Yes"
					proposals[a,propose]=1
					col=a

					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
					temp3[change]=col

					tempfavors[a,]=0
					tempfavors[,propose]=0

					Switches[a]="Closed"
				}
			}
			#Unnecessary?
			#if(length(which(tempfavors[a,]==1))>1){
			#	proposes=which(tempfavors[a,]==1)
			#	matched[a]="Left"
			#	proposed[proposes]="Yes"
			#	proposals[a,proposes]=1
			#	col=a
			#
			#	change=which(xaxis %in% fan[[proposes]])
			#	temp3[change]=col
			#
			#	tempfavors[a,]=0
			#	tempfavors[,proposes]=0
			#
			#	Switches[a]="Closed"
			#
			#	OneLeftR=TRUE
			#}
		}
		Continue=TRUE
		if(length(which(matched=="Free")) == 0){
			Continue=FALSE
		}

		while(length(which(matched=="Free")) != 0 | !(matequal(proposals[which(matched=="Free"),], matrix(1, length(which(matched=="Free")), nrclusters))) | Continue!=FALSE){
			#for(a in which(matched=="Free")){
			#if(length(which(tempfavors[a,]!=0))==1){
			#	propose=which.max(tempfavors[a,])
			#
			#	matched[a]=propose
			#	proposed[propose]="Yes"
			#	proposals[a,propose]=1
			#	col=a
			#
			#	change=which(xaxis %in% fan[[propose]])
			#	temp3[change]=col
			#
			#	tempfavors[a,propose]=0
			#}

			#else{
			a=which(matched=="Free")[1]
			propose=which.max(tempfavors[a,])
			if(tempfavors[a,propose]==0){
				if(length(which(matched=="Free"))==1){
					Continue=FALSE
				}
				matched[a]="Left"
			}
			else{
				if(proposed[propose]=="No"){
					proposed[propose]="Yes"
					matched[a]=propose
					proposals[a,propose]=1
					col=a

					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
					temp3[change]=col

					tempfavors[a,propose]=0

					if(length(which(tempfavors[a,]==0))==dim(tempfavors)[2]){
						Switches[a]="Closed"
						tempfavors[,propose]=0

						c=1
						while(c < a){
							if(Switches[c] != "Closed" & length(which(tempfavors[c,]==0))==dim(tempfavors)[2]){
								Switches[c]="Closed"
								if(matched[c]=="Left"){
									tempfavors[c,]=0
								}
								else{
									tempfavors[,matched[c]]=0
								}
								c=1
							}
							else{
								c=c+1
							}
						}
					}
				}
				else if(proposed[propose]=="Yes"){
					if(favors[a,propose] > max(favors[which(matched==propose),propose]) & Switches[which(matched==propose)]=="Open"){


						#first undo then replace
						#tempfavors[which(matched==propose),propose]=favors[which(matched==propose),propose]

						changeback=which(xaxis.names %in%  rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[changeback]=match(xaxis.names,rownames(List[[j]]$Clust$DistM))[changeback]
						matched[which(matched==propose)]="Free"

						matched[a]=propose
						proposals[a,propose]=1
						col=a
						change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[change]=col

						tempfavors[a,propose]=0
					}
					else if(length(which(tempfavors[a,]!=0))==1){
						#if only 1 remains, these MUST BE matched
						changeback=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[changeback]=match(xaxis.names,rownames(List[[j]]$Clust$DistM))[changeback]
						matched[which(matched==propose)]="Free"

						matched[a]=propose
						proposals[a,propose]=1
						col=a

						change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[propose]]])
						temp3[change]=col

						tempfavors[a,propose]=0

					}
					else{
						proposals[a,propose]=1
						tempfavors[a,propose]=0
					}

				}
			}
			if(length(which(matched=="Free"))==0){
				Continue=FALSE
			}
		}
		fusions=0
		for( i in unique(matched)){
			if(length(which(!(seq(1,nrclusters) %in% matched)))>=1){
				fusions=length(which(!(seq(1,nrclusters) %in% matched)))
			}
		}

		if(fusions != 0 & fusionsLog==FALSE){
			message(paste("specify",fusions,"more color(s) and put fusionsLog equal to TRUE",sep=" "))
		}
		premiumcol=c()
		for (i in 1:(fusions)){
			premiumcol=c(premiumcol,length(matched)+i)
		}

		if((length(which(matched=="Left"))!=0) | (length(which(proposed=="No"))!=0)){
			if(length(which(proposed=="No"))!=0){
				for(i in 1:length(which(proposed=="No"))){
					Left=which(proposed=="No")[1]
					maxLeft=which(favors[,Left]==max(favors[,Left]))

					proposed[Left]="Yes"
					proposals[maxLeft,Left]=1
					col=premiumcol[i]

					change=which(xaxis.names %in% rownames(List[[j]]$Clust$DistM)[fan[[Left]]])
					temp3[change]=col

					tempfavors[,Left]=0
					tempfavors[maxLeft,]=0
				}

			}
			if(length(which(matched=="Left"))!=0){
				for (i in 1:length(which(matched=="Left")))
					Left=which(matched=="Left")[1]
				message(paste("Cluster",Left,"of the reference has found no suitable match.",sep=" "))
				#maxLeft=which(favors[Left,]==max(favors[Left,]))

			}
		}

		Ordered[[j]]=temp3

		if(max(ordercolorsj)>nrclusters){
			Ordered[[j]][which(ordercolorsj>nrclusters)]=ordercolorsj[which(ordercolorsj>nrclusters)]
		}
	}

	Matrix=c()
	for(j in 1:length(Ordered)){
		Matrix=rbind(Matrix,Ordered[[j]])
	}
	colnames(Matrix)=List[[1]]$Clust$Clust$order.lab
	rownames(Matrix)=names
	return(Matrix)

}


#' @title Find a selection of objects in the output of \code{ReorderToReference}
#'
#' @description \code{FindCluster} selects the objects belonging to a cluster
#' after the results of the methods have been rearranged by
#' \code{ReorderToReference}.
#'
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param select The row (the method) and the number of the cluster to select. Default is c(1,1).
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}:
#' indicator for the fusion of clusters. Default is TRUE.
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be
#' used for the outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only
#' the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return A character vector containing the names of the objects in the
#' selected cluster.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res <- list(Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
#'             Cluster(L[[2]], type = "data", distmeasure = "euclidean"))
#' Comps <- FindCluster(List = res, nrclusters = 5, select = c(1, 2),
#'                      weightclust = FALSE)
#' }
#' @export FindCluster
FindCluster<-function(List,nrclusters=NULL,select=c(1,1),fusionsLog=TRUE, weightclust=TRUE,names=NULL){
	if(length(List)==1 & attributes(List[[1]])$method == "Weighted" & weightclust==TRUE){
		T=List[[1]]$Clust
		attr(T,"method")="Single Clustering"
		List=list(T)
	}

	if(length(List)==1){
		Matrix=stats::cutree(List[[1]]$Clust,nrclusters)
		names(Matrix)=rownames(List[[1]]$DistM)
		clusternr=select[2]
		Comps=names(which(Matrix==clusternr))
	}

	else{
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		methodnr=select[1]
		clusternr=select[2]
		Comps=names(which(Matrix[methodnr,]==clusternr))
	}
	return(Comps)
}


#' @title Find an element in a data structure
#'
#' @description The function \code{FindElement} is used internally in the
#' \code{PreparePathway} function but might come in handy for other uses as well.
#' Given the name of an object, the function searches for that object in the data
#' structure and extracts it. When multiple objects have the same name, all are
#' extracted.
#'
#' @param what A character string indicating which object to look for. Default is NULL.
#' @param object The data structure to look into. Only the classes data frame
#' and list are supported. Default is NULL.
#' @param element Not to be specified by the user.
#' @return The returned value is a list with an element for each object found.
#' The element contains everything the object contained in the original data
#' structure.
#' @examples
#' \dontrun{
#' obj <- list(a = 1, b = list(TopDE = 42, c = list(TopDE = 7)))
#' FindElement(what = "TopDE", object = obj)
#' }
#' @export FindElement
FindElement<-function(what=NULL,object=NULL,element=list()){
	#str(Object)
	if(is.data.frame(object)){
		#search in columns
		if(what %in% colnames(object)){
			element[[length(element)+1]]<-object[,what]
			names(element)[length(element)]=paste(what,"_",length(element),sep="")
		}
		else if(what %in% rownames(object)){
			element[[length(element)+1]]<-object[what,]
			names(element)[length(element)]=paste(what,"_",length(element),sep="")
		}
	}
	if(is.list(object)){
		#Element=list()

		for(i in 0:length(object)){
			if(i==0){
				Names=names(object)
				if(what%in%Names){
					for(j in which(what==Names)){
						element[length(element)+1]=object[j]
						names(element)[length(element)]=paste(what,"_",length(element),sep="")
						return(element)
					}
				}
			}
			else if(class(object[[i]])[1]=="list"){
				#Names=names(Object[[i]])
				#if(What%in%Names){
				#	for(j in which(What==Names)){
				#			Element[length(Element)+1]=Object[[i]][j]
				#		names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
				#
				#	}
				#}
				element=FindElement(what,object[[i]],element=element)
				#for(j in 1:length(temp)){
				#	Element[length(Element)+1]=temp[j]
				#	names(Element)[length(Element)]=paste(What,"_",length(Element),sep="")
				#}

			}
			else if(class(object[[i]])[1]=="data.frame"){
				element=FindElement(what,object[[i]],element=element)

			}
		}
	}
	return(element)
}


#' @title Objects shared across clusterings
#'
#' @description \code{SharedComps} determines, per cluster, which objects are
#' shared across all listed clustering results. The clusterings are first
#' rearranged to a common reference with \code{ReorderToReference} and the
#' intersection of the lead objects of every method is taken per cluster.
#'
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param nrclusters If List is the output of several clustering methods, it has
#' to be provided in how many clusters to cut the dendrograms in. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}:
#' indicator for the fusion of clusters. Default is FALSE.
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be
#' used for the outputs of CEC, WeightedClust or WeightedSimClust. If TRUE, only
#' the result of the Clust element is considered. Default is FALSE.
#' @param names Names of the methods or clusters. Default is NULL.
#' @return A list containing the shared objects of all listed elements per cluster.
#' @examples
#' \dontrun{
#' data(mosaic_toy)
#' L <- mosaic_toy$List
#' res <- list(Cluster(L[[1]], type = "data", distmeasure = "euclidean"),
#'             Cluster(L[[2]], type = "data", distmeasure = "euclidean"))
#' Comps <- SharedComps(List = res, nrclusters = 5, fusionsLog = FALSE,
#'                      weightclust = FALSE, names = c("S1", "S2"))
#' }
#' @export SharedComps
SharedComps<-function(List,nrclusters=NULL,fusionsLog=FALSE,weightclust=FALSE,names=NULL){
	FoundComps=NULL

	FoundComps=FindElement("objects",List)

	if(is.null(FoundComps)|(is.list(FoundComps) & length(FoundComps) == 0)){
		ListNew=list()
		element=0
		for(i in 1:length(List)){
			if(attributes(List[[i]])$method != "CEC" & attributes(List[[i]])$method != "Weighted" & attributes(List[[i]])$method!= "WeightedSim"){
				ResultsClust=list()
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1
				ListNew[[element]]=ResultsClust[[1]]
				#attr(ListNew[element],"method")="Weights"
			}
			else if(attributes(List[[i]])$method=="CEC" | attributes(List[[i]])$method=="Weighted" | attributes(List[[i]])$method == "WeightedSim"){
				ResultsClust=list()
				if(weightclust==TRUE){
					ResultsClust[[1]]=list()
					if(attributes(List[[i]])$method != "WeightedSim"){
						ResultsClust[[1]][[1]]=List[[i]]$Clust
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1
						ListNew[[element]]=ResultsClust[[1]]
						attr(ListNew[element],"method")="Weights"
					}
					else{
						ResultsClust[[1]]=list()
						ResultsClust[[1]][[1]]=List[[i]]
						names(ResultsClust[[1]])[1]="Clust"
						element=element+1
						ListNew[[element]]=ResultsClust[[1]]
					}
				}
				else{
					for (j in 1:length(List[[i]]$Results)){
						ResultsClust[[j]]=list()
						ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
						names(ResultsClust[[j]])[1]="Clust"
						element=element+1
						ListNew[[element]]=ResultsClust[[j]]
						attr(ListNew[element],"method")="Weights"
					}
				}
			}
		}

		if(is.null(names)){
			names=seq(1,length(ListNew),1)
			for(i in 1:length(ListNew)){
				names[i]=paste("Method",i,sep=" ")
			}
		}
		names(ListNew)=names


		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		List=ListNew
		Comps=list()
		for (k in 1:dim(MatrixClusters)[1]){
			clusters=MatrixClusters[k,]

			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			clusters=MatrixClusters[k,]
			Method=list()
			for (i in 1:nrclusters){

				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)]
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds")
				names(temp)="objects"
				Method[[i]]=temp
				names(Method)[i]=paste("Cluster",i,sep=" ")
			}
			Comps[[k]]=Method

		}
		names(Comps)=names
		List=Comps
	}

	for(a in 1:length(List)){
		if(a==1){
			Comps1=FindElement("LeadCpds",List[[a]])
		}
		else if(is.list(Comps1) & length(Comps1) != 0){
			Comps2=FindElement("LeadCpds",List[[a]])
			for(b in 1:length(Comps1)){
				Comps1[[b]]=intersect(Comps1[[b]],Comps2[[b]])
			}
		}

	}

	namesCl=c()
	for(c in 1:length(Comps1)){
		namesCl=c(namesCl,paste("Cluster",c,sep=" "))
	}

	names(Comps1)=namesCl


	return(Comps1)
}
