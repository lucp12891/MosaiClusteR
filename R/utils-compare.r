## Cluster-comparison and weight-determination utilities.
## Ported verbatim (with R 4.x safety fixes and suggested-package guards) from the
## original MoSaiClusteR project. These functions rely on the package-provided
## Distance(), Normalization(), Cluster(), SimilarityMeasure(),
## ReorderToReference(), ColorsNames() and FindElement().

#' @title Interactive comparison of single and multiple source clustering results
#'
#' @description The function \code{CompareInteractive} produces an interactive plot of the
#' multiple source clustering results in \code{ListM}. By identifying a cluster
#' or a method, a comparison with the single source clustering results in
#' \code{ListS} is shown in a new plotting window.
#' @export CompareInteractive
#' @param ListM A list of the outputs from the multiple source clusterings to be compared.
#' @param ListS A list of the outputs from the single source clusterings to be compared.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the names of the colours. Default is NULL.
#' @param fusionsLogM The fusionsLog parameter for the elements in ListM. Default is FALSE.
#' @param fusionsLogS The fusionsLog parameter for the elements in ListS. Default is FALSE.
#' @param weightclustM The weightclust parameter for the elements in ListM. Default is FALSE.
#' @param weightclustS The weightclust parameter for the elements in ListS. Default is FALSE.
#' @param namesM Optional. Names of the multiple source clusterings. Default is NULL.
#' @param namesS Optional. Names of the single source clusterings. Default is NULL.
#' @param marginsM Optional. Margins to be used for the ListM plot. Default is c(2,2.5,2,2.5).
#' @param marginsS Optional. Margins to be used for the ListS plot. Default is c(8,2.5,2,2.5).
#' @param Interactive Optional. Do you want an interactive plot? Default is TRUE.
#' @param n The number of methods/clusters you want to identify. Default is 1.
#' @return A plot of the comparison of the elements of ListM, on which multiple clusters
#' and/or methods can be identified to compare them to the elements in ListS.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(fingerprintMat,targetMat)
#'
#' MCF7_W=WeightedClust(List=L,type="data",distmeasure=c("tanimoto","tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(1,0,-0.1),weightclust=0.5,
#' clust="agnes",linkage="ward",StopRange=FALSE)
#'
#' ListM=list(MCF7_W)
#' namesM=c(seq(1.0,0.0,-0.1))
#'
#' ListS=list(MCF7_F,MCF7_T)
#' namesS=c("FP","TP")
#'
#' CompareInteractive(ListM,ListS,nrclusters=7,cols=Colors1,fusionsLogM=FALSE,
#' fusionsLogS=FALSE,weightclustM=FALSE,weightclustS=TRUE,namesM,namesS,
#' marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,n=1)
#' }
CompareInteractive<-function(ListM,ListS,nrclusters=NULL,cols=NULL,fusionsLogM=FALSE,fusionsLogS=FALSE,weightclustM=FALSE,weightclustS=FALSE,namesM=NULL,namesS=NULL,marginsM=c(2,2.5,2,2.5),marginsS=c(8,2.5,2,2.5),Interactive=TRUE,n=1){
	if (!requireNamespace("plotrix", quietly=TRUE)) stop("CompareInteractive() requires the suggested package 'plotrix'.")

	MatrixColorsM=ReorderToReference(ListM,nrclusters,fusionsLogM,weightclustM,namesM)

	NamesM=ColorsNames(MatrixColorsM,cols)

	nobsM=dim(MatrixColorsM)[2]
	nmethodsM=dim(MatrixColorsM)[1]

	if(is.null(namesM)){
		for(j in 1:nmethodsM){
			namesM[j]=paste("Method",j,sep=" ")
		}
	}

	similarM=round(SimilarityMeasure(MatrixColorsM),2)

	grDevices::dev.new()
	graphics::par(mar=marginsM)
	plotrix::color2D.matplot(MatrixColorsM,cellcolors=NamesM,show.values=FALSE,axes=FALSE,xlab="",ylab="")
	graphics::axis(1,at=seq(0.5,(nobsM-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsM-0.5)),labels=rev(namesM),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsM-0.5)),labels=rev(similarM),cex.axis=0.65,las=2)

	if(Interactive==TRUE){
		yseq=c(seq(dim(MatrixColorsM)[1]-0.5,0.5,-1))
		for(i in seq(dim(MatrixColorsM)[1]-0.5,0.5,-1)){
			yseq=c(yseq,rep(i,dim(MatrixColorsM)[2]))
		}
		ids=graphics::identify(x=c(rep(-1,dim(MatrixColorsM)[1]),rep(seq(0.5,dim(MatrixColorsM)[2]-0.5),dim(MatrixColorsM)[1])),y=yseq,n=n,plot=FALSE)

		comparison<-function(id){
			if(id%in%seq(dim(MatrixColorsM)[1])){
				grDevices::dev.new()
				graphics::layout(matrix(c(1,2),nrow=2), heights=c(1,2))
				NamesMSel=NamesM[id,]
				namesMSel=namesM[id]

				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(t(as.matrix(MatrixColorsM[id,])),cellcolors=NamesMSel,show.values=FALSE,axes=FALSE,xlab="",ylab="")
				graphics::axis(2,at=c(0.5),labels=rev(namesMSel),cex.axis=0.65,las=2)
				#axis(4,at=c(0.5),labels=rev(similarSel),cex.axis=0.65,las=2)

				#Find reference for MatrixColorsM
				if(weightclustM==FALSE){
					temp=FindElement("Results",ListM[1])
					if(!(is.null(temp))&length(temp)!=0){
						Ref=list(Clust=temp$Results_1[[1]])
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else if(length(temp)==0){
						temp=FindElement("Clust",ListM[1])
						Ref=list(Clust=temp$Clust_1)
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else{
						message('Cannot find a reference for the second plot, try: weightclust=TRUE')
					}
				}
				else{
					Ref=ListM[[1]]
					attr(Ref,'method')<-attributes(ListM[[1]])$method
				}

				L=c(Ref,ListS)
				for(i in 1:length(L)){
					if(i==1){
						attr(L[[1]],'method')<-"Ref"
					}
					else{
						attr(L[[i]],"method")<-attributes(ListS[[i-1]])$method
					}
				}
				MatrixColorsS=ReorderToReference(L,nrclusters,fusionsLogS,weightclustS,names=c("Ref",namesS))
				MatrixColorsS=MatrixColorsS[-1,]
				NamesS=ColorsNames(MatrixColorsS,cols)

				nobs=dim(MatrixColorsS)[2]
				nmethodS=dim(MatrixColorsS)[1]

				if(is.null(namesS)){
					for(j in 1:nmethodS){
						namesS[j]=paste("Method",j,sep=" ")
					}
				}

				similarS=round(SimilarityMeasure(MatrixColorsS),2)

				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")
				graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
				graphics::axis(2,at=c(seq(0.5,nmethodS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
				graphics::axis(4,at=c(seq(0.5,nmethodS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)
			}
			else{
				SelCluster=t(MatrixColorsM)[id-nrow(MatrixColorsM)]
				Temp=sapply(seq(nrow(MatrixColorsM)),function(i) ncol(MatrixColorsM)*i)
				Row=which(Temp>(id-nrow(MatrixColorsM)))[1]
				Index=which(MatrixColorsM[Row,]!=SelCluster)

				grDevices::dev.new()
				graphics::layout(matrix(c(1,2),nrow=2), heights=c(1,2))
				NamesMSel=NamesM[Row,]
				NamesMSel[Index]="white"
				namesMSel=namesM[Row]

				graphics::par(mar=marginsM)
				plotrix::color2D.matplot(t(as.matrix(MatrixColorsM[Row,])),cellcolors=NamesMSel,show.values=FALSE,axes=FALSE,xlab="",ylab="")
				graphics::axis(2,at=c(0.5),labels=rev(namesMSel),cex.axis=0.65,las=2)

				#Find reference for MatrixColorsM
				if(weightclustM==FALSE){
					temp=FindElement("Results",ListM[1])
					if(!(is.null(temp))&length(temp)!=0){
						Ref=list(Clust=temp$Results_1[[1]])
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else if(length(temp)==0){
						temp=FindElement("Clust",ListM[1])
						Ref=list(Clust=temp$Clust_1)
						attr(Ref,'method')<-attributes(ListM[[1]])$method
					}
					else{
						message('Cannot find a reference for the second plot, try: weightclust=TRUE')
					}
				}
				else{
					Ref=ListM[[1]]
					attr(Ref,'method')<-attributes(ListM[[1]])$method
				}

				L=c(Ref,ListS)
				for(i in 1:length(L)){
					if(i==1){
						attr(L[[1]],'method')<-"Ref"
					}
					else{
						attr(L[[i]],"method")<-attributes(ListS[[i-1]])$method
					}
				}

				MatrixColorsS=ReorderToReference(L,nrclusters,fusionsLogS,weightclustS,names=c("Ref",namesS))
				MatrixColorsS=MatrixColorsS[-1,]

				IndexS=lapply(seq(nrow(MatrixColorsS)),function(i) which(MatrixColorsS[i,]!=SelCluster))

				NamesS=ColorsNames(MatrixColorsS,cols)
				for(i in 1:nrow(NamesS)){
					NamesS[i,IndexS[[i]]]="white"
				}

				nobs=dim(MatrixColorsS)[2]
				nmethodsS=dim(MatrixColorsS)[1]

				if(is.null(namesS)){
					for(j in 1:nmethodsS){
						namesS[j]=paste("Method",j,sep=" ")
					}
				}

				similarS=round(SimilarityMeasure(MatrixColorsS),2)

				graphics::par(mar=marginsS)
				plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")
				graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
				graphics::axis(2,at=c(seq(0.5,nmethodsS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
				graphics::axis(4,at=c(seq(0.5,nmethodsS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)


			}
		}

		plots=sapply(seq(length(ids)),function(i) comparison(ids[i]))

	}

}

#' @title Compares medoid clustering results based on silhouette widths
#'
#' @description The function \code{CompareSilCluster} compares the results of two medoid
#' clusterings. The null hypothesis is that the clustering is identical. A test
#' statistic is calculated and a p-value obtained with bootstrapping.
#' @export CompareSilCluster
#' @param List A list of data matrices. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices or distance
#' matrices. Type should be one of "data" or "dist".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices. Default is c(FALSE, FALSE).
#' @param method A method of normalization. Default is c(NULL,NULL).
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param nboot Number of bootstraps to be run. Default is 100.
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "new".
#' @param location If plottype is "pdf", a location should be provided here. Default is NULL.
#' @return A plot of the density of the statistic under the null hypothesis. Further a list with two
#' elements: \item{Observed Statistic}{The observed statistical value} \item{P-Value}{The P-value of the
#' obtained statistic retrieved after bootstrapping}.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#'
#' List=list(fingerprintMat,targetMat)
#'
#' Comparison=CompareSilCluster(List=List,type="data",
#' distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' nrclusters=7,names=NULL,nboot=100,plottype="new",location=NULL)
#'
#' Comparison
#' }
CompareSilCluster<-function(List,type=c("data","dist"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),nrclusters=NULL,names=NULL,nboot=100,plottype="new",location=NULL){

	type=match.arg(type)

	if(is.null(names)){
		names=c()
		for(i in 1:length(List)){
			names=c(names,paste("Method",i,sep=" "))
		}
	}

	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names

	}
	else{
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters)$silinfo$widths)
		names(silwidth)=names
	}


	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){

		}
	}

	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}

	regressioncomb=gtools::permutations(n=length(List),r=2,repeats.allowed=T)

	StatRSq<-function(regressioncomb,silwidth,ordernames,names){

		regressRSq<-function(x,silwidth,ordernames,names){

			i1=x[1]
			i2=x[2]
			L1=silwidth[[i1]][,3][ordernames]
			L2=silwidth[[i2]][,1][ordernames]

			regress<-stats::lm(L1~L2)
			Rsq<-summary(regress)$r.squared
			#names(Rsq)=paste("RSquared_",names[i1],names[i2],sep="_")
			return(Rsq)

			#paste names on this object!!!
		}

		RSqs=apply(regressioncomb,1,function(x) regressRSq(x,silwidth,ordernames,names))

		#for (i in 1:nrow(regressioncomb)){
		#	names(RSqs)[i]=paste("RSquared_",names[regressioncomb[i,1]],names[regressioncomb[i,2]],sep="_")
		#}

		stat=0
		xx=0
		xy=0
		for(i in 1:nrow(regressioncomb)){
			if(regressioncomb[i,1]==regressioncomb[i,2]){
				xx=xx+RSqs[i]
			}
			else{
				xy=xy+RSqs[i]
			}
		}
		stat=abs(xx-xy)  #check this formula with Nolen
		names(stat)=NULL
		return(stat)

	}

	StatRSqObs=StatRSq(regressioncomb,silwidth,ordernames=rownames(Dist[[1]]),names)


	#bootstrapping
	statNULL=c(1:nboot)
	perm.rowscols <- function (D, n)
	{
		s <- sample(1:n)
		D=D[s, s]
		return(D)
	}

	for(i in 1:nboot){
		set.seed(i)
		DistNULL=Dist
		DistNULL[[1]] <- perm.rowscols(DistNULL[[1]],nrow(DistNULL[[1]]))

		silwidthNULL=lapply(DistNULL,function(x) cluster::pam(x,nrclusters)$silinfo$widths)

		statNULL[i]=StatRSq(regressioncomb,silwidthNULL,ordernames=rownames(DistNULL[[1]]),names)

	}

	pval=(sum(abs(statNULL)<=abs(StatRSqObs))+1)/(nboot+1)

	plottypein(plottype,location)
	graphics::plot(stats::density(statNULL),type="l",main="The Density of the Statistic under the H0")
	graphics::abline(v=StatRSqObs)

	out=list()
	out[[1]]=StatRSqObs
	#out[[2]]=statNULL
	out[[2]]=pval
	names(out)=c("Observed Statistic","P-Value")

	return(out)
}

#' @title Comparison of clustering results for the single and multiple source clustering.
#'
#' @description The function \code{CompareSvsM} plots the comparison of the single source
#' clustering results on the left and that of the multiple source clustering results on the
#' right such that a visual comparison is possible.
#' @export CompareSvsM
#' @param ListS A list of the outputs from the single source clusterings to be compared.
#' @param ListM A list of the outputs from the multiple source clusterings to be compared.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the names of the colours. Default is NULL.
#' @param fusionsLogS The fusionslog parameter for the elements in ListS. Default is FALSE.
#' @param fusionsLogM The fusionsLog parameter for the elements in ListM. Default is FALSE.
#' @param weightclustS The weightclust parameter for the elements in ListS. Default is FALSE.
#' @param weightclustM The weightclust parameter for the elements in ListM. Default is FALSE.
#' @param namesS Optional. Names of the single source clusterings. Default is NULL.
#' @param namesM Optional. Names of the multiple source clusterings. Default is NULL.
#' @param margins Optional. Margins to be used for the plot. Default is c(8.1,3.1,3.1,4.1).
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "new".
#' @param location If plottype is "pdf", a location should be provided here. Default is NULL.
#' @return A plot with on the left the comparison over the objects in ListS and on the right a
#' comparison over the objects in ListM.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(fingerprintMat,targetMat)
#'
#' MCF7_W=WeightedClust(List=L,type="data", distmeasure=c("tanimoto","tanimoto"),
#' normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(1,0,-0.1),weightclust=0.5
#' ,clust="agnes",linkage="ward",StopRange=FALSE)
#'
#' ListM=list(MCF7_W)
#' namesM=seq(1.0,0.0,-0.1)
#'
#' ListS=list(MCF7_F,MCF7_T)
#' namesS=c("FP","TP")
#'
#' CompareSvsM(ListS,ListM,nrclusters=7,cols=Colors1,fusionsLogS=FALSE,
#' fusionsLogM=FALSE,weightclustS=FALSE,weightclustM=FALSE,namesS,
#' namesM,plottype="new",location=NULL)
#' }
CompareSvsM<-function(ListS,ListM,nrclusters=NULL,cols=NULL,fusionsLogS=FALSE,fusionsLogM=FALSE,weightclustS=FALSE,weightclustM=FALSE,namesS=NULL,namesM=NULL,margins=c(8.1,3.1,3.1,4.1),plottype="new",location=NULL){
	if (!requireNamespace("plotrix", quietly=TRUE)) stop("CompareSvsM() requires the suggested package 'plotrix'.")
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new(wdith=14,height=7)
		}
		if(plottype=="sweave"){

		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}
	nmethodsS=0
	nmethodsM=0

	MatrixColorsS=ReorderToReference(ListS,nrclusters,fusionsLogS,weightclustS,namesS)
	MatrixColorsM=ReorderToReference(c(ListS[1],ListM),nrclusters,fusionsLogM,weightclustM,c("ref",namesM))

	similarS=round(SimilarityMeasure(MatrixColorsS),2)
	similarM=round(SimilarityMeasure(MatrixColorsM),2)

	MatrixColorsM=MatrixColorsM[-c(1),]

	NamesM=ColorsNames(MatrixColorsM,cols)
	NamesS=ColorsNames(MatrixColorsS,cols)

	nobsM=dim(MatrixColorsM)[2]
	nmethodsM=dim(MatrixColorsM)[1]

	nobsS=dim(MatrixColorsS)[2]
	nmethodsS=dim(MatrixColorsS)[1]

	if(is.null(namesS)){
		for(j in 1:nmethodsS){
			namesS[j]=paste("Method",j,sep=" ")
		}
	}

	if(is.null(namesM)){
		for(j in 1:nmethodsM){
			namesM[j]=paste("Method",j,sep=" ")
		}
	}


	plottypein(plottype,location)
	graphics::par(mfrow=c(1,2),mar=margins)
	plotrix::color2D.matplot(MatrixColorsS,cellcolors=NamesS,show.values=FALSE,axes=FALSE,xlab="",ylab="")
	graphics::axis(1,at=seq(0.5,(nobsS-0.5)),labels=colnames(MatrixColorsS),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsS-0.5)),labels=rev(namesS),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsS-0.5)),labels=rev(similarS),cex.axis=0.65,las=2)

	plotrix::color2D.matplot(MatrixColorsM,cellcolors=NamesM,show.values=FALSE,axes=FALSE,xlab="",ylab="")
	graphics::axis(1,at=seq(0.5,(nobsM-0.5)),labels=colnames(MatrixColorsM),las=2,cex.axis=0.70)
	graphics::axis(2,at=seq(0.5,(nmethodsM-0.5)),labels=rev(namesM),cex.axis=0.65,las=2)
	graphics::axis(4,at=seq(0.5,(nmethodsM-0.5)),labels=rev(similarM[-1]),cex.axis=0.65,las=2)
	plottypeout(plottype)

}

#' @title Determines an optimal weight for weighted clustering by silhouettes widths.
#'
#' @description The function \code{DetermineWeight_SilClust} determines an optimal weight
#' for weighted similarity clustering by calculating silhouettes widths. For each given
#' weight, a linear combination of the distance matrices of the single data sources is
#' obtained, medoid clustering is performed and the silhouette widths are regressed against
#' the cluster memberships. A statistic is derived and a p-value obtained via bootstrapping.
#' @export DetermineWeight_SilClust
#' @param List A list of matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. Type should be one of "data", "dist" or "clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices. Default is c(FALSE, FALSE).
#' @param method A method of normalization. Default is c(NULL,NULL).
#' @param weight Optional. A list of different weight combinations for the data sets in List. Defaults to seq(0,1,by=0.01).
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param nboot Number of bootstraps to be run. Default is 10.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one
#' should be standardized to have so. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "new".
#' @param location If plottype is "pdf", a location should be provided here. Default is NULL.
#' @return Two plots and a list with two elements: \item{Result}{A data frame with the statistic for each
#' weight combination} \item{Weight}{The optimal weight}.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#'
#' MC7_Weight=DetermineWeight_SilClust(List=L,type="clusters",distmeasure=
#' c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),
#' weight=seq(0,1,by=0.01),nrclusters=c(7,7),names=c("FP","TP"),nboot=10,
#' StopRange=FALSE,plottype="new",location=NULL)
#' }
DetermineWeight_SilClust<-function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(0,1,by=0.01),nrclusters=NULL,names=NULL,nboot=10,StopRange=FALSE,plottype="new",location=NULL){

	type=match.arg(type)

	if(is.null(names)){
		names=c()
		for(i in 1:length(List)){
			names=c(names,paste("Method",i,sep=" "))
		}
	}

	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE  & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}


	if(type=="data"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}
		Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i],normalize[i],method[i]))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))

		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters[i])$silinfo$widths)
		names(silwidth)=names

	}
	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}
		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters[i])$silinfo$widths)
		names(silwidth)=names
	}
	else{
		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))
		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}
		silwidth=lapply(Dist,function(x) cluster::pam(x,nrclusters[i])$silinfo$widths)
		names(silwidth)=names
	}


	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){

		}
	}

	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}


	namesw=c()
	for(i in 1:length(names)){
		namesw=c(namesw,paste("w_",names[i],sep=""))
	}

	labels<-c(namesw,"Observed Statistic","P-Value")

	ResultsWeight<-matrix(0,ncol=length(labels),nrow=length(weight))
	colnames(ResultsWeight)=labels

	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))
	}
	else if(inherits(weight,"list") & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}

	if(!inherits(weight,"list")){
		condition<-function(l){
			l=as.numeric(l)
			if( sum(l)==1 ){  #working wit characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		weight=t2[which(t3!=0)]
	}

	if(inherits(weight,"list") & "x" %in% weight[[1]]){ #x indicates a free weight
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
	}

	DistW=lapply(weight,weightedcomb,Dist)

	nrclus=ceiling(mean(nrclusters))

	silwidthW=lapply(DistW,function(x) cluster::pam(x,nrclus)$silinfo$widths)


	StatRSq<-function(silwidthW,silwidth,ordernames,names){
		n=length(silwidth)

		regressRSq<-function(silwidthweight,silwidth,ordernames,names){

			L1=silwidthweight[,3][ordernames]
			L2W=silwidthweight[,1][ordernames]

			regressWW<-stats::lm(L1~L2W)
			RsqWW<-summary(regressWW)$r.squared

			regressWX<-function(silw,L1,ordernames){
				L2<-silw[,1][ordernames]

				regresWX<-stats::lm(L1~L2)
				RsqWX<-summary(regresWX)$r.squared
				return(RsqWX)
			}

			RsqWX<-sapply(silwidth,regressWX,L1=L1,ordernames=ordernames)

			Rsq=c(RsqWW,RsqWX)
			return(Rsq)
		}


		RSqs=lapply(c(1:length(silwidthW)),function(x) regressRSq(silwidthW[[x]],silwidth,ordernames,names))

		statfunction<-function(RS){
			stat=0
			xx=0
			xy=0
			for(i in 1:length(RS)){
				if(i==1){
					xx=xx+RS[i]
				}
				else{
					xy=xy+RS[i]
				}
			}
			stat=abs(n*xx-xy)  #check this formula with Nolen
			names(stat)=NULL
			return(stat)
		}

		Stats=sapply(RSqs,statfunction)

	}

	StatRSqObs=StatRSq(silwidthW,silwidth,ordernames=rownames(Dist[[1]]),names)

	#bootstrapping
	statNULL=matrix(0,nrow=length(weight),ncol=nboot)
	perm.rowscols <- function (D, n)
	{
		s <- sample(1:n)
		D=D[s, s]
		return(D)
	}

	for(i in 1:nboot){
		set.seed(i)
		DistNULL=Dist
		DistNULL[[1]] <- perm.rowscols(DistNULL[[1]],nrow(DistNULL[[1]]))
		silwidthNULL=lapply(1:length(DistNULL),function(x) cluster::pam(DistNULL[[x]],nrclusters[x])$silinfo$widths)


		DistWNULL=lapply(weight,weightedcomb,DistNULL)
		silwidthWNULL=lapply(DistWNULL,function(x) cluster::pam(x,nrclus)$silinfo$widths)

		statNULL[,i]=StatRSq(silwidthWNULL,silwidthNULL,ordernames=rownames(DistNULL[[1]]),names)
	}

	PVals=lapply(c(1:nrow(statNULL)),function(x) (1+sum(abs(statNULL[x,])<=abs(StatRSqObs[x])))/(nboot+1))

	ResultsWeight=t(mapply(c,weight,StatRSqObs,PVals))
	colnames(ResultsWeight)=labels

	#Choose weight with smallest observed test statistic

	Weight=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),c(1:length(List))]

	plottypein(plottype,location)
	graphics::plot(x=ResultsWeight[,1],y=ResultsWeight[,"Observed Statistic"],xlim=c(0,max(ResultsWeight[,1])),ylim=c(min(ResultsWeight[,"Observed Statistic"]),max(ResultsWeight[,"Observed Statistic"])),xlab="",ylab="Observed Statistic",pch=19,col="black")
	graphics::points(ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),"Observed Statistic"],pch=19,col="red")
	graphics::mtext("Weight Combinations", side=1, line=4)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],line=2)
	plottypeout(plottype)

	plottypein(plottype,location)
	graphics::plot(x=ResultsWeight[,1],y=ResultsWeight[,"P-Value"],xlim=c(0,max(ResultsWeight[,1])),ylim=c(min(ResultsWeight[,"P-Value"]),max(ResultsWeight[,"P-Value"])),xlab="",ylab="P-Value",pch=19,col="black")
	graphics::points(ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),"P-Value"],pch=19,col="red")
	graphics::mtext("Weight Combinations", side=1, line=4)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=ResultsWeight[which.min(abs(ResultsWeight[,3]-0)),1],line=2)
	plottypeout(plottype)


	out=list()
	out[[1]]=ResultsWeight
	out[[2]]=Weight
	names(out)=c("Result","Weight")

	return(out)

}

#' @title Determines an optimal weight for weighted clustering by similarity weighted clustering.
#'
#' @description The function \code{DetermineWeight_SimClust} determines an optimal weight
#' for performing weighted similarity clustering. For each given weight, each separate clustering
#' is compared to the clustering on a weighted dissimilarity matrix and a Jaccard coefficient is
#' calculated. The ratio of the Jaccard coefficients closest to one indicates an optimal weight.
#' @export DetermineWeight_SimClust
#' @param List A list of matrices of the same type. It is assumed the rows are corresponding with the objects.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. Type should be one of "data", "dist" or "clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Defaults to c("tanimoto","tanimoto").
#' @param normalize Logical. Indicates whether to normalize the distance matrices. Default is c(FALSE, FALSE).
#' @param method A method of normalization. Default is c(NULL,NULL).
#' @param weight Optional. A list of different weight combinations for the data sets in List. Defaults to seq(0,1,by=0.01).
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param clust Choice of clustering function (character). Defaults to "agnes".
#' @param linkage Choice of inter group dissimilarity for the individual clusterings. Defaults to c("flexible","flexible").
#' @param linkageF Choice of inter group dissimilarity for the final clustering. Defaults to "ward".
#' @param alpha The parameter alpha to be used in the "flexible" linkage of the agnes function. Defaults to 0.625.
#' @param gap Logical. Whether or not to calculate the gap statistic. Only if type="data". Default is FALSE.
#' @param maxK The maximal number of clusters to consider in calculating the gap statistic. Default is 15.
#' @param names The labels to give to the elements in List. Default is NULL.
#' @param StopRange Logical. Indicates whether the distance matrices with values not between zero and one
#' should be standardized to have so. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "new".
#' @param location If plottype is "pdf", a location should be provided here. Default is NULL.
#' @return The returned value is a list with three elements: \item{ClustSep}{The result of \code{Cluster}
#' for each single element of List} \item{Result}{A data frame with the Jaccard coefficients and their
#' ratios for each weight} \item{Weight}{The optimal weight}.
#' @references Perualila-Tan, N. et al. (2016). Weighted similarity-based clustering of chemical
#' structures and bioactivity data in early drug discovery. Journal of Bioinformatics and Computational Biology.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",alpha=0.625,gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#'
#' MCF7_Weight=DetermineWeight_SimClust(List=L,type="clusters",weight=seq(0,1,by=0.01),
#' nrclusters=c(7,7),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),
#' method=c(NULL,NULL),clust="agnes",linkage=c("flexible","flexible"),linkageF="ward",
#' alpha=0.625,gap=FALSE,maxK=50,names=c("FP","TP"),StopRange=FALSE,plottype="new",location=NULL)
#' }
DetermineWeight_SimClust<-function(List,type=c("data","dist","clusters"),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),weight=seq(0,1,by=0.01),nrclusters=NULL,clust="agnes",linkage=c("flexible","flexible"),linkageF="ward",alpha=0.625,gap=FALSE,maxK=15,names=NULL,StopRange=FALSE,plottype="new",location=NULL){

	CheckDist<-function(Dist,StopRange){
		if(StopRange==FALSE & !(0<=min(Dist) & max(Dist)<=1)){
			message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
			Dist=Normalization(Dist,method="Range")
		}
		else{
			Dist=Dist
		}
	}


	type<-match.arg(type)
	if(type=="data"){

		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange))

		Dist=lapply(seq(length(List)),function(i) Clusterings[[i]]$DistM)
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))

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
		out<-list(Clusterings)
		names(out)="ClusterSep"
	}

	else if(type=="dist"){
		OrderNames=rownames(List[[1]])
		for(i in 1:length(List)){
			List[[i]]=List[[i]][OrderNames,OrderNames]
		}

		Clusterings=lapply(seq(length(List)),function(i) Cluster(List[[i]],type,distmeasure[i],normalize[i],method[i],clust,linkage[i],alpha,gap,maxK,StopRange))

		Dist=List
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))

		if(is.null(nrclusters)){
			if(gap==FALSE){
				stop("Specify a number of clusters or put gap to TRUE")
			}
			else{
				clusters=sapply(seq(length(List)),function(i) Clusterings[[i]]$k$Tibs2001SEmax)
				nrclusters=ceiling(mean(clusters))
			}
		}

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"

	}
	else{

		Dist=lapply(seq(length(List)),function(i) return(List[[i]]$DistM))
		Dist=lapply(seq(length(Dist)),function(i) CheckDist(Dist[[i]],StopRange))

		OrderNames=rownames(Dist[[1]])
		for(i in 1:length(Dist)){
			Dist[[i]]=Dist[[i]][OrderNames,OrderNames]
		}

		Clusterings=List

		for(i in 1:length(Clusterings)){
			names(Clusterings)[i]=paste("Clust",i,sep=' ')
		}
		out<-list(Clusterings)
		names(out)="ClusterSep"

	}


	namesw=c()
	for(i in 1:length(names)){
		namesw=c(namesw,paste("w_",names[i],sep=""))
	}
	namesJ=c()
	for(i in 1:length(names)){
		namesJ=c(namesJ,paste("J(sim",names[i],",simW)",sep=""))
	}
	namesR=c()
	combs=utils::combn(seq(length(List)),m=2,simplify=FALSE)
	for(i in 1:length(combs)){
		namesR=c(namesR,paste("J_",names[combs[[i]][1]],"/J_",names[combs[[i]]][2],sep=""))
	}

	labels<-c(namesw,namesJ,namesR)

	ResultsWeight<-matrix(0,ncol=length(labels),nrow=length(weight))
	#data.frame(col1=numeric(),col2=numeric(),col3=numeric(),col4=numeric())
	colnames(ResultsWeight)=labels

	if(is.null(weight)){
		equalweights=1/length(List)
		weight=list(rep(equalweights,length(List)))
	}
	else if(inherits(weight,"list") & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}

	if(!inherits(weight,"list")){
		condition<-function(l){
			l=as.numeric(l)
			if( sum(l)==1 ){  #working wit characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		t1=gtools::permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		weight=t2[which(t3!=0)]
	}

	if(inherits(weight,"list") & "x" %in% weight[[1]]){ #x indicates a free weight
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
	}
	DistM=lapply(weight,weightedcomb,Dist)

	hclustOr=lapply(seq(length(List)),function(i) stats::cutree(Clusterings[[i]]$Clust,nrclusters[i]))
	nrclus=ceiling(mean(nrclusters))
	hclustW=lapply(seq(length(weight)),function(i) stats::cutree(cluster::agnes(DistM[[i]],diss=TRUE,method=linkageF),nrclus))

	Counts=function(clusterlabs1,clusterlabs2){
		index=c(1:length(clusterlabs1))
		allpairs=utils::combn(index,2,simplify=FALSE)  #all pairs of indices: now check clutserlabels for every pair==> only 1 for loop
		n11=n10=n01=n00=0

		counts<-function(pair){
			if(clusterlabs1[pair[1]]==clusterlabs1[pair[2]]){
				if(clusterlabs2[pair[1]]==clusterlabs2[pair[2]]){
					n11=n11+1
				}
				else{
					n10=n10+1
				}
			}
			else{
				if(clusterlabs2[pair[1]]==clusterlabs2[pair[2]]){
					n01=n01+1
				}
				else{
					n00=n00+1
				}

			}
			return(c(n11,n10,n01,n00))
		}

		n=lapply(seq(length(allpairs)),function(i) counts(allpairs[[i]]))
		nn=Reduce("+",n)
		#2: compute jaccard coefficient
		Jac=nn[1]/(nn[1]+nn[2]+nn[3])
		return(Jac)
	}

	Jaccards<-function(hclust){
		jacs=lapply(seq(length(hclustOr)),function(i) Counts(clusterlabs1=hclustOr[[i]],clusterlabs2=hclust))
		return(unlist(jacs))
	}

	AllJacs=lapply(hclustW,Jaccards)  #make this faster:lapply + transfrom to data frame with plyr package

	Ratios<-function(Jacs){
		combs=utils::combn(seq(length(List)),m=2,simplify=FALSE)
		ratio<-function(v,Jacs){
			return(Jacs[v[1]]/Jacs[v[2]])
		}

		ratios=lapply(seq(length(combs)),function(i) ratio(v=combs[[i]],Jacs=Jacs))

	}

	AllRatios=lapply(seq(length(AllJacs)),function(i) unlist(Ratios(AllJacs[[i]])))


	ResultsWeight=t(mapply(c,weight,AllJacs,AllRatios))
	colnames(ResultsWeight)=labels

	#Choose weight with ratio closest to one==> smallest where this happens: ##### START HERE WITH OPTIMIZATION #####
	ResultsWeight=cbind(ResultsWeight,rep(0,nrow(ResultsWeight)))
	colnames(ResultsWeight)[ncol(ResultsWeight)]="trick"
	Weight=ResultsWeight[which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),c(1:length(List))]

	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
		if(plottype=="sweave"){

		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}

	plottypein(plottype,location)
	graphics::plot(x=0,y=0,type="n",xlim=c(0,dim(ResultsWeight)[1]),ylim=c(min(ResultsWeight[,namesR]),max(ResultsWeight[,namesR])),xlab="",ylab="Ratios")
	if(is.null(ncol(ResultsWeight[,namesR]))){
		L=1
	}
	else{
		L=ncol(ResultsWeight[,namesR])
	}
	for(i in 1:L){
		graphics::lines(1:dim(ResultsWeight)[1],y=ResultsWeight[,namesR[i]],col=i)
	}
	graphics::abline(h=0,v=which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),col="black",lwd=2)
	graphics::mtext("Weight Combinations", side=1, line=3)
	graphics::axis(1,labels=paste("Optimal weights:", paste(Weight,collapse=", "),sep=" "), at=which.min(rowSums(abs(ResultsWeight[,c(namesR,"trick")]-1))),line=1,tck=1,lwd=2)
	plottypeout(plottype)

	ResultsWeight=ResultsWeight[,-ncol(ResultsWeight)]
	out[[2]]=ResultsWeight
	out[[3]]=Weight
	names(out)=c("ClusterSep","Result","Weight")



	return(out)

}
