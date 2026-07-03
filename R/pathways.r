#' @title Find differentially expressed genes
#'
#' @description The \code{DiffGenes} function looks for differentially expressed
#' genes for each cluster of a clustering result by comparing the objects of a
#' cluster against all other objects. The limma method is applied to find the
#' differentially expressed genes.
#' @export
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param Selection If differential gene expression should be investigated for
#' a specific selection of objects, this selection can be provided here.
#' Selection can be of the type "character" (names of the objects) or
#' "numeric" (the number of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' Default is NULL.
#' @param method The method to applied to look for DE genes. For now, only the
#' limma method is available. Default is "limma".
#' @param sign The significance level to be handled. Default is 0.05.
#' @param topG Overrules sign. The number of top genes to be shown. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The result is a list with an element per method. Each element is
#' again a list with for each cluster the found differentially expressed genes.
#' @details This function relies on the suggested Bioconductor package 'limma'.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#'
#' MCF7_FT_DE = DiffGenes(List=L,geneExpr=geneMat,nrclusters=7,method="limma",sign=0.05,
#' topG=10,fusionsLog=TRUE,weightclust=TRUE,names=names)
#' }
DiffGenes=function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method="limma",sign=0.05,topG=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if (!requireNamespace("limma", quietly = TRUE)) stop("DiffGenes() requires the suggested package 'limma' (Bioconductor).")
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	}
	if(!is.null(Selection)){
		ResultLimma=DiffGenesSelection(List,Selection,geneExpr,nrclusters,method,sign,topG,fusionsLog,weightclust,names)
	}
	else{

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

		if(is.null(names)){
			for(j in 1:length(List)){
				names[j]=paste("Method",j,sep=" ")
			}
		}

		names(ListNew)=names

		if(is.null(topG)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}


		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
		List=ListNew
		ResultLimma=list()
		maxclus=0
		for (k in 1:dim(MatrixClusters)[1]){
			clusters=MatrixClusters[k,]
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			Genes=list()
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			for (i in clust){

				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)]
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds") #names of the objects

				label = rep(0,length(names(clusters)))
				label[which(clusters==i)] = 1

				label.factor = factor(label)
				GeneExpr.2=geneExpr[,names(clusters)]

				if(class(GeneExpr.2)[1]=="ExpressionSet"){

					if (!requireNamespace("a4Base", quietly = TRUE)) {
						stop("a4Base needed for this function to work. Please install it.",
								call. = FALSE)
					}


					GeneExpr.2$LeadCmpds<-label.factor
					DElead <- a4Base::limmaTwoLevels(GeneExpr.2,"LeadCpds")

					allDE <-a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")

					if(is.null(allDE$ID)){
						allDE$ID<- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}

					if(top1==TRUE){
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")

					}
					else if(top1==FALSE){
						topG=length(which(allDE$adj.P.Val<=sign))
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")

					}

				}
				else{

					design = stats::model.matrix(~label.factor)
					fit = limma::lmFit(GeneExpr.2,design=design)
					fit = limma::eBayes(fit)
					allDE=limma::topTable(fit,n=dim(geneExpr)[1],coef=2,adjust="fdr",sort.by="P")

					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}

					if(top1==TRUE){
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
					}
					else if(top1==FALSE){
						topG=length(which(allDE$adj.P.Val<=sign))
						result = list(allDE[1:topG,],allDE)
						names(result)=c("TopDE","AllDE")
					}
				}

				temp[[2]]=result

				names(temp)=c("objects","Genes")

				Genes[[i]]=temp

				names(Genes)[i]=paste("Cluster",i,sep=" ")
			}
			ResultLimma[[k]]=Genes

		}
		names(ResultLimma)=names
		for(i in 1:length(ResultLimma)){
			for(k in 1:length(ResultLimma[[i]])){
				if(is.null(ResultLimma[[i]][[k]])[1]){
					ResultLimma[[i]][[k]]=NA
					names(ResultLimma[[i]])[k]=paste("Cluster",k,sep=" ")
				}
			}
			if(length(ResultLimma[[i]]) != maxclus){
				extra=maxclus-length(ResultLimma[[i]])
				#temp=length(ResultLimma[[i]])
				for(j in 1:extra){
					ResultLimma[[i]][[length(ResultLimma[[i]])+1]]=NA
					names(ResultLimma[[i]])[length(ResultLimma[[i]])]=paste("Cluster",length(ResultLimma[[i]]),sep=" ")
				}
			}
		}

	}
	return(ResultLimma)
}

#' @title Differential expression for a selection of objects
#' @param List A list of the clustering outputs to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param Selection If differential gene expression should be investigated for
#' a specific selection of objects, this selection can be provided here.
#' Selection can be of the type "character" (names of the objects) or
#' "numeric" (the number of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for DE genes. For now, only the
#' limma method is available. Default is "limma".
#' @param sign The significance level to be handled. Default is 0.05.
#' @param topG Overrules sign. The number of top genes to be shown. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @description Internal function of \code{DiffGenes}.
#' @export
#' @examples
#' \dontrun{
#' MCF7_FT_DE = DiffGenes(List=L,Selection=c("Cpd1","Cpd2"),geneExpr=geneMat,method="limma")
#' }
DiffGenesSelection=function(List,Selection,geneExpr=NULL,nrclusters=NULL,method="limma",sign=0.05,topG=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if (!requireNamespace("limma", quietly = TRUE)) stop("DiffGenesSelection() requires the suggested package 'limma' (Bioconductor).")
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	}

	if(is.null(topG)){
		top1=FALSE
	}
	else{
		top1=TRUE
	}

	if(inherits(Selection,"character")){
		ResultLimma=list()
		Genes=list()
		temp=list()

		LeadCpds=Selection #names of the objects
		OrderedCpds=colnames(geneExpr)
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")

		label = rep(0,dim(geneExpr)[2])
		label[which(colnames(geneExpr)%in%Selection)] = 1
		label.factor = factor(label)


		if(class(geneExpr)[1]=="ExpressionSet"){

			if (!requireNamespace("a4Base", quietly = TRUE)) {
				stop("a4Base needed for this function to work. Please install it.",
						call. = FALSE)
			}

			geneExpr$LeadCmpds<-label.factor
			DElead <- a4Base::limmaTwoLevels(geneExpr,"LeadCpds")

			allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
			if(is.null(allDE$ID)){
				allDE$ID <- rownames(allDE)
			}
			else
			{
				allDE$ID=allDE$ID
			}
			if(top1==TRUE){
				result = list(allDE[1:topG,],allDE)
				names(result)=c("TopDE","AllDE")

			}
			else if(top1==FALSE){
				topG=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")

			}

		}
		else{

			design = stats::model.matrix(~label.factor)
			fit = limma::lmFit(geneExpr,design=design)
			fit = limma::eBayes(fit)

			allDE=limma::topTable(fit,coef=2,n=dim(geneExpr)[1],adjust="fdr",sort.by="P")
			if(is.null(allDE$ID)){
				allDE$ID <- rownames(allDE)
			}
			else
			{
				allDE$ID=allDE$ID
			}
			if(top1==TRUE){
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")

			}
			else if(top1==FALSE){
				topG=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[0:topG,],allDE)
				names(result)=c("TopDE","AllDE")

			}

		}
		temp[[2]]=result

		names(temp)=c("objects","Genes")
		ResultLimma[[1]]=temp
		names(ResultLimma)="Selection"

	}
	else if(inherits(Selection,"numeric") & !(is.null(List))){

		ListNew=list()
		element=0

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

		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)

		List=ListNew
		ResultLimma=list()
		for(k in 1:dim(Matrix)[1]){
			cluster=Selection

			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]

			Genes=list()
			temp=list()
			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the objects
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")

			label = rep(0,dim(Matrix)[2])
			label[which(Matrix[k,]==cluster)] = 1
			label.factor = factor(label)

			GeneExpr.2=geneExpr[,colnames(Matrix)]

			if(class(GeneExpr.2)[1]=="ExpressionSet"){

				if (!requireNamespace("a4Base", quietly = TRUE)) {
					stop("a4Base needed for this function to work. Please install it.",
							call. = FALSE)
				}

				GeneExpr.2$LeadCmpds<-label.factor
				DElead <- a4Base::limmaTwoLevels(GeneExpr.2,"LeadCpds")

				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(top1==TRUE){
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")

				}
				else if(top1==FALSE){
					topG=length(which(allDE$adj.P.Val<=sign))
					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")

				}

			}
			else{


				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(GeneExpr.2,design=design)
				fit = limma::eBayes(fit)

				allDE=limma::topTable(fit,coef=2,n=dim(geneExpr)[1],adjust="fdr",sort.by="P")
				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}
				else
				{
					allDE$ID=allDE$ID
				}
				if(top1==TRUE){
					result = list(allDE[1:topG,],allDE)
					names(result)=c("TopDE","AllDE")

				}
				else if(top1==FALSE){
					topG=length(which(allDE$adj.P.Val<=sign))

					result = list(allDE[0:topG,],allDE)
					names(result)=c("TopDE","AllDE")

				}

			}
			temp[[2]]=result

			names(temp)=c("objects","Genes")
			ResultLimma[[k]]=temp
			names(ResultLimma)[k]=paste(names[k],": Cluster", cluster, sep="")
		}
	}

	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultLimma)

}

#' @title Find the genes shared across methods
#'
#' @description The \code{FindGenes} function collects the differentially
#' expressed genes per cluster per method (the output of \code{DiffGenes}) and
#' determines which genes are found across the methods.
#' @export
#' @param dataLimma The output of the \code{DiffGenes} function.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The result is a list with for each cluster a list of the found genes
#' and the methods that found them.
#' @examples
#' \dontrun{
#' MCF7_SharedGenes=FindGenes(dataLimma=MCF7_FT_DE,names=c("FP","TP"))
#' }
FindGenes<-function(dataLimma,names=NULL){

	FoundGenes=list()

	if(is.null(names)){
		for(j in 1:length(dataLimma)){
			names[j]=paste("Method",j,sep=" ")
		}
	}

	nrclusters=length(dataLimma[[1]])


	for(j in 1:nrclusters){
		FoundGenes[[j]]=list()
	}

	for(i in 1:length(dataLimma)){ #i == method

		for(j in 1:nrclusters){ #j == cluster
			if(!(is.na(dataLimma[[i]][[j]]))[1]){

				tempgenes=dataLimma[[i]][[j]]$Genes$TopDE$ID


				if(!(is.null(tempgenes))){
					tempgenes=tempgenes[!(is.na(tempgenes))]
					for(k in 1:length(tempgenes)){
						if(!(tempgenes[k] %in% names(FoundGenes[[j]]))){
							#FoundGenes[[j]][[length(FoundGenes[[j]])+1]]=c()
							FoundGenes[[j]][[length(FoundGenes[[j]])+1]]=names[i]
							names(FoundGenes[[j]])[length(FoundGenes[[j]])]=tempgenes[k]

						}
						else if (tempgenes[k] %in% names(FoundGenes[[j]])){
							found=which(names(FoundGenes[[j]])==tempgenes[k])
							FoundGenes[[j]][[found]]=c(FoundGenes[[j]][[found]],names[i])
						}
					}
				}
			}
		}
	}

	for(l in 1:nrclusters){
		namesl=names(FoundGenes[[l]])
		if(!(is.null(namesl))){

			for(m in l:nrclusters){
				namesm=names(FoundGenes[[m]])
				Templist=FoundGenes[[m]]
				if(length(namesm)!=0){
					if(l != m){

						for(k in 1:length(namesm)){

							if (namesm[k] %in% namesl){
								found=which(namesl==namesm[k])
								methods=Templist[[k]]
								del=which(names(FoundGenes[[m]])==namesm[k])
								FoundGenes[[m]][[del]]=c()
								for(a in 1:length(methods)){
									methods[a]=paste(methods[a],"_",m,sep="")
								}
								FoundGenes[[l]][[found]]=c(FoundGenes[[l]][[found]],methods)
							}

						}
					}
				}
			}
		}
	}


	for(i in 1:length(FoundGenes)){
		names(FoundGenes)[i]=paste("Cluster",i,sep=" ")
	}

	return(FoundGenes)
}

#' @title Intersection over resulting gene sets of \code{PathwaysIter} function
#'
#' @description The function \code{Geneset.intersect} collects the results of the
#' \code{PathwaysIter} function per method for each cluster and takes the
#' intersection over the iterations per cluster per method. This is to see if
#' over the different resamplings of the data, similar pathways were
#' discovered.
#' @export Geneset.intersect
#'
#' @param PathwaysOutput The output of the \code{PathwaysIter} function.
#' @param Selection Logical. Indicates whether or not the output of the
#' pathways function were concentrated on a specific selection of objects. If
#' this was the case, Selection should be put to TRUE. Otherwise, it should be
#' put to FALSE. Default is TRUE.
#' @param sign The significance level to be handled for cutting of the
#' pathways. Default is 0.05.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param seperatetables Logical. If TRUE, a separate element is created per
#' cluster containing the pathways for each iteration. Default is FALSE.
#' @param separatepvals Logical. If TRUE, the p-values of the each iteration of
#' each pathway in the intersection is given. If FALSE, only the mean p-value
#' is provided. Default is FALSE.
#' @return The output is a list with an element per method. For each method, it
#' is portrayed per cluster which pathways belong to the intersection over all
#' iterations and their corresponding mean p-values.
#' @examples
#'
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#'
#' MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr = geneMat, nrclusters = 7, method =
#' c("limma", "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE,
#' weightclust = TRUE, names =names)
#'
#' MCF7_Paths_intersection=Geneset.intersect(PathwaysOutput=MCF7_Paths_FandT,
#' sign=0.05,names=c("FP","TP"),seperatetables=FALSE,separatepvals=FALSE)
#'
#' str(MCF7_Paths_intersection)
#' }
Geneset.intersect<-function(PathwaysOutput,Selection=FALSE,sign=0.05,names=NULL,seperatetables=FALSE,separatepvals=FALSE){

	if (!requireNamespace("plyr", quietly = TRUE)) stop("Geneset.intersect() requires the suggested package 'plyr'.")

	if(Selection==TRUE){
		if(length(PathwaysOutput$'Iteration 1')==1){
			names="Selection"
		}
		Intersect=Geneset.intersectSelection(PathwaysOutput,sign,names,seperatetables,separatepvals)
	}

	else{

		if(is.null(names)){
			for(j in 1:length(PathwaysOutput$"Iteration 1")){
				names[j]=paste("Method",j,sep=" ")
			}
		}


		#put all of same method together:preparation of lists
		subsets=list()
		nmethods=length(PathwaysOutput$"Iteration 1")
		for(i in 1:nmethods){
			subsets[[i]]=list()

		}
		names(subsets)=names

		#put all of same method together: go through PathwaysOutput
		for(j in 1:length(PathwaysOutput)){
			name1=names(PathwaysOutput)[j]
			for(k in 1:nmethods){
				name2=names[k]
				subsets[[name2]][[name1]]=PathwaysOutput[[name1]][[name2]]

			}

		}

		#for every subset (= every method) take intersection over the interations per cluster
		Intersect=list()

		for(i in 1:length(subsets)){
			Method=subsets[[i]]
			Clusters=list()
			nclus=length(Method[[1]])
			for(j in 1:length(Method)){
				name3=paste("Iteration",j,sep=" ")
				for(k in 1:nclus){
					name4=paste("Cluster",k,sep=" ")
					if(!(is.na(Method[[name3]][[name4]])[1])){
						Clusters[[name4]][[name3]]=Method[[name3]][[name4]]
					}
					else{
						Clusters[[name4]][[name3]]=NA
					}
				}

			}

			IntersectM=list()

			for(a in 1:length(Clusters)){ #per cluster
				if(!(is.na(Clusters[[a]])[1])){
					result.out=list()
					result.name = c()
					for(b in 1:length(Clusters[[a]])){#per iteration
						if(b==1){
							objects=Clusters[[a]][[1]]$objects
							Genes=Clusters[[a]][[1]]$Genes
							Names=data.frame("description"=Clusters[[a]][[1]]$Pathways$AllPaths$geneSetDescription,"genesetcode"=rownames(Clusters[[a]][[1]]$Pathways$AllPaths))
							Names$description=as.character(Names$description)
							Names$genesetcode=as.character(Names$genesetcode)
						}
						cut = Clusters[[a]][[b]]$Pathways$AllPaths[Clusters[[a]][[b]]$Pathways$AllPaths$geneSetPValue<=sign,]
						colnames(cut)[4] = paste("pvalues.",b,sep="")
						colnames(cut)[2] = paste("testedgenesetsize.",b,sep="")
						colnames(cut)[3] = paste("genesetstatistic.",b,sep="")
						cut=cut[,c(1,5,2,3,4)]
						result.out[[b]] = cut
						result.name = c(result.name,paste("genesettable",b,sep=""))

					}



					names(result.out) = result.name

					genesets.table.intersect = plyr::join_all(result.out,by=c("totalGeneSetSize","geneSetDescription"),type="inner")
					genesets.table.intersect$mean_testedGeneSetSize=round(apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='testedgenesetsize')],1,mean),1)
					genesets.table.intersect$mean_geneSetStatistic=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='genesetstatistic')],1,mean)
					genesets.table.intersect$mean_geneSetPValue=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='pvalues')],1,mean)

					rownames(genesets.table.intersect)=as.character(Names[which(genesets.table.intersect$geneSetDescription%in%Names[,1]),2])

					class(genesets.table.intersect)=c("MLP","data.frame")
					attr(genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]][[1]]$Pathways$AllPaths)$geneSetSource


					result.out$genesets.table.intersect = genesets.table.intersect



					if(separatepvals==FALSE){
						result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,(ncol(genesets.table.intersect)-2):ncol(genesets.table.intersect))]
						class(result.out$genesets.table.intersect)=c("MLP","data.frame")
						attr(result.out$genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]][[1]]$Pathways$AllPaths)$geneSetSource
					}


					if(seperatetables==FALSE){
						result.out=result.out$genesets.table.intersect
						class(result.out)=c("MLP","data.frame")
						attr(result.out,'geneSetSource')=attributes(Clusters[[1]][[1]]$Pathways$AllPaths)$geneSetSource
					}

					newresult=list(objects=objects,Genes=Genes,Pathways=result.out)


					IntersectM[[a]]=newresult
					names(IntersectM)[a]=names(Clusters)[[a]]
				}
				else{
					IntersectM[[a]]=NA
					names(IntersectM)[a]=names(Clusters)[[a]]
				}
			}

			Intersect[[i]]=IntersectM

		}
	}
	names(Intersect)=names
	return(Intersect)
}

#' @title Intersection over resulting gene sets of \code{PathwaysIter} function for a selection of objects
#' @param list.output The output of the \code{PathwaysIter} function.
#' @param sign The significance level to be handled for cutting of the
#' pathways. Default is 0.05.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param seperatetables Logical. If TRUE, a separate element is created per
#' cluster containing the pathways for each iteration. Default is FALSE.
#' @param separatepvals Logical. If TRUE, the p-values of the each iteration of
#' each pathway in the intersection is given. If FALSE, only the mean p-value
#' is provided. Default is FALSE.
#' @description Internal function of \code{Geneset.intersect}.
#' @export
#' @examples
#' \dontrun{
#' Geneset.intersectSelection(MCF7_Paths_FandT,sign=0.05)
#' }
Geneset.intersectSelection<-function(list.output,sign=0.05,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	if (!requireNamespace("plyr", quietly = TRUE)) stop("Geneset.intersectSelection() requires the suggested package 'plyr'.")
	if(is.null(names)){
		for(j in 1:length(list.output$"Iteration 1")){
			names[j]=paste("Method",j,sep=" ")
		}
	}

	#put all of same method together:preparation of lists
	subsets=list()
	nmethods=length(list.output$"Iteration 1")
	for(i in 1:nmethods){
		subsets[[i]]=list()
	}
	names(subsets)=names

	#put all of same method together: go through list.output
	for(j in 1:length(list.output)){
		name1=names(list.output)[j]
		for(k in 1:nmethods){
			name2=k
			subsets[[name2]][[name1]]=list.output[[name1]][[name2]]

		}

	}

	#for every subset (= every method) take intersection over the interations per cluster
	Intersect=list()

	for(i in 1:length(subsets)){
		Method=subsets[[i]]
		Clusters=list()
		nclus=1
		for(j in 1:length(Method)){
			name3=paste("Iteration",j,sep=" ")
			Clusters[[name3]]=Method[[name3]]
		}

		IntersectM=list()

		result.out=list()
		result.name = c()
		for(a in 1:length(Clusters)){ #per cluster
			if(a==1){
				objects=Clusters[[a]]$objects
				Genes=Clusters[[a]]$Genes
				Names=data.frame("description"=Clusters[[a]]$Pathways$AllPaths$geneSetDescription,"genesetcode"=rownames(Clusters[[a]]$Pathways$AllPaths))
				Names$description=as.character(Names$description)
				Names$genesetcode=as.character(Names$genesetcode)
			}
			cut = Clusters[[a]]$Pathways$AllPaths[Clusters[[a]]$Pathways$AllPaths$geneSetPValue<=sign,]
			colnames(cut)[4] = paste("pvalues.",a,sep="")
			colnames(cut)[2] = paste("testedgenesetsize.",a,sep="")
			colnames(cut)[3] = paste("genesetstatistic.",a,sep="")
			cut=cut[,c(1,5,2,3,4)]
			result.out[[a]] = cut
			result.name = c(result.name,paste("genesettable",a,sep=""))
		}

		names(result.out) = result.name

		genesets.table.intersect = plyr::join_all(result.out,by=c("totalGeneSetSize","geneSetDescription"),type="inner")
		genesets.table.intersect$mean_testedGeneSetSize=round(apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='testedgenesetsize')],1,mean),1)
		genesets.table.intersect$mean_geneSetStatistic=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='genesetstatistic')],1,mean)
		genesets.table.intersect$mean_geneSetPValue=apply(genesets.table.intersect[,which(substring(colnames(genesets.table.intersect),1,nchar(colnames(genesets.table.intersect))-nchar(".1"))=='pvalues')],1,mean)

		rownames(genesets.table.intersect)=as.character(Names[which(genesets.table.intersect$geneSetDescription%in%Names[,1]),2])

		class(genesets.table.intersect)=c("MLP","data.frame")
		attr(genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource


		result.out$genesets.table.intersect = genesets.table.intersect

		if(separatepvals==FALSE){
			result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,(ncol(genesets.table.intersect)-2):ncol(genesets.table.intersect))]
			class(result.out$genesets.table.intersect)=c("MLP","data.frame")
			attr(result.out$genesets.table.intersect,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		}


		if(seperatetables==FALSE){
			result.out=result.out$genesets.table.intersect
			class(result.out)=c("MLP","data.frame")
			attr(result.out,'geneSetSource')=attributes(Clusters[[1]]$Pathways$AllPaths)$geneSetSource
		}


		newresult=list(objects=objects,Genes=Genes,Pathways=result.out)


		#IntersectM[[a]]=
		#names(IntersectM)[a]=names(Clusters)[[a]]
		Intersect[[i]]=newresult

	}
	names(Intersect)=names
	return(Intersect)
}

#' @title Pathway analysis with intersection over iterations
#'
#' @description The \code{PathwayAnalysis} function performs pathway analysis
#' multiple times via \code{PathwaysIter} and takes the intersection over the
#' iterations via \code{Geneset.intersect}.
#' @export
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in. Default is NULL.
#' @param method The methods for gene and pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function, defaults to "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each cluster. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned. Default is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param niter The number of times to perform pathway analysis. Default is 10.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param seperatetables Logical. Default is FALSE.
#' @param separatepvals Logical. Default is FALSE.
#' @return The intersection over the iterations of the pathway analysis.
#' @details This function relies on the suggested Bioconductor packages 'MLP', 'biomaRt' and 'org.Hs.eg.db'.
#' @examples
#' \dontrun{
#' MCF7_PathsFandT=PathwayAnalysis(List=L, geneExpr = geneMat, nrclusters = 7,
#' method = c("limma","MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP",niter=2)
#' }
PathwayAnalysis<-function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,niter=10,fusionsLog=TRUE,weightclust=TRUE,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}

	Pathways=PathwaysIter(List,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign,niter,fusionsLog,weightclust,names)

	if(is.null(Selection)){
		Selection=FALSE
	}
	else{
		Selection=TRUE
	}

	Intersection=Geneset.intersect(PathwaysOutput=Pathways,Selection,sign,names,seperatetables,separatepvals)

	return(Intersection)

}

#' @title Pathway analysis for multiple clustering results
#'
#' @description A pathway analysis per cluster per method is conducted.
#' @export Pathways
#' @details After finding differently expressed genes, it can be investigated whether
#' pathways are related to those genes. This can be done with the help of the
#' function \code{Pathways} which makes use of the \code{MLP} function of the
#' MLP package. Given the output of a method, the cutree function is performed
#' which results into a specific number of clusters. For each cluster, the
#' limma method is performed comparing this cluster to the other clusters. This
#' to obtain the necessary p-values of the genes. These are used as the input
#' for the \code{MLP} function to find interesting pathways. By default the
#' candidate gene sets are determined by the \code{AnnotateEntrezIDtoGO}
#' function. The default source will be GOBP, but this can be altered.
#' Further, it is also possible to provide own candidate gene sets in the form
#' of a list of pathway categories in which each component contains a vector of
#' Entrez Gene identifiers related to that particular pathway. The default
#' values for the minimum and maximum number of genes in a gene set for it to
#' be considered were used. For MLP this is respectively 5 and 100. If a list
#' of outputs of several methods is provided as data input, the cluster numbers
#' are rearranged according to a reference method. The first method is taken as
#' the reference and ReorderToReference is applied to get the correct ordering.
#' When the clusters haven been re-appointed, the pathway analysis as described
#' above is performed for each cluster of each method.
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function, defaults to
#' "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Defaults is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return The returned value is a list with an element per cluster per method.
#' @details This function relies on the suggested Bioconductor packages 'MLP', 'biomaRt' and 'org.Hs.eg.db'.
#' @examples
#'
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#'
#' MCF7_PathsFandT=Pathways(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma",
#' "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
#' topG = NULL, GENESET = NULL, sign = 0.05,fusionsLog = TRUE, weightclust = TRUE,
#'  names =names)
#'  }
Pathways<-function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=NULL){

	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}



	if(!(is.null(Selection))){
		ResultMLP=PathwaysSelection(List,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign,fusionsLog,weightclust,names)

	}
	else if(inherits(List,"ChosenClusters")){
		ResultMLP=list()
		for(i in 1:length(List)){
			Selection=List[[i]]$objects$LeadCpds
			L=List[i]
			ResultMLP[[i]]=PathwaysSelection(List=L,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign,fusionsLog,weightclust,names)
			names(ResultMLP)=paste("Choice",i,sep=' ')
		}
	}
	else{

		#Check for gene expression data: if not, reordering to ListNew not necessary
		DataPrepared<-plyr::try_default(PreparePathway(List[[1]],geneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(DataPrepared)){

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

			maxclus=0
			DataPrepared=list()
			for (k in 1:dim(MatrixClusters)[1]){
				message(k)
				clusters=MatrixClusters[k,]

				if(max(clusters)>maxclus){
					maxclus=max(clusters)
				}

				check<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(check)){
					Temp=List[[k]]

					for(i in unique(clusters)){
						objects=list()
						objects$LeadCpds=names(clusters)[which(clusters==i)]
						objects$OrderedCpds=stats::as.hclust(List[[k]]$Clust$Clust)$labels[stats::as.hclust(List[[k]]$Clust$Clust)$order]

						Temp[[i+1]]=list(objects=objects)
						names(Temp)[i+1]=paste("Cluster",i,sep=" ")

					}
					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)

				}
			}
		}

		else{
			for(k in 1:length(List)){
				DataPrepared[[k]]<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(DataPrepared[[k]])){
					Temp=List[[k]]

					for(i in unique(clusters)){
						objects=list()
						objects$LeadCpds=List[[k]]$objects$LeadCpds
						objects$OrderedCpds=List[[k]]$objects$OrderedCpds

						Temp[[i+1]]=list(objects=objects)
						names(Temp)[i+1]=paste("Cluster",i,sep=" ")

					}
					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)

				}
			}


		}



		method.test = function(sign.method,path.method){
			method.choice = FALSE

			if( sign.method=="limma"  & path.method=="MLP"  ){
				method.choice = TRUE
			}
			if(method.choice==TRUE){
				return(list(sign.method=sign.method,path.method=path.method))
			}
			else{
				stop("Incorrect choice of method.")
			}

		}

		method.out = method.test(method[1],method[2])

		sign.method = method.out$sign.method
		path.method = method.out$path.method

		if(length(geneInfo$ENTREZID)==1){
			geneInfo$ENTREZID = colnames(geneExpr)
		}

		# Determining the genesets if they were not given with the function input
		if((inherits(GENESET,"geneSetMLP"))[1] ){
			geneSet <- GENESET
		}
		else{
			geneSet <- MLP::getGeneSets(species = "Human",geneSetSource = geneSetSource,entrezIdentifiers = geneInfo$ENTREZID)
		}

		if(is.null(topP)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}


		ResultMLP=list()

		for (k in 1:length(DataPrepared)){
			message(k)

			PathwaysResults=list()

			for (i in 1:length(DataPrepared[[k]]$pvalsgenes)){
				message(paste(k,i,sep='.'))
				temp=list()
				temp[[1]]=DataPrepared[[k]]$objects[[i]] # the objects
				temp[[2]]=DataPrepared[[k]]$Genes[[i]]		# the genes
				pvalscluster=DataPrepared[[k]]$pvalsgenes[[i]]

				Entrezs=sapply(names(pvalscluster),function(x) return(geneInfo$ENTREZID[which(geneInfo$SYMBOL==x)]))

				if(path.method=="MLP"){
					## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR

					names(pvalscluster) = Entrezs

					out.mlp <- MLP::MLP(
							geneSet = geneSet,
							geneStatistic = pvalscluster,
							minGenes = 5,
							maxGenes = 100,
							rowPermutations = TRUE,
							nPermutations = 100,
							smoothPValues = TRUE,
							probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9,addGeneSetDescription=TRUE)

					output = list()
					#output$gene.p.values = p.adjust(p.values,method="fdr")

					#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
					#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
					#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)

					#if(is.null(topP)){
					#		topP=length(ranked.genesets.table$p.values<=sign)
					#}

					if(is.null(topP)){
						topP=length(which(out.mlp$geneSetPValue<=sign))
					}

					#TopPaths=ranked.genesets.table[1:topP,]
					#AllPaths=ranked.genesets.table

					TopPaths=out.mlp[1:topP,]
					AllPaths=out.mlp

					output$TopPaths=TopPaths
					attr(output$TopPaths,'geneSetSource')=geneSetSource
					output$AllPaths=AllPaths
					attr(output$AllPaths,'geneSetSource')=geneSetSource
					#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]


					#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) 	)
					#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
					#output$nr.genesets = nr.genesets

					#output$object = out.mlp
					#output$method = "MLP"

					temp[[3]]=output
				}
				names(temp)=c("objects","Genes","Pathways")
				PathwaysResults[[i]]=temp
				names(PathwaysResults)[i]=paste("Cluster",i,sep=" ")

			}

			ResultMLP[[k]]=PathwaysResults
		}
		names(ResultMLP)=names
		for(i in 1:length(ResultMLP)){
			for(k in 1:length(ResultMLP[[i]])){
				if(is.null(ResultMLP[[i]][[k]])[1]){
					ResultMLP[[i]][[k]]=NA
					names(ResultMLP[[i]])[k]=paste("Cluster",k,sep=" ")
				}
			}
			if(length(ResultMLP[[i]]) != maxclus){
				extra=maxclus-length(ResultMLP[[i]])
				for(j in 1:extra){
					ResultMLP[[i]][[length(ResultMLP[[i]])+j]]=NA
					names(ResultMLP[[i]])[length(ResultMLP[[i]])]=paste("Cluster",length(ResultMLP[[i]]),sep=" ")
				}
			}
		}
	}
	return(ResultMLP)
}

#' @title Iterations of the pathway analysis
#'
#' @description The MLP method to perform pathway analysis is based on resampling of the
#' data. Therefore it is recommended to perform the pathway analysis multiple
#' times to observe how much the results are influenced by a different
#' resample. The function \code{PathwaysIter} performs the pathway analysis as
#' described in \code{Pathways} a specified number of times. The input can be
#' one data set or a list as in \code{Pathway.2} and \code{Pathways}.
#' @export PathwaysIter
#'
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix of the objects. The rows should
#' correspond with the genes.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function ("GOBP",
#' "GOMF","GOCC", "KEGG" or "REACTOME"). Default is "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Default is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param niter The number of times to perform pathway analysis. Default is 10.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return A list with an element per iteration, each being the output of \code{Pathways}.
#' @details This function relies on the suggested Bioconductor packages 'MLP', 'biomaRt' and 'org.Hs.eg.db'.
#' @examples
#'
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#' data(GeneInfo)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#'
#' MCF7_Paths_FandT=PathwaysIter(List=L, geneExpr = geneMat, nrclusters = 7, method =
#' c("limma", "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE,
#' weightclust = TRUE, names =names)
#' }
PathwaysIter<-function(List,Selection=NULL,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,niter=10,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}

	PathwaysOutput = list()
	for (i in 1:niter){
		message(paste("Iteration",i,sep=" "))
		mlp = Pathways(List,Selection,geneExpr,nrclusters,method,geneInfo,geneSetSource,topP,topG,GENESET,sign=sign,fusionsLog,weightclust,names)
		PathwaysOutput [[length(PathwaysOutput )+1]] = mlp
		names(PathwaysOutput )[i]=paste("Iteration",i,sep=" ")
	}


	return(PathwaysOutput)
}

#' @title Pathway analysis for a selection of objects
#' @param List A list of clustering outputs or output of the\code{DiffGenes}
#' function. The first element of the list will be used as the reference in
#' \code{ReorderToReference}. The output of \code{ChooseFeatures} is also
#' accepted.
#' @param Selection If pathway analysis should be conducted for a specific
#' selection of objects, this selection can be provided here. Selection can
#' be of the type "character" (names of the objects) or "numeric" (the number
#' of specific cluster). Default is NULL.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' The rows should correspond with the genes.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' The number of clusters should not be specified if the interest lies only in
#' a specific selection of objects which is known by name.  Otherwise, it is
#' required. Default is NULL.
#' @param method The method to applied to look for differentially expressed genes and related pathways. For now, only the
#' limma method is available for gene analysis and the MLP method for pathway analysis. Default is c("limma","MLP").
#' @param geneInfo A data frame with at least the columns ENTREZID and SYMBOL.
#' This is necessary to connect the symbolic names of the genes with their
#' EntrezID in the correct order. The order of the gene is here not in the
#' order of the rownames of the gene expression matrix but in the order of
#' their significance. Default is NULL.
#' @param geneSetSource The source for the getGeneSets function, defaults to
#' "GOBP".
#' @param topP Overrules sign. The number of pathways to display for each
#' cluster. If not specified, only the significant genes are shown. Default is NULL.
#' @param topG Overrules sign. The number of top genes to be returned in the
#' result. If not specified, only the significant genes are shown. Defaults is NULL.
#' @param GENESET Optional. Can provide own candidate gene sets. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @description Internal function of \code{Pathways}.
#' @export
#' @details This function relies on the suggested Bioconductor packages 'MLP', 'biomaRt' and 'org.Hs.eg.db'.
#' @examples
#' \dontrun{
#' PathwaysSelection(List=L,Selection=c("Cpd1","Cpd2"),geneExpr=geneMat,geneInfo=GeneInfo)
#' }
PathwaysSelection<-function(List=NULL,Selection,geneExpr=NULL,nrclusters=NULL,method=c("limma", "MLP"),geneInfo=NULL,geneSetSource = "GOBP",topP=NULL,topG=NULL,GENESET=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if (!requireNamespace("MLP", quietly = TRUE)) {
		stop("MLP needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("biomaRt", quietly = TRUE)) {
		stop("biomaRt needed for this function to work. Please install it.",
				call. = FALSE)
	}

	if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
		stop("org.Hs.eg.db needed for this function to work. Please install it.",
				call. = FALSE)
	}


	method.test = function(sign.method,path.method){
		method.choice = FALSE

		if( sign.method=="limma"  & path.method=="MLP"  ){
			method.choice = TRUE
		}
		if(method.choice==TRUE){
			return(list(sign.method=sign.method,path.method=path.method))
		}
		else{
			stop("Incorrect choice of method.")
		}

	}

	method.out = method.test(method[1],method[2])

	sign.method = method.out$sign.method
	path.method = method.out$path.method

	if(length(geneInfo$ENTREZID)==1){
		geneInfo$ENTREZID = colnames(geneExpr)
	}

	# Determining the genesets if they were not given with the function input
	if((inherits(GENESET,"geneSetMLP"))[1] ){
		geneSet <- GENESET
	}
	else{
		geneSet <- MLP::getGeneSets(species = "Human",geneSetSource = geneSetSource,entrezIdentifiers = geneInfo$ENTREZID)
	}


	if(inherits(Selection,"character")){
		ResultMLP=list()

		DataPrepared<-plyr::try_default(PreparePathway(List[[1]],geneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(DataPrepared)){
			Temp=List[[1]]
			objects=list()
			objects$LeadCpds=Selection
			objects$OrderedCpds=colnames(geneExpr)
			Temp[[length(Temp)+1]]=list(objects=objects)
			names(Temp)[length(Temp)]=paste("Cluster")

			DataPrepared<-PreparePathway(Temp,geneExpr,topG,sign)
		}


		temp=list()
		temp[[1]]=DataPrepared$objects[[1]] #names of the objects
		temp[[2]]=DataPrepared$Genes[[1]]




		if(path.method=="MLP"){
			## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
			p.values=DataPrepared$pvalsgenes[[1]]

			Entrezs=sapply(names(p.values),function(x) return(geneInfo$ENTREZID[which(geneInfo$SYMBOL==x)]))

			names(p.values) = Entrezs
			out.mlp <- MLP::MLP(
					geneSet = geneSet,
					geneStatistic = p.values,
					minGenes = 5,
					maxGenes = 100,
					rowPermutations = TRUE,
					nPermutations = 100,
					smoothPValues = TRUE,
					probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)

			output = list()
			#output$gene.p.values = p.adjust(p.values,method="fdr")

			#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
			#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
			#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)

			#if(is.null(topP)){
			#		topP=length(ranked.genesets.table$p.values<=sign)
			#}

			if(is.null(topP)){
				topP=length(which(out.mlp$geneSetPValue<=sign))
			}

			#TopPaths=ranked.genesets.table[1:topP,]
			#AllPaths=ranked.genesets.table

			TopPaths=out.mlp[1:topP,]
			AllPaths=out.mlp

			output$TopPaths=TopPaths
			attr(output$TopPaths,'geneSetSource')=geneSetSource
			output$AllPaths=AllPaths
			attr(output$AllPaths,'geneSetSource')=geneSetSource
			#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]

			#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
			#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
			#output$nr.genesets = nr.genesets

			#output$object = out.mlp
			#output$method = "MLP"

			temp[[3]]=output


		}
		names(temp)=c("objects","Genes","Pathways")
		ResultMLP[[1]]=temp
		names(ResultMLP)="Selection"
	}



	else if(inherits(Selection,"numeric") & !(is.null(List))){
		check<-plyr::try_default(PreparePathway(List[[1]],geneExpr,topG,sign),NULL,quiet=TRUE)
		if(is.null(check)){

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
			Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
			List=ListNew

			DataPrepared=list()
			for (k in 1:dim(Matrix)[1]){

				cluster=Selection

				check<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(check)){
					Temp=List[[k]]
					objects=list()
					objects$LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)]
					objects$OrderedCpds=stats::as.hclust(List[[k]]$Clust$Clust)$labels[stats::as.hclust(List[[k]]$Clust$Clust)$order]
					Temp[[length(Temp)+1]]=list(objects=objects)
					names(Temp)[length(Temp)]=paste("Cluster")

					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)
				}
			}
			names(DataPrepared)=names
		}

		else{
			for(k in 1:length(List)){
				DataPrepared[[k]]<-plyr::try_default(PreparePathway(List[[k]],geneExpr,topG,sign),NULL,quiet=TRUE)
				if(is.null(DataPrepared[[k]])){
					Temp=List[[k]]


					objects=list()
					objects$LeadCpds=List[[k]]$objects$LeadCpds
					objects$OrderedCpds=List[[k]]$objects$OrderedCpds

					Temp[[length(Temp)+1]]=list(objects=objects)
					names(Temp)[length(Temp)]=paste("Cluster")

					DataPrepared[[k]]<-PreparePathway(Temp,geneExpr,topG,sign)

				}
			}


		}

		ResultMLP=list()
		for (k in 1:length(DataPrepared)){
			message(k)
			cluster=Selection


			temp=list()
			temp[[1]]=DataPrepared[[k]]$objects[[1]] #names of the objects
			temp[[2]]=DataPrepared[[k]]$Genes[[1]]


			if(path.method=="MLP"){
				## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
				p.values=DataPrepared[[k]]$pvalsgenes[[1]]
				Entrezs=sapply(names(p.values),function(x) return(geneInfo$ENTREZID[which(geneInfo$SYMBOL==x)]))

				names(p.values) = Entrezs

				out.mlp <- MLP::MLP(
						geneSet = geneSet,
						geneStatistic = p.values,
						minGenes = 5,
						maxGenes = 100,
						rowPermutations = TRUE,
						nPermutations = 100,
						smoothPValues = TRUE,
						probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)

				output = list()
				#output$gene.p.values = p.adjust(p.values,method="fdr")

				#ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
				#ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
				#ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)

				#if(is.null(topP)){
				#		topP=length(ranked.genesets.table$p.values<=sign)
				#}

				if(is.null(topP)){
					topP=length(which(out.mlp$geneSetPValue<=sign))
				}

				#TopPaths=ranked.genesets.table[1:topP,]
				#AllPaths=ranked.genesets.table

				TopPaths=out.mlp[1:topP,]
				AllPaths=out.mlp

				output$TopPaths=TopPaths
				attr(output$TopPaths,'geneSetSource')=geneSetSource
				output$AllPaths=AllPaths
				attr(output$AllPaths,'geneSetSource')=geneSetSource
				#output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=sign,]

				#nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
				#names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
				#output$nr.genesets = nr.genesets

				#output$object = out.mlp
				#output$method = "MLP"

				temp[[3]]=output
			}
			names(temp)=c("objects","Genes","Pathways")
			ResultMLP[[k]]=temp

		}
		names(ResultMLP)=names

	}
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}

	return(ResultMLP)
}

#' @title Prepare data for pathway analysis
#'
#' @description The \code{PreparePathway} function prepares the input for the
#' pathway analysis by computing differentially expressed genes (via limma) for
#' a selection of objects if these were not provided.
#' @export
#' @param Object A clustering output or DiffGenes output element.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' @param topG Overrules sign. The number of top genes. Default is NULL.
#' @param sign The significance level to be handled. Default is 0.05.
#' @return A list with the p-values of the genes, the objects and the genes.
#' @details This function may rely on the suggested Bioconductor package 'limma'.
#' @examples
#' \dontrun{
#' PreparePathway(Object=MCF7_F,geneExpr=geneMat,topG=NULL,sign=0.05)
#' }
PreparePathway<-function(Object,geneExpr,topG,sign){
	FoundGenes=NULL
	FoundComps=NULL

	FoundGenes=FindElement("Genes",Object)

	if(is.null(FoundGenes)|(is.list(FoundGenes) & length(FoundGenes) == 0)){
		FoundComps=FindElement("objects",Object)
		if(is.null(FoundComps)|(is.list(FoundComps) & length(FoundComps) == 0)){
			stop("Specify either the p-values of the genes or a selection of objects to test for DE genes.")
		}

		pvalsgenes=list()
		FoundGenes=list()
		CompsP=list()
		TopDEP=list()
		for(i in 1:length(FoundComps)){
			LeadCpds=FoundComps[[i]]$LeadCpds
			CompsP[[i]]=FoundComps[[i]]
			names(CompsP)[[i]]=paste("objects_",i,sep="")
			if(is.null(LeadCpds)){
				stop("In the objects element, specify an element LeadCpds")
			}

			group <- factor(ifelse(colnames(geneExpr) %in% LeadCpds, 1, 0))

			if(class(geneExpr)[1]=="ExpressionSet"){
				geneExpr$LeadCmpds<-group
				if (!requireNamespace("a4Base", quietly = TRUE)) {
					stop("a4Base needed for this function to work. Please install it.",
							call. = FALSE)
				}

				DElead <- a4Base::limmaTwoLevels(geneExpr,"LeadCmpds")

				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")

				if(is.null(allDE$ID)){
					allDE$Genes <- rownames(allDE)
				}
				else
				{
					allDE$Genes=allDE$ID
				}

				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE <- allDE[1:topG, ]

				Genes <- list(TopDE,allDE)
				names(Genes)<-c("TopDE","AllDE")
			}
			else{
				label.factor = factor(group)
				design = stats::model.matrix(~label.factor)
				fit = limma::lmFit(geneExpr,design=design)
				fit = limma::eBayes(fit)
				allDE = limma::topTable(fit,coef=2,adjust="fdr",n=dim(geneExpr)[1], sort.by="p")

				if(is.null(allDE$ID)){
					allDE$ID <- rownames(allDE)
				}

				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}


				TopDE<-allDE[1:topG,]

				Genes <- list(TopDE,allDE)
				names(Genes)<-c("TopDE","AllDE")

			}
			FoundGenes[[i]]=Genes
			TopDEP[[i]]=FoundGenes[[i]]
			names(TopDEP)[i]=paste("genes_",i,sep="")
			names(FoundGenes)=paste("Genes_",i,sep="")
			pvalsgenes[[i]]=Genes$AllDE$P.Value
			names(pvalsgenes[[i]])=Genes$AllDE$ID
			names(pvalsgenes)[i]=paste("pvals_",i,sep="")
		}
	}

	else{
		pvalsgenes=list()
		TopDEP=list()
		for(i in 1:length(FoundGenes)){
			#names(FoundGenes)[i]=paste("Genes_",i,sep="")
			TopDEP[[i]]=FoundGenes[[i]]
			names(TopDEP)[i]=paste("genes_",i,sep="")
			pvalsgenes[[i]]=FoundGenes[[i]]$AllDE$P.Value
			names(pvalsgenes[[i]])=FoundGenes[[i]]$AllDE$ID
			names(pvalsgenes)[i]=paste("pvals_",i,sep="")
		}

		FoundComps=FindElement("objects",Object)
		CompsP=list()
		for(i in 1:length(FoundComps)){
			if(is.null(FoundComps[[i]]$LeadCpds)){
				CompsP[[i]]="No LeadCpds specified"
			}
			else{
				CompsP[[i]]=FoundComps[[i]]
				names(CompsP)[[i]]=paste("objects_",i,sep="")
			}
		}
	}


	return(list(pvalsgenes=pvalsgenes,objects=CompsP,Genes=TopDEP))
}

#' @title Shared genes, pathways and features across methods
#'
#' @description The \code{SharedGenesPathsFeat} function takes the intersection
#' of the differentially expressed genes, the pathways and the characteristic
#' features over the methods per cluster.
#' @export
#' @param DataLimma Optional. The output of a \code{DiffGenes} function. Default is NULL.
#' @param DataMLP Optional. The output of \code{Geneset.intersect} function. Default is NULL.
#' @param DataFeat Optional. The output of \code{CharacteristicFeatures} function. Default is NULL.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param Selection Logical. Whether a selection of objects was considered. Default is FALSE.
#' @return A list with a Table summarizing the shared elements and a Which
#' element detailing the shared genes, pathways and features per cluster.
#' @examples
#' \dontrun{
#' SharedGenesPathsFeat(DataLimma=MCF7_FT_DE,DataMLP=MCF7_Paths_intersection,names=c("FP","TP"))
#' }
SharedGenesPathsFeat<-function(DataLimma=NULL,DataMLP=NULL,DataFeat=NULL,names=NULL,Selection=FALSE){  #Input=result of DiffGenes and Geneset.intersect
	#Include sharedLimma and SharedMLP inside the function
	if(Selection==TRUE){
		ResultShared=SharedSelection(DataLimma,DataMLP,DataFeat,names)
	}

	else if(is.null(DataLimma) & is.null(DataMLP) & is.null(DataFeat)){
		stop("At least one Data set should be specified")
	}

	else{

		List=list(DataLimma,DataMLP,DataFeat)
		AvailableData=sapply(seq(length(List)),function(i) if(!(is.null(List[[i]]))) return(i))
		AvailableData=unlist(AvailableData)
		len=c()
		for(i in AvailableData){
			len=c(len,length(List[[i]]))
		}
		if(length(unique(len))!=1){
			stop("Unequal number of methods for limma and MLP")
		}
		else{
			DataSets=lapply(AvailableData,function(i)  return(List[[i]]))
			nmethods=length(DataSets[[1]])
			nclusters=length(DataSets[[1]][[1]])
		}

		if(is.null(names)){
			for(j in 1:length(DataSets[[1]])){
				names[j]=paste("Method",j,sep=" ")
			}
		}

		which=list()
		table=c()


		for (i in 1:nclusters){

			name=paste("Cluster",i,sep=" ")

			comps=c()

			temp1g=c()
			temp1p=c()
			temp1f=list()


			for(j in 1:nmethods){

				if(!(is.na(DataSets[[1]][[j]][[i]])[1])){
					comps=c(comps,length(DataSets[[1]][[j]][[i]]$objects$LeadCpds))
					names(comps)[j]=paste("Ncomps", names[j],sep=" ")
				}
				else{
					comps=c(comps,"-")
				}

				if(!(is.null(DataLimma))){
					if(!(is.na(DataLimma[[j]][[i]])[1])){
						temp1g=c(temp1g,length(DataLimma[[j]][[i]]$Genes$TopDE$ID))
					}
					else{
						temp1g=c(temp1g,"-")
					}
					names(temp1g)[j]=names[j]
				}
				else{
					temp1g=NULL
				}

				if(!(is.null(DataMLP))){
					if(!(is.na(DataMLP[[j]][[i]])[1])){
						temp1p=c(temp1p,length(DataMLP[[j]][[i]][[3]]$geneSetDescription))
					}
					else{
						temp1p=c(temp1g,"-")
					}
					names(temp1p)[j]=names[j]
				}
				else{
					temp1p=NULL
				}

				if(!(is.null(DataFeat))){
					temp=c()
					for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
						if(!(is.na(DataFeat[[j]][[i]])[1])){
							temp=c(temp,length(DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names))
						}
						else{
							temp=c(temp,"-")
						}

						names(temp)[f]=names(DataFeat[[j]][[i]]$Characteristics)[f]
					}
					temp1f[[j]]=temp
					names(temp1f)[j]=names[j]
				}
				else{
					temp1f=NULL
				}
			}

			j=1
			Continue=TRUE
			while (Continue==TRUE){
				cont=c()
				for(d in 1:length(DataSets)){
					cont=c(cont,!(is.na(DataSets[[d]][[j]][[i]])[1]))
				}

				if(any(cont)){

					sharedcomps=DataSets[[1]][[j]][[i]]$objects$LeadCpds
					nsharedcomps=length(sharedcomps)
					names(nsharedcomps)="nsharedcomps"


					if(!(is.null(DataLimma))){
						sharedgenes=DataLimma[[j]][[i]]$Genes$TopDE$ID
						nsharedgenes=length(sharedgenes)
						names(nsharedgenes)="Nshared"
						#pvalsg=DataLimma[[j]][[i]]$Genes$TopDE$adj.P.Val
					}
					else{
						sharedgenes=NULL
						nsharedgenes=NULL
					}


					if(!(is.null(DataMLP))){
						sharedpaths=DataMLP[[j]][[i]][[3]]$geneSetDescription
						nsharedpaths=length(sharedpaths)
						names(nsharedpaths)="Nshared"
						#pvalsp=DataMLP[[j]][[i]][[3]]$mean_geneSetPValue

					}
					else{
						sharedpaths=NULL
						nsharedpaths=NULL
					}
					if(!(is.null(DataFeat))){
						sharedfeat=list()
						nsharedfeat=list()
						#pvalsf=list()
						for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
							sharedfeat[[f]]=DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names
							names(sharedfeat)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")

							nsharedfeat[[f]]=length(sharedfeat[[f]])
							names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")

							#pvalsf[[f]]=DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$adj.P.Val
							#names(pvalsf)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")

						}

					}
					else{
						sharedfeat=NULL
						nsharedfeat=NULL
					}

					Continue=FALSE
				}
				j=j+1
			}

			if(nmethods>=2){
				for (j in 2:nmethods){
					cont=c()
					for(d in 1:length(DataSets)){
						cont=c(cont,!(is.na(DataSets[[d]][[j]][[i]])[1]))
					}
					if(any(cont)){

						sharedcomps=intersect(sharedcomps,DataSets[[1]][[j]][[i]]$objects$LeadCpds)
						nsharedcomps=length(sharedcomps)
						names(nsharedcomps)="Nsharedcomps"


						if(!(is.null(DataLimma))){
							sharedgenes=intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$TopDE$ID)
							nsharedgenes=length(sharedgenes)
							names(nsharedgenes)="Nshared"

						}
						if(!(is.null(DataMLP))){
							sharedpaths=intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$geneSetDescription)
							nsharedpaths=length(sharedpaths)
							names(nsharedpaths)="Nshared"

						}
						if(!(is.null(DataFeat))){

							for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
								sharedfeat[[f]]=intersect(sharedfeat[[f]],DataFeat[[j]][[i]]$Characteristics[[f]]$TopFeat$Names)
								names(sharedfeat)[f]=paste("shared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")

								nsharedfeat[[f]]=length(sharedfeat[[f]])
								names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[j]][[i]]$Characteristics)[f],sep="")

							}

						}

					}
				}
			}

			pvalsgenes=list()
			meanpvalsgenes=c()

			pvalspaths=list()
			meanpvalspaths=c()

			pvalsfeat=list()
			meanpvalsfeat=c()


			if(!(is.null(sharedgenes))&length(sharedgenes)!=0){
				for(c in 1:nmethods){
					pvalsg=c()
					for(g in sharedgenes){
						if(!(is.na(DataLimma[[c]][[i]])[1])){
							pvalsg=c(pvalsg,DataLimma[[c]][[i]]$Genes$TopDE$adj.P.Val[DataLimma[[c]][[i]]$Genes$TopDE$ID==g])
						}
					}

					pvalsgenes[[c]]=pvalsg
					names(pvalsgenes)[c]=paste("P.Val.",names[c],sep="")
				}

				for(g1 in 1:length(sharedgenes)){
					pvalstemp=c()
					for(c in 1:nmethods){
						if(!(is.na(DataLimma[[c]][[i]])[1])){
							pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
						}
					}
					meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))
				}
				pvalsgenes[[nmethods+1]]=meanpvalsgenes
				names(pvalsgenes)[nmethods+1]="Mean pvals genes"
			}
			else{pvalsgenes=NULL}

			if(!(is.null(sharedpaths))&length(sharedpaths)!=0){
				for(c in 1:nmethods){
					pvalsp=c()
					if(!(is.na(DataMLP[[c]][[i]])[1])){
						for(p in sharedpaths){
							pvalsp=c(pvalsp,DataMLP[[c]][[i]][[3]][DataMLP[[c]][[i]][[3]]$geneSetDescriptions==p,5][1])
						}
					}

					pvalspaths[[c]]=pvalsp
					names(pvalspaths)[c]=paste("P.Val.",names[c],sep="")
				}


				for(p1 in 1:length(sharedpaths)){
					pvalstemp1=c()
					for(c in 1:nmethods){
						if(!(is.na(DataMLP[[c]][[i]])[1])){
							pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])

						}

					}

					meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
				}
				pvalspaths[[nmethods+1]]=meanpvalspaths
				names(pvalspaths)[nmethods+1]="Mean pvals paths"
			}
			else{pvalpaths=NULL}


			if(!(is.null(sharedfeat))){
				for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){

					if(length(sharedfeat[[f]])!=0){
						pvalschar=list()
						for(c in 1:nmethods){
							pvalsf=c()
							if(!(is.na(DataFeat[[c]][[i]])[1])){
								for(s in sharedfeat[[f]]){
									pvalsf=c(pvalsf,DataFeat[[c]][[i]]$Characteristics[[f]]$TopFeat$adj.P.Val[DataFeat[[c]][[i]]$Characteristics[[f]]$TopFeat$Names==s])
								}
							}
							pvalschar[[c]]=pvalsf
							names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
						}
					}
					else{
						pvalschar=list()
						pvalschar[1:nmethods]=0
						for(c in 1:nmethods){
							names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
						}
					}
					pvalsfeat[[f]]=pvalschar
					names(pvalsfeat)[f]=names(DataFeat[[j]][[i]]$Characteristics)[f]

				}

				for(f in 1:length(DataFeat[[j]][[i]]$Characteristics)){
					meanpvalsfeat=c()

					if((length(sharedfeat[[f]])!=0)){
						for(f1 in 1:length(sharedfeat[[f]])) {

							pvalstemp=c()
							for(c in 1:nmethods){
								if(!(is.na(DataFeat[[c]][[i]])[1])){
									pvalstemp=c(pvalstemp,pvalsfeat[[f]][[c]][[f1]])
								}
							}
							meanpvalsfeat=c(meanpvalsfeat,mean(pvalstemp))
						}
					}
					else{
						meanpvalsfeat=0
					}

					pvalsfeat[[f]][[nmethods+1]]=meanpvalsfeat
					names(pvalsfeat[[f]])[nmethods+1]="Mean pvals feat"
				}
				lenchar=length(DataFeat[[j]][[i]]$Characteristics)

			}
			else{
				pvalsfeat=NULL
				lenchar=0
			}


			if(!(is.null(temp1f))){
				temp1f=do.call(rbind.data.frame, temp1f)
				colnames(temp1f)=names(DataFeat[[j]][[i]]$Characteristics)
				temp1f=as.matrix(temp1f)
				nsharedfeat=do.call(cbind.data.frame, nsharedfeat)
				nsharedfeat=as.matrix(nsharedfeat)
			}
			part1=cbind(cbind(temp1g,temp1p),temp1f)
			part1=as.matrix(part1)
			if(is.null(nsharedgenes) & is.null(nsharedpaths) &is.null(nsharedfeat)){
				if(!(is.null(temp1g))){
					nsharedgenes=0
				}
				if(!(is.null(temp1p))){
					nsharedpaths=0
				}
				if(!(is.null(temp1f))){
					nsharedfeat=rep(0,length(temp1f))
				}
			}
			part2=cbind(cbind(nsharedgenes,nsharedpaths),nsharedfeat)
			part2=as.matrix(part2)
			colnames(part1)=NULL
			colnames(part2)=NULL
			rownames(part2)="NShared"
			part3=c()
			for(r in 1:length(comps)){
				part3=rbind(part3,rep(comps[r],dim(part1)[2]))
			}
			colnames(part3)=NULL
			rownames(part3)=names(comps)
			part4=rep(nsharedcomps,dim(part1)[2])
			names(part4)=NULL
			temp=rbind(part1,part2,part3,part4)
			rownames(temp)[nrow(temp)]="Nsharedcomps"


			table=cbind(table,temp)
			colnames(table)=seq(1,ncol(table))

			if(!(is.null(pvalsgenes))){
				SharedGenes=cbind(sharedgenes,do.call(cbind.data.frame, pvalsgenes))
			}
			else{
				SharedGenes=NULL
			}
			if(!(is.null(pvalspaths))){
				SharedPaths=cbind(sharedpaths,do.call(cbind.data.frame, pvalspaths))
			}
			else{
				SharedPaths=NULL
			}
			if(!(is.null(pvalsfeat))){
				SharedFeat=list()
				for(f in 1:lenchar){
					if(length(sharedfeat[[f]])==0){
						SharedFeat[[f]]=NULL
					}
					else{
						SharedFeat[[f]]=cbind(sharedfeat[[f]],do.call(cbind.data.frame, pvalsfeat[[f]]))
						names(SharedFeat)[f]=names(pvalsfeat[[1]])[f]
					}

				}
			}
			else{
				SharedFeat=NULL
			}
			which[[i]]=list(SharedComps=sharedcomps,SharedGenes=SharedGenes,SharedPaths=SharedPaths,SharedFeat=SharedFeat)
			names(which)[i]=paste("Cluster",i,sep=" ")

		}
		#Sep for all situations?
		if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar+2))){
				number=seq(1,dim(table)[2],2+lenchar)[i]
				colnames(table)[number]=paste("G.Cluster",i,sep=" ")
				colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")
				for(u in seq(1:lenchar)){
					colnames(table)[number+1+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}

			}
		}

		else if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),is.null(DataFeat))){
			for (i in 1:length(seq(1,dim(table)[2],2))){
				number=seq(1,dim(table)[2],2)[i]
				colnames(table)[number]=paste("G.Cluster",i,sep=" ")
				colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")

			}
		}

		else if(all(!(is.null(DataLimma)),is.null(DataMLP),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
				number=seq(1,dim(table)[2],1+lenchar)[i]
				colnames(table)[number]=paste("G.Cluster",i,sep=" ")
				for(u in seq(1:lenchar)){
					colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}

			}
		}

		else if(all((is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
				number=seq(1,dim(table)[2],lenchar+1)[i]
				colnames(table)[number]=paste("P.Cluster",i,sep=" ")
				for(u in seq(1:lenchar)){
					colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}

			}
		}

		else if(all(!(is.null(DataLimma)),(is.null(DataMLP)),(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],1))){
				colnames(table)[i]=paste("G.Cluster",i,sep=" ")
			}
		}

		else if(all((is.null(DataLimma)),!(is.null(DataMLP)),(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],1))){
				colnames(table)[i]=paste("P.Cluster",i,sep=" ")
			}
		}

		else if(all((is.null(DataLimma)),(is.null(DataMLP)),!(is.null(DataFeat)))){
			for (i in 1:length(seq(1,dim(table)[2],lenchar))){
				number=seq(1,dim(table)[2],2)[i]
				for(u in c(0,1)){
					colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]][[1]]$Characteristics)[u+1],sep=""),paste(".Cluster",i,sep=" "),sep="")
				}

			}
		}

		ResultShared=list(Table=table,Which=which)

	}
	return(ResultShared)

}

#' @title Intersection of genes and pathways over multiple methods for a selection of objects.
#' @param DataLimma Optional. The output of a \code{DiffGenes} function. Default is NULL.
#' @param DataMLP Optional. The output of \code{Geneset.intersect} function. Default is NULL.
#' @param DataFeat Optional. The output of \code{CharacteristicFeatures}
#' function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @description Internal function of \code{SharedGenesPathsFeat}.
#' @export
#' @return A list with a Table and a Which element for the selection.
#' @examples
#' \dontrun{
#' SharedSelection(DataLimma=MCF7_FT_DE,names="Selection")
#' }
SharedSelection<-function(DataLimma=NULL,DataMLP=NULL,DataFeat=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	if(is.null(DataLimma) & is.null(DataMLP) & is.null(DataFeat)){
		stop("At least one Data set should be specified")
	}

	List=list(DataLimma,DataMLP,DataFeat)
	AvailableData=sapply(seq(length(List)),function(i) if(!(is.null(List[[i]]))) return(i))
	AvailableData=unlist(AvailableData)
	len=c()
	for(i in AvailableData){
		len=c(len,length(List[[i]]))
	}
	if(length(unique(len))!=1){
		stop("Unequal number of methods for limma and MLP")
	}
	else{
		DataSets=lapply(AvailableData,function(i)  return(List[[i]]))
		nmethods=length(DataSets[[1]])
		#nclusters=length(DataSets[[1]][[1]])
	}

	if(is.null(names)){
		for(j in 1:nmethods){
			names[j]=paste("Method",j,sep=" ")
		}
	}

	which=list()
	table=c()

	comps=c()

	temp1g=c()
	temp1p=c()
	temp1f=list()

	for (i in 1:nmethods){

		if(!(is.na(DataSets[[1]][[i]])[1])){
			comps=c(comps,length(DataSets[[1]][[i]]$objects$LeadCpds))
			names(comps)[i]=paste("Ncomps", names[i],sep=" ")
		}
		else{
			comps=c(comps,"-")
		}

		if(!(is.null(DataLimma))){
			if(!(is.na(DataLimma[[i]])[1])){
				temp1g=c(temp1g,length(DataLimma[[i]]$Genes$TopDE$ID))
			}
			else{
				temp1g=c(temp1g,"-")
			}
			names(temp1g)[i]=names[i]
		}
		else{
			temp1g=NULL
		}

		if(!(is.null(DataMLP))){
			if(!(is.na(DataMLP[[i]])[1])){
				temp1p=c(temp1p,length(DataMLP[[i]][[3]]$geneSetDescription))
			}
			else{
				temp1p=c(temp1g,"-")
			}
			names(temp1p)[i]=names[i]
		}
		else{
			temp1p=NULL
		}

		if(!(is.null(DataFeat))){
			temp=c()
			for(f in 1:length(DataFeat[[i]]$Characteristics)){
				if(!(is.na(DataFeat[[i]])[1])){
					temp=c(temp,length(DataFeat[[i]]$Characteristics[[f]]$TopFeat$Names))
				}
				else{
					temp=c(temp,"-")
				}

				names(temp)[f]=names(DataFeat[[i]]$Characteristics)[f]
			}
			temp1f[[i]]=temp
			names(temp1f)[i]=names[i]
		}
		else{
			temp1f=NULL
		}

		if (i==1){
			cont=c()
			for(d in 1:length(DataSets)){
				cont=c(cont,!(is.na(DataSets[[d]][[i]])[1]))
			}

			if(any(cont)){

				sharedcomps=DataSets[[1]][[i]]$objects$LeadCpds
				nsharedcomps=length(sharedcomps)
				names(nsharedcomps)="Nsharedcomps"


				if(!(is.null(DataLimma))){
					sharedgenes=DataLimma[[i]]$Genes$TopDE$ID
					nsharedgenes=length(sharedgenes)
					names(nsharedgenes)="Nshared"
				}
				else{
					sharedgenes=NULL
					nsharedgenes=0
				}


				if(!(is.null(DataMLP))){
					sharedpaths=DataMLP[[i]][[3]]$geneSetDescription
					nsharedpaths=length(sharedpaths)
					names(nsharedpaths)="Nshared"

				}
				else{
					sharedpaths=NULL
					nsharedpaths=0
				}
				if(!(is.null(DataFeat))){
					sharedfeat=list()
					nsharedfeat=list()
					for(f in 1:length(DataFeat[[i]]$Characteristics)){
						sharedfeat[[f]]=DataFeat[[i]]$Characteristics[[f]]$TopFeat$Names
						names(sharedfeat)[f]=paste("Shared: ",names(DataFeat[[1]]$Characteristics)[f],sep="")

						nsharedfeat[[f]]=length(sharedfeat[[f]])
						names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[1]]$Characteristics)[f],sep="")

					}

				}
				else{
					sharedfeat=NULL
					nsharedfeat=0
				}
			}
		}
		else{

			sharedcomps=intersect(sharedcomps,DataSets[[1]][[i]]$objects$LeadCpds)
			nsharedcomps=length(sharedcomps)
			names(nsharedcomps)="Nsharedcomps"


			if(!(is.null(DataLimma))){
				sharedgenes=intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$ID)
				nsharedgenes=length(sharedgenes)
				names(nsharedgenes)="Nshared"

			}
			if(!(is.null(DataMLP))){
				sharedpaths=intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription)
				nsharedpaths=length(sharedpaths)
				names(nsharedpaths)="Nshared"

			}
			if(!(is.null(DataFeat))){

				for(f in 1:length(DataFeat[[i]]$Characteristics)){
					sharedfeat[[f]]=intersect(sharedfeat[[f]],DataFeat[[i]]$Characteristics[[f]]$TopFeat$Names)
					names(sharedfeat)[f]=paste("Shared: ",names(DataFeat[[1]]$Characteristics)[f],sep="")

					nsharedfeat[[f]]=length(sharedfeat[[f]])
					names(nsharedfeat)[f]=paste("Nshared: ",names(DataFeat[[i]]$Characteristics)[f],sep="")


				}

			}
		}
	}

	pvalsgenes=list()
	meanpvalsgenes=c()

	pvalspaths=list()
	meanpvalspaths=c()

	pvalsfeat=list()
	meanpvalsfeat=c()


	if(!(is.null(sharedgenes)) & nsharedgenes != 0 ){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalsg=c(pvalsg,DataLimma[[c]]$Genes$TopDE$adj.P.Val[DataLimma[[c]]$Genes$TopDE$ID==g])
				}
			}

			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("P.Val.",names[c],sep="")
		}

		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()
			for(c in 1:nmethods){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{
		pvalsgenes=NULL
		nsharedgenes=NULL
	}

	if(!(is.null(sharedpaths)) & nsharedpaths != 0){
		for(c in 1:nmethods){
			pvalsp=c()
			if(!(is.na(DataMLP[[c]])[1])){
				for(p in sharedpaths){
					pvalsp=c(pvalsp,DataMLP[[c]][[3]][DataMLP[[c]][[3]]$geneSetDescription==p,5][1])
				}
			}

			pvalspaths[[c]]=pvalsp
			names(pvalspaths)[c]=paste("P.Val.",names[c],sep="")
		}


		for(p1 in 1:length(sharedpaths)){
			pvalstemp1=c()
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])

				}

			}

			meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
		}
		pvalspaths[[nmethods+1]]=meanpvalspaths
		names(pvalspaths)[nmethods+1]="Mean pvals paths"
	}
	else{
		pvalpaths=NULL
		nsharedpaths=NULL
	}

	if(!(is.null(sharedfeat)) & !(is.null(Reduce("+",nsharedfeat)))){
		for(f in 1:length(DataFeat[[1]]$Characteristics)){
			pvalschar=list()
			for(c in 1:nmethods){
				pvalsf=c()
				if(!(is.na(DataFeat[[c]])[1])){
					for(s in sharedfeat[[f]]){
						pvalsf=c(pvalsf,DataFeat[[c]]$Characteristics[[f]]$TopFeat$adj.P.Val[DataFeat[[c]]$Characteristics[[f]]$TopFeat$Names==s])
					}
				}

				pvalschar[[c]]=pvalsf
				names(pvalschar)[c]=paste("P.Val.",names[c],sep="")
			}
			pvalsfeat[[f]]=pvalschar
			names(pvalsfeat)[f]=names(DataFeat[[1]]$Characteristics)[f]

		}

		for(f in 1:length(DataFeat[[1]]$Characteristics)){
			meanpvalsfeat=c()
			for(f1 in 1:length(sharedfeat[[f]])){
				pvalstemp=c()
				for(c in 1:nmethods){
					if(!(is.na(DataFeat[[c]])[1]) & pvalsfeat[[f]]){
						pvalstemp=c(pvalstemp,pvalsfeat[[f]][[c]][[f1]])
					}
				}
				meanpvalsfeat=c(meanpvalsfeat,mean(pvalstemp))
			}

			pvalsfeat[[f]][[nmethods+1]]=meanpvalsfeat
			names(pvalsfeat[[f]])[nmethods+1]="Mean pvals feat"
		}
		lenchar=length(DataFeat[[1]]$Characteristics)

	}
	else{
		if(is.null(sharedfeat)){
			pvalsfeat=NULL
			nsharedfeat=NULL
			lenchar=0
		}
		else if(!(is.null(sharedfeat))){
			pvalsfeat=NULL
			lenchar=lenchar=length(DataFeat[[1]]$Characteristics)
		}
	}

	if(!(is.null(temp1f))){
		temp1f=do.call(rbind.data.frame, temp1f)
		colnames(temp1f)=names(DataFeat[[1]]$Characteristics)
		temp1f=as.matrix(temp1f)
		nsharedfeat=do.call(cbind.data.frame, nsharedfeat)
		nsharedfeat=as.matrix(nsharedfeat)
	}
	part1=cbind(cbind(temp1g,temp1p),temp1f)
	part1=as.matrix(part1)
	if(is.null(nsharedgenes) | is.null(nsharedpaths)  | is.null(nsharedfeat)){
		if(!(is.null(temp1g))){
			nsharedgenes=0
		}
		if(!(is.null(temp1p))){
			nsharedpaths=0
		}
		if(!(is.null(temp1f))){
			nsharedfeat=rep(0,length(temp1f))
		}
	}
	part2=cbind(cbind(nsharedgenes,nsharedpaths),nsharedfeat)
	part2=as.matrix(part2)
	rownames(part2)="NShared"
	#print(str(part2))
	colnames(part1)=NULL
	colnames(part2)=NULL

	part3=c()
	for(r in 1:length(comps)){
		part3=rbind(part3,rep(comps[r],dim(part1)[2]))
	}
	colnames(part3)=NULL
	rownames(part3)=names(comps)
	part4=rep(nsharedcomps,dim(part1)[2])
	names(part4)=NULL
	temp=rbind(part1,part2,part3,part4)
	rownames(temp)[nrow(temp)]="Nsharedcomps"


	table=cbind(table,temp)
	colnames(table)=seq(1,ncol(table))

	if(!(is.null(pvalsgenes))){
		SharedGenes=cbind(SharedGenes=sharedgenes,do.call(cbind.data.frame, pvalsgenes))
	}
	else{
		SharedGenes=NULL
	}
	if(!(is.null(pvalspaths))){
		SharedPaths=cbind(SharedPaths=sharedpaths,do.call(cbind.data.frame, pvalspaths))
	}
	else{
		SharedPaths=NULL
	}
	if(!(is.null(pvalsfeat))){
		SharedFeat=list()
		for(f in 1:lenchar){
			SharedFeat[[f]]=cbind(SharedFeat=sharedfeat[[f]],do.call(cbind.data.frame, pvalsfeat[[f]]))
			names(SharedFeat)[f]=names(pvalsfeat)[f]
		}
	}
	else{
		SharedFeat=NULL
	}

	which[[1]]=list(SharedComps=sharedcomps,SharedGenes=SharedGenes,SharedPaths=SharedPaths,SharedFeat=SharedFeat)
	names(which)[1]="Selection"

	#Sep for all situations?
	if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+2))){
			number=seq(1,dim(table)[2],2+lenchar)[i]
			colnames(table)[number]=paste("G.Cluster",sep=" ")
			colnames(table)[number+1]=paste("P.Cluster",sep=" ")
			for(u in seq(1:lenchar)){
				colnames(table)[number+1+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u],sep=""),paste(".Cluster",sep=" "),sep="")
			}

		}
	}

	else if(all(!(is.null(DataLimma)),!(is.null(DataMLP)),is.null(DataFeat))){
		for (i in 1:length(seq(1,dim(table)[2],2))){
			number=seq(1,dim(table)[2],2)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")
			colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")

		}
	}

	else if(all(!(is.null(DataLimma)),is.null(DataMLP),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
			number=seq(1,dim(table)[2],1+lenchar)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")
			for(u in seq(1:lenchar)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}

		}
	}

	else if(all((is.null(DataLimma)),!(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar+1))){
			number=seq(1,dim(table)[2],lenchar+1)[i]
			colnames(table)[number]=paste("P.Cluster",i,sep=" ")
			for(u in seq(1:lenchar)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u],sep=""),paste(".Cluster",i,sep=" "),sep="")
			}

		}
	}

	else if(all(!(is.null(DataLimma)),(is.null(DataMLP)),(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],1))){
			colnames(table)[i]=paste("G.Cluster",sep=" ")
		}
	}

	else if(all((is.null(DataLimma)),!(is.null(DataMLP)),(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],1))){
			colnames(table)[i]=paste("P.Cluster",sep=" ")
		}
	}

	else if(all((is.null(DataLimma)),(is.null(DataMLP)),!(is.null(DataFeat)))){
		for (i in 1:length(seq(1,dim(table)[2],lenchar))){
			number=seq(1,dim(table)[2],2)[i]
			for(u in c(0,1)){
				colnames(table)[number+u]=paste(paste('Feat.',names(DataFeat[[1]]$Characteristics)[u+1],sep=""),paste(".Cluster",sep=" "),sep="")
			}

		}
	}



	ResultShared=list(Table=table,Which=which)
	return(ResultShared)


}

#' @title Intersection of genes over multiple methods for a selection of objects.
#' @param DataLimma Optional. The output of a \code{DiffGenes} function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @description Internal function of \code{SharedGenesPathsFeat}.
#' @export
#' @return A list with a Table and a Which element for the shared genes.
#' @examples
#' \dontrun{
#' SharedSelectionLimma(DataLimma=MCF7_FT_DE)
#' }
SharedSelectionLimma<-function(DataLimma=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect

	which=list()
	table=c()

	if(is.null(names)){
		for(j in 1:length(DataLimma)){
			names[j]=paste("Method",j,sep=" ")
		}
	}

	nmethods=length(DataLimma)

	temp1g=c()
	comps=c()

	pvalsg=c()
	for (i in 1:nmethods){

		temp1g=c(temp1g,length(DataLimma[[i]]$Genes$TopDE$Genes))
		comps=c(comps,length(DataLimma[[i]]$objects$LeadCpds))


		names(temp1g)[i]=names[i]
		names(comps)[i]=paste("Ncomps",names[i],i,sep=" ")




		if (i==1){
			if(!(is.na(DataLimma[[i]])[1])){
				sharedcomps=DataLimma[[i]]$objects$LeadCpds
				sharedgenes=DataLimma[[i]]$Genes$TopDE$Genes


				pvalsg=c(pvalsg,DataLimma[[i]]$Genes$TopDE$adj.P.Val)


				nsharedcomps=length(DataLimma[[i]]$objects$LeadCpds)
				nsharedgenes=length(DataLimma[[i]]$Genes$TopDE$Genes)

				names(nsharedgenes)="nshared"

				names(nsharedcomps)="nsharedcomps"
			}

		}
		else{
			sharedcomps=intersect(sharedcomps,DataLimma[[i]]$objects$LeadCpds)
			sharedgenes=intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$Genes)


			nsharedcomps=length(intersect(sharedcomps,DataLimma[[i]]$objects$LeadCpds))
			nsharedgenes=length(intersect(sharedgenes,DataLimma[[i]]$Genes$TopDE$Genes))

			names(nsharedgenes)="nshared"

			names(nsharedcomps)="nsharedcomps"
		}

	}
	pvalsgenes=list()
	meanpvalsgenes=c()

	if(nsharedgenes != 0){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalsg=c(pvalsg,DataLimma[[c]]$Genes$TopDE$adj.P.Val[DataLimma[[c]]$Genes$TopDE$Genes==g])
				}
			}

			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("Method",c,sep=" ")
		}

		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()
			for(c in 1:nmethods){
				if(!(is.na(DataLimma[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{pvalsgenes=0}


	temp=rbind(temp1g,nsharedgenes,comps,nsharedcomps)

	table=cbind(table,temp)

	which[[1]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes)
	#names(which)[1]=paste("Cluster",i,sep=" ")


	ResultShared=list(Table=table,Which=which)
	return(ResultShared)
}

#' @title Intersection of pathways over multiple methods for a selection of objects.
#' @param DataMLP Optional. The output of \code{Geneset.intersect} function. Default is NULL.
#' @param names Optional. Names of the methods or "Selection" if it only
#' considers a selection of objects. Default is NULL.
#' @description Internal function of \code{SharedGenesPathsFeat}.
#' @export
#' @return A list with a Table and a Which element for the shared pathways.
#' @examples
#' \dontrun{
#' SharedSelectionMLP(DataMLP=MCF7_Paths_intersection)
#' }
SharedSelectionMLP<-function(DataMLP=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect

	which=list()
	table=c()


	nmethods=length(DataMLP)

	if(is.null(names)){
		for(j in 1:length(DataMLP)){
			names[j]=paste("Method",j,sep=" ")
		}
	}

	temp1g=c()
	temp1p=c()
	comps=c()

	pvalsg=c()
	pvalsp=c()
	for (i in 1:nmethods){

		temp1g=c(temp1g,length(DataMLP[[i]]$Genes$TopDE$Genes))
		temp1p=c(temp1p,length(DataMLP[[i]][[3]]$geneSetDescription))
		comps=c(comps,length(DataMLP[[i]]$objects$LeadCpds))


		names(temp1g)[i]=names[i]
		names(temp1p)[i]=names[i]
		names(comps)[i]=paste("Ncomps",names[i],i,sep=" ")


		if (i==1){
			if(!(is.na(DataMLP[[i]])[1]) | !(is.na(DataMLP[[i]])[1])){
				sharedcomps=DataMLP[[i]]$objects$LeadCpds
				sharedgenes=DataMLP[[i]]$Genes$TopDE$Genes
				sharedpaths=DataMLP[[i]][[3]]$geneSetDescription

				pvalsg=c(pvalsg,DataMLP[[i]]$Genes$TopDE$adj.P.Val)
				pvalsp=c(pvalsp,DataMLP[[i]]$mean_geneSetPValue)

				nsharedcomps=length(DataMLP[[i]]$objects$LeadCpds)
				nsharedgenes=length(DataMLP[[i]]$Genes$TopDE$Genes)
				nsharedpaths=length(DataMLP[[i]][[3]]$geneSetDescription)
				names(nsharedgenes)="nshared"
				names(nsharedpaths)="nshared"
				names(nsharedcomps)="nsharedcomps"
			}

		}
		else{
			sharedcomps=intersect(sharedcomps,DataMLP[[i]]$objects$LeadCpds)
			sharedgenes=intersect(sharedgenes,DataMLP[[i]]$Genes$TopDE$Genes)
			sharedpaths=intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription)

			nsharedcomps=length(intersect(sharedcomps,DataMLP[[i]]$objects$LeadCpds))
			nsharedgenes=length(intersect(sharedgenes,DataMLP[[i]]$Genes$TopDE$Genes))
			nsharedpaths=length(intersect(sharedpaths,DataMLP[[i]][[3]]$geneSetDescription))
			names(nsharedgenes)="nshared"
			names(nsharedpaths)="nshared"
			names(nsharedcomps)="nsharedcomps"
		}

	}
	pvalsgenes=list()
	meanpvalsgenes=c()
	meanpvalspaths=c()
	pvalspaths=list()

	if(nsharedgenes != 0){
		for(c in 1:nmethods){
			pvalsg=c()
			for(g in sharedgenes){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalsg=c(pvalsg,DataMLP[[c]]$Genes$TopDE$adj.P.Val[DataMLP[[c]]$Genes$TopDE$Genes==g])
				}
			}

			pvalsgenes[[c]]=pvalsg
			names(pvalsgenes)[c]=paste("Method",c,sep=" ")
		}

		for(g1 in 1:length(sharedgenes)){
			pvalstemp=c()
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
				}
			}
			meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))
		}
		pvalsgenes[[nmethods+1]]=meanpvalsgenes
		names(pvalsgenes)[nmethods+1]="Mean pvals genes"
	}
	else{pvalsgenes=0}

	if(nsharedpaths!=0){
		for(c in 1:nmethods){
			pvalsp=c()
			if(!(is.na(DataMLP[[c]])[1])){
				for(p in sharedpaths){
					pvalsp=c(pvalsp,DataMLP[[c]][[3]][DataMLP[[c]][[3]]$geneSetDescription==p,5])
				}
			}

			pvalspaths[[c]]=pvalsp


			names(pvalspaths)[c]=paste("Method",c,sep=" ")
		}


		for(p1 in 1:length(sharedpaths)){
			pvalstemp1=c()
			for(c in 1:nmethods){
				if(!(is.na(DataMLP[[c]])[1])){
					pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])

				}

			}

			meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))

		}
		pvalspaths[[nmethods+1]]=meanpvalspaths
		names(pvalspaths)[nmethods+1]="Mean pvals paths"
	}
	else{pvalpaths=0}

	temp=rbind(cbind(temp1g,temp1p),cbind(nsharedgenes,nsharedpaths),cbind(comps,comps),cbind(nsharedcomps,nsharedcomps))

	table=cbind(table,temp)

	which[[1]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes,sharedpaths=sharedpaths,pvalspaths=pvalspaths)
	#names(which)[1]=paste("Cluster",i,sep=" ")


	for (i in 1:length(seq(1,dim(table)[2],2))){
		number=seq(1,dim(table)[2],2)[i]
		colnames(table)[number]=c("G.Cluster")
		colnames(table)[number+1]=c("P.Cluster")

	}
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)


}
