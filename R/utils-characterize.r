#' @title Determine the characteristic features of clusters
#'
#' @description \code{CharacteristicFeatures} identifies the features that are characteristic
#' for the objects of each cluster of one or more clustering results. Binary
#' features are evaluated with Fisher's exact test, continuous features with the
#' t-test. Multiplicity correction is included.
#' @export
#' @param List A list of clustering outputs to be compared. The first element of
#' the list is used as the reference in \code{ReorderToReference}.
#' @param Selection If a selection of objects is provided, the function only
#' investigates the features of this selection. Default is NULL.
#' @param binData A list of the binary feature data matrices. Default is NULL.
#' @param contData A list of continuous feature data matrices. Default is NULL.
#' @param datanames A vector with the names of the data matrices. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param sign The significance level. Default is 0.05.
#' @param topChar The number of top characteristics to return. If NULL, only the
#' significant characteristics are saved. Default is NULL.
#' @param fusionsLog Logical. Indicator for the fusion of clusters. Default is TRUE.
#' @param weightclust Logical. For the outputs of CEC, WeightedClust or
#' WeightedSimClust, if TRUE only the result of the Clust element is considered.
#' Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @return A list with for each method the found (top) characteristics of the
#' feature data per cluster.
#' @references Perualila-Tan N. et al. (2016). Weighted similarity-based clustering
#' of chemical structures and bioactivity data in early drug discovery.
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
#' MCF7_Char=CharacteristicFeatures(List=L,Selection=NULL,binData=list(fingerprintMat),
#' contData=NULL,datanames=c("FP"),nrclusters=7,topChar=20)
#' }
CharacteristicFeatures<-function(List,Selection=NULL,binData=NULL,contData=NULL,datanames=NULL,nrclusters=NULL,sign=0.05,topChar=NULL,fusionsLog=TRUE,weightclust=TRUE,names=NULL){
	if(is.null(datanames)){
		for(j in 1:(length(binData)+length(contData))){
			datanames[j]=paste("Data",j,sep=" ")
		}
	}

	if(!(is.null(Selection))){
		ResultFeat=FeatSelection(List,Selection,binData,contData,datanames,nrclusters,topChar,sign,fusionsLog,weightclust)
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
		names(ListNew)=names
		MatrixClusters=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)

		List=ListNew

		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else if(!is.null(contData)){
			cpdSet <- rownames(contData[[1]])
		}
		else{
			stop("Specify a data set in binData and/or in contData")
		}

		ResultFeat=list()
		maxclus=0
		for (k in 1:dim(MatrixClusters)[1]){

			clusters=MatrixClusters[k,]
			if(max(clusters)>maxclus){
				maxclus=max(clusters)
			}
			Characteristics=list()
			clust=sort(unique(clusters)) #does not matter: Genes[i] puts right elements on right places
			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]
			for (i in clust){

				temp=list()
				LeadCpds=names(clusters)[which(clusters==i)]
				temp[[1]]=list(LeadCpds,OrderedCpds)
				names(temp[[1]])=c("LeadCpds","OrderedCpds") #names of the objects

				group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest

				#Determine characteristic features for the objects: fishers exact test
				result=list()
				if(!is.null(binData)){
					for(i in 1:length(binData)){
						binData[[i]]=binData[[i]]+0
						binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
					}


					for(j in 1: length(binData)){
						binMat=binData[[j]]

						pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)
						pFish <- sort(pFish)
						adjpFish<-stats::p.adjust(pFish, method = "fdr")

						AllFeat=data.frame(Names=names(pFish),P.Value=pFish,adj.P.Val=adjpFish)
						AllFeat$Names=as.character(AllFeat$Names)

						if(is.null(topChar)){
							topChar=length(which(pFish<sign))
						}

						TopFeat=AllFeat[0:topChar,]
						TopFeat$Names=as.character(TopFeat$Names)
						temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
						result[[j]]<-temp1
						names(result)[j]=datanames[length(binData)+j]

					}
				}

				resultC=list()
				if(!is.null(contData)){
					for(j in 1:length(contData)){
						contMat=contData[[j]]

						group1=which(group==1)
						group2=which(group==0)


						pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)

						pTTest <- sort(pTTest)
						adjpTTest<-stats::p.adjust(pTTest, method = "fdr")

						AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
						AllFeat$Names=as.character(AllFeat$Names)
						if(is.null(topChar)){
							topChar=length(which(pTTest<sign))
						}

						TopFeat=data.frame(Names=as.character(names(pTTest[0:topChar])),P.Value=pTTest[0:topChar],adj.P.Val=adjpTTest[0:topChar])
						TopFeat$Names=as.character(TopFeat$Names)
						temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
						resultC[[j]]<-temp1
						names(resultC)[j]=datanames[j]

					}
				}

				temp[[2]]=c(result,resultC)

				names(temp)=c("objects","Characteristics")

				Characteristics[[i]]=temp

				names(Characteristics)[i]=paste("Cluster",i,sep=" ")
			}
			ResultFeat[[k]]=Characteristics

		}
		names(ResultFeat)=names
		for(i in 1:length(ResultFeat)){
			for(k in 1:length(ResultFeat[[i]])){
				if(is.null(ResultFeat[[i]][[k]])[1]){
					ResultFeat[[i]][[k]]=NA
					names(ResultFeat[[i]])[k]=paste("Cluster",k,sep=" ")
				}
			}
			if(length(ResultFeat[[i]]) != maxclus){
				extra=maxclus-length(ResultFeat[[i]])
				for(j in 1:extra){
					ResultFeat[[i]][[length(ResultFeat[[i]])+1]]=NA
					names(ResultFeat[[i]])[length(ResultFeat[[i]])]=paste("Cluster",length(ResultFeat[[i]]),sep=" ")
				}
			}
		}

	}

	return(ResultFeat)
}


#' @title Interactive plot to determine DE Genes and DE features for a specific cluster
#'
#' @description If desired, the function produces a dendrogram of a clustering result. One
#' or multiple clusters can be indicated by a mouse click. From these clusters
#' DE genes and characteristic features are determined. It is also possible to
#' provide the objects of interest without producing the plot.
#' @export
#' @param Interactive Logical. Whether an interactive plot should be made. Defaults to TRUE.
#' @param leadCpds A list of the objects of the clusters of interest. Defaults to NULL.
#' @param clusterResult The output of one of the aggregated cluster functions. Default is NULL.
#' @param colorLab The clustering result the dendrogram should be colored after. Default is NULL.
#' @param binData A list of the binary feature data matrices. Default is NULL.
#' @param contData A list of continuous data sets of the objects. Default is NULL.
#' @param datanames A vector with the names of the data matrices. Default is c("FP").
#' @param geneExpr A gene expression matrix, may also be an ExpressionSet. Default is NULL.
#' @param topChar The number of top characteristics to return. Default is 20.
#' @param topG The number of top genes to return. Default is 20.
#' @param sign The significance level. Default is 0.05.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols The colors to use in the dendrogram. Default is NULL.
#' @param n The number of clusters one wants to identify by a mouse click. Default is 1.
#' @return A list with one element per cluster of interest with elements
#' objects, Characteristics and (optionally) Genes.
#' @details This function may require the suggested packages 'a4Base', 'a4Core'
#' and 'limma' when a gene expression matrix is supplied.
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
#' MCF7_Interactive=ChooseCluster(Interactive=TRUE,leadCpds=NULL,clusterResult=MCF7_T,
#' colorLab=MCF7_F,binData=list(fingerprintMat),datanames=c("FP"),geneExpr=geneMat,
#' topChar = 20, topG = 20,nrclusters=7,n=1)
#' }
ChooseCluster=function(Interactive=TRUE,leadCpds=NULL,clusterResult=NULL,colorLab=NULL,binData=NULL,contData=NULL,datanames=c("FP"),geneExpr=NULL,topChar = 20, topG = 20,sign=0.05,nrclusters=NULL,cols=NULL,n=1){
	if(is.null(datanames)){
		for(j in 1:(length(binData)+length(contData))){
			datanames[j]=paste("Data",j,sep=" ")
		}
	}
	OrInteractive=Interactive

	if(Interactive==TRUE){
		#windows()
		ClusterPlot(clusterResult,colorLab,nrclusters,cols)
		hc1<-stats::as.hclust(clusterResult$Clust)
		ClusterSpecs<-list()
		ClusterSpecs=graphics::identify(hc1, N=n, MAXCLUSTER = nrow(binData[[1]]), function(j) ChooseCluster(Interactive=FALSE,leadCpds=rownames(binData[[1]][j,]),clusterResult,colorLab=NULL,binData=binData,contData=contData,datanames=datanames,geneExpr=geneExpr,topChar=topChar,topG=topG,sign=sign,nrclusters=nrclusters,cols=cols))

		names(ClusterSpecs)<-sapply(seq(1,n),FUN=function(x) paste("Choice",x,sep=" "))

	}
	else{

		DistM<-clusterResult$DistM
		Clust<-clusterResult$Clust
		if(is.null(Clust)){
			clusterResult$Clust=clusterResult
			Clust<-clusterResult$Clust
		}

		hc <- stats::as.hclust(Clust)
		OrderedCpds <- hc$labels[hc$order]

		if(inherits(leadCpds,"character")){
			leadCpds=list(leadCpds)
		}

		Specs=list()

		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else if(!is.null(contData)){
			cpdSet <- rownames(contData[[1]])
		}
		else if(!is.null(geneExpr)){
			cpdSet <- colnames(geneExpr)
		}
		else{
			stop("Specify a data set in binData, contData and/or geneExpr")
		}

		for(a in 1:length(leadCpds)){
			objects=list(leadCpds[[a]],OrderedCpds)
			names(objects)=c("LeadCpds","OrderedCpds")
			group <- factor(ifelse(cpdSet %in% leadCpds[[a]], 1, 0)) #identify the group of interest

			#Determine characteristic features for the objects: fishers exact test
			Characteristics=list()

			resultB=list()
			if(!is.null(binData)){

				if(!inherits(binData,"list")){
					stop("The binary data matrices must be put into a list")
				}
				for(i in 1:length(binData)){
					binData[[i]]=binData[[i]]+0
					binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
				}
				for(j in 1: length(binData)){

					binMat=binData[[j]]

					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)

					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")

					AllFeat=data.frame(Names=names(pFish),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)

					if(is.null(topChar)){
						topChar=length(which(pFish<sign))
					}

					TopFeat=AllFeat[0:topChar,]
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultB[[j]]<-temp1
					names(resultB)[j]=datanames[length(binData)+j]

				}
			}

			resultC=list()
			if(!is.null(contData)){

				if(!inherits(contData,"list")){
					stop("The continuous data matrices must be put into a list")
				}
				for(j in 1:length(contData)){
					contMat=contData[[j]]

					group1=which(group==1)
					group2=which(group==0)

					pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)

					pTTest <- sort(pTTest)
					adjpTTest<-stats::p.adjust(pTTest, method = "fdr")

					AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topChar)){
						topChar=length(which(pTTest<sign))
					}

					TopFeat=data.frame(Names=as.character(names(pTTest[0:topChar])),P.Value=pTTest[0:topChar],adj.P.Val=adjpTTest[0:topChar])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultC[[j]]<-temp1
					names(resultC)[j]=datanames[j]

				}

			}

			Characteristics=c(resultB,resultC)
			names(Characteristics)=datanames

			#Determine DE Genes with limma --> make difference between "regular" data matrix and "expression set"
			#GeneExpr.2=GeneExpr[,colnames(Matrix)]
			if(!is.null(geneExpr)){
				cpdSetG <-colnames(geneExpr)
				groupG <- factor(ifelse(cpdSetG %in% leadCpds[[a]], 1, 0)) #identify the group of interest
				if(class(geneExpr)[1]=="ExpressionSet"){
					if (!requireNamespace("limma", quietly = TRUE)) stop("ChooseCluster() requires the suggested package 'limma'.")
					if (!requireNamespace("a4Core", quietly = TRUE)) stop("ChooseCluster() requires the suggested package 'a4Core'.")
					geneExpr$LeadCmpds<-groupG


					if (!requireNamespace("a4Base", quietly = TRUE)) {
						stop("a4Base needed for this function to work. Please install it.",
								call. = FALSE)
					}

					DElead <- a4Base::limmaTwoLevels(geneExpr,"LeadCpds")

					#allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), resort.by = "logFC",sort.by="p")
					allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					if(is.null(topG)){
						topG=length(which(allDE$adj.P.Val<=sign))
					}
					TopDE <- allDE[0:topG, ]
					#TopAdjPval<-TopDE$adj.P.Val
					#TopRawPval<-TopDE$P.Value

					#RawpVal<-allDE$P.Value
					#AdjpVal <- allDE$adj.P.Val
					#genesEntrezId <- allDE$ENTREZID

					Genes<-list(TopDE,allDE)
					names(Genes)<-c("TopDE","AllDE")
					#Genes <- list(TopDE$SYMBOL,TopAdjPval,TopRawPval,genesEntrezId,RawpVal,AdjpVal)
					#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
				}
				else{
					if (!requireNamespace("limma", quietly = TRUE)) stop("ChooseCluster() requires the suggested package 'limma'.")

					label.factor = factor(groupG)
					design = stats::model.matrix(~label.factor)
					fit = limma::lmFit(geneExpr,design=design)
					fit = limma::eBayes(fit)

					#allDE = topTable(fit,coef=2,adjust="fdr",n=nrow(GeneExpr),resort.by = "logFC", sort.by="p")
					allDE = limma::topTable(fit,coef=2,adjust="fdr",n=nrow(geneExpr), sort.by="p")
					if(is.null(allDE$ID)){
						allDE$ID <- rownames(allDE)
					}
					else
					{
						allDE$ID=allDE$ID
					}
					if(is.null(topG)){
						topG=length(which(allDE$adj.P.Val<=sign))
					}
					TopDE=allDE[0:topG,]
					#TopAdjPval<-TopDE$adj.P.Val
					#TopRawPval<-TopDE$P.Value

					#RawpVal<-allDE$P.Value
					#AdjpVal <- allDE$adj.P.Val

					Genes<-list(TopDE,allDE)
					names(Genes)<-c("TopDE","AllDE")
					#Genes <- list(TopDE[,1],TopAdjPval,TopRawPval,allDE[,1],RawpVal,AdjpVal)
					#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")

				}

				out=list(objects,Characteristics,Genes)
				names(out)=c("objects","Characteristics","Genes")
				Specs[[a]]=out
				names(Specs)[a]=paste("Choice",a,sep=" ")
			}
			else{
				out=list(objects,Characteristics)
				names(out)=c("objects","Characteristics")
				Specs[[a]]=out
				names(Specs)[a]=paste("Choice",a,sep=" ")
			}
		}

		if(OrInteractive==TRUE|length(Specs)==1){
			return(out)
		}
		else{
			return(Specs)
		}
	}
	class(ClusterSpecs)="ChosenClusters"
	return(ClusterSpecs)
}


#' @title Determine the characteristic features of a cluster or selection
#'
#' @description \code{FeatSelection} determines the characteristic features of a selection of
#' objects or of a specific cluster across the methods. Binary features are
#' evaluated with Fisher's exact test, continuous features with the t-test.
#' @export
#' @param List A list of clustering outputs. Default is NULL.
#' @param Selection Either a character vector of objects, or a numeric cluster
#' number. Default is NULL.
#' @param binData A list of the binary feature data matrices. Default is NULL.
#' @param contData A list of continuous feature data matrices. Default is NULL.
#' @param datanames A vector with the names of the data matrices. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param topChar The number of top characteristics to return. Default is NULL.
#' @param sign The significance level. Default is 0.05.
#' @param fusionsLog Logical. Indicator for the fusion of clusters. Default is TRUE.
#' @param weightclust Logical. For the outputs of CEC, WeightedClust or
#' WeightedSimClust, if TRUE only the result of the Clust element is considered.
#' Default is TRUE.
#' @return A list with the found (top) characteristics of the feature data.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F)
#' MCF7_FeatSel=FeatSelection(List=L,Selection=4,binData=list(fingerprintMat),
#' datanames=c("FP"),nrclusters=7,topChar=20)
#' }
FeatSelection<-function(List,Selection=NULL,binData=NULL,contData=NULL,datanames=NULL,nrclusters=NULL,topChar=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE){


	if(is.null(datanames)){
		for(j in 1:(length(binData)+length(contData))){
			datanames[j]=paste("Data",j,sep=" ")
		}
	}

	if(inherits(Selection,"character")){


		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else if(!is.null(contData)){
			cpdSet <- rownames(contData[[1]])
		}
		else{
			stop("Specify a data set in binData and/or in contData")
		}
		ResultFeat=list()
		Characteristics=list()
		temp=list()

		LeadCpds=Selection #names of the objects
		OrderedCpds=cpdSet
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")

		group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest

		#Determine characteristic features for the objects: fishers exact test
		result=list()

		if(!is.null(binData)){
			for(i in 1:length(binData)){
				binData[[i]]=binData[[i]]+0
				binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
			}
			for(j in 1: length(binData)){
				binMat=binData[[j]]


				if(length(LeadCpds)==1){
					FP=which(binMat[LeadCpds,]==1)
					Ranks=sort(colSums(binMat[,FP]),)

					N=c(names(Ranks),colnames(binMat)[which(!colnames(binMat)%in%names(Ranks))])

					AllFeat=data.frame(Names=as.character(N))
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topChar)){
						topChar=length(Ranks)
					}
					else{
						Ranks=Ranks[c(1:topChar)]
					}
					TopFeat=data.frame(Names=as.character(names(Ranks)))
					TopFeat$Names=as.character(TopFeat$Names)


					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(result)[j]=datanames[j]


				}
				else{
					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)

					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")

					AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topChar)){
						topChar=length(which(pFish<sign))
					}

					TopFeat=data.frame(Names=as.character(names(pFish[0:topChar])),P.Value=pFish[0:topChar],adj.P.Val=adjpFish[0:topChar])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(result)[j]=datanames[j]
				}
			}

		}

		resultC=list()
		if(!is.null(contData)){
			for(j in 1:length(contData)){
				contMat=contData[[j]]

				group1=which(group==1)
				group2=which(group==0)


				pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)

				pTTest <- sort(pTTest)
				adjpTTest<-stats::p.adjust(pTTest, method = "fdr")

				AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
				AllFeat$Names=as.character(AllFeat$Names)
				if(is.null(topC)){
					topC=length(which(pTTest<sign))
				}

				TopFeat=data.frame(Names=as.character(names(pTTest[0:topC])),P.Value=pTTest[0:topC],adj.P.Val=adjpTTest[0:topC])
				TopFeat$Names=as.character(TopFeat$Names)
				temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
				resultC[[j]]<-temp1
				names(resultC)[j]=datanames[length(binData)+j]

			}
		}

		temp[[2]]=c(result,resultC)


		names(temp)=c("objects","Characteristics")

		ResultFeat[[1]]=temp
		names(ResultFeat)="Selection"



	}

	else if(inherits(Selection,"numeric") & !(is.null(List))){

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


		if(!is.null(binData)){
			cpdSet <- rownames(binData[[1]])
		}
		else{
			cpdSet <- rownames(contData[[1]])
		}

		ResultFeat=list()
		for(k in 1:dim(Matrix)[1]){
			cluster=Selection

			hc<-stats::as.hclust(List[[k]]$Clust$Clust)
			OrderedCpds <- hc$labels[hc$order]

			Genes=list()
			temp=list()

			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the objects
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")

			group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest

			#Determine characteristic features for the objects: fishers exact test
			result=list()
			if(!is.null(binData)){
				for(i in 1:length(binData)){
					binData[[i]]=binData[[i]]+0
					binData[[i]]<-binData[[i]][,which(colSums(binData[[i]]) != 0 & colSums(binData[[i]]) != nrow(binData[[i]]))]
				}

				for(j in 1: length(binData)){
					binMat=binData[[j]]
					pFish <- apply(binMat, 2, function(x) stats::fisher.test(table(x, group))$p.value)

					pFish <- sort(pFish)
					adjpFish<-stats::p.adjust(pFish, method = "fdr")

					AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topC)){
						topC=length(which(pFish<0.05))
					}

					TopFeat=AllFeat[0:topC,]
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					result[[j]]<-temp1
					names(resultC)[j]=datanames[length(binData)+j]

				}
			}

			resultC=list()
			if(!is.null(contData)){
				for(j in 1:length(contData)){
					contMat=contData[[j]]

					group1=which(group==1)
					group2=which(group==0)


					pTTest <- apply(contMat, 2, function(x) stats::t.test(x[group1],x[group2])$p.value)

					pTTest <- sort(pTTest)
					adjpTTest<-stats::p.adjust(pTTest, method = "fdr")

					AllFeat=data.frame(Names=as.character(names(pTTest)),P.Value=pTTest,adj.P.Val=adjpTTest)
					AllFeat$Names=as.character(AllFeat$Names)
					if(is.null(topC)){
						topC=length(which(pTTest<sign))
					}

					TopFeat=data.frame(Names=as.character(names(pTTest[0:topC])),P.Value=pTTest[0:topC],adj.P.Val=adjpTTest[0:topC])
					TopFeat$Names=as.character(TopFeat$Names)
					temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
					resultC[[j]]<-temp1
					names(resultC)[j]=datanames[j]

				}
			}

			temp[[2]]=c(result,resultC)

			names(temp)=c("objects","Characteristics")
			ResultFeat[[k]]=temp
			names(ResultFeat)[k]=paste(names[k],": Cluster",cluster, sep=" ")
		}
	}

	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultFeat)

}


#' @title List all features present in a selected cluster of objects
#'
#' @description The function \code{FeaturesOfCluster} lists the number of features objects
#' of the cluster have in common. A threshold can be set selecting among how
#' many objects of the cluster the features should be shared. An optional
#' plot of the features is available.
#' @export
#' @param leadCpds A character vector containing the objects one wants to
#' investigate in terms of features.
#' @param data The data matrix.
#' @param threshold The number of objects the features at least should be
#' shared amongst. Default is set to 1.
#' @param plot Logical. Indicates whether or not a BinFeaturesPlot should be
#' set up. Default is TRUE.
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "new".
#' @param location If plottype is "pdf", a location should be provided. Default is NULL.
#' @return A list with 2 items: the number of shared features among the objects,
#' and a character vector of the plotted features.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#'
#' Lead=rownames(fingerprintMat)[1:5]
#'
#' FeaturesOfCluster(leadCpds=Lead,data=fingerprintMat,
#' threshold=1,plot=TRUE,plottype="new",location=NULL)
#' }
FeaturesOfCluster<-function(leadCpds,data,threshold=1,plot=TRUE,plottype="new",location=NULL){
	SubsetData=as.matrix(data[which(rownames(data)%in%leadCpds),])

	Common=SubsetData%*%t(SubsetData)

	if(threshold>length(leadCpds)){
		stop("threshold is larger than the number of LeadCpds. This number can maximally be the number of LeadCpds.")
	}


	Features=colnames(SubsetData[,which(apply(SubsetData,2,sum)>=threshold)])

	if(is.null(Features)){
		Features=""
		Plot=FALSE
	}

	#SharedAmongLeadCpds=SubsetData[,Features]

	if(plot==TRUE){
		BinFeaturesPlot_SingleData(leadCpds=leadCpds,orderLab=rownames(data),
				features=as.character(Features),data=data,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Shared Features Among Selected objects"),
				margins=c(8.5,2.0,0.5,9.5),plottype=plottype,location=location)
	}

	Out=list("Number_of_Common_Features"=Common,"SharedFeatures"=Features)

	return(Out)
}


#' @title Track a cluster or a selection of objects across multiple methods
#'
#' @description \code{TrackCluster} follows a selection of objects or a specific cluster
#' across a sequence of clustering results, recording how the objects are
#' distributed over the clusters of each method and optionally producing
#' tracking plots.
#' @export
#' @param List A list of clustering outputs to be compared. The first element of
#' the list is used as the reference in \code{ReorderToReference}.
#' @param Selection Either a character vector of objects, or a numeric cluster number.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param followMaxComps Logical. Whether to follow the cluster with the maximum
#' number of objects. Default is FALSE.
#' @param followClust Logical. Whether to follow the cluster of interest. Default is TRUE.
#' @param fusionsLog Logical. Indicator for the fusion of clusters. Default is TRUE.
#' @param weightclust Logical. For the outputs of CEC, WeightedClust or
#' WeightedSimClust, if TRUE only the result of the Clust element is considered.
#' Default is TRUE.
#' @param names Optional. Names of the methods. Default is NULL.
#' @param selectionPlot Logical. Whether to produce the selection plot. Default is FALSE.
#' @param table Logical. Whether to produce a table of the tracking. Default is FALSE.
#' @param completeSelectionPlot Logical. Whether to produce the complete selection
#' plot. Default is FALSE.
#' @param ClusterPlot Logical. Whether to produce the cluster plot. Default is FALSE.
#' @param cols The colors to use in the plots. Default is NULL.
#' @param legendposx The x-coordinate of the legend. Default is 0.5.
#' @param legendposy The y-coordinate of the legend. Default is 2.4.
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "sweave".
#' @param location If plottype is "pdf", a location should be provided. Default is NULL.
#' @return A list with one element per method describing how the selection is
#' distributed across the clusters, and optionally a tracking table.
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
#' MCF7_Track=TrackCluster(List=L,Selection=4,nrclusters=7,followMaxComps=FALSE,
#' followClust=TRUE,selectionPlot=TRUE,table=TRUE,cols=NULL)
#' }
TrackCluster <- function(List,Selection,nrclusters=NULL,followMaxComps=FALSE,followClust=TRUE,fusionsLog=TRUE,weightclust=TRUE,names=NULL,selectionPlot=FALSE,table=FALSE,completeSelectionPlot=FALSE,ClusterPlot=FALSE,cols=NULL,legendposx=0.5,legendposy=2.4,plottype="sweave",location=NULL){

	ClusterDistribution.2<-function(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,weightclust,names){
		FoundClusters=list()
		FoundCl=list()
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)

		ClusterNumber=NULL
		if(inherits(Selection,"numeric")){
			ClusterNumber=Selection
			Selection=colnames(Matrix)[which(Matrix[1,]==Selection)]
		}


		if(followClust==TRUE){
			cluster.interest=NULL
			m=1
			while(m<=dim(Matrix)[1] & is.null(cluster.interest)){
				clustersfound=unique(Matrix[m,which(colnames(Matrix)%in%Selection)])
				if(length(clustersfound)==1){
					cluster.interest=clustersfound
				}
				else{
					m=m+1
				}
			}
			if(is.null(cluster.interest)){
				message("This selection is not found to be part of a cluster. FollowMaxComps will be put to true and the function will proceed")
				followMaxComps=TRUE
				followClust=FALSE
			}

		}

		for(i in 1:dim(Matrix)[1]){


			temp=list()
			temp[[1]]=Selection
			names(temp)[1]="Selection"

			clusternumbers=unique(Matrix[i,which(colnames(Matrix)%in%Selection)])

			nr.clusters=length(clusternumbers)
			temp[[2]]=nr.clusters
			names(temp)[2]="nr.clusters"

			min.together=min(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			max.together=max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))

			nr.min.max.together=c(min.together,max.together)
			temp[[3]]=nr.min.max.together
			names(temp)[3]="nr.min.max.together"

			min.perc.together <- min.together/length(Selection) *100
			max.perc.together <- max.together/length(Selection) *100

			perc.min.max.together =c(min.perc.together,max.perc.together)
			temp[[4]]=perc.min.max.together
			names(temp)[4]="perc.min.max.together"

			temp[[5]]=list()
			names(temp)[5]="AllClusters"

			for(a in 1:length(clusternumbers)){
				temp[[5]][[a]]=list()
				names(temp[[5]])[a]=paste("Cluster",clusternumbers[a],sep=" ")

				temp[[5]][[a]][[1]]=clusternumbers[a]
				temp[[5]][[a]][[2]]=names(which(Matrix[i,]==clusternumbers[a])) #complete cluster
				temp[[5]][[a]][[3]]=intersect(Selection,temp[[5]][[a]][[2]]) #Objects from original selection in this cluster
				temp[[5]][[a]][[4]]=temp[[5]][[a]][[2]][which(!(temp[[5]][[a]][[2]] %in% Selection))] #Objects extra to this cluster
				names(temp[[5]][[a]])=c("clusternumber","Complete cluster","Objects from original selection in this cluster","Objects extra to this cluster")
			}

			if(followMaxComps==TRUE){

				maxcluster=names(which(table(Matrix[i,which(colnames(Matrix)%in%Selection)])==max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))))
				temp[[6]]=maxcluster
				names(temp)[6]="Cluster with max Objects"
				complabels=rownames(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(maxcluster)),drop=FALSE])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(maxcluster))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"

			}

			if(followClust==TRUE){

				temp[[6]]=cluster.interest
				names(temp)[6]="Cluster"
				complabels=rownames(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(cluster.interest)),drop=FALSE])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(cluster.interest))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"

			}

			FoundClusters[[i]]=temp
		}

		if(ClusterPlot==TRUE & !(is.null(ClusterNumber))){

			for(i in 1:dim(Matrix)[1]){
				temp=list()
				if(i==1){
					temp[[1]]=names(Matrix[i,which(Matrix[i,]==ClusterNumber)])
					names(temp)[1]=paste("Cluster ", ClusterNumber,sep="")
					PrevCluster=temp[[1]]
				}
				else{
					temp[[1]]=names(Matrix[i,which(Matrix[i,]==ClusterNumber)])
					names(temp)[1]=paste("Cluster ",  ClusterNumber,sep="")

					Diss=list()
					DissComps=NULL
					if(length((which(!(PrevCluster%in%temp[[1]]))!=0))){
						DissComps=PrevCluster[(which(!(PrevCluster%in%temp[[1]])))]
					}
					if(!(is.null(DissComps))){ #Dissapeared comps: what is missing form PrevClust, where to did they move?
						cl=i
						for(t in 1:nrclusters){
							TempComps=which(DissComps%in%names(which(Matrix[cl,]==t)))
							disscl=NULL
							if(length(TempComps)!=0){
								disscl=list()
								disscl[[1]]=DissComps[TempComps]
								disscl[[2]]=t
								names(disscl)=c("Cpds","Cluster")
							}
							Diss[[t]]=disscl
						}
					}
					if(length(Diss)!=0){
						r=c()
						for(l in 1:length(Diss)){
							if(is.null(Diss[[l]])){
								r=c(r,l)
							}
						}
						if(!is.null(r)){
							Diss=Diss[-r]
						}
						for(l in 1:length(Diss)){
							names(Diss)[l]=paste("Comps Diss To Cluster " ,Diss[[l]][[2]],sep="")

						}
					}
					temp[[2]]=Diss
					names(temp)[2]="Dissapeared"

					Joined=list()
					JoinComps=NULL
					if(length((which(!(temp[[1]]%in%PrevCluster))!=0))){
						JoinComps=temp[[1]][(which(!(temp[[1]]%in%PrevCluster)))]
					}
					if(!(is.null(JoinComps))){ #Dissapeared comps: what is missing form PrevClust, where to did they move?
						prevcl=i-1
						for(t in 1:nrclusters){
							TempComps=which(JoinComps%in%names(which(Matrix[prevcl,]==t)))
							joincl=NULL
							if(length(TempComps)!=0){
								joincl=list()
								joincl[[1]]=JoinComps[TempComps]
								joincl[[2]]=t
								names(joincl)=c("Cpds","Cluster")

							}
							Joined[[t]]=joincl
						}
					}
					if(length(Joined)!=0){
						r=c()
						for(l in 1:length(Joined)){
							if(is.null(Joined[[l]])){
								r=c(r,l)
							}
						}
						if(!is.null(r)){
							Joined=Joined[-r]
						}

						for(l in 1:length(Joined)){
							names(Joined)[l]=paste("Comps Joined From Cluster " ,Joined[[l]][[2]],sep="")

						}
					}
					temp[[3]]=Joined
					names(temp)[3]="Joined"
					PrevCluster=temp[[1]]
				}
				FoundCl[[i]]=temp
			}

		}

		FoundClusters[[length(FoundClusters)+1]]=FoundCl

		FoundClusters[[length(FoundClusters)+1]]=Matrix

		return(FoundClusters)
	}

	Found=ClusterDistribution.2(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,weightclust,names=names)

	Matrix=Found[[length(Found)]]
	Found=Found[-length(Found)]

	FoundCl=Found[[length(Found)]]
	Found=Found[-length(Found)]

	if(is.null(names)){
		for(j in 1:dim(Matrix)[1]){
			names[j]=paste("Method",j,1)
		}
	}

	ClusterNumber=NULL
	if(inherits(Selection,"numeric")){
		ClusterNumber=Selection
		Selection=colnames(Matrix)[which(Matrix[1,]==Selection)]

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


	if(selectionPlot==TRUE){


		if(followMaxComps==TRUE){
			lab1=c("Maximum of objects of original cluster together")
			labelcluster=c()
			for(z in 1:length(Found)){
				labelcluster=c(labelcluster,Found[[z]]$Cluster)
			}
		}
		else{lab1=c("Number of objects still in original cluster")}

		nrcluster=c()
		nrcomps=c()
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}

		if(is.null(ClusterNumber)){
			xl=c(0,length(Found)+0.5)
		}
		else{
			xl=c(1,length(Found)+0.5)
		}

		if(!(is.null(location))){
			location=paste(location,"_SelectionPlot.pdf",sep="")

		}
		plottypein(plottype,location)

		graphics::plot(type="n",x=0,y=0,xlim=xl,ylim=c(0,max(nrcluster,nrcomps)+legendposy+0),xlab="",ylab="",xaxt="n",yaxt="n",cex.lab=1.25)
		graphics::lines(x=seq(1,length(Found)),y=nrcluster,lty=1,col="red",lwd=1.5)
		graphics::points(x=seq(1,length(Found)),y=nrcluster,pch=19,col="red",cex=1.5)

		if(is.null(ClusterNumber)){
			graphics::lines(x=c(0,1),y=c(length(Selection),nrcomps[1]),lty=1,col="black",lwd=1.5)
			graphics::points(x=0,y=length(Selection),pch=19,col="black",cex=1.5)
		}

		graphics::lines(x=seq(1,length(Found)),y=nrcomps,lty=1,col="blue",lwd=1.5)
		graphics::points(x=seq(1,length(Found)),y=nrcomps,pch=19,col="blue",cex=1.5)

		if(is.null(ClusterNumber)){
			graphics::text(0,length(Selection), "S",cex=1.5,pos=1,col="black",font=2)
		}
		if(selectionPlot==TRUE & followMaxComps==TRUE)
			graphics::text(seq(1,length(Found)),nrcomps, labelcluster,cex=1.5,pos=1,col="black",font=2)

		if(is.null(ClusterNumber)){
			graphics::axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1)
		}
		else{
			graphics::axis(1,at=seq(1,length(Found)),labels=c(names),las=2,cex=1)
		}
		graphics::axis(2,at=seq(0,max(nrcluster,nrcomps)+2.0),labels=seq(0,max(nrcluster,nrcomps)+2.0),cex=1)


		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,lty=c(1,1,0),pch=c(19,19,0),col=c("blue","red","black"),legend=c(lab1,"Number of clusters original cluster divided amongst","Cluster number"),bty="n",cex=1.2)

		plottypeout(plottype)
	}
	if(completeSelectionPlot==TRUE){
		if(!(is.null(location))){
			location=paste(location,"_CompleteSelectionPlot.pdf",sep="")

		}
		plottypein(plottype,location)
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}



		graphics::plot(type="n",x=0,y=0,xlim=c(0,length(Found)+0.5),ylim=c(0,max(nrcluster,nrcomps)+legendposy+0),xlab="",ylab="Number of objects",xaxt="n",yaxt="n")

		xnext=c()
		ynext=c()
		colorsp=c()

		if(is.null(ClusterNumber)){
			p=seq(0,length(Found))
		}
		else{
			p=seq(0,length(Found))
		}

		for(m in p){
			if(m==0){
				if(!(is.null(ClusterNumber))){
					howmany=length(Found[[1]]$AllClusters)
				}
				else{
					howmany=1
				}
			}
			else{
				howmany=length(Found[[m]]$AllClusters)
			}

			if(m==0){
				xnext=0.1
				ynext=c(length(Selection))
				if(!(is.null(ClusterNumber))){
					colorsp=cols[ClusterNumber]
				}
				else{
					colorsp=c("black")
				}


				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				if(!(is.null(ClusterNumber))){
					labelcluster=ClusterNumber
				}
				else{
					labelcluster="S"
				}
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)

				L1=list()
				for(n in 1:howmany){
					if(is.null(ClusterNumber)){
						L1[[n]]=Selection
					}
					else{
						L1[[n]]=Found[[1]]$AllClusters[[n]][[3]]
					}
				}

			}
			else{
				xprev=xnext
				yprev=ynext
				ynext=c()
				colorsp=c()
				xnext=rep(seq(1,length(Found))[m],howmany)
				for(n in 1:howmany){
					ynext=c(ynext,length(Found[[m]]$AllClusters[[n]][[3]]))

					colorsp=c(colorsp,cols[Found[[m]]$AllClusters[[n]]$clusternumber])

					if(length(ynext)>1){
						for(t in 1:(length(ynext)-1)){
							if(ynext[t]==ynext[length(ynext)]){
								ynext[length(ynext)]=ynext[length(ynext)]-0.3
							}
						}
					}

					graphics::points(x=xnext[n],y=ynext[length(ynext)],col=colorsp[length(colorsp)],pch=19,cex=1.25)
					labelcluster=Found[[m]]$AllClusters[[n]]$clusternumber
					position=3
					if(!(is.integer(ynext[length(ynext)]))){
						position=1
					}
					graphics::text(xnext[n],ynext[length(ynext)],labelcluster,cex=1.5,pos=position,col="black",font=2)
				}


				L2=L1
				L1=list()
				for(n in 1:howmany){
					L1[[n]]=Found[[m]]$AllClusters[[n]][[3]]
				}

				for(q in 1:length(L1)){
					for(p in 1:length(L2)){
						if(length(which(L2[[p]] %in% L1[[q]])) != 0){
							graphics::segments(x0=xprev[p],y0=yprev[p],x1=xnext[q],y1=ynext[q],col=colorsp[q],lwd=2)
						}
					}
				}
			}
		}
		if(is.null(ClusterNumber)){
			graphics::axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1.5)
		}
		else{
			graphics::axis(1,at=seq(0,length(Found)-1),labels=c(names),las=2,cex=1.5)
		}
		graphics::axis(2,at=seq(0,max(nrcluster,nrcomps)),labels=seq(0,max(nrcluster,nrcomps)),cex=1.5)
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)


		plottypeout(plottype)
	}

	if(ClusterPlot==TRUE & !(is.null(ClusterNumber))){

		#FoundCl=Found[[length(Found)]]

		if(!(is.null(location))){
			location=paste(location,"_SelectionPlot.pdf",sep="")

		}
		plottypein(plottype,location)
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(FoundCl)){
			nrcomps=c(nrcomps,length(FoundCl[[j]][[1]]))
		}


		graphics::plot(type="n",x=0,y=0,xlim=c(0,length(FoundCl)-0.5),ylim=c(0,max(nrcluster,nrcomps)+legendposy+0),xlab="",ylab="Number of objects",xaxt="n",yaxt="n")

		xnext=c()
		ynext=c()
		colorsp=c()


		p=seq(1,length(FoundCl))

		for(m in p){

			if(m==1){
				xnext=0.1
				ynext=length(FoundCl[[1]][[1]])
				colorsp=cols[ClusterNumber]
				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)

				labelcluster=ClusterNumber

				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)


			}

			else{
				xprev=xnext
				yprev=ynext

				ynext=c()
				colorsp=c()

				xnext=m-1
				ynext=length(FoundCl[[m]][[1]])
				colorsp=cols[ClusterNumber]

				graphics::points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				graphics::segments(x0=xprev,y0=yprev,x1=xnext,y1=ynext,col=colorsp)
				labelcluster=ClusterNumber

				position=3
				if(!(is.integer(ynext[length(ynext)]))){
					position=1
				}
				graphics::text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)

				#Dissapeared objects
				if(length(FoundCl[[m]]$Dissapeared)!=0){

					xdiss=rep(xnext,length(FoundCl[[m]]$Dissapeared))
					ydiss=c()
					colorsd=c()
					labs=c()
					for(d in 1:length(FoundCl[[m]]$Dissapeared)){
						ydiss=c(ydiss,length(FoundCl[[m]]$Dissapeared[[d]][[1]]))
						colorsd=c(colorsd,cols[FoundCl[[m]]$Dissapeared[[d]][[2]]])
						labs=c(labs,FoundCl[[m]]$Dissapeared[[d]][[2]])
					}

					if(length(ydiss)>1){
						for(t in 1:(length(ydiss)-1)){
							if(ydiss[t]==ydiss[length(ydiss)]){
								ydiss[length(ydiss)]=ydiss[length(ydiss)]-0.3
							}
						}
					}

					graphics::points(x=xdiss,y=ydiss,col=colorsd,pch=19,cex=1.25)
					labelcluster=labs
					position=3
					if(!(is.integer(ydiss[length(ydiss)]))){
						position=1
					}
					graphics::text(xdiss,ydiss,labelcluster,cex=1.5,pos=position,col="black",font=2)

					for(p in 1:length(ydiss)){
						graphics::segments(x0=xprev,y0=yprev,x1=xdiss[p],y1=ydiss[p],col=colorsd[p])
					}

				}


				if(length(FoundCl[[m]]$Joined)!=0){

					xjoin=rep(xprev,length(FoundCl[[m]]$Joined))
					yjoin=c()
					colorsj=c()
					labs=c()
					for(d in 1:length(FoundCl[[m]]$Joined)){
						yjoin=c(yjoin,length(FoundCl[[m]]$Joined[[d]][[1]]))
						colorsj=c(colorsj,cols[FoundCl[[m]]$Joined[[d]][[2]]])
						labs=c(labs,FoundCl[[m]]$Joined[[d]][[2]])
					}

					if(length(yjoin)>1){
						for(t in 1:(length(yjoin-1))){
							if(yjoin[t]==yjoin[length(yjoin)]){
								yjoin[length(yjoin)]=yjoin[length(yjoin)]-0.3
							}
						}
					}

					graphics::points(x=xjoin,y=yjoin,col=colorsj,pch=19,cex=1.25)
					labelcluster=labs
					position=3
					if(!(is.integer(yjoin[length(yjoin)]))){
						position=1
					}
					graphics::text(xjoin,yjoin,labelcluster,cex=1.5,pos=position,col="black",font=2)

					for(p in 1:length(yjoin)){
						graphics::segments(x0=xjoin[p],y0=yjoin[p],x1=xnext,y1=ynext,col=colorsj[p])
					}

				}

			}

		}
		graphics::axis(1,at=c(0.1,seq(1,length(FoundCl)-1)),labels=c(names),las=2,cex=1.5)
		graphics::axis(2,at=seq(0,max(nrcomps)),labels=seq(0,max(nrcomps)),cex=1.5)
		graphics::legend(legendposx,max(nrcluster,nrcomps)+legendposy,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)

		plottypeout(plottype)

	}

	if(table==TRUE & selectionPlot==TRUE){
		SharedComps=Selection
		Extra=list()
		temp=c()
		for(a in 1:length(Found)){
			SharedComps=intersect(SharedComps,Found[[a]]$Complabels)
		}
		for(a in 1:length(Found)){
			Extra[[a]]=Found[[a]]$Complete.new.cluster[which(!(Found[[a]]$Complete.new.cluster%in%SharedComps))]
			names(Extra)[a]=names[a]
			if(followMaxComps==TRUE){
				names(Extra)[a]=paste(names(Extra)[a],"_",labelcluster[a],sep="")
			}
			temp=c(temp,length(Extra[[a]]))
		}

		ExtraOr=Selection[which(!(Selection%in%SharedComps))]

		collength=max(length(SharedComps),length(ExtraOr),temp)

		if(length(SharedComps)<collength){
			spaces=collength-length(SharedComps)
			SharedComps=c(SharedComps,rep(" ",spaces))
		}

		if(length(ExtraOr)<collength){
			spaces=collength-length(ExtraOr)
			ExtraOr=c(ExtraOr,rep(" ",spaces))
		}

		for(b in 1:length(Extra)){
			if(length(Extra[[b]])<collength){
				spaces=collength-length(Extra[[b]])
				Extra[[b]]=c(Extra[[b]],rep(" ",spaces))
			}

		}
		table1=data.frame(ExtraOr=ExtraOr,SharedComps=SharedComps)
		for(b in 2:length(Extra)){
			table1[1+b]=Extra[[b]]
			colnames(table1)[1+b]=names(Extra)[b]
		}

	}
	else{table1=list()}

	names(Found)=names
	Found[[length(Found)+1]]=table1
	names(Found)[length(Found)]="Table"
	return(Found)
}
