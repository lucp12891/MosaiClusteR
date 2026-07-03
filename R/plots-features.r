
#' @title Visualization of characteristic binary features of a single data set
#'
#' @description  A tool to visualize characteristic binary features of a set of objects in comparison with the remaining objects for a single data set. The result is a matrix with coloured cells. Columns represent
#' objects and rows represent the specified features. A feature which is present is give a coloured cell while an absent feature is represented by a grey cell. The labels on the right indicate the names of the features while the labels on the bottom are the names of the objects.
#' @param leadCpds A character vector with the names of the objects in a first group, i.e., the group for which the specified features are characteristic. Default is NULL.
#' @param orderLab A character vector with the order of the objects. Default is NULL.
#' @param features A character vector with the names of the features to be visualized. Default is NULL.
#' @param data The data matrix. Default is NULL.
#' @param colorLab Optional. A clustering object if the objects are to be coloured accoring to their clustering order. Default is NULL.
#' @param nrclusters Optional. The number of clusters to divide the dendrogram of ColorLab. Default is NULL.
#' @param cols Optional. A character vector with the colours of the different clusters. Default is NULL.
#' @param name A character string with the name of the data. Default is "Data".
#' @param colors1 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in LeadCpds. Default is c('gray90','blue').
#' @param colors2 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in the remaining objects. Default is c('gray90','green').
#' @param highlightFeat Optional. A character vector with names of features to be highlighted. The names of the features are coloured purple. Default is NULL.
#' @param margins A vector with the margings of the plot. Default is c(5.5,3.5,0.5,5.5).
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot visualizing the presence and absence of the characteristic binary features.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' Comps=FindCluster(list(MCF7_F),nrclusters=10,select=c(1,8))
#'
#' MCF7_Char=CharacteristicFeatures(List=list(fingerprintMat),Selection=Comps,
#' binData=list(fingerprintMat),datanames=c("FP"),nrclusters=NULL,topC=NULL,
#' sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=c("FP"))
#' Feat=MCF7_Char$Selection$Characteristics$FP$TopFeat$Names[c(1:10)]
#'
#' BinFeaturesPlot_SingleData(leadCpds=Comps,orderLab=MCF7_Char$Selection$
#' objects$OrderedCpds,features=Feat,data=fingerprintMat,colorLab=NULL,
#' nrclusters=NULL,cols=NULL,name=c("FP"),colors1=c('gray90','blue'),colors2=
#' c('gray90','green'),highlightFeat=NULL,margins=c(5.5,3.5,0.5,5.5),
#' plottype="new",location=NULL)
#' }
#' @export BinFeaturesPlot_SingleData
BinFeaturesPlot_SingleData<-function(leadCpds=c(),orderLab=c(),features=c(),data=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Data"),colors1=c('gray90','blue'),
		colors2=c('gray90','green'),highlightFeat=NULL,margins=c(5.5,3.5,0.5,5.5),plottype="new",location=NULL){

	if(all(leadCpds%in%rownames(data))){
		data=t(data)
	}

	if(!is.null(orderLab)){
		if(inherits(orderLab,"character")){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			data=data[,match(orderlabs,colnames(data))]

		}
	}
	else{
		orderlabs=colnames(data)
	}

	temp=orderlabs[which(!(orderlabs%in%leadCpds))]
	AllCpds=c(leadCpds,temp)


	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			grDevices::pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			grDevices::dev.new()
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			grDevices::dev.off()
		}
	}

	x<-c(1:length(AllCpds)) #x=comps
	y<-c(1:length(features)) #y=feat
	PlotData<-t(data[as.character(features),AllCpds,drop=FALSE])
	plottypein(plottype,location)
	graphics::par(mar=margins)
	graphics::image(x,y,PlotData,col=colors1,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	if(length(unique(as.vector(PlotData[1:length(leadCpds),])))==1){
		if(unique(as.vector(PlotData[1:length(leadCpds),]))==1){
			PlotData[1:length(leadCpds),]=2
			colors2=c("gray90","blue","green")
			graphics::image(x,y,PlotData,col=colors2,xlab="",axes=FALSE,add=TRUE,ann=FALSE,xaxt='n')
		}
		else{
			colors2=c("gray90","blue")
			graphics::image(x,y,PlotData,col=colors2,xlab="",axes=FALSE,add=TRUE,ann=FALSE,xaxt='n')
		}
	}
	else{
		graphics::image(x[1:length(leadCpds)],y,PlotData[1:length(leadCpds),,drop=FALSE],col=colors2,add=TRUE,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	}

	if(!(is.null(colorLab)) & !is.null(nrclusters)){
		Data1 <- colorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters)

		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab

		ClustData1=ClustData1[Data1$order]


		order=seq(1,nrclusters)

		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}

		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)
	}
	else{
		colors1<-rep("green",length(leadCpds))
		colors2<-rep("black",length(temp))
		colors=c(colors1,colors2)
		names(colors)=AllCpds
	}
	if(!is.null(highlightFeat)){
		colfeat=rep("black",ncol(PlotData))
		colfeat[which(colnames(PlotData)%in%highlightFeat)]="purple"
		graphics::mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0.2, las=2,cex=0.8,col=colfeat)
	}
	else{
		graphics::mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0.2, las=2,cex=0.8)
	}
	graphics::mtext(name, side = 2,  line=1, las=0, cex=1)
	graphics::mtext(rownames(PlotData), side = 1, at= c(1:nrow(PlotData)), line=0.2, las=2, cex=0.8,col=colors[AllCpds])
	plottypeout(plottype)
}

#' @title Visualization of characteristic binary features of multiple data sets
#'
#' @description A tool to visualize characteristic binary features of a set of objects in comparison with the remaining objects for multiple data sets. The result is a matrix with coloured cells. Columns represent
#' objects and rows represent the specified features. A feature which is present is give a coloured cell while an absent feature is represented by a grey cell. The labels on the right indicate the names of the features while the labels on the bottom are the names of the objects.
#' @param leadCpds A character vector with the names of the objects in a first group, i.e., the group for which the specified features are characteristic. Default is NULL.
#' @param orderLab A character vector with the order of the objects. Default is NULL.
#' @param features A list with as elements character vectors with the names of the features to be visualized for each data set. Default is NULL.
#' @param data A list with the different data sets. Default is NULL.
#' @param validate Optional. A list with validation data sets. If a feature has a validation reference, these are added in a red colour. Default is NULL.
#' @param colorLab Optional. A clustering object if the objects are to be coloured accoring to their clustering order. Default is NULL.
#' @param nrclusters Optional. The number of clusters to divide the dendrogram of ColorLab. Default is NULL.
#' @param cols Optional. A character vector with the colours of the different clusters. Default is NULL.
#' @param name A character string with the names of the data sets. Default is c("Data1" ,"Data2") for two data sets.
#' @param colors1 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in LeadCpds. Default is c('gray90','blue').
#' @param colors2 A character vector with the colours to indicate the presence (first element) or the absence of the features for the objects in the remaining objects. Default is c('gray90','green').
#' @param margins A vector with the margings of the plot. Default is c(5.5,3.5,0.5,5.5).
#' @param cexB The font size of the labels on the bottom: the object labels. Default is 0.80.
#' @param cexL The font size of the labels on the left: the data labels. Default is 0.80.
#' @param cexR The font size of the labels on the right: the feature labels. Default is 0.80.
#' @param spaceNames A percentage of the height of the figure to be reserved for the names of the objects. Default is 0.20.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot visualizing the presence and absence of the characteristic binary features across multiple data sets.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' Comps=FindCluster(list(MCF7_F),nrclusters=10,select=c(1,8))
#'
#' MCF7_Char=CharacteristicFeatures(List=NULL,Selection=Comps,binData=
#' list(fingerprintMat,targetMat),datanames=c("FP","TP"),nrclusters=NULL,
#' topC=NULL,sign=0.05,fusionsLog=TRUE,weightclust=TRUE,names=c("FP","TP"))
#'
#' FeatFP=MCF7_Char$Selection$Characteristics$FP$TopFeat$Names[c(1:10)]
#' FeatTP=MCF7_Char$Selection$Characteristics$TP$TopFeat$Names[c(1:10)]
#'
#' BinFeaturesPlot_MultipleData(leadCpds=Comps,orderLab=MCF7_Char$Selection$
#' objects$OrderedCpds,features=list(FeatFP,FeatTP),data=list(fingerprintMat,targetMat),
#' validate=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("FP","TP"),colors1=
#' c('gray90','blue'),colors2=c('gray90','green'),margins=c(5.5,3.5,0.5,5.5),cexB=0.80,
#' cexL=0.80,cexR=0.80,spaceNames=0.20,plottype="new",location=NULL)
#' }
#' @export BinFeaturesPlot_MultipleData
BinFeaturesPlot_MultipleData<-function(leadCpds,orderLab,features=list(),data=list(),validate=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,name=c("Data1" ,"Data2"),colors1=c('gray90','blue'),colors2=c('gray90','green'),margins=c(5.5,3.5,0.5,5.5),cexB=0.80,cexL=0.80,cexR=0.80,spaceNames=0.20,plottype="new",location=NULL){

	if (!requireNamespace("plotrix", quietly = TRUE)) {
		stop("BinFeaturesPlot_MultipleData() requires the suggested package 'plotrix'.")
	}

	if(all(leadCpds%in%rownames(data))){
		Data=t(Data)
	}

	if(!is.null(orderLab)){
		if(inherits(orderLab,"character")){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			data=data[,match(orderlabs,colnames(data))]
		}
	}
	else{
		orderlabs=rownames(data[[1]])
	}



	temp=orderlabs[which(!(orderlabs%in%leadCpds))]
	AllCpds=c(leadCpds,temp)


	if(!is.null(validate)){
		for(v in 1:length(validate)){
			validate[[v]]=validate[[v]][AllCpds,]
		}
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

	x<-c(1:length(AllCpds)) #x=comps

	Data_new=list()
	for(j in 1:length(data)){
		tempD=data[[j]]
		Data_temp=tempD[AllCpds,as.character(features[[j]]),drop=FALSE]
		Data_new[[j]]=as.matrix(Data_temp)
#		 if(ncol(Data_new[[j]])==1){
#			 colnames(Data_new[[j]])=Features[[j]][which(as.character(Features[[j]])%in%colnames(tempD))]
#			 rownames(Data_new[[j]])=rownames(Data[[1]])
#		 }
	}
	names(Data_new)=name

	Draw=matrix(0,ncol=(length(name)+2),nrow=length(unique(unlist(features))))
	i=0
	for(f in unique(unlist(features))){
		Row=c(f,rep(0,length(name)+1))
		for(j in 1:length(Data_new)){
			if(f%in%colnames(Data_new[[j]])){
				Row[j+1]=name[j]
			}
		}
		Nr=length(which((Row[c(2:(ncol(Draw)-1))])!="0"))
		Row[ncol(Draw)]=Nr
		i=i+1
		Draw[i,]=Row
	}
	Draw=Draw[order(as.numeric(Draw[,ncol(Draw)]),decreasing=TRUE),,drop=FALSE]


	Temp_Draw=Draw[,-c(1,ncol(Draw)),drop=FALSE]
	Fin_Draw=c()
	for(a in 1:nrow(Temp_Draw[!duplicated(Temp_Draw),,drop=FALSE])){
		example=Temp_Draw[!duplicated(Temp_Draw),,drop=FALSE][a,,drop=FALSE]
		join=c()
		for(b in 1:nrow(Temp_Draw)){
			if(all(Temp_Draw[b,]==example)){
				join=c(join,b)
			}
		}

		t=paste(Draw[join,1],collapse=",")
		Fin_Draw=rbind(Fin_Draw,cbind(t,paste(example,collapse=",")))
	}

	plottypein(plottype,location)
	Mat <- matrix(c(1:(nrow(Fin_Draw)+1)),nrow = (nrow(Fin_Draw)+1),ncol = 1,byrow = TRUE)
	H=rep(0,nrow(Fin_Draw))
	Nrtotal=length(unlist(features))
	for(m in 1:nrow(Fin_Draw)){
		Methods=unlist(strsplit(Fin_Draw[m,2],","))
		NrofMethods=length(which(Methods!="0"))
		if(NrofMethods==1){
			Nrtotal=Nrtotal+1
		}
	}

	for(m in 1:nrow(Fin_Draw)){
		NrofFeat=length(unlist(strsplit(Fin_Draw[m,1],",")))
		Methods=unlist(strsplit(Fin_Draw[m,2],","))
		NrofMethods=length(which(Methods!="0"))
		if(!is.null(validate)){
			NrofMethods=NrofMethods+length(validate)
		}
		if(NrofFeat*NrofMethods==1){
			NrofMethods=1.5
		}
		H[m]=((NrofFeat*NrofMethods)/Nrtotal)*(1-spaceNames)

	}
	graphics::layout(mat = Mat,heights=c(H,spaceNames))

	if(length(validate)>0){
		Data_new=c(Data_new,validate)
	}

	for(n in 1:nrow(Fin_Draw)){
		print(n)
		Shared=unlist(strsplit(Fin_Draw[n,1],","))
		Datasets=unlist(strsplit(Fin_Draw[n,2],","))
		V=0
		if(!is.null(validate)){
			V=length(validate)
			Datasets=c(Datasets,names(validate))
		}

		if(any("0"%in%Datasets)){
			Datasets=Datasets[-which(Datasets==0)]
		}

		ImageData=c()
		print(Datasets)
		for(d in Datasets){
			if(all(Shared%in%colnames(Data_new[[d]]))){
				ImageData=cbind(ImageData,Data_new[[d]][,Shared,drop=FALSE])
			}
			else{
				ImageData=cbind(ImageData,Data_new[[d]][,Shared[which(Shared%in%colnames(Data_new[[d]]))],drop=FALSE])
				Fill=Shared[which(!Shared%in%colnames(Data_new[[d]]))]
				for(f in 1:length(Fill)){
					Data_new[[d]]=cbind(Data_new[[d]],rep(0,nrow(Data_new[[d]])))
					colnames(Data_new[[d]])[ncol(Data_new[[d]])]=Fill[f]
				}
				ImageData=cbind(ImageData,Data_new[[d]][,Shared,drop=FALSE])
			}
		}


		if(V>0){
			for(i in 1:(length(Shared)*length(validate))){
				ImageData[,(ncol(ImageData)-(i-1))][which(ImageData[,(ncol(ImageData)-(i-1))]==1)]=3
			}
		}

		ImageData=as.matrix(ImageData)
		ImageData=ImageData[,order(colnames(ImageData)),drop=FALSE]
		colnames(ImageData)=paste(colnames(ImageData),rep(Datasets,length(Shared)),sep="_")
		rownames(ImageData)=rownames(Data_new[[1]])
		#ImageData=ImageData[,c(ncol(ImageData):1),drop=FALSE]
		Colors=matrix(0,nrow(ImageData),ncol(ImageData))
		for(nr in 1:nrow(ImageData)){
			for(nc in 1:ncol(ImageData)){
				if(rownames(ImageData)[nr]%in%leadCpds){
					if(ImageData[nr,nc]==1){
						Colors[nr,nc]="green"
					}
					else if(ImageData[nr,nc]==3){
						Colors[nr,nc]=grDevices::adjustcolor("red", alpha.f = 0.3)
					}
					else{
						Colors[nr,nc]="grey90"
					}
				}
				else{
					if(ImageData[nr,nc]==1){
						Colors[nr,nc]="blue"
					}
					else if(ImageData[nr,nc]==3){
						Colors[nr,nc]=grDevices::adjustcolor("red", alpha.f = 0.3)
					}
					else{
						Colors[nr,nc]="grey90"
					}
				}
			}
		}

		graphics::par(mar=margins)


		plotrix::color2D.matplot(t(ImageData),cellcolors=t(Colors),show.values=FALSE,axes=FALSE,xlab="",ylab="",border=NA)

		ColorsF=rep("black",ncol(ImageData))
		names(ColorsF)=colnames(ImageData)

		if(!(is.null(colorLab)) & !is.null(nrclusters)){
			Data1 <- colorLab$Clust
			ClustData1=stats::cutree(Data1,nrclusters)

			ordercolors=ClustData1[Data1$order]
			names(ordercolors)=Data1$order.lab

			ClustData1=ClustData1[Data1$order]


			order=seq(1,nrclusters)

			for (k in 1:length(unique(ClustData1))){
				select=which(ClustData1==unique(ClustData1)[k])
				ordercolors[select]=order[k]
			}

			colors<- cols[ordercolors]
			names(colors) <-names(ordercolors)
		}
		else{
			Colors1<-rep("green",length(leadCpds))
			Colors2<-rep("black",length(temp))
			Colors=c(Colors1,Colors2)
			names(Colors)=AllCpds

		}

		graphics::mtext(colnames(ImageData), side = 4, at= c(ncol(ImageData):1), line=0.2, las=2,cex=cexR,col=ColorsF)
		if(length(Datasets)==length(Data_new)){
			graphics::mtext("All", side = 2,  line=1, las=1, cex=2)
		}
		else{
			if(length(Datasets)>2){
				N=c()
				for(w in 1:length(Datasets)){
					if(w%%2==0&w!=length(Datasets)){
						N=c(N," & ",Datasets[w],"\n")
					}
					else if(w==1){

						N=c(N,Datasets[w])
					}
					else{
						N=c(N," & ",Datasets[w])
					}
				}
				graphics::mtext(paste(N,collapse=""), side = 2,  line=1, las=1, cex=cexL)
			}
			else{
				graphics::mtext(paste(Datasets,collapse=" & "), side = 2,  line=1, las=1, cex=cexL)
			}
		}

		if(n==nrow(Fin_Draw)){
			graphics::mtext(rownames(ImageData), side = 1, at= c(1:nrow(ImageData)), line=0.2, las=2, cex=cexB,col=Colors[AllCpds])
		}

	}

	plottypeout(plottype)
}

#' @title Plot of continuous features
#'
#' @description The function \code{ContFeaturesPlot} plots the values of continuous features.
#' It is possible to separate between objects of interest and the
#' other objects.
#' @param leadCpds A character vector containing the objects one wants to
#' separate from the others.
#' @param data The data matrix.
#' @param nrclusters Optional. The number of clusters to consider if colorLab
#' is specified. Default is NULL.
#' @param orderLab Optional. If the objects are to set in a specific order of
#' a specific method. Default is NULL.
#' @param colorLab The clustering result that determines the color of the
#' labels of the objects in the plot. If NULL, the labels are black. Default is NULL.
#' @param cols The colors for the labels of the objects. Default is NULL.
#' @param ylab The lable of the y-axis. Default is "features".
#' @param addLegend Logical. Indicates whether a legend should be added to the
#' plot. Default is TRUE.
#' @param margins Optional. Margins to be used for the plot. Default is c(5.5,3.5,0.5,8.7).
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot in which the values of the features of the leadCpds are
#' separeted from the others.
#' @examples
#' \dontrun{
#' data(Colors1)
#' Comps=c("Cpd1", "Cpd2", "Cpd3", "Cpd4", "Cpd5")
#'
#' Data=matrix(sample(15, size = 50*5, replace = TRUE), nrow = 50, ncol = 5)
#' colnames(Data)=colnames(Data, do.NULL = FALSE, prefix = "col")
#' rownames(Data)=rownames(Data, do.NULL = FALSE, prefix = "row")
#' for(i in 1:50){
#' 	rownames(Data)[i]=paste("Cpd",i,sep="")
#' }
#'
#' ContFeaturesPlot(leadCpds=Comps,orderLab=rownames(Data),colorLab=NULL,data=Data,
#' nrclusters=7,cols=Colors1,ylab="features",addLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),
#' plottype="new",location=NULL)
#' }
#' @export ContFeaturesPlot
ContFeaturesPlot<-function(leadCpds,data,nrclusters=NULL,orderLab=NULL,colorLab=NULL,cols=NULL,ylab="features",addLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),plottype="new",location=NULL){

	if(all(leadCpds%in%rownames(data))){
		data=t(data)
	}

	if(!is.null(orderLab)){
		if(inherits(orderLab,"character")){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			data=data[,match(orderlabs,colnames(data))]
		}
	}
	else{
		orderlabs=colnames(data)
	}

	temp=orderlabs[which(!(orderlabs%in%leadCpds))]
	AllCpds=c(leadCpds,temp)

	if(is.null(dim(data))){
		data=t(as.matrix(data))
		rownames(data)="Feature"
	}
	data=data[,AllCpds,drop=FALSE]

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
	graphics::par(mar=margins)
	graphics::plot(x=0,y=0,xlim=c(0,(ncol(data)+3)),ylim=c(min(data)-0.5,max(data)+0.5),type="n",ylab=ylab,xlab='',xaxt='n')
	for(i in c(1:nrow(data))){
		graphics::lines(x=seq(1,length(leadCpds)),y=data[i,which(colnames(data)%in%leadCpds)],col=i)
	}
	for(i in c(1:nrow(data))){
		graphics::lines(x=seq(length(leadCpds)+4,(ncol(data)+3)),y=data[i,which(!(colnames(data)%in%leadCpds))],col=i)
	}



	if(!(is.null(colorLab)) | is.null(nrclusters)){
		Data1 <- colorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters)

		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab

		ClustData1=ClustData1[Data1$order]


		order=seq(1,nrclusters)

		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}

		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)
	}
	else{
		colors1<-rep("green",length(leadCpds))
		colors2<-rep("black",length(temp))
		colors=c(colors1,colors2)
		names(colors)=AllCpds
	}

	graphics::mtext(leadCpds,side=1,at=seq(1,length(leadCpds)),line=0.6,las=2,cex=0.70,col=colors[leadCpds])
	graphics::mtext(temp, side = 1, at=c(seq(length(leadCpds)+4,(ncol(data)+3))), line=0.5, las=2, cex=0.70,col=colors[temp])
	if(addLegend==TRUE){

		labels=rownames(data)
		colslegend=seq(1,length(rownames(data)))

		graphics::par(xpd=T,mar=margins)
		graphics::legend(ncol(data)+5,mean(c(min(c(min(data)-0.5,max(data)+0.5)),max(c(min(data)-0.5,max(data)+0.5)))),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=0.8)

	}
	plottypeout(plottype)
}

#' @title A GO plot of a pathway analysis output.
#'
#' @description The \code{PlotPathways} function takes an output of the
#' \code{PathwayAnalysis} function and plots a GO graph with the help of the
#' \code{plotGOgraph} function of the MLP package.
#' @param Pathways One element of the output list returned by
#' \code{PathwayAnalysis} or \code{Geneset.intersect}.
#' @param nRow Number of GO IDs for which to produce the plot. Default is 5.
#' @param main Title of the plot. Default is NULL.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return The output is a GO graph.
#' @examples
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
#' MCF7_PathsFandT=PathwayAnalysis(List=L, geneExpr = geneMat, nrclusters = 7, method = c("limma",
#' "MLP"), geneInfo = GeneInfo, geneSetSource = "GOBP", topP = NULL,
#' topG = NULL, GENESET = NULL, sign = 0.05,niter=2,fusionsLog = TRUE, weightclust = TRUE,
#'  names =names,seperatetables=FALSE,separatepvals=FALSE)
#'
#' PlotPathways(MCF7_PathsFandT$FP$"Cluster 1"$Pathways,nRow=5,main=NULL)
#' }
#' @export PlotPathways
PlotPathways<-function(Pathways,nRow=5,main=NULL,plottype="new",location=NULL){

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
	#preparing data structure for the plotGOGraph
	colnames(Pathways)[3:ncol(Pathways)]=sub("mean_","",colnames(Pathways)[3:ncol(Pathways)])
	#plot GOgraph
	plottypein(plottype,location)
	MLP::plotGOgraph(Pathways,nRow=nRow,main=main)
	plottypeout(plottype)

}

#' @title Plotting gene profiles
#'
#' @description In \code{ProfilePlot}, the gene profiles of the significant genes for a
#' specific cluster are shown on 1 plot. Therefore, each gene is normalized by
#' subtracting its the mean.
#'
#' @param Genes The genes to be plotted.
#' @param Comps The objects to be plotted or to be separated from the other
#' objects.
#' @param geneExpr The gene expression matrix or ExpressionSet of the objects.
#' @param raw Logical. Should raw p-values be plotted? Default is FALSE.
#' @param orderLab Optional. If the objects are to set in a specific order of
#' a specific method. Default is NULL.
#' @param colorLab The clustering result that determines the color of the
#' labels of the objects in the plot. Default is NULL.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' @param cols Optional. The color to use for the objects in Clusters for each
#' method.
#' @param addLegend Optional. Whether a legend of the colors should be added to
#' the plot.
#' @param margins Optional. Margins to be used for the plot. Default is margins=c(8.1,4.1,1.1,6.5).
#' @param extra The space between the plot and the legend. Default is 5.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot which contains multiple gene profiles. A distinction is made
#' between the values for the objects in Comps and the others.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(geneMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#'		method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#'		method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' L=list(MCF7_F,MCF7_T)
#' names=c('FP','TP')
#'
#' MCF7_FT_DE = DiffGenes(List=L,geneExpr=geneMat,nrclusters=7,method="limma",sign=0.05,topG=10,
#' fusionsLog=TRUE,weightclust=TRUE)
#'
#' Comps=SharedComps(list(MCF7_FT_DE$`Method 1`$"Cluster 1",MCF7_FT_DE$`Method 2`$"Cluster 1"))[[1]]
#'
#' MCF7_SharedGenes=FindGenes(dataLimma=MCF7_FT_DE,names=c("FP","TP"))
#'
#' Genes=names(MCF7_SharedGenes[[1]])[-c(2,4,5)]
#'
#' colscl=ColorPalette(colors=c("red","green","purple","brown","blue","orange"),ncols=9)
#'
#' ProfilePlot(Genes=Genes,Comps=Comps,geneExpr=geneMat,raw=FALSE,orderLab=MCF7_F,
#' colorLab=NULL,nrclusters=7,cols=colscl,addLegend=TRUE,margins=c(16.1,6.1,1.1,13.5),
#' extra=4,plottype="sweave",location=NULL)
#' }
#' @export ProfilePlot
ProfilePlot<-function(Genes,Comps,geneExpr=NULL,raw=FALSE,orderLab=NULL,colorLab=NULL,nrclusters=NULL,cols=NULL,addLegend=TRUE,margins=c(8.1,4.1,1.1,6.5),extra=5,plottype="new",location=NULL){
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



	if(inherits(geneExpr,"ExpressionSet")){
		if (!requireNamespace("Biobase", quietly = TRUE)) {
			stop("ProfilePlot() requires the suggested package 'Biobase'.")
		}
		geneExpr <- Biobase::exprs(geneExpr)

	}

	if(!is.null(orderLab)){
		if(inherits(orderLab,"character")){
			orderlabs=orderLab
		}
		else{
			orderlabs=orderLab$Clust$order.lab
			geneExpr=geneExpr[,match(orderlabs,colnames(geneExpr))]
		}
	}
	else{
		orderlabs=colnames(geneExpr)
	}

	if(!is.null(colorLab)){
		Data1 <- colorLab$Clust
		ClustData1=stats::cutree(Data1,nrclusters)

		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab

		ClustData1=ClustData1[Data1$order]

		order=seq(1,nrclusters)

		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}

		if(!is.null(orderLab)){
			if(inherits(orderLab,"character")){
				ordernames=orderLab
			}
			else{
				ordernames=orderLab$Clust$order.lab
			}
			ordercolors=ordercolors[ordernames]
		}

		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)

	}
	else{
		#colors1<-rep("green",length(Comps))
		#colors2<-rep("black",length(orderlabs[which(!(orderlabs%in%Comps))]))
		colors1<-rep("blue",length(Comps))
		colors2<-rep("red",length(orderlabs[which(!(orderlabs%in%Comps))]))
		colors=c(colors1,colors2)
		AllCpds=c(Comps,orderlabs[which(!(orderlabs%in%Comps))])
		names(colors)=AllCpds
	}

	#yvalues=c()
	#allvalues=c()
	#for(i in 1:length(Genes)){
	#	yvalues=as.vector(GeneExpr[which(rownames(GeneExpr)==Genes[i]),])
	#	allvalues=c(allvalues,yvalues-mean(yvalues))
	#}
	yvalues=geneExpr[Genes,]
	if(raw==FALSE & !inherits(yvalues,"numeric")){
		allvalues=as.vector(apply(yvalues,1,function(c) c-mean(c)))
	}
	else if(raw==FALSE & inherits(yvalues,"numeric")){
		allvalues=as.vector(sapply(yvalues,function(c) c-mean(yvalues)))
	}
	else{
		allvalues=as.vector(yvalues)
	}
	ylims=c(min(allvalues)-0.1,max(allvalues)+0.1)

	plottypein(plottype,location)
	graphics::par(mar=margins,xpd=TRUE)
	graphics::plot(type="n",x=0,y=0,xlim=c(0,ncol(geneExpr)),ylim=ylims,ylab=expression(log[2] ~ paste("fold ", "change")),xlab=" ",xaxt="n",cex.axis=1.5,cex.lab=2)
	#ylims=c()
	Indices=c(colnames(geneExpr)[which(colnames(geneExpr)%in%Comps)],colnames(geneExpr)[which(!colnames(geneExpr)%in%Comps)])

	for(i in 1:length(Genes)){
		GenesComps=as.numeric(geneExpr[which(rownames(geneExpr)==Genes[i]),colnames(geneExpr)%in%Comps])
		Others=as.numeric(geneExpr[which(rownames(geneExpr)==Genes[i]),!(colnames(geneExpr)%in%Comps)])
		if(length(Others)==0){
			Continue=FALSE
		}else{Continue=TRUE}

		#ylims=c(ylims,c(GenesComps,Others)-mean(c(GenesComps,Others)))
		if(raw==FALSE){
			yvalues1=GenesComps-mean(c(GenesComps,Others))
		}
		else{
			yvalues1=GenesComps
		}

		graphics::lines(x=seq(1,length(GenesComps)),y=yvalues1,lty=1,col=i,lwd=1.6)
		#points(x=seq(1,length(GenesComps)),y=yvalues1,pch=19,col=i)
		graphics::segments(x0=1,y0=mean(yvalues1[1:length(GenesComps)]),x1=length(GenesComps),y1=mean(yvalues1[1:length(GenesComps)]),lwd=1.5,col=i)


		if(Continue==TRUE){
			if(raw==FALSE){
				yvalues2=Others-mean(c(GenesComps,Others))
			}
			else{
				yvalues2=Others
			}


			graphics::lines(x=seq(length(GenesComps)+1,ncol(geneExpr)),y=yvalues2,lty=1,col=i,lwd=1.6)
			graphics::segments(x0=length(GenesComps)+1,y0=mean(yvalues2[1:length(Others)]),x1=ncol(geneExpr),y1=mean(yvalues2[1:length(Others)]),lwd=1.5,col=i)

		}


	}
	#Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	if(!is.null(colorLab)){
		graphics::axis(1, labels=FALSE)
		#box("outer")
		graphics::mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(geneExpr)-0.5)), line=0.2, las=2, cex=1.5,col=colors[Indices])
		#mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(GeneExpr)-0.5)), line=0.2, las=2, cex=0.70,col=c(rep("blue",7),rep("black",(56-7))))
	}
	else{
		#axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)),labels=Indices,las=2,cex.axis=0.70,xlab=" ",col=colors[Indices])
		#axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)), labels=FALSE)
		graphics::mtext(Indices, side = 1,  at=seq(0.5,(ncol(geneExpr)-0.5)),line=0.2, las=2, cex=1.5,col=colors[Indices])
	}
	graphics::axis(2,ylab=expression(log[2] ~ paste("fold ", "change")),cex=2,cex.axis=1.5)

	if(addLegend==TRUE){

		labels=Genes
		colslegend=seq(1,length(Genes))

		graphics::par(xpd=T,mar=margins)
		graphics::legend(ncol(geneExpr)+extra,max(ylims),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=1.5)

	}
	plottypeout(plottype)
}
