# =============================================================================
# plots-heatmap.r  - heatmap / box-plot comparison utilities
# =============================================================================

#' @title Box plots of one distance matrix categorized against another distance
#' matrix.
#'
#' @description Given two distance matrices, the function categorizes one distance matrix
#' and produces a box plot from the other distance matrix against the created
#' categories. The option is available to choose one of the plots or to have
#' both plots. The function also works on outputs from ADEC and CEC functions
#' which do not have distance matrices but incidence matrices.
#'
#'
#' @param Data1 The first data matrix, cluster outcome or distance matrix to be
#' plotted.
#' @param Data2 The second data matrix, cluster outcome or distance matrix to
#' be plotted.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".
#' @param distmeasure A vector of the distance measures to be used on each data matrix. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to c("tanimoto","tanimoto").
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is c(NULL,NULL) for two data sets.
#' @param lab1 The label to plot for Data1.
#' @param lab2 The label to plot for Data2.
#' @param limits1 The limits for the categories of Data1.
#' @param limits2 The limits for the categories of Data2.
#' @param plot The type of plots: 1 - Plot the values of Data1 versus the
#' categories of Data2. 2 - Plot the values of Data2 versus the categories of
#' Data1. 3 - Plot both types 1 and 2.
#' @param StopRange Logical. Indicates whether the distance matrices with
#' values not between zero and one should be standardized to have so. If FALSE
#' the range normalization is performed. See \code{Normalization}. If TRUE, the
#' distance matrices are not changed. This is recommended if different types of
#' data are used such that these are comparable. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return One/multiple box plots.
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
#'
#' BoxPlotDistance(MCF7_F,MCF7_T,type="cluster",lab1="FP",lab2="TP",limits1=c(0.3,0.7),
#' limits2=c(0.3,0.7),plot=1,StopRange=FALSE,plottype="new", location=NULL)
#' }
#' @export BoxPlotDistance
BoxPlotDistance<-function(Data1,Data2,type=c('data','dist','clusters'),distmeasure=c("tanimoto","tanimoto"),normalize=c(FALSE,FALSE),method=c(NULL,NULL),lab1,lab2,limits1=NULL,limits2=NULL,plot=1,StopRange=FALSE,plottype="new",location=NULL){
	if (!requireNamespace("ggplot2", quietly=TRUE)) stop("BoxPlotDistance() requires the suggested package 'ggplot2'.")
	if (!requireNamespace("gridExtra", quietly=TRUE)) stop("BoxPlotDistance() requires the suggested package 'gridExtra'.")
	type<-match.arg(type)
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

	C1=D1=C2=D2=NULL
	if(type=='clusters'){

		Dist1<-Data1$DistM
		Dist2<-Data2$DistM


	}
	else if(type=='data'){
		Dist1<-Distance(Data1,distmeasure[1],normalize[1],method[1])
		Dist2<-Distance(Data2,distmeasure[2],normalize[2],method[2])
		DistL=list(Dist1,Dist2)
		for(i in 1:2){
			if(StopRange==FALSE & !(0<=min(DistL[[i]]) & max(DistL[[i]])<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				DistL[[i]]=Normalization(DistL[[i]],method="Range")
			}
		}

	}
	else if(type=='dist'){
		Dist1=Data1
		Dist2=Data2
		DistL=list(Dist1,Dist2)
		for(i in 1:2){
			if(StopRange==FALSE &  !(0<=min(DistL[[i]]) & max(DistL[[i]])<=1)){
				message("It was detected that a distance matrix had values not between zero and one. Range Normalization was performed to secure this. Put StopRange=TRUE if this was not necessary")
				DistL[[i]]=Normalization(DistL[[i]],method="Range")
			}
		}
	}

	OrderNames=rownames(Dist1)
	Dist2=Dist2[OrderNames,OrderNames]

	Dist1lower <- Dist1[lower.tri(Dist1)]
	Dist2lower <- Dist2[lower.tri(Dist2)]

	Categorize<-function(Distlower,limits){
		Cat=c(rep(0,length(Distlower)))
		for(j in 1:(length(limits)+1)){
			if(j==1){
				Cat[Distlower<=limits[j]]=j
			}
			else if(j<=length(limits)){
				Cat[Distlower>limits[j-1] & Distlower<=limits[j]]=j
			}
			else{
				Cat[Distlower>limits[j-1]]=j
			}
		}
		Cat<-factor(Cat)
		return(Cat)

	}

	#plot2
	if(!(is.null(limits1))){
		Dist1cat<-Categorize(Dist1lower,limits1)

		dataBox2<-data.frame(D2=Dist2lower,C1=Dist1cat)
		p2<-ggplot2::ggplot(dataBox2,ggplot2::aes(factor(C1),D2)) #x,y
		p2<-p2+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+ggplot2::xlab(lab1)+ggplot2::ylab(lab2)
	}
	#plot1
	if(!(is.null(limits2))){
		Dist2cat<-Categorize(Dist2lower,limits2)
		dataBox1<-data.frame(D1=Dist1lower,C2=Dist2cat)
		p1<-ggplot2::ggplot(dataBox1,ggplot2::aes(factor(C2),D1)) #x,y
		p1<-p1+ggplot2::geom_boxplot(outlier.shape=NA)+ggplot2::geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+ggplot2::xlab(lab2)+ggplot2::ylab(lab1)

	}
	if(plot==3){
		if(plottype=="pdf"){
			location=paste(location,'_type3.pdf',sep="")
		}
		plottypein(plottype,location)
		gridExtra::grid.arrange(p1, p2, ncol=2,nrow=1)

	}
	else if(plot==1){
		if(plottype=="pdf"){
			location=paste(location,'_type1.pdf',sep="")
		}
		plottypein(plottype,location)
		print(p2)

	}
	else if(plot==2){
		if(plottype=="pdf"){
			location=paste(location,'_type2.pdf',sep="")
		}
		plottypein(plottype,location)
		print(p2)
	}

}


#' @title Determine the distance in a heatmap
#' @param Data1 The resulting clustering of method 1.
#' @param Data2 The resulting clustering of method 2.
#' @param names The names of the objects in the data sets. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @description Internal function of \code{HeatmapPlot}
distanceheatmaps<-function(Data1,Data2,names=NULL,nrclusters=7){
	ClustData1=stats::cutree(Data1,nrclusters) #original clusters (aggregated data clustering)
	ClustData2=stats::cutree(Data2,nrclusters) #clusters of changed method

	ClustData1=ClustData1[Data1$order]
	ClustData2=ClustData2[Data2$order]

	trueorder1=sort(Data1$order,index.return = TRUE)
	trueorder2=sort(Data2$order,index.return = TRUE)

	ordercolors=ClustData1
	order=seq(1,nrclusters)

	for (k in 1:length(unique(ClustData1))){
		select=which(ClustData1==unique(ClustData1)[k])
		ordercolors[select]=order[k]
	}

	ordercolors2=ClustData2


	for (k in 1:length(unique(ClustData2))){
		select=which(ClustData2==unique(ClustData2)[k])
		ordercolors2[select]=order[k]
	}

	ClustData1=ordercolors[trueorder1$ix]
	ClustData2=ordercolors2[trueorder2$ix]

	out=matrix(0,length(ClustData1),length(ClustData2)) #want the rows to be the other method and the columns to be the aggregated data clustering

	#names=names[Data2$order]

	rownames(out)=names
	colnames(out)=names

	for(i in 1:length(names)){
		focus=names[i] #defines the column

		label=ClustData2[i] #color of the cluster is defined by the original cluster that contains focus (1 to 7)

		for(j in 1:length(names)){ #go over the rows
			other=names[j]
			found=FALSE  #find cluster of other
			k=1
			while(found==FALSE & k<=nrclusters){
				label2=k
				if(other %in% names[ClustData1==label2]){
					othercluster=names[ClustData1==label2]
					found=TRUE
				}
				k=k+1
			}

			if(focus %in% othercluster){ #if other and focus still together: give it color of cluster defined by focus
				out[j,i]=label
			}
		}
	}
	return(out)
}


#' @title A heatmap of the comparison of two clustering results.
#'
#' @description A function to plot a heatmap of the agreement between two clustering
#' results.
#'
#' These are the clusters to which a comparison is made. A matrix is set up of
#' which the columns are determined by the ordering of clustering of method 2
#' and the rows by the ordering of method 1. Every column represent one object
#' just as every row and every column represent the color of its cluster. A
#' function visits every cell of the matrix. If the objects represented by the
#' cell are still together in a cluster, the color of the column is passed to
#' the cell. This creates the distance matrix which can be given to the
#' HeatmapCols function to create the heatmap.
#'
#' @param Data1 The resulting clustering of method 1.
#' @param Data2 The resulting clustering of method 2.
#' @param names The names of the objects in the data sets. Default is NULL.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the colours for the clusters. Default is NULL.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A heatmap based on the distance matrix created by the function with
#' the dendrogram of method 2 on top of the plot and the one from method 1 on
#' the left. The names of the objects are depicted on the bottom in the order
#' of clustering of method 2 and on the right by the ordering of method 1.
#' Vertically the cluster of method 2 can be seen while horizontally those of
#' method 1 are portrayed.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' clust="agnes",linkage="flexible",gap=FALSE,maxK=15)
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' clust="agnes",linkage="flexible",gap=FALSE,maxK=15)
#'
#'
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#'
#' HeatmapPlot(Data1=MCF7_T,Data2=MCF7_F,names=rownames(fingerprintMat)
#' ,nrclusters=7,cols=Colors1,plottype="new", location=NULL)
#' }
#'
#' @export HeatmapPlot
HeatmapPlot<-function(Data1,Data2,names=NULL,nrclusters=NULL,cols=NULL,plottype="new",location=NULL){
	if (!requireNamespace("gplots", quietly=TRUE)) stop("HeatmapPlot() requires the suggested package 'gplots'.")
	data1=Data1$Clust
	data2=Data2$Clust
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
	DistM=distanceheatmaps(data1,data2,names,nrclusters)
	plottypein(plottype,location)
	gplots::heatmap.2(DistM,Rowv =stats::as.dendrogram(data1), Colv=stats::as.dendrogram(data2),trace="none",col=cols,key=FALSE)
	plottypeout(plottype)
}

#' @title A function to select a group of objects via the similarity heatmap.
#'
#' @description The function \code{HeatmapSelection} plots the similarity values between
#' objects. The plot is similar to the one produced by
#' \code{SimilarityHeatmap} but without the dendrograms on the sides. The
#' function is rather explorative and experimental and is to be used with some
#' caution. By clicking in the plot, the user can select a group of objects
#' of interest. See more in \code{Details}.
#'
#' A similarity heatmap is created in the same way as in
#' \code{SimilarityHeatmap}. The user is now free to select two points on the
#' heatmap. It is advised that these two points are in opposite corners of a
#' square that indicates a high similarity among the objects. The points do
#' not have to be the exact corners of the group of interest, a little
#' deviation is allowed as rows and columns of the selected subset of the
#' matrix with sum equal to 1 are filtered out. A sum equal to one, implies
#' that the compound is only similar to itself.
#'
#' The function is meant to be explorative but is experimental. The goal was to
#' make the selection of interesting objects easier as sometimes the labels
#' of the dendrograms are too distorted to be read. If the figure is exported
#' to a pdf file with an appropriate width and height, the labels can be become
#' readable again.
#'
#' @param Data The data of which a heatmap should be drawn.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" or "clusters".
#' @param distmeasure The distance measure. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param cutoff Optional. If a cutoff value is specified, all values lower are
#' put to zero while all other values are kept. This helps to highlight the
#' most similar objects. Default is NULL.
#' @param percentile Logical. The cutoff value can be a percentile. If one want
#' the cutoff value to be the 90th percentile of the data, one should specify
#' cutoff = 0.90 and percentile = TRUE. Default is FALSE.
#' @param dendrogram Optional. If the clustering results of the data is already
#' available and should not be recalculated, this results can be provided here.
#' Otherwise, it will be calculated given the data. This is necessary to have
#' the objects in their order of clustering on the plot. Default is NULL.
#' @param width The width of the plot to be made. This can be adjusted since
#' the default size might not show a clear picture. Default is 7.
#' @param height The height of the plot to be made. This can be adjusted since
#' the default size might not show a clear picture. Default is 7.
#' @return A heatmap with the names of the objects on the right and bottom.
#' Once points are selected, it will return the names of the objects that are
#' in the selected square provided that these show similarity among each other.
#' @examples
#'
#' \dontrun{
#' data(fingerprintMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55)
#'
#' HeatmapSelection(Data=MCF7_F$DistM,type="dist",cutoff=0.90,percentile=TRUE,
#' dendrogram=MCF7_F,width=7,height=7)
#' }
#'
#' @export HeatmapSelection
HeatmapSelection<-function(Data,type=c("data","dist","clust","sim"),distmeasure="tanimoto",normalize=FALSE,method=NULL,linkage="flexible",cutoff=NULL,percentile=FALSE,dendrogram=NULL,width=7,height=7){

	#create binary similarity heatmap first
	if(type=="data"){
		ClustData<-Cluster(Data=Data,distmeasure=distmeasure,normalize=normalize,method=method,clust="agnes",linkage=linkage,gap=FALSE,maxK=55,StopRange=FALSE)
		Data=ClustData$DistM
		type="dist"
	}


	if(type=="clust"){
		Dist=Data$DistM
		if(0<=min(Dist) & max(Dist)<=1){
			SimData=1-Dist
		}
		else{
			NormData=Normalization(Dist,method="Range")
			SimData=1-NormData
		}
		if(is.null(dendrogram)){
			dendrogram=Data
		}
	}

	else if(type=="dist"){
		if(0<=min(Data) & max(Data)<=1){
			SimData=1-Data
			if(is.null(dendrogram)){
				dendrogram=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		else{
			NormData=Normalization(Data,method="Range")
			SimData=1-NormData
			if(is.null(dendrogram)){
				dendrogram=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=TRUE,method="Q",clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}


	}
	else if(type=="sim"){
		SimData=Data
		if(0<=min(SimData) & max(SimData)<=1){
			if(is.null(dendrogram)){
				DistData=1-Data
				ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
		else{
			if(is.null(dendrogram)){
				NormData=Normalization(Dist,method="Range")
				DistData=1-Data
				ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
			}
		}
	}


	if(!is.null(cutoff)){
		if(percentile==TRUE){
			cutoff=stats::quantile(SimData[lower.tri(SimData)], cutoff)
		}

		SimData_bin <- ifelse(SimData<=cutoff,0,SimData) # Every value higher than the 90ieth percentile is kept, all other are put to zero
	}

	else{
		SimData_bin=SimData
	}



	dend <- stats::as.dendrogram(dendrogram$Clust)
	Ind <- stats::order.dendrogram(dend)

	SimData_bin=SimData_bin[Ind,Ind]

	#Layout<-rbind(4:3, 2:1)
	#lhei <- c(0.4, 4)
	#lwid <- c(0.4, 4)
	#layout(Layout, widths = lwid, heights = lhei, respect = FALSE)
	grDevices::dev.new(width=width,height=height)
	graphics::par(mar = c(9,7, 7, 9))
	graphics::image(x=1:nrow(SimData_bin),y=1:ncol(SimData_bin),z=t(SimData_bin),col=(grDevices::gray(seq(0.9,0,len=1000))),axes=FALSE,xlab="",ylab="")
	graphics::axis(1, 1:ncol(SimData_bin), labels = colnames(SimData_bin), las = 2, line =0, tick = 0, cex.axis = 0.6)
	graphics::axis(4, 1:nrow(SimData_bin), labels = rownames(SimData_bin), las = 2, line = 0, tick = 0, cex.axis = 0.6)

	points=graphics::locator(n=2,type="l")
	cols=c(floor(points$x[1]),ceiling(points$x[2]))
	rows=c(floor(points$y[1]),ceiling(points$y[2]))

	if(cols[1]>cols[2]){
		colseq=seq(cols[2],cols[1],1)
	}
	else{
		colseq=seq(cols[1],cols[2],1)
	}

	if(rows[1]>rows[2]){
		rowseq=seq(rows[2],rows[1],1)
	}
	else{
		rowseq=seq(rows[1],rows[2],1)
	}

	print(rowseq)
	print(colseq)
	SubsetData=SimData_bin[rowseq,colseq]
#	DelRows=rownames(SubsetData)[which(rowSums(SubsetData)==1)]
#	DelCols=colnames(SubsetData)[which(colSums(SubsetData)==1)]
#
#	if(length(DelRows)!=0 & length(DelCols)!=0){
#		Subset=SubsetData[-which(rownames(SubsetData)%in%c(DelRows,DelCols)),-which(colnames(SubsetData)%in%c(DelRows,DelCols))]
#	}
#	else if(length(DelRows)!=0 & length(DelCols)==0){
#		Subset=SubsetData[-which(rownames(SubsetData)%in%c(DelRows)),]
#
#	}
#	else if(length(DelRows)==0 & length(DelCols)!=0){
#		Subset=SubsetData[,-which(colnames(SubsetData)%in%c(DelCols))]
#
#	}
#	else if(length(DelRows)==0 & length(DelCols)==0){
#		Subset=SubsetData
#
#	}
#	SelComps=colnames(Subset)
	SelComps=colnames(SubsetData)

	return(SelComps)
}

#' @title A heatmap of similarity values between objects
#'
#' @description The function \code{SimilarityHeatmap} plots the similarity values between
#' objects. The darker the shade, the more similar objects are. The option
#' is available to set a cutoff value to highlight the most similar objects.
#'
#' @details If data is of type "clust", the distance matrix is extracted from the result
#' and transformed to a similarity matrix. Possibly a range normalization is
#' performed. If data is of type "dist", it is also transformed to a similarity
#' matrix and cluster is performed on the distances. If data is of type "sim",
#' the data is tranformed to a distance matrix on which clustering is
#' performed. Once the similarity mattrix is obtained, the cutoff value is
#' applied and a heatmap is drawn. If no cutoff value is desired, one can leave
#' the default NULL specification.
#'
#' @param Data The data of which a heatmap should be drawn.
#' @param type indicates whether the provided matrices in "List" are either data matrices, distance
#' matrices or clustering results obtained from the data. If type="dist" the calculation of the distance
#' matrices is skipped and if type="clusters" the single source clustering is skipped.
#' Type should be one of "data", "dist" ,"sim" or "clusters".
#' @param distmeasure The distance measure. Should be one of "tanimoto", "euclidean", "jaccard", "hamming". Defaults to "tanimoto".
#' @param normalize	Logical. Indicates whether to normalize the distance matrices or not, defaults to c(FALSE, FALSE) for two data sets. This is recommended if different distance types are used. More details on normalization in \code{Normalization}.
#' @param method A method of normalization. Should be one of "Quantile","Fisher-Yates", "standardize","Range" or any of the first letters of these names. Default is NULL.
#' @param linkage Choice of inter group dissimilarity (character). Defaults to "flexible".
#' @param cutoff Optional. If a cutoff value is specified, all values lower are
#' put to zero while all other values are kept. This helps to highlight the
#' most similar objects. Default is NULL.
#' @param percentile Logical. The cutoff value can be a percentile. If one want
#' the cutoff value to be the 90th percentile of the data, one should specify
#' cutoff = 0.90 and percentile = TRUE. Default is FALSE.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A heatmap with the names of the objects on the right and bottom
#' and a dendrogram of the clustering at the left and top.
#' @examples
#'
#' \dontrun{
#' data(fingerprintMat)
#'
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55)
#'
#' SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.90,percentile=TRUE)
#' SimilarityHeatmap(Data=MCF7_F,type="clust",cutoff=0.75,percentile=FALSE)
#'
#' }
#'
#' @export SimilarityHeatmap
SimilarityHeatmap<-function(Data,type=c("data","clust","sim","dist"),distmeasure="tanimoto",normalize=FALSE,method=NULL,linkage="flexible",cutoff=NULL,percentile=FALSE,plottype="new",location=NULL){
	if (!requireNamespace("gplots", quietly=TRUE)) stop("SimilarityHeatmap() requires the suggested package 'gplots'.")
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

	if(type=="data"){
		ClustData<-Cluster(Data=Data,distmeasure=distmeasure,normalize=normalize,method=method,clust="agnes",linkage=linkage,gap=FALSE,maxK=55,StopRange=FALSE)
		Data=ClustData$DistM
		type="dist"
	}


	if(type=="clust"){
		Dist=Data$DistM
		if(0<=min(Dist) & max(Dist)<=1){
			SimData=1-Dist
		}
		else{
			NormData=Normalization(Dist,method="Range")
			SimData=1-NormData
		}

		ClustData=Data
	}

	else if(type=="dist"){
		if(0<=min(Data) & max(Data)<=1){
			SimData=1-Data
			ClustData=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

		}
		else{
			NormData=Normalization(Data,method="Range")
			SimData=1-NormData
			ClustData=Cluster(Data=Data,type="dist",distmeasure="tanimoto",normalize=TRUE,method="Q",clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

		}


	}
	else if(type=="sim"){
		SimData=Data
		if(0<=min(SimData) & max(SimData)<=1){
			DistData=1-Data
			ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

		}
		else{
			NormData=Normalization(Dist,method="Range")
			DistData=1-Data
			ClustData=Cluster(Data=DistData,type="dist",distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
		}
	}


	if(!is.null(cutoff)){
		if(percentile==TRUE){
			cutoff=stats::quantile(SimData[lower.tri(SimData)], cutoff)
		}

		SimData_bin <- ifelse(SimData<=cutoff,0,SimData) # Every value higher than the 90ieth percentile is kept, all other are put to zero
	}

	else{
		SimData_bin=SimData
	}

	plottypein(plottype,location)
	gplots::heatmap.2(SimData_bin,
			Rowv = stats::as.dendrogram(stats::as.hclust(ClustData$Clust)), Colv=stats::as.dendrogram(stats::as.hclust(ClustData$Clust)),trace="none",
			col=(grDevices::gray(seq(0.9,0,len=1000))),
			cexRow=0.6, cexCol=0.6,
			margins=c(9,9),
			key=FALSE,
			keysize=0.4,
			symkey=FALSE,
			sepwidth=c(0.01,0.01),
			sepcolor="black",
			colsep=c(0,ncol(SimData_bin)),
			rowsep=c(0,nrow(SimData_bin))
	)
	plottypeout(plottype)

}
