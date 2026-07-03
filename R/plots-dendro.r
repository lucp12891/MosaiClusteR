#' @title Comparison of clustering results over multiple methods
#'
#' @description A visual comparison of the clustering results of several methods.
#' The function relies on \code{ReorderToReference} and \code{ColorsNames} and
#' renders the resulting matrix either as a rectangular heatmap (using the
#' \code{plotrix} package) or in a circular format (using the \code{circlize}
#' package).
#' @param List A list of the outputs from the methods to be compared. The first
#' element of the list is used as the reference in \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the colours to be used. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}. Default is FALSE.
#' @param weightclust Logical. To be handed to \code{ReorderToReference}. Default is FALSE.
#' @param names Optional. Names of the methods to be used as labels. Default is NULL.
#' @param margins Optional. Margins for the plot. Default is c(8.1,3.1,3.1,4.1).
#' @param circle Logical. Whether the figure should be circular (TRUE) or a rectangle (FALSE). Default is FALSE.
#' @param canvaslims The limits for the circular dendrogram. Default is c(-1.0,1.0,-1.0,1.0).
#' @param Highlight Optional. A list of character vectors of objects to be highlighted. Default is NULL.
#' @param substr Optional. A vector of length two giving start and stop positions to shorten labels. Default is NULL.
#' @param cex.highlight Magnification for the highlight text. Default is 1.
#' @param cex.labels Magnification for the labels. Default is 0.7.
#' @param trackheightdend The height of the dendrogram track in the circular plot. Default is 0.5.
#' @param plottype Should be one of "pdf","new" or "sweave". Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided. Default is NULL.
#' @return A plot which translates the matrix output of \code{ReorderToReference}.
#' @export
ComparePlot<-function(List,nrclusters=NULL,cols=NULL,fusionsLog=FALSE,weightclust=FALSE,names=NULL,margins=c(8.1,3.1,3.1,4.1),circle=FALSE,canvaslims=c(-1.0,1.0,-1.0,1.0),Highlight=NULL,substr=NULL,cex.highlight=1,cex.labels=0.7,trackheightdend=0.5,plottype="new",location=NULL){
	if (!requireNamespace("circlize", quietly = TRUE)) stop("ComparePlot() requires the suggested package 'circlize'.")
	if (!requireNamespace("plotrix", quietly = TRUE)) stop("ComparePlot() requires the suggested package 'plotrix'.")
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

	for(i in 1:length(List)){
		if(attributes(List[[i]])$method == "Weighted" & weightclust==TRUE){
			T=List[[i]]$Clust
			attr(T,"method")="Single Clustering"
			List[[i]]=T
		}
	}

	MatrixColors=ReorderToReference(List,nrclusters,fusionsLog,weightclust,names)
	#capture singletons
	for(i in 1:nrow(MatrixColors)){
		if(any(table(MatrixColors[i,])==1)){
			clusters=names(table(MatrixColors[i,]))[which(table(MatrixColors[i,])==1)]
			MatrixColors[i,][which(MatrixColors[i,]%in%as.numeric(clusters))]=5000
		}
	}

	Names=ColorsNames(MatrixColors,cols)
	colnames(Names)=colnames(MatrixColors)
	nobs=dim(MatrixColors)[2]
	nmethods=dim(MatrixColors)[1]

	if(is.null(names)){
		for(j in 1:nmethods){
			names[j]=paste("Method",j,sep=" ")
		}
	}

	if(circle){
		plottypein(plottype, location)
		graphics::par(mar=margins)
		circlize::circos.initialize(factors =c(1:ncol(MatrixColors)) , xlim = c(0, ncol(MatrixColors)))
		for(i in 1:length(List)){
			circlize::circos.trackPlotRegion(factors = c(1:ncol(MatrixColors)), ylim = c(0,1),track.height = 0.05)
		}
		track=c(length(List):1)
		for(i in 1:nrow(MatrixColors)){
			for(j in 1:ncol(MatrixColors)){
				circlize::highlight.sector(sector.index=j, track.index = track[i],col=Names[i,j])
			}
		}
		hc=List[[1]]$Clust
		#labels_short=substr(hc$order.lab,1,7)
		labels=hc$order.lab

		if(!is.null(Highlight)){
			for(h in 1:length(Highlight)){
				Name=names(Highlight)[h]
				HL=which(labels%in%Highlight[[h]])

				Sims=c()
				for(i in 1:length(List)){
					Values=List[[i]]$DistM[lower.tri(List[[i]]$DistM)]
					Sims=c(Sims,as.numeric(1-List[[i]]$DistM[Highlight[[h]],Highlight[[h]]][lower.tri(List[[i]]$DistM[Highlight[[h]],Highlight[[h]]])]))
				}
				MedSim=round(median(Sims),2)

				circlize::draw.sector(circlize::get.cell.meta.data("cell.start.degree", sector.index = min(HL)),
						circlize::get.cell.meta.data("cell.end.degree", sector.index = max(HL)),
						rou1 = 1, col = "#00000020")

				circlize::highlight.sector(sector.index=c(min(HL):max(HL)), track.index = 1, text = paste(Name,": ",MedSim,sep=""),
						facing = "bending.inside", niceFacing = TRUE, text.vjust = -1.5,,cex=cex.highlight)

			}
		}
		circlize::circos.clear()
		graphics::par(new = TRUE)
		hc=List[[1]]$Clust
		max_height=max(hc$height)
		dend=as.dendrogram(hc)

		if(!is.null(substr)){
			labelsshort=substr(hc$order.lab,substr[1],substr[2])
		}
		else{
			labelsshort=hc$order.lab
		}

		#ct=cutree(dend,6)
		circlize::circos.par("canvas.xlim" = c(canvaslims[1], canvaslims[2]), "canvas.ylim" = c(canvaslims[3], canvaslims[4]))
		circlize::circos.initialize(factors =1 , xlim = c(0, ncol(MatrixColors)))
		circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = trackheightdend,bg.border=NA,
				panel.fun = function(x, y) {
					for(i in seq_len(ncol(MatrixColors))) {
						circlize::circos.text(i-0.5, 0, labelsshort[i], adj = c(0, 0.5),
								facing = "clockwise", niceFacing = TRUE,
								col = Names[1,colnames(MatrixColors)[i]],
								cex = cex.labels,font=2)
					}
				})
		circlize::circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,track.height = 0.4, panel.fun = function(x, y) {
					circlize::circos.dendrogram(dend, max_height = max_height)})
		circlize::circos.clear()




		plottypeout(plottype)


	}

	else{
		#similar=round(SimilarityMeasure(MatrixColors),2)
		plottypein(plottype,location)
		graphics::par(mar=margins)
		plotrix::color2D.matplot(MatrixColors,cellcolors=Names,show.values=FALSE,axes=FALSE,xlab="",ylab="")
		graphics::axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColors),las=2,cex.axis=1.5)
		graphics::axis(2,at=seq(0.5,(nmethods-0.5)),labels=rev(names),cex.axis=1.5,las=2)
		#axis(4,at=seq(0.5,(nmethods-0.5)),labels=rev(similar),cex.axis=0.65,las=2)
		plottypeout(plottype)
	}
}


#' @title Comparison of clustering results over multiple results in circular format
#'
#' @description A visual comparison of all methods is handy to see which objects will
#' always cluster together independent of the applied methods. To this aid the
#' function \code{Cyclogram} has been written. The function relies on methods
#' of the \code{circlize} package.
#' @details This function makes use of the functions \code{ReorderToReference} and
#' \code{Colorsnames}. Given a list with the outputs of several methods, the
#' first step is to call upon \code{ReorderToReference} and to produce a
#' matrix of which the columns are ordered according to the ordering of the
#' objects of the first method in the list. Each cell represent the number of
#' the cluster the object belongs to for a specific method indicated by the
#' rows. The clusters are arranged in such a way that these correspond to that
#' one cluster of the referenced method that they have the most in common with.
#' The `circlize` package is used to visualize the matrix is a circular format.
#' The inner element of the circle is the first element of the List parameter, the second inner
#' element is the second element of the list and so on. The object names are written on the
#' circular dendrogram portayed in the center of the circle.
#' @param List A list of the outputs from the methods to be compared. The first
#' element of the list will be used as the reference in
#' \code{ReorderToReference}.
#' @param nrclusters The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols A character vector with the colours to be used. Default is NULL.
#' @param fusionsLog Logical. To be handed to \code{ReorderToReference}: indicator for the fusion of clusters. Default is TRUE
#' @param weightclust Logical. To be handed to \code{ReorderToReference}: to be used for the outputs of CEC,
#' WeightedClust or WeightedSimClust. If TRUE, only the result of the Clust element is considered. Default is TRUE.
#' @param names Optional. Names of the methods to be used as labels for the
#' columns. Default is NULL.
#' @param canvaslims The limits for the circular dendrogam. Default is c(-1.0,1.0,-1.0,1.0).
#' @param margins Optional. Margins to be used for the plot. Default is c(8.1,3.1,3.1,4.1).
#' @param Highlight Optional. A list of character vectors of objects to be highlighted. Default is NULL.
#' @param cex.highlight Magnification for the highlight text. Default is 1.
#' @param substr Optional. A vector of length two giving start and stop positions to shorten labels. Default is NULL.
#' @param cex.labels Magnification for the labels. Default is 0.7.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document, i.e. no new device is
#' opened and the plot appears in the current device or document. Default is "new".
#' @param location Optional. If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @return A plot which translates the matrix output of the function
#' \code{ReorderToReference} into a circular format.
#' @export Cyclogram
Cyclogram<-function (List, nrclusters = NULL, cols = NULL, fusionsLog = FALSE,
		weightclust = FALSE, names = NULL, canvaslims=c(-1.0,1.0,-1.0,1.0), margins = c(8.1, 3.1,
				3.1, 4.1), Highlight=NULL,cex.highlight=1,substr=NULL,cex.labels=0.7, plottype = "new", location = NULL)
{
	if (!requireNamespace("circlize", quietly = TRUE)) stop("Cyclogram() requires the suggested package 'circlize'.")
	plottypein <- function(plottype, location) {
		if (plottype == "pdf" & !(is.null(location))) {
			grDevices::pdf(paste(location, ".pdf", sep = ""))
		}
		if (plottype == "new") {
			grDevices::dev.new()
		}
		if (plottype == "sweave") {
		}
	}
	plottypeout <- function(plottype) {
		if (plottype == "pdf") {
			grDevices::dev.off()
		}
	}
	for (i in 1:length(List)) {
		if (attributes(List[[i]])$method == "Weighted" & weightclust ==
				TRUE) {
			T = List[[i]]$Clust
			attr(T, "method") = "Single Clustering"
			List[[i]] = T
		}
	}
	MatrixColors = ReorderToReference(List, nrclusters, fusionsLog,
			weightclust, names)
	Names = ColorsNames(MatrixColors, cols)
	colnames(Names)=colnames(MatrixColors)
	nobs = dim(MatrixColors)[2]
	nmethods = dim(MatrixColors)[1]
	if (is.null(names)) {
		for (j in 1:nmethods) {
			names[j] = paste("Method", j, sep = " ")
		}
	}
	plottypein(plottype, location)

	circlize::circos.initialize(factors =c(1:ncol(MatrixColors)) , xlim = c(0, ncol(MatrixColors)))
	for(i in 1:length(List)){
		circlize::circos.trackPlotRegion(factors = c(1:ncol(MatrixColors)), ylim = c(0,1),track.height = 0.05)
	}
	track=c(length(List):1)
	for(i in 1:nrow(MatrixColors)){
		for(j in 1:ncol(MatrixColors)){
			circlize::highlight.sector(sector.index=j, track.index = track[i],col=Names[i,j])
		}
	}
	hc=List[[1]]$Clust
	#labels=substr(hc$order.lab,1,5)
	labels=hc$order.lab
	if(!is.null(Highlight)){
		for(h in 1:length(Highlight)){
			Name=names(Highlight)[h]
			#HL=which(labels%in%substr(Highlight[[h]],1,5))
			HL=which(labels%in%Highlight[[h]])

			Sims=c()
			for(i in 1:length(List)){
				Values=List[[i]]$DistM[lower.tri(List[[i]]$DistM)]
				Sims=c(Sims,as.numeric(1-List[[i]]$DistM[Highlight[[h]],Highlight[[h]]][lower.tri(List[[i]]$DistM[Highlight[[h]],Highlight[[h]]])]))
			}
			MedSim=round(stats::median(Sims),2)

			circlize::draw.sector(circlize::get.cell.meta.data("cell.start.degree", sector.index = min(HL)),
					circlize::get.cell.meta.data("cell.end.degree", sector.index = max(HL)),
					rou1 = 1, col = "#00000020")

			circlize::highlight.sector(sector.index=c(min(HL):max(HL)), track.index = 1, text = paste(Name,": ",MedSim,sep=""),
					facing = "bending.inside", niceFacing = TRUE, text.vjust = -1.5,cex=cex.highlight)

		}
	}
	circlize::circos.clear()
	graphics::par(new = TRUE)
	hc=List[[1]]$Clust
	max_height=max(hc$height)
	dend=stats::as.dendrogram(hc)
	if(!is.null(substr)){
		labels=substr(hc$order.lab,substr[1],substr[2])
	}
	else{
		labels=hc$order.lab
	}
	circlize::circos.par("canvas.xlim" = c(canvaslims[1], canvaslims[2]), "canvas.ylim" = c(canvaslims[3], canvaslims[4]))
	circlize::circos.initialize(factors =1 , xlim = c(0, ncol(MatrixColors)))
	circlize::circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.4,bg.border=NA,
			panel.fun = function(x, y) {
				for(i in seq_len(ncol(MatrixColors))) {
					circlize::circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5),
							facing = "clockwise", niceFacing = TRUE,
							col = Names[1,colnames(MatrixColors)[i]],
							cex = cex.labels,font=2)
				}
			})
	circlize::circos.trackPlotRegion(ylim = c(0, max_height), bg.border = NA,track.height = 0.4, panel.fun = function(x, y) {
				circlize::circos.dendrogram(dend, max_height = max_height)})
	circlize::circos.clear()




	plottypeout(plottype)
}



#' @title Function that annotates colors to their names
#'
#' @description The \code{ColorsNames} function is used on the output of the
#' \code{ReorderToReference} and matches the cluster numbers indicated by the
#' cell with the names of the colors.  This is necessary to produce the plot of
#' the \code{ComparePlot} function and is therefore an internal function of
#' this function but can also be applied separately.
#'
#' @param matrixColors The output of the ReorderToReference function.
#' @param cols A character vector with the names of the colours to be used. Default is NULL.
#' @return A matrix containing the hex code of the color that corresponds to
#' each cell of the matrix to be colored. This function is called upon by the
#' \code{ComparePlot} function.
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
#' L=list(MCF7_F,MCF7_T)
#' names=c("FP","TP")
#'
#' MatrixColors=ReorderToReference(List=L,nrclusters=7,fusionsLog=TRUE,weightclust=TRUE,
#' names=names)
#'
#' Names=ColorsNames(matrixColors=MatrixColors,cols=Colors1)
#' }
#' @export ColorsNames
ColorsNames<-function(matrixColors,cols=NULL){
	Names=matrix(0,nrow=dim(matrixColors)[1],ncol=dim(matrixColors)[2])
	for (i in 1:dim(matrixColors)[1]){
		for(j in 1:dim(matrixColors)[2]){
			if(is.na(matrixColors[i,j])){
				Names[i,j]="grey"

			}
			else if(length(unique(matrixColors[i,]))==1){
				Names[i,j]="grey"

			}
			else if(matrixColors[i,j]==5000){
				Names[i,j]="white"

			}
			else{
				temp=matrixColors[i,j]
				Color=cols[temp]
				Names[i,j]=Color
			}
		}
	}
	return(Names)
}


#' @title Helper that colours dendrogram leaves by cluster membership
#'
#' @description Internal helper applied via \code{dendrapply} that colours the
#' leaves and edges of a dendrogram according to the cluster a leaf belongs to,
#' or to a specific selection of objects. Used by \code{ClusterPlot}.
#' @param x A node of a dendrogram (passed by \code{dendrapply}).
#' @param Data The clustering result (an \code{agnes}/\code{hclust} object) used to determine cluster membership.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in. Default is NULL.
#' @param cols The colours for the clusters if nrclusters is specified. Default is NULL.
#' @param colorComps Optional. A character vector of objects to highlight in red. Default is NULL.
#' @return The dendrogram node \code{x} with coloured \code{nodePar} and \code{edgePar} attributes.
#' @export
ClusterCols <- function(x,Data,nrclusters=NULL,cols=NULL,colorComps=NULL) {

	if(is.null(nrclusters) & is.null(colorComps)){
		return(x)
	}
	else if(!is.null(nrclusters)){
		if(length(cols)<nrclusters){
			stop("Not for every cluster a color is specified")
		}
	}

	if(!is.null(nrclusters)){
		Clustdata=stats::cutree(Data,nrclusters)
		Clustdata=Clustdata[Data$order]

		ordercolors=Clustdata
		order=seq(1,nrclusters)

		for (k in 1:length(unique(Clustdata))){
			select=which(Clustdata==unique(Clustdata)[k])
			ordercolors[select]=cols[order[k]]
		}
		names(ordercolors)=Data$order.lab

	}
	else{
		cols=rep("black",length(Data$order.lab))
		names(cols)=Data$order.lab
		cols[which(names(cols)%in%colorComps)]="red"
		ordercolors=cols

	}

	colfunc=function(x,cols,colorComps){
#		if(is.null(colorComps)){
#			color=cols[which(names(cols)==x)]
#			indextemp=which(attr(Data$diss,"Labels")==x)
#			index1=which(Data$order==indextemp)
#
#			index2=ordercolors[index1]
#
#			color=cols[index2]
#
#		}
#		else{
#			color=cols[which(names(cols)==x)]
#		}
		color=cols[which(names(cols)==x)]
		return(color)
	}

	if (stats::is.leaf(x)) {
		## fetch label
		label <- attr(x, "label")
		## set label color to clustercolor
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,ordercolors,colorComps),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,ordercolors,colorComps))
	}
	return(x)
}


#' @title Colouring clusters in a dendrogram
#'
#' @description Plot a dendrogram with leaves colored by a result of choice.
#' @param Data1 The resulting clustering of a method which contains the
#' dendrogram to be colored.
#' @param Data2 Optional. The resulting clustering of another method , i.e. the
#' resulting clustering on which the colors should be based. Default is NULL.
#' @param nrclusters Optional. The number of clusters to cut the dendrogram in.
#' If not specified the dendrogram will be drawn without colours to discern the
#' different clusters. Default is NULL.
#' @param cols The colours for the clusters if nrclusters is specified. Default is NULL.
#' @param colorComps If only a specific set of objects needs to be
#' highlighted, this can be specified here. The objects should be given in a
#' character vector. If specified, all other compound labels will be colored
#' black. Default is NULL.
#' @param hangdend A specification for the length of the brances of the dendrogram. Default is 0.02.
#' @param plottype Should be one of "pdf","new" or "sweave". If "pdf", a
#' location should be provided in "location" and the figure is saved there. If
#' "new" a new graphic device is opened and if "sweave", the figure is made
#' compatible to appear in a sweave or knitr document. Default is "new".
#' @param location If plottype is "pdf", a location should be provided in
#' "location" and the figure is saved there. Default is NULL.
#' @param \dots Other options which can be given to the plot function.
#' @return A plot of the dendrogram of the first clustering result with colored
#' leaves. If a second clustering result is given in Data2, the colors are
#' based on this clustering result.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' data(targetMat)
#' data(Colors1)
#'
#' MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' ClusterPlot(MCF7_T ,nrclusters=7,cols=Colors1,plottype="new",location=NULL,
#' main="Clustering on Target Predictions: Dendrogram",ylim=c(-0.1,1.8))
#' }
#' @export ClusterPlot
ClusterPlot<-function(Data1,Data2=NULL,nrclusters=NULL,cols=NULL,colorComps=NULL,hangdend=0.02,plottype="new",location=NULL,...){
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

	cx=Data1$Clust
	if(is.null(Data2)){
		Data=Data1$Clust
	}
	else{
		Data=Data2$Clust
	}

	d_temp<- stats::dendrapply(stats::as.dendrogram(stats::as.hclust(cx),hang=hangdend),ClusterCols,Data,nrclusters,cols,colorComps)
	plottypein(plottype,location)
	graphics::plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2)
	graphics::axis(side = 2, lwd = 2)
	plottypeout(plottype)
}


#' @title Create a color palette to be used in the plots
#'
#' @description In order to facilitate the visualization of the influence of the different
#' methods on the clustering of the objects, colours can be used. The function
#' \code{ColorPalette} is able to pick out as many colours as there are
#' clusters. This is done with the help of the \code{ColorRampPalette} function
#' of the grDevices package
#' @param colors A vector containing the colors of choice
#' @param ncols The number of colors to be specified. If higher than the number
#' of colors, it specifies colors in the region between the given colors.
#' @return A vector containing the hex codes of the chosen colors.
#' @examples
#'
#' Colors1<-ColorPalette(c("cadetblue2","chocolate","firebrick2",
#' "darkgoldenrod2", "darkgreen","blue2","darkorchid3","deeppink2"), ncols=8)
#' @export ColorPalette
ColorPalette<-function(colors=c("red","green"),ncols=5){
	my_palette=grDevices::colorRampPalette(colors)(ncols)

	return(my_palette)

}


#' @title Helper that colours specific dendrogram leaves
#'
#' @description Internal helper applied via \code{dendrapply} that colours the
#' leaves and edges of a dendrogram for one or two selections of objects. Used
#' by \code{LabelPlot}.
#' @param x A node of a dendrogram (passed by \code{dendrapply}).
#' @param Sel1 The selection of objects to be colored.
#' @param Sel2 An optional second selection to be colored. Default is NULL.
#' @param col1 The color for the first selection. Default is NULL.
#' @param col2 The color for the optional second selection. Default is NULL.
#' @return The dendrogram node \code{x} with coloured \code{nodePar} and \code{edgePar} attributes.
#' @export
LabelCols <- function(x,Sel1,Sel2=NULL,col1=NULL,col2=NULL) {
	colfunc=function(x,Sel1,Sel2,col1,col2){
		if (x %in% Sel1){
			return(col1)
		}
		else if(x %in% Sel2){
			return(col2)
		}
		else{
			return("black")
		}
	}

	if (stats::is.leaf(x)) {
		## fetch label
		label <- attr(x, "label")
		## set label color to red for SelF, to black otherwise
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,Sel1,Sel2,col1,col2),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,Sel1,Sel2,col1,col2))
	}
	return(x)
}

#' @title Coloring specific leaves of a dendrogram
#'
#' @description The function plots a dendrogrmam of which specific leaves are coloured.
#'
#' @param Data The result of a method which contains the dendrogram to be
#' colored.
#' @param sel1 The selection of objects to be colored. Default is NULL.
#' @param sel2 An optional second selection to be colored. Default is NULL.
#' @param col1 The color for the first selection. Default is NULL.
#' @param col2 The color for the optional second selection. Default is NULL.
#' @return A plot of the dendrogram of which the leaves of the selection(s) are
#' colored.
#' @examples
#' \dontrun{
#' data(fingerprintMat)
#' MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
#' method=NULL,clust="agnes",linkage="flexible",gap=FALSE,maxK=55,StopRange=FALSE)
#'
#' ClustF_6=cutree(MCF7_F$Clust,6)
#'
#' SelF=rownames(fingerprintMat)[ClustF_6==6]
#' SelF
#'
#' LabelPlot(Data=MCF7_F,sel1=SelF,sel2=NULL,col1='darkorchid')
#' }
#'
#' @export LabelPlot
LabelPlot<-function(Data,sel1,sel2=NULL,col1=NULL,col2=NULL){
	x=Data$Clust

	d_temp <- stats::dendrapply(stats::as.dendrogram(x,hang=0.02),LabelCols,sel1,sel2,col1,col2)

	graphics::plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2)
	graphics::axis(side = 2, lwd = 2)
}
