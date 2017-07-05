plotTreeDataMDS <-
function(dataList, main="Tree MDS Comparisons", calcMLE=TRUE, mleTitles=NULL, dotColors=NULL, 
		dotSizes=NULL, showNames=FALSE, returnCoords=FALSE, data=NULL){
	### Fix any old arguments
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList))
		stop("At least 1 valid data set is required.")
	
	if(is.null(dotColors))
		dotColors <- c("red", "orange", "blue", "green", "yellow", "purple")
	
	if(is.null(dotSizes))
		dotSizes <- c(1, 2)
	if(length(dotSizes) != 2)
		stop("dotSizes must contain 2 integers.")
	
	numGrps <- length(dataList)
	
	### Add group names if we dont have them
	if(is.null(names(dataList))){
		grpNames <- paste("Data Set", 1:numGrps)
	}else{
		grpNames <- names(dataList)
	}
	
	### If we don't have enough colors, make a new set of colors
	if(length(dotColors) < numGrps)
		dotColors <- rainbow(numGrps)
	
	sizes <- NULL
	colors <- NULL
	if(calcMLE){
		mleDotSizes <- rep(0, numGrps)
		mleColors <- rep(0, numGrps)
	}
	### Figure out colors and sizes for all points
	for(i in 1:numGrps){
		tempData <- dataList[[i]]
		colors <- c(colors, rep(dotColors[i], ncol(tempData)))
		sizes <- c(sizes, rep(dotSizes[1], ncol(tempData)))
		
		if(calcMLE){
			mleColors[i] <- dotColors[i]
			mleDotSizes[i] <- dotSizes[2]	
		}
	}
	if(calcMLE){
		colors <- c(colors, mleColors, "black")
		sizes <- c(sizes, mleDotSizes, dotSizes[2])
	}
	
	### Get point positions and plot
	tData <- t(mergeDataSets(dataList, calcMLE))
	loc <- cmdscale(dist(tData), k=2)
	x <- loc[,1]
	y <- -loc[,2]
	plot(x, y, pch=19, xlab="MDS 1", ylab="MDS 2", pty="s", col=colors, cex=sizes, main=main)
	
	### Add mle titles if we aren't showing all the names
	if(calcMLE && !showNames){
		grpNames <- c(grpNames, "Combined MLE")
		text(x[(nrow(tData) - numGrps):nrow(tData)], y[(nrow(tData) - numGrps):nrow(tData)], grpNames, pos=3, cex=.75)
	}
	
	### Add the names of the samples
	if(showNames) 
		text(x, y, rownames(tData), pos=3, cex=.75)
	
	if(returnCoords)
		return(list(x=x, y=y))
}
