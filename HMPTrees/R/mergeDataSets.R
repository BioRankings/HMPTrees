mergeDataSets <-
function(dataList, calcMLE=FALSE, uniqueNames=FALSE, data=NULL){
	### Fix any old arguments
	if(!is.null(data)){
		warning("'data' is deprecated. It has been replaced with dataList. View the help files for details.")
		dataList <- data
	}
	
	if(missing(dataList))
		stop("dataList is missing.")
	
	if(length(dataList) < 2)
		stop("At least 2 data sets are needed to merge")
	
	### Merge all the data sets into 1
	mles <- vector("list", length(dataList))
	newData <- NULL
	for(i in 1:length(dataList)){
		if(uniqueNames)
			colnames(dataList[[i]]) <- paste("Data", i, "-", colnames(dataList[[i]]))
		newData <- merge(newData, dataList[[i]], by=0, all=TRUE)
		rownames(newData) <- newData[,1]
		newData <- newData[,-1]
		
		if(calcMLE)
			mles[[i]] <- getMLEandLoglike(dataList[[i]])$mleTree
	}
	newData[is.na(newData)] <- 0 #set any NA's to 0
	
	### Merge all the MLEs into the data sets
	if(calcMLE){
		mleAll <- getMLEandLoglike(newData)$mleTree
		colnames(mleAll) <- "Combined MLE"
		
		for(i in 1:length(mles)){
			newData <- merge(newData, mles[[i]], by=0, all=TRUE)
			rownames(newData) <- newData[,1]
			newData <- newData[,-1]
			
			colnames(newData)[ncol(newData)] <- paste("data", i, "mle", sep="")
		}	
		
		newData <- cbind(newData, mleAll)
		newData[is.na(newData)] <- 0 #set any NA's to 0
	}
	
	return(newData)
}
