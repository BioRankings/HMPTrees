checkTreeValidity <-
function(data, samples=NULL, epsilon=0.0001, split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Check the validity of each tree
	overAllValid <- rep(FALSE, ncol(data))
	for(i in 1:ncol(data)){
		tempData <- data[,i, drop=FALSE]
		nameSplit <- strsplit(rownames(tempData), split, fixed=TRUE)
		
		### Find the starting point of the tree
		startingPoint <- which(length(nameSplit[[i]]) == 1)
		if(length(startingPoint) > 1){
			next
		}else{
			overAllValid[i] <- checkTreeValidHelp(tempData, startingPoint, epsilon, split)
		}
	}
	
	return(overAllValid)
}
