removeUnclass <-
function(data, remove=TRUE){
	falsePos <- grep(".U", rownames(data), fixed=TRUE)
	unclass <- grep("U", rownames(data), fixed=TRUE)
	
	if(length(unclass) == 0) # Nothing to remove
		return(data)
	
	if(length(falsePos) != 0) # We have false positives
		unclass <- setdiff(unclass, falsePos)
	
	if(remove){
		data <- data[-unclass,,drop=FALSE]
	}else{
		data[unclass,] <- 0
	}
	
	return(data)		
}
