transformHMPTreetoHMP <-
function(data, elimZero=FALSE, zeroValue=.00001){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Find 0 taxa
	dataSum <- rowSums(data)
	zeroRows <- which(dataSum==0)
	
	### Increase the first sample by zero value so taxa aren't all 0
	data[zeroRows, 1] <- zeroValue
	
	if(elimZero)
		data <- data[-zeroRows,]
	
	return(t(data))
}
