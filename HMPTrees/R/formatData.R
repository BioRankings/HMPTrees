formatData <-
function(data, countThreshold=1000, normalizeThreshold=10000){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Order the data and turn any NAs to 0
	data <- data[order(rownames(data)),, drop=FALSE]
	data[is.na(data)] <- 0
	
	### Keep only samples where the top level is above the count Threshold
	data <- data[, data[1,] >= countThreshold, drop=FALSE]
	
	### Make sure we havent removed everything
	if(ncol(data) == 0)
		stop("'countThreshold' is too high.")
	
	# Normalize the read counts
	if(normalizeThreshold > 0){
		for(i in ncol(data):1)
			data[,i] <- data[,i] * (normalizeThreshold/data[1, i])
	}
	
	return(data)
}
