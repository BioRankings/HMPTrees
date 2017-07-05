buildTree <-
function(data, split="."){
	### Make a copy to attach our new levels too
	retData <- data
	
	### Go through every taxa
	for(i in 1:nrow(data)){ 
		fullNameSplit <- strsplit(rownames(data)[i], split, fixed=TRUE)[[1]]
		
		### Skip top level taxa
		if(length(fullNameSplit) == 1)
			next
		
		### Build a full branch from the name
		for(j in 1:(length(fullNameSplit)-1)){
			name <- paste(fullNameSplit[1:j], collapse=split)
			
			if(name %in% rownames(retData)){
				loc <- which(rownames(retData) %in% name)
				retData[loc,] <- retData[loc,] + data[i,]
			}else{
				retData <- rbind(temp=unlist(data[i,]), retData)
				rownames(retData)[1] <- name
			}
		}
	}
	
	### Reorder the taxa names
	retData <- retData[order(rownames(retData)),]
	
	return(retData)
}
