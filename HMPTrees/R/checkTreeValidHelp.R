checkTreeValidHelp <-
function(data, place, epsilon, split){
	myName <- rownames(data)[place]	
	myVal <- data[place, 1]
	
	### Pull apart the starting point name
	data <- data[-place,, drop=FALSE] 
	childRows <- grep(paste(myName, split, sep=""), rownames(data), fixed=TRUE)
	
	### Go through all child branches
	if(length(childRows) != 0){ # we know we have children
		newData <- data[childRows,, drop=FALSE]
		
		### Pull apart child names
		rownames(newData) <- substring(rownames(newData), nchar(myName)+2, nchar(rownames(newData)))	
		splitLength <- unlist(lapply(strsplit(rownames(newData), split, fixed=TRUE), length))
		startPlace <- which(splitLength == 1)
		
		### Tranverse all child branches
		childReadTot <- 0
		for(t in startPlace){ 	
			childReadTot <- childReadTot + newData[t, 1]
			if((myVal+epsilon) < childReadTot) # the child is greater so bad tree
				return(FALSE)
			
			valid <- checkTreeValidHelp(newData, t, epsilon, split)
			if(!valid)
				return(FALSE)
		}
	}
	
	return(myVal >= 0) # check bottom nodes arent negative
}
