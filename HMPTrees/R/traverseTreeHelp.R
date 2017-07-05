traverseTreeHelp <-
function(data, place, treeLvl, maxTaxaDepth, split){
	REMOVELBL <- TRUE  	# set to false to show all labels, otherwise only non-zero labels will be shown	
	childStr <- ""				
	myName <- rownames(data)[place]		
	myVal <- data[place, 1]
	stopTree <- FALSE
	noSiblings <- TRUE
	
	### Make sure we haven't reached our max depth yet
	if(treeLvl >= maxTaxaDepth)
		stopTree <- TRUE
	
	### Find all the children from our starting point
	data <- data[-place,, drop=FALSE] 
	childRows <- grep(paste(myName, split, sep=""), rownames(data), fixed=TRUE)
	
	### Go through all child branches
	if(length(childRows) != 0 && !stopTree){ # we know we have children
		newData <- data[childRows,, drop=FALSE]
		
		### Pull apart child branch names
		rownames(newData) <- substring(rownames(newData), nchar(myName)+2, nchar(rownames(newData)))
		splitLength <- unlist(lapply(strsplit(rownames(newData), split, fixed=TRUE), length))
		startPlace <- which(splitLength == 1)
		
		### Tranverse any child branches
		for(t in startPlace){ 	
			temp_hstr <- traverseTreeHelp(newData, t, treeLvl+1, maxTaxaDepth, split)
			if(childStr != ""){
				childStr <- paste(childStr, ",", temp_hstr, sep="")
			}else{
				childStr <- temp_hstr
			}
			
			noSiblings <- FALSE
		}
		
		### Add child branches to our newick format to return
		if(myVal == 0){ # set up string to be returned
			retStr <- paste("(", childStr, "):", myVal,  sep="")
		}else{
			retStr <- paste("(", childStr, ")", myName, ":", myVal,  sep="")
		}
		
		if(noSiblings) # adds an extra 0 value branch if there are no siblings
			retStr <- paste(retStr, ",:0.0", sep="")
	}else{ # no children found
		if(myVal == 0 && REMOVELBL == TRUE){
			retStr <- paste(":", myVal, sep="")
		}else{
			retStr <- paste(myName, ":", myVal, sep="")
		}
		
		retStr <- paste(retStr, ",:0.0", sep="")
	}
	
	return(retStr)
}
