createTrees <-
function(data, samples=NULL, level="genus", split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	if(any(grepl(")", rownames(data), fixed=TRUE)) || any(grepl("(", rownames(data), fixed=TRUE)) || any(grepl(":", rownames(data), fixed=TRUE)))
		stop("Using parentheses and/or colons in the taxa names is not allowed.")
	
	### Sort the data based on taxa names
	data <- data[order(rownames(data)),, drop=FALSE]
	
	### Create a newick format for each tree
	allTrees <- vector("list", length(samples))
	for(i in 1:ncol(data)){
		oneSamp <- data[,i, drop=FALSE]
		if(sum(oneSamp) <= 0) # skips entries without data
			next
		
		allTrees[[i]] <- traverseTree(oneSamp, level, split)
	}
	names(allTrees) <- colnames(data)
	
	return(allTrees)
}
