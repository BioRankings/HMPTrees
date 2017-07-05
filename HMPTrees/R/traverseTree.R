traverseTree <-
function(data, level, split){	
	maxTaxaDepth <- getTaxaDepth(level)
	
	### Pull apart our starting name
	splitLength <- unlist(lapply(strsplit(rownames(data), split, fixed=TRUE), length))
	startPlace <- which(splitLength == 1)
	
	### Go through every starting branch
	myStr <- ""
	for(i in startPlace){ 	
		tempStr <- traverseTreeHelp(data[, 1, drop=FALSE], i, 1, maxTaxaDepth, split)
		if(tempStr == "") 
			next 
		
		if(myStr != ""){
			myStr <- paste(myStr, ",", tempStr, sep="")
		}else{
			myStr <- tempStr
		}
	}
	### Add the final newick touches to the string
	myTree <- paste("", myStr, ";", sep="")
	
	### Turn the newick format into a 'phylo' tree
	retTree <- ape::read.tree(text=myTree) 
	retTree <- ape::collapse.singles(retTree)
	
	return(retTree)
}
