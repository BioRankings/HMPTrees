
library(ape)			# plotting and creating tree objects
library(HMP)			# creating new trees
library(dirmult)		# creating new trees
library(doParallel)		# parallelizing pvalue calc


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### External
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### pvalue functions
### ~~~~~~~~~~~~~~~~~~~~~
compareTwoDataSets <- function(data1, data2, numPerms=1000, parallel=FALSE, cores=3, maxSteps=50, delta=10^(-6)){
	if(missing(data1) || missing(data2))
		stop("Two valid data sets are required.")
	
	if(numPerms <= 0)
		stop("The number of boostraps must be an integer greater than 0.")
	
	### get the subject numbers
	numSub1 <- ncol(data1)
	numSub2 <- ncol(data2)
	numSubC <- numSub1 + numSub2
	
	### Merge data1 and data2 together
	if(any(rownames(data1) != rownames(data2))){
		dataComb <- merge(data1, data2, by=0, all=TRUE)
		rownames(dataComb) <- dataComb[,1]
		dataComb <- dataComb[,-1]
		dataComb[is.na(dataComb)] <- 0
		
		### Pull data1 and data2 back out if we had to merge
		data1 <- dataComb[,1:numSub1, drop=FALSE]
		data2 <- dataComb[,-c(1:numSub1), drop=FALSE]
	}else{
		dataComb <- cbind(data1, data2)
	}
	
	### Get our starting Loglik
	ll1 <- getMLEandLoglike(data1, maxSteps, delta=delta)$Loglik
	ll2 <- getMLEandLoglike(data2, maxSteps, delta=delta)$Loglik
	llC <- getMLEandLoglike(dataComb, maxSteps, delta=delta)$Loglik
	LRTobs <- -2*(llC-ll1-ll2)
	if(LRTobs == 0) # Exactly the same so no need to test farther
		return(1)
	
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numPerms)) 
		doParallel::registerDoParallel(cl)
		tryCatch({ 
					LRTbootstrap <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .export=c("getMLEandLoglike")) %dopar%{
						### Get a random sampling of the data and break into 2 groups
						samps <- sample(numSubC, numSub1, replace=FALSE)
						temp1 <- dataComb[,samps, drop=FALSE]
						temp2 <- dataComb[,-samps, drop=FALSE]
						
						### Get logliks for the sampled data
						tempLL1 <- getMLEandLoglike(temp1, maxSteps, delta=delta)$Loglik
						tempLL2 <- getMLEandLoglike(temp2, maxSteps, delta=delta)$Loglik
						LRT <- -2*(llC-tempLL1-tempLL2)
						
						return(LRT)
					}
				}, finally = {
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{ 	
		LRTbootstrap <- rep(0, numPerms)
		for(i in 1:numPerms){
			### Get a random sampling of the data and break into 2 groups
			samps <- sample(numSubC, numSub1, replace=FALSE)
			temp1 <- dataComb[,samps, drop=FALSE]
			temp2 <- dataComb[,-samps, drop=FALSE]
			
			### Get logliks for the sampled data
			tempLL1 <- getMLEandLoglike(temp1, maxSteps, delta=delta)$Loglik
			tempLL2 <- getMLEandLoglike(temp2, maxSteps, delta=delta)$Loglik
			LRTbootstrap[i] <- -2*(llC-tempLL1-tempLL2)
		}
	}
	
	pValue <- (sum(LRTbootstrap >= LRTobs)+1)/(numPerms+1)
	return(pValue)
}

getMLEandLoglike <- function(data, maxSteps=50, weightCols=NULL, delta=10^(-6)){
	if(missing(data))
		stop("A valid data set is required.")
		
	numSubs <- ncol(data)
	numEdges <- nrow(data)
	
	### If subject weighting we need to redefine numSubs
	if(!is.null(weightCols)) 
		numSubs <- sum(weightCols)
	
	### Generate inital g* from mean
	if(is.null(weightCols)){ 	# No weighting
		gstar <- rowSums(data)/ncol(data)
	}else{ 						# Subject weighting
		gstar <- apply(data, 1, function(x){
					(x%*%weightCols)/sum(weightCols)
				})
	}
	
	### Storage matrix for outputs
	fvals <- matrix(0, maxSteps+1, 4)
	colnames(fvals) <- c("f", "deltaf", "tau", "LL")
	
	### Set up starting points for our searching
	# f - Sum of squares of each observed connectome to g*
	calcfBase <- sqrt(colSums((data-gstar)^2))
	calcf <- sum(calcfBase)
	deltaf <- 1
	tau <- numSubs/calcf
	logLikBase <- lgamma(numEdges/2) - lgamma(numEdges) - numEdges*log(2) - numEdges/2*log(base::pi)
	logLik <- numSubs*(logLikBase + log(tau)) - numSubs	
	
	### Search for the best gstar
	count <- 1
	fvals[count,] <- c(calcf, deltaf, tau, logLik)
	while((count < maxSteps) && (deltaf > delta)){ 
		### Get an updated distance to every graph from the g*
		mse <- calcfBase
		
		### Get the num and den for g* calculation
		tempMSE <- matrix(rep(mse, nrow(data)), nrow(data), byrow=TRUE)
		numBase <- data[,mse!=0]/tempMSE[,mse!=0]
		denBase <- 1/mse[mse!=0]
		if(is.null(weightCols)){	# No weighting
			num <- rowSums(numBase) 
			den <- sum(denBase)
		}else{ 						# Subject weighting
			num <- numBase %*% weightCols
			den <- as.vector(denBase %*% weightCols)
		}
		
		### Update g* and calcf
		gstar <- num/den
		calcfBase <- sqrt(colSums((data-gstar)^2))
		if(is.null(weightCols)){	# No weighting
			calcf <- sum(calcfBase)
		}else{ 						# Subject weighting
			calcf <- as.vector(calcfBase %*% weightCols)
		}
		
		### Recalc new stat values
		deltaf <- abs(fvals[count, 1] - calcf)
		tau <- numSubs/calcf
		logLik <- numSubs*(logLikBase + log(tau)) - numSubs	
		
		### Save our values
		count <- count + 1
		fvals[count,] <- c(calcf, deltaf, tau, logLik)
	}
	
	### Trim our return stat values file before returning
	fvals <- fvals[1:count,]
	
	ret <- list(iters=count, Loglik=logLik, tau=tau, mleTree=as.data.frame(gstar), fvals=fvals)
	return(ret)
}

pairedCompareTwoDataSets  <- function(data1, data2, numPerms=1000, parallel=FALSE, cores=3, maxSteps=50, delta=10^(-6)){
	if(missing(data1) || missing(data2))
		stop("data is missing.")
	
	if(numPerms <= 0)
		stop("The number of permutations must be an integer greater than 0.")
	
	numSub <- ncol(data1)
	if(numSub != ncol(data2))
		stop("Groups must have the same number of subjects.")
	
	### Merge data1 and data2 together
	if(any(rownames(data1) != rownames(data2))){
		dataComb <- merge(data1, data2, by=0, all=TRUE)
		rownames(dataComb) <- dataComb[,1]
		dataComb <- dataComb[,-1]
		dataComb[is.na(dataComb)] <- 0
		
		### Pull data1 and data2 back out if we had to merge
		data1 <- dataComb[,1:numSub, drop=FALSE]
		data2 <- dataComb[,-c(1:numSub), drop=FALSE]
	}else{
		dataComb <- cbind(data1, data2)
	}
	
	### Get the starting gstar distance
	gstar1 <- getMLEandLoglike(data1, maxSteps, delta=delta)$mleTree
	gstar2 <- getMLEandLoglike(data2, maxSteps, delta=delta)$mleTree
	gstarDistance <- sqrt(sum((gstar1-gstar2)^2))
	
	if(parallel){
		cl <- parallel::makeCluster(min(cores, numPerms)) 
		doParallel::registerDoParallel(cl)
		tryCatch({ 
					permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .export=c("getMLEandLoglike")) %dopar%{
						### Randomly split each pair into seperate groups
						samps <- sample(0:1, numSub, replace=TRUE)
						samps <- samps*numSub + 1:numSub
						
						### Get gstar distance
						gstar1 <- getMLEandLoglike(dataComb[,samps], maxSteps, delta=delta)$mleTree
						gstar2 <- getMLEandLoglike(dataComb[,-samps], maxSteps, delta=delta)$mleTree
						return(sqrt(sum((gstar1-gstar2)^2)))
					}	
				}, finally = {				
					parallel::stopCluster(cl) # Close the parallel connections
				}
		)
	}else{
		permDistances <- rep(0, numPerms)
		for(i in 1:numPerms){ 	
			### Randomly split each pair into seperate groups
			samps <- sample(0:1, numSub, replace=TRUE)
			samps <- samps*numSub + 1:numSub
			
			### Get gstar distance
			gstar1 <- getMLEandLoglike(dataComb[,samps], maxSteps, delta=delta)$mleTree
			gstar2 <- getMLEandLoglike(dataComb[,-samps], maxSteps, delta=delta)$mleTree
			permDistances[i] <- sqrt(sum((gstar1-gstar2)^2))
		}
	}
	
	pvalue <- (sum(permDistances >= gstarDistance)+1)/(numPerms+1) 
	return(pvalue)
}



### ~~~~~~~~~~~~~~~~~~~~~
### transformation functions
### ~~~~~~~~~~~~~~~~~~~~~
trimToTaxaLevel <- function(data, level="genus", trimNodes="all", trimNames=TRUE, split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	# Get depth
	maxLevel <- getTaxaDepth(level)
	
	# Pull out row names
	rNames <- rownames(data)
	rownames(data) <- 1:nrow(data)
	nameSplit <- strsplit(rNames, split, fixed=TRUE)
	nameLength <- sapply(nameSplit, length)
	
	# Check we will still have data
	if(maxLevel > max(nameLength))
		stop("maxLevel is greater than the number of levels.")
	
	# Find the rows we want based on depth
	if(trimNodes == tolower("all")){
		sel <- nameLength == maxLevel
	}else if(trimNodes == tolower("above")){
		sel <- nameLength >= maxLevel
	}else if(trimNodes == tolower("below")){
		sel <- nameLength <= maxLevel
	}else{
		stop("trimNodes must be 'all', 'above', or 'below'.")
	}
	
	# Pull out the rows
	dataTemp <- data[sel,, drop=FALSE]
	
	# Return early if we are only trimming lower levels
	if(trimNodes == tolower("below") || (trimNodes == tolower("all") && !trimNames)){
		rownames(dataTemp) <- rNames[sel]
		return(dataTemp)
	}
	
	# Pull out row names again
	rNamesTemp <- rNames[sel]
	nameSplit <- stringr::str_split_fixed(rNamesTemp, stringr::fixed(split), maxLevel)
	
	# Check for duplicates and fix them
	dups <- duplicated(nameSplit[,maxLevel])
	dupCount <- 1
	if(any(dups))
		warning("Duplicated names detected and fixed.")
	while(any(dups)){
		# Find duplicated taxa
		dupTaxa <- unique(nameSplit[dups, maxLevel])
		
		# Go through each duplicated taxa
		for(t in dupTaxa){
			loc <- which(nameSplit[,maxLevel] == t)
			
			# Go through each duplicated row
			for(i in 1:length(loc))
				nameSplit[loc[i], maxLevel] <- paste(nameSplit[loc[i], maxLevel-dupCount], nameSplit[loc[i], maxLevel], sep="-")
		}
		
		# Ready for next loop
		dupCount <- dupCount + 1
		dups <- duplicated(nameSplit[,maxLevel])
	}
	
	# Add rownames back to data
	rownames(dataTemp) <- nameSplit[,maxLevel]
	
	return(dataTemp)
}

formatData <- function(data, countThreshold=1000, normalizeThreshold=10000){
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

mergeDataSets <- function(dataList, calcMLE=FALSE, uniqueNames=FALSE){
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



### ~~~~~~~~~~~~~~~~~~~~~
### plotting functions
### ~~~~~~~~~~~~~~~~~~~~~
plotTree <- function(treeList, colors=NULL, divisions=NULL, main=NULL, sub="", showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE){	
	if(missing(treeList))
		stop("At least one valid tree of type 'phylo' is required inside a list.")
	
	if(displayLegend)
		displayLegend(colors, divisions)
	
	### Plot every tree passed
	par(cex=.75)
	for(i in 1:length(treeList)){
		tree <- treeList[[i]]
		
		### Determine what color and size each branch should be
		branchData <- getBranchSizes(tree$edge.length, colors, divisions)
		
		### Figure out what we are titleing the plot
		if(length(main) == length(treeList)){
			main2 <- main[i]
		}else if(is.null(main)){
			main2 <- names(treeList)[i]
		}else{
			main2 <- main
		}	
		
		### Use ape to plot the tree
		ape::plot.phylo(tree, type="r", root.edge=FALSE, edge.color=branchData$edgecol, edge.width=branchData$edgewid, 
				show.tip.label=showTipLabel, show.node.label=showNodeLabel, main=main2, sub=sub)
	}
}

plotTreeDataMDS <- function(dataList, main="Tree MDS Comparisons", calcMLE=TRUE, mleTitles=NULL, dotColors=NULL, 
		dotSizes=NULL, showNames=FALSE, returnCoords=FALSE){
	if(missing(dataList))
		stop("At least 1 valid data set is required.")
	
	if(is.null(dotColors))
		dotColors <- c("red", "orange", "blue", "green", "yellow", "purple")
	
	if(is.null(dotSizes))
		dotSizes <- c(1, 2)
	if(length(dotSizes) != 2)
		stop("dotSizes must contain 2 integers.")
	
	numGrps <- length(dataList)
	
	### Add group names if we dont have them
	if(is.null(names(dataList))){
		grpNames <- paste("Data Set", 1:numGrps)
	}else{
		grpNames <- names(dataList)
	}
	
	### If we don't have enough colors, make a new set of colors
	if(length(dotColors) < numGrps)
		dotColors <- rainbow(numGrps)
	
	sizes <- NULL
	colors <- NULL
	if(calcMLE){
		mleDotSizes <- rep(0, numGrps)
		mleColors <- rep(0, numGrps)
	}
	### Figure out colors and sizes for all points
	for(i in 1:numGrps){
		tempData <- dataList[[i]]
		colors <- c(colors, rep(dotColors[i], ncol(tempData)))
		sizes <- c(sizes, rep(dotSizes[1], ncol(tempData)))
		
		if(calcMLE){
			mleColors[i] <- dotColors[i]
			mleDotSizes[i] <- dotSizes[2]	
		}
	}
	if(calcMLE){
		colors <- c(colors, mleColors, "black")
		sizes <- c(sizes, mleDotSizes, dotSizes[2])
	}
	
	### Get point positions and plot
	tData <- t(mergeDataSets(dataList, calcMLE))
	loc <- cmdscale(dist(tData), k=2)
	x <- loc[,1]
	y <- -loc[,2]
	plot(x, y, pch=19, xlab="MDS 1", ylab="MDS 2", pty="s", col=colors, cex=sizes, main=main)
	
	### Add mle titles if we aren't showing all the names
	if(calcMLE && !showNames){
		grpNames <- c(grpNames, "Combined MLE")
		text(x[(nrow(tData) - numGrps):nrow(tData)], y[(nrow(tData) - numGrps):nrow(tData)], grpNames, pos=3, cex=.75)
	}
	
	### Add the names of the samples
	if(showNames) 
		text(x, y, rownames(tData), pos=3, cex=.75)
	
	if(returnCoords)
		return(list(x=x, y=y))
}

createAndPlot <- function(data, samples=NULL, level="genus", colors=NULL, divisions=NULL, main=NULL, sub="", 
		showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE, onePerPage=FALSE, split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Turn the trees into neuwick format
	trees <- createTrees(data, level=level)
	
	### Change layout to plot 1 or 4 trees per page
	if(onePerPage){
		par(layout(1))
	}else{ #4 per page
		par(layout(matrix(c(1,3,2,4), 2, 2)))
	}
	
	### Plot the trees
	plotTree(trees, colors, divisions, main, sub, showTipLabel, showNodeLabel, displayLegend)
}

displayLegend <- function(colors=NULL, divisions=NULL, title="Confidence Value"){
	### Get colors and divisions to use
	csi <- getColorSizeInfo(colors, divisions)
	colors <- csi$col
	divisions <- csi$div
	
	### Get the legend wording
	lgd <- rep(NA, length(divisions))
	for(i in length(divisions):2)
		lgd[i] <- paste(divisions[i], "-", divisions[i-1], sep="")
	lgd[1] <- paste(divisions[1], "-0", sep="")
	
	### Plot the legend on a new page
	grDevices::palette(colors)
	plot(0, 0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
	legend(-.75, .75, legend=rev(lgd), col=rev(colors), pch=19, title=title)
}



### ~~~~~~~~~~~~~~~~~~~~~
### other functions
### ~~~~~~~~~~~~~~~~~~~~~
createTrees <- function(data, samples=NULL, level="genus", split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	if(any(grepl(")", rownames(data), fixed=TRUE)) || any(grepl("(", rownames(data), fixed=TRUE)) || 
			any(grepl(":", rownames(data), fixed=TRUE)) || any(grepl(",", rownames(data), fixed=TRUE)))
		stop("Using parentheses, commas, and/or colons in the taxa names is not allowed.")
	
	### Sort the data based on taxa names
	data <- data[order(rownames(data)),, drop=FALSE]
	
	### Create a newick format for each tree
	allTrees <- vector("list", ncol(data))
	for(i in 1:ncol(data)){
		oneSamp <- data[,i, drop=FALSE]
		if(sum(oneSamp) <= 0) # skips entries without data
			next
		
		allTrees[[i]] <- traverseTree(oneSamp, level, split)
	}
	names(allTrees) <- colnames(data)
	
	return(allTrees)
}

checkTreeValidity <- function(data, samples=NULL, epsilon=0.0001, split="."){
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

generateTree <- function(data, numReadsPerSamp, numSamps, est="mom", level="genus", split=".", theta=NULL){
	if(missing(data) || missing(numReadsPerSamp))
		stop("data and/or numReadsPerSamp missing.")
	
	### Take a full tree and pull out a single level
	tempdata <- trimToTaxaLevel(data, level, "all", FALSE, split)
	tempdata <- transformHMPTreetoHMP(tempdata, TRUE)
	
	### Get our starting shape
	if(tolower(est) == "mle"){
		dirfit <- dirmult::dirmult(tempdata)
	}else if(tolower(est) == "mom"){
		dirfit <- HMP::DM.MoM(tempdata)
	}else{
		stop("'est' must be either mom or mle.")
	}
	if(is.null(theta)){
		dirgamma <- dirfit$gamma
	}else{
		dirgamma <- dirfit$pi * ((1 - theta)/theta)
	}
	
	### Generate the data using HMP
	gendata <- HMP::Dirichlet.multinomial(rep(numReadsPerSamp, numSamps), dirgamma)
	colnames(gendata) <- colnames(tempdata)
	
	### Rotate back into HMPTree format
	gendata <- as.data.frame(t(gendata))
	
	### Build a full tree back out
	gendata <- buildTree(gendata, split)
	
	return(gendata)
}



### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Internal
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### ~~~~~~~~~~~~~~~~~~~~~
### recursion help functions
### ~~~~~~~~~~~~~~~~~~~~~
traverseTreeHelp <- function(data, place, treeLvl, maxTaxaDepth, split){
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
	}else{ # no children found
		if(myVal == 0 && REMOVELBL == TRUE){
			retStr <- paste(":", myVal, sep="")
		}else{
			retStr <- paste(myName, ":", myVal, sep="")
		}
	}
	
	return(retStr)
}

checkTreeValidHelp <- function(data, place, epsilon, split){
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



### ~~~~~~~~~~~~~~~~~~~~~
### other functions
### ~~~~~~~~~~~~~~~~~~~~~
traverseTree <- function(data, level, split){	
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
	
	return(retTree)
}

buildTree <- function(data, split="."){
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



### ~~~~~~~~~~~~~~~~~~~~~
### get info functions
### ~~~~~~~~~~~~~~~~~~~~~
getBranchSizes <- function(edgeLength, colors, divisions){
	numEdges <- length(edgeLength)
	
	edgeColor <- rep("white", numEdges)
	edgeWidth <- rep(0, numEdges)
	
	### Catch an all 0 tree
	if(max(edgeLength) == 0)
		return(list(edgecol=edgeColor, edgewid=edgeWidth))
	
	### Get colors and divisions to use
	csi <- getColorSizeInfo(colors, divisions)
	colors <- csi$col
	divisions <- csi$div
	
	### Get color and size for current tree
	for(i in 1:numEdges){ 	
		if(edgeLength[i] != 0){
			divLoc <- min(which(edgeLength[i] <= divisions))
			
			edgeColor[i] <- colors[divLoc]
			edgeWidth[i] <- floor(4*(divLoc+1)/length(divisions))
		}
	}
	
	retData <- list(edgecol=edgeColor, edgewid=edgeWidth)
	return(retData)
}

getColorSizeInfo <- function(colors, divisions){
	if(is.null(divisions))
		divisions <- c(.1, 1, 10, 100, 1000, 10000, Inf)
	divisions <- sort(divisions)
	
	if(is.null(colors))
		colors <- c("red", "orange", "yellow", "green" , "cyan", "blue", "black")
	
	numCols <- length(colors)
	numDivs <- length(divisions)
	
	### Check if we need more colors and add them
	if(numDivs > (numCols+1)) 
		colors <- rainbow(numDivs)
	
	return(list(col=colors, div=divisions))
}

getTaxaDepth <- function(level){
	if(tolower(level) == "kingdom" || tolower(level) == "k"){
		return(1)
	}else if(tolower(level) == "phylum" || tolower(level) == "p"){
		return(2)
	}else if(tolower(level) == "class" || tolower(level) == "c"){
		return(3)
	}else if(tolower(level) == "order" || tolower(level) == "o"){
		return(4)
	}else if(tolower(level) == "family" || tolower(level) == "f"){
		return(5)
	}else if(tolower(level) == "genus" || tolower(level) == "g"){
		return(6)
	}else if(tolower(level) == "species" || tolower(level) == "s"){
		return(7)
	}else if(tolower(level) == "subspecies" || tolower(level) == "ss"){
		return(8)
	}else if(tolower(level) == "variety" || tolower(level) == "v"){
		return(9)
	}else if(tolower(level) == "form" || tolower(level) == "fm"){
		return(10)
	}else{
		if(is.na(suppressWarnings(as.numeric(level))))
			stop(sprintf("%s isn't recognized.", as.character(level)))
		
		lvl <- as.numeric(level)
		if(lvl <= 0)
			stop("'level' must be a positive integer.")
		
		return(lvl)
	}
}



### ~~~~~~~~~~~~~~~~~~~~~
### hidden data functions
### ~~~~~~~~~~~~~~~~~~~~~
removeUnclass <- function(data, remove=TRUE){
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

subsetData <- function(data, site, region){
	if(!missing(site)){
		data <- subset(data, data[,1] == site)
		data <- data[,-1, drop=FALSE]
	}
	
	if(!missing(region)){
		data <- subset(data, data[,1] == region) 
		data <- data[,-1, drop=FALSE]
	}
	
	return(data)
}



### ~~~~~~~~~~~~~~~~~~~~~
### Only used in tree generation
### ~~~~~~~~~~~~~~~~~~~~~
transformHMPTreetoHMP <- function(data, elimZero=FALSE, zeroValue=.00001){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Find 0 taxa
	dataSum <- rowSums(data)
	zeroRows <- which(dataSum == 0)
	
	if(length(zeroRows) != 0){
		### Increase the first sample by zero value so taxa aren't all 0
		data[zeroRows, 1] <- zeroValue
		
		if(elimZero)
			data <- data[-zeroRows,]
	}
	
	return(t(data))
}

transformHMPtoHMPTree <- function(data){
	if(missing(data))
		stop("A valid data set is required.")
	
	data <- as.data.frame(t(data))
	
	return(data)
}











