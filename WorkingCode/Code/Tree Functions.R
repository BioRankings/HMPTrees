
library(ape)			# plotting and creating tree objects
library(HMP)			# creating new trees
library(dirmult)		# creating new trees
library(doParallel)		# parallelizing pvalue calc

traverseTreeHelp <- function(data, place, treeLvl=1, maxTaxaDepth, split){
	# Set some base
	REMOVELBL <- TRUE  	#set to false to show all labels, otherwise only non-zero labels will be shown	
	childStr <- ""				
	myName <- rownames(data)[place]		
	myVal <- data[place, 1]
	stopTree <- FALSE
	
	if(treeLvl >= maxTaxaDepth)
		stopTree <- TRUE
	
	data <- data[-place, , drop=FALSE] 
	childRows <- grep(paste(myName, split, sep=""), rownames(data), fixed=TRUE)
	
	noSiblings <- TRUE
	if(length(childRows) != 0 && !stopTree){ #we know we have children
		newData <- data[childRows,, drop=FALSE]
		
		rownames(newData) <- substring(rownames(newData), nchar(myName)+2, nchar(rownames(newData)))
		nameSplit <- strsplit(rownames(newData), split, fixed=TRUE)
		
		for(t in 1:nrow(newData)){ 	
			if(length(nameSplit[[t]]) == 1){ 
				temp_hstr <- traverseTreeHelp(newData, t, treeLvl+1, maxTaxaDepth, split)
				if(childStr != ""){
					childStr <- paste(childStr, ",", temp_hstr, sep="")
				}else{
					childStr <- temp_hstr
				}
				
				noSiblings <- FALSE
			}
		}
		
		if(myVal == 0){ #set up string to be returned
			retStr <- paste("(", childStr, "):", myVal,  sep="")
		}else{
			retStr <- paste("(", childStr, ")", myName, ":", myVal,  sep="")
		}
		
		if(noSiblings) #adds an extra 0 value branch if there are no siblings
			retStr <- paste(retStr, ",:0.0", sep="")
	}else{ #no childern found
		if(myVal == 0 && REMOVELBL == TRUE){
			retStr <- paste(":", myVal, sep="")
		}else{
			retStr <- paste(myName, ":", myVal, sep="")
		}
		
		retStr <- paste(retStr, ",:0.0", sep="")
	}
	
	return(retStr)
}

traverseTree <- function(data, level="genus", split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	if(any(grepl(")", rownames(data), fixed=TRUE)) || any(grepl("(", rownames(data), fixed=TRUE)) || any(grepl(":", rownames(data), fixed=TRUE)))
		stop("Using parentheses and/or colons in the taxa names is not allowed.")
	
	myStr <- ""
	splitLength <- unlist(lapply(strsplit(rownames(data), split, fixed=TRUE), length))
	startPlace <- which(splitLength == 1)
	maxTaxaDepth <- getTaxaDepth(level)
	
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
	
	myTree <- paste("", myStr, ";", sep="")
	
	retTree <- ape::read.tree(text=myTree) #turns the newick format into a tree
	retTree <- ape::collapse.singles(retTree)
	
	return(retTree)
}

getTaxaDepth <- function(level="genus"){	
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
	}else{
		if(is.na(suppressWarnings(as.numeric(level))))
			stop(sprintf("%s isn't recognized.", as.character(level)))
		
		lvl <- as.numeric(level)
		if(lvl > 0)
			return(lvl)
		
		stop("'level' cannot be negative.")
	}
}

formatData <- function(data, countThreshold=1000, normalizeThreshold=10000){
	if(missing(data))
		stop("A valid data set is required.")
	
	data <- data[order(rownames(data)),,drop=FALSE]
	data[is.na(data)] <- 0
	data <- data[, data[1,] >= countThreshold, drop=FALSE]
	
	if(ncol(data) == 0)
		stop("'countThreshold' is too high.")
	
	if(normalizeThreshold > 0){
		for(i in ncol(data):1)
			data[,i] <- data[,i] * (normalizeThreshold/data[1,i])
	}
	
	return(data)
}

createTrees <- function(data, samples=1, level="genus"){
	if(missing(data))
		stop("A valid data set is required.")
	
	if(samples[1] == 0)
		samples <- 1:ncol(data)
	
	data <- data[order(rownames(data)),,drop=FALSE]
	
	allTrees <- list()
	i <- 1
	for(sample in samples){
		if(is.na(as.numeric(sample)))
			stop(sprintf("%s must be a number", as.character(sample)))
		if(sample > ncol(data)) 
			stop(sprintf("%s is larger than the bounds of the data set", as.character(sample)))
		if(sample < 1)
			stop(sprintf("%s is smaller than the bounds of the data set", as.character(sample)))
		
		oneSamp <- data[, sample, drop=FALSE]
		
		if(sum(oneSamp) <= 0)#skips entries without data
			next
		tempTree <- traverseTree(oneSamp, level)
		allTrees[[i]] <- tempTree
		i <- i + 1
	}
	
	names(allTrees) <- colnames(data)[samples]
	
	return(allTrees)
}

plotTree <- function(trees, colors, divisions, main, sub, 
		showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE){
	
	if(missing(trees))
		stop("At least one valid tree of type 'phylo' is required inside a list.")
	
	if(missing(sub))
		sub <- ""
	
	genMain <- missing(main)
	if(genMain)
		main <- NULL
	passMain <- length(main) == length(trees)
	
	if(displayLegend)
		displayLegend(colors, divisions)
	
	par(cex=.75)
	
	for(i in 1:length(trees)){
		tr <- trees[[i]]
		bs <- getBranchSizes(tr, colors, divisions)
		
		if(passMain){
			main2 <- main[i]
		}else if(genMain){
			main2 <- names(trees)[i]
		}else{
			main2 <- main[1]
		}	
		
		ape::plot.phylo(tr, type="r", root.edge=FALSE, edge.color=bs$edgecol, edge.width=bs$edgewid, 
				show.tip.label=showTipLabel, show.node.label=showNodeLabel, main=main2, sub=sub)
	}
}

plotTreeDataMDS <- function(data, main="Tree MDS Comparisons", calcMLE=TRUE, mleTitles, 
		dotColors=c("red", "orange", "blue", "green", "yellow", "purple"), dotSizes=c(1, 2), 
		showNames=FALSE, returnCoords=FALSE){
	normalDotSizes <- NULL
	mleDotSizes <- NULL
	normalColors <- NULL
	mleColors <- NULL
	mleLabels <- NULL
	
	if(missing(data))
		stop("At least 1 valid data set is required.")
	
	if(class(data) == "data.frame" || class(data) == "matrix"){ #turn a single dataset into a list
		tempData <- data
		data <- NULL
		data[[1]] <- tempData
	}
	
	if(length(dotColors) < length(data))
		dotColors <- rainbow(length(data))
	
	twoColors <- length(dotColors) >= length(data)*2 #2 colors per data set
	dataCount <- length(data)
	titles <- !missing(mleTitles)
	if(titles) #make sure we have the same number of titles as data sets
		titles <- length(mleTitles) == length(data)
	
	for(i in 1:dataCount){
		tempData <- data[[i]]
		if(twoColors){
			colorLoc <- i*2-1
			colorLoc2 <- colorLoc+1
		}else{
			colorLoc <- i
			colorLoc2 <- i
		}
		
		normalColors <- c(normalColors, rep(dotColors[colorLoc], ncol(tempData)))
		mleColors <- c(mleColors, dotColors[colorLoc2])
		normalDotSizes <- c(normalDotSizes, rep(dotSizes[1], ncol(tempData)))
		mleDotSizes <- c(mleDotSizes, dotSizes[2])
		if(calcMLE && titles)
			mleLabels <- c(mleLabels, mleTitles[i])
	}
	
	mData <- mergeDataSets(data, calcMLE)
	tData <- t(mData)
	loc <- cmdscale(dist(tData), k=2)
	x <- loc[,1]
	y <- -loc[,2]
	
	colors <- c(normalColors, mleColors, "black")
	sizes <- c(normalDotSizes, mleDotSizes, dotSizes[2])
	plot(x, y, pch=19, xlab="MDS 1", ylab="MDS 2", pty="s", col=colors, cex=sizes, main=main)
	
	if(calcMLE && !showNames){ #only plot mle titles if we are calculating them
		if(dataCount > 1){ #dont add a combined mle if we dont make one
			if(titles){
				mleTitles <- c(mleTitles, "Combined MLE")
			}else{
				mleTitles <- c(paste("Data", 1:length(data), "MLE"),  "Combined MLE")
			}
			dataCount <- dataCount + 1
		}
		text(x[(nrow(tData)-dataCount+1):nrow(tData)], y[(nrow(tData)-dataCount+1):nrow(tData)], mleTitles, pos=3, cex=.75)
	}
	
	if(showNames) #Plot the names of the samples
		text(x, y, colnames(mData), pos=3, cex=.75)
	
	if(returnCoords)
		return(list(x=x, y=y))
}

mergeDataSets <- function(data, calcMLE=TRUE, uniqueNames=FALSE){
	newData <- NULL
	mles <- NULL
	
	if(missing(data))
		stop("At least 1 valid data set is required.")
	
	if(class(data) == "data.frame" || class(data) == "matrix"){ #turn a single dataset into a list
		temp <- data
		data <- NULL
		data[[1]] <- temp
	}
	
	for(i in 1:length(data)){
		if(uniqueNames)
			colnames(data[[i]]) <- paste("Data", i, "-", colnames(data[[i]]))
		newData <- merge(newData, data[[i]], by=0, all=TRUE)
		rownames(newData) <- newData[,1]
		newData <- newData[,-1]
		
		if(calcMLE)
			mles[[i]] <- getMLEandLoglike(data[[i]])$mleTree
	}
	newData[is.na(newData)] <- 0 #set any NA's to 0
	
	if(calcMLE){
		mleAll <- getMLEandLoglike(newData)$mleTree
		colnames(mleAll) <- "combined mle"
		
		if(length(data) > 1){
			for(i in 1:length(mles)){
				newData <- merge(newData, mles[[i]], by=0, all=TRUE)
				rownames(newData) <- newData[,1]
				newData <- newData[,-1]
				
				colnames(newData)[ncol(newData)] <- paste("data", i, "mle", sep="")
			}	
			newData[is.na(newData)] <- 0 #set any NA's to 0
		}
		
		newData <- cbind(newData, mleAll)
	}
	
	return(newData)
}

createAndPlot <- function(data, samples=1, level="genus", colors, divisions, main, 
		sub, showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE, onePerPage=FALSE){
	
	if(missing(data))
		stop("A valid data set is required.")
	
	trees <- createTrees(data, samples, level)
	
	if(onePerPage){
		par(layout(1))
	}else{ #4 per page
		par(layout(matrix(c(1,3,2,4), 2, 2)))
	}
	plotTree(trees, colors, divisions, main, sub, 
			showTipLabel, showNodeLabel, displayLegend)
	
	par(layout(1)) #reset layout back
}

displayLegend <- function(colors, divisions, title="Confidence Value"){
	if(missing(divisions)) 
		divisions <- c(0, .1, 1, 10, 100, 1000, 10000)
	
	divisions <- sort(divisions)
	if(missing(colors)) 
		colors <- c("red", "orange", "yellow", "green" , "cyan", "blue")
	
	if(length(divisions) > (length(colors)+1)) #need more colors, dont care if more colors than divisons
		colors <- c(colors, rep(colors[length(colors)], (length(divisions) - length(colors)-1)))
	
	lgd <- NULL
	for(num in length(divisions):2)
		lgd <- c(lgd, paste(divisions[num], "-", divisions[num-1], sep=""))
	
	palette(colors)
	plot(0, 0, type="n", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
	legend(-.75, .75, legend=lgd, col=rev(palette()), pch=19, title=title)
}

getBranchSizes <- function(tree, colors, divisions){
	if(missing(tree) || class(tree) != "phylo")
		stop("A valid tree is required.")
	
	edgeColor <- NULL
	edgeWidth <- NULL
	edgeLength <- tree$edge.length
	
	if(max(edgeLength) == 0){ #catch an all 0 tree
		retData <- list(edgecol=rep(0, length(edgeLength)), edgewid=rep(0, length(edgeLength)))
		return(retData)
	}
	
	if(missing(divisions))
		divisions <- c(.1, 1, 10, 100, 1000, 10000)
	
	divisions <- sort(divisions)
	
	if(missing(colors))
		colors <- c("red", "orange", "yellow", "green" , "cyan", "blue")
	
	if(length(divisions) > (length(colors)+1)) #need more colors, dont care if more colors than divisons
		colors <- c(colors, rep(colors[length(colors)], (length(divisions) - length(colors)-1)))
	
	palette(colors)
	
	for(i in 1:length(edgeLength)){ 	
		if(edgeLength[i] == 0){ #0 value so make it white
			edgeColor <- append(edgeColor, 0)
			edgeWidth <- append(edgeWidth, 0)
		}else{
			for(j in 1:length(divisions)){
				if(edgeLength[i] <= divisions[j]){
					cval <- j
					len <- floor(4*(j+1)/length(divisions)) 
					edgeColor <- append(edgeColor, cval)
					edgeWidth <- append(edgeWidth, len)
					break
				}
			}
		}
	}
	
	retData <- list(edgecol=edgeColor, edgewid=edgeWidth)
	return(retData)
}

compareTwoDataSets <- function(data1, data2, numBootStraps=1000, enableMC=FALSE, cores=3, maxSteps=50, delta=10^(-6)){
	if(missing(data1) || missing(data2))
		stop("Two valid data sets are required.")
	
	if(numBootStraps <= 0)
		stop("The number of boostraps must be an integer greater than 0.")
	
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
	
	if(enableMC){
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		LRTbootstrap <- foreach::foreach(i=1:numBootStraps, .combine=c, .inorder=FALSE, .multicombine=TRUE, .export=c("getMLEandLoglike")) %dopar%{
			### Get a random sampling of the data and break into 2 groups
			samps <- sample(numSubC, replace=TRUE)
			tempC <- dataComb[,samps, drop=FALSE]
			temp1 <- tempC[,1:numSub1, drop=FALSE]
			temp2 <- tempC[,-c(1:numSub1), drop=FALSE]
			
			### Get logliks for the sampled data
			tempLL1 <- getMLEandLoglike(temp1, maxSteps, delta=delta)$Loglik
			tempLL2 <- getMLEandLoglike(temp2, maxSteps, delta=delta)$Loglik
			tempLLC <- getMLEandLoglike(tempC, maxSteps, delta=delta)$Loglik
			LRT <- -2*(tempLLC-tempLL1-tempLL2)
			
			return(LRT)
		}
		parallel::stopCluster(cl) 
	}else{ 	
		LRTbootstrap <- rep(0, numBootStraps)
		for(i in 1:numBootStraps){
			### Get a random sampling of the data and break into 2 groups
			samps <- sample(numSubC, replace=TRUE)
			tempC <- dataComb[,samps, drop=FALSE]
			temp1 <- tempC[,1:numSub1, drop=FALSE]
			temp2 <- tempC[,-c(1:numSub1), drop=FALSE]
			
			### Get logliks for the sampled data
			tempLL1 <- getMLEandLoglike(temp1, maxSteps, delta=delta)$Loglik
			tempLL2 <- getMLEandLoglike(temp2, maxSteps, delta=delta)$Loglik
			tempLLC <- getMLEandLoglike(tempC, maxSteps, delta=delta)$Loglik
			LRTbootstrap[i] <- -2*(tempLLC-tempLL1-tempLL2)
		}
	}
	
	pValue <- sum(LRTbootstrap > LRTobs)/numBootStraps	
	return(pValue)
}

getMLEandLoglike <- function(data, maxSteps=50, weightCols=NULL, delta=10^(-6), weight=NULL){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Fix for anyone using the old weighting argument
	if(is.null(weightCols) && !is.null(weight))
		weightCols <- weight
	
	numSubs <- ncol(data)
	numEdges <- nrow(data)
	
	### If subject weighting we need to redefine numSubs
	if(!is.null(weightCols)) 
		numSubs <- sum(weightCols)
	
	### Work around for adding mse into applys
	dataTemp <- as.data.frame(data)
	
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
	
	### f - Sum of squares of each observed connectome to g*
	calcf <- sum(apply(data, 2, function(x, gs) {sqrt(sum((x - gs)^2))}, gs=gstar))
	deltaf <- 1
	tau <- numSubs/calcf
	logLikBase <- lgamma(numEdges/2) - lgamma(numEdges) - numEdges*log(2) - numEdges/2*log(base::pi)
	logLik <- numSubs*(logLikBase + log(tau)) - numSubs	
	
	count <- 1
	fvals[count,] <- c(calcf, deltaf, tau, logLik)
	while((count < maxSteps) && (deltaf > delta)){ 
		### Get an updated distance to every graph from the g*
		mse <- apply(data, 2, function(x, gs) {sqrt(sum((x - gs)^2))}, gs=gstar)
		dataTemp[nrow(data)+1,] <- mse
		
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
		calcfBase <- apply(data, 2, function(x, gs) {sqrt(sum((x - gs)^2))}, gs=gstar)
		if(is.null(weightCols)){	# No weighting
			calcf <- sum(calcfBase)
		}else{ 												# Subject weighting
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

checkTreeValidity <- function(data, samples=1, epsilon=0.0001, split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	overAllValid <- NULL
	if(samples[1] == 0)
		samples <- 1:ncol(data)
	
	for(sample in samples){
		if(is.na(as.numeric(sample)))
			stop(sprintf("%s must be a number", as.character(sample)))
		if(sample > ncol(data))
			stop(sprintf("%s is larger than the bounds of the data set", as.character(sample)))
		if(sample < 1)
			stop(sprintf("%s is smaller than the bounds of the data set", as.character(sample)))
		
		tempData <- data[, sample, drop=FALSE]
		nameSplit <- strsplit(rownames(tempData), split, fixed=TRUE)
		
		for(i in 1:nrow(tempData)){ 	
			if(length(nameSplit[[i]]) == 1)	
				valid <- checkTreeValidHelp(tempData, i, epsilon, split)
		}
		overAllValid <- c(overAllValid, valid)
	}
	
	return(overAllValid)
}

checkTreeValidHelp <- function(data, place, epsilon, split){
	myName <- rownames(data)[place]	
	myVal <- data[place, 1]
	
	data <- data[-place,, drop=FALSE] 
	childRows <- grep(paste(myName, split, sep=""), rownames(data), fixed=TRUE)
	
	if(length(childRows) != 0){ #we know we have children
		newData <- data[childRows,, drop=FALSE]
		rownames(newData) <- substring(rownames(newData), nchar(myName)+2, nchar(rownames(newData)))	
		splitLength <- unlist(lapply(strsplit(rownames(newData), split, fixed=TRUE), length))
		startPlace <- which(splitLength == 1)
		
		childCount <- 0
		for(t in startPlace){ 	
			childCount <- childCount + newData[t, 1]
			if((myVal+epsilon) < childCount) #the child is greater so bad tree
				return(FALSE)
			valid <- checkTreeValidHelp(newData, t, epsilon, split)
			if(!valid)
				return(FALSE)
		}
	}
	
	return(myVal >= 0) #check bottom nodes arent negative
}

generateTree <- function(data, nreads=10000, nsamps=50, theta=0, level="genus", split="."){
	if(missing(data))
		stop("A valid data set is required.")
	if(nreads <= 0)
		stop("'nreads' must be positive and greater than 0.")
	if(nsamps <= 0)
		stop("'nsamps' must be positive and greater than 0.")
	
	tempdata <- trimToTaxaLevel(data, level, FALSE, split=split)
	tempdata <- transformHMPTreetoHMP(tempdata, TRUE)
	
	if(theta > 0 && theta < 1){
		dirfit <- dirmult::dirmult(tempdata) 
		dirgamma <- dirfit$pi * ((1 - theta)/theta)
	}else{
		dirfit <- HMP::DM.MoM(tempdata)
		dirgamma <- dirfit$gamma
	}
	
	gendata <- HMP::Dirichlet.multinomial(rep(nreads, nsamps), dirgamma)
	colnames(gendata) <- colnames(tempdata)
	gendata <- transformHMPtoHMPTree(gendata)
	
	gendata <- buildTree(gendata, split)
	
	return(gendata)
}

trimToTaxaLevel <- function(data, level="genus", eliminateParentNodes=FALSE, trimBelow=NULL, split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
	maxLevel <- getTaxaDepth(level)
	
	nameSplit <- strsplit(as.character(rownames(data)), split, fixed=TRUE)
	lowerLevels <- NULL
	for(l in 1:nrow(data)){ 
		if(length(nameSplit[[l]]) == maxLevel){
			lowerLevels <- c(lowerLevels, l)
			if(eliminateParentNodes && is.null(trimBelow))
				rownames(data)[l] <- nameSplit[[l]][maxLevel]
		}else if(!is.null(trimBelow) && is.logical(trimBelow)){
			if(length(nameSplit[[l]]) < maxLevel && trimBelow){
				lowerLevels <- c(lowerLevels, l)
			}else if(length(nameSplit[[l]]) > maxLevel && !trimBelow){
				lowerLevels <- c(lowerLevels, l)
			}
		}
	}
	data <- data[lowerLevels, , drop=FALSE]
	
	return(data)
}

transformHMPTreetoHMP <- function(data, elimZero=FALSE, zeroValue=.00001){
	if(missing(data))
		stop("A valid data set is required.")
	
	dataSum <- rowSums(data)
	loc <- NULL
	for(r in 1:nrow(data)){
		if(dataSum[r] == 0){
			loc <- c(loc, r)
			data[r, 2] <- zeroValue
		}
	}
	if(elimZero)
		data <- data[-loc,]
	
	return(t(data))
}

transformHMPtoHMPTree <- function(data){
	if(missing(data))
		stop("A valid data set is required.")
	
	data <- as.data.frame(t(data))
	
	return(data)
}

#Works on our data only
removeUnclass <- function(data, remove=TRUE){
	falsePos <- grep(".U", rownames(data), fixed=TRUE)
	unclass <- grep("U", rownames(data), fixed=TRUE)
	
	if(length(unclass) == 0) #nothing to remove
		return(data)
	
	if(length(falsePos) != 0) #we have false positives
		unclass <- setdiff(unclass, falsePos)
	
	if(remove){
		data <- data[-unclass,,drop=FALSE]
	}else{
		data[unclass,] <- 0
	}
	
	return(data)		
}

#Works on our data only
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
		cl <- parallel::makeCluster(cores) 
		doParallel::registerDoParallel(cl)
		
		permDistances <- foreach::foreach(i=1:numPerms, .combine=c, .inorder=FALSE, .multicombine=TRUE, .export=c("getMLEandLoglike")) %dopar%{
			### Randomly split each pair into seperate groups
			samps <- sample(0:1, numSub, replace=TRUE)
			samps <- samps*numSub + 1:numSub
			
			### Get gstar distance
			gstar1 <- getMLEandLoglike(dataComb[,samps], maxSteps, delta=delta)$mleTree
			gstar2 <- getMLEandLoglike(dataComb[,-samps], maxSteps, delta=delta)$mleTree
			return(sqrt(sum((gstar1-gstar2)^2)))
		}	
		parallel::stopCluster(cl) 
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
	
	Pvalue <- sum(permDistances >= gstarDistance)/numPerms 
	return(Pvalue)
}

buildTree <- function(data, split="."){
	if(missing(data))
		stop("A valid data set is required.")
	
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
