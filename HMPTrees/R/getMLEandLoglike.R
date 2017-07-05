getMLEandLoglike <-
function(data, maxSteps=50, weightCols=NULL, delta=10^(-6), weight=NULL){
	if(missing(data))
		stop("A valid data set is required.")
	
	### Fix for anyone using the old weighting argument
	if(!is.null(weight)){
		warning("'weight' is deprecated. It has been replaced with weightCols. View the help files for details.")
		weightCols <- weight
	}
	
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
	
	### Set up starting points for our searching
	# f - Sum of squares of each observed connectome to g*
	calcf <- sum(apply(data, 2, function(x, gs) {sqrt(sum((x - gs)^2))}, gs=gstar))
	deltaf <- 1
	tau <- numSubs/calcf
	logLikBase <- lgamma(numEdges/2) - lgamma(numEdges) - numEdges*log(2) - numEdges/2*log(base::pi)
	logLik <- numSubs*(logLikBase + log(tau)) - numSubs	
	
	### Search for the best gstar
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
