pairedCompareTwoDataSets <-
function(data1, data2, numPerms=1000, parallel=FALSE, cores=3, maxSteps=50, delta=10^(-6)){
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
