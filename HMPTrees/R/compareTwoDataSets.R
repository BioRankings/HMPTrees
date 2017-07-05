compareTwoDataSets <-
function(data1, data2, numPerms=1000, parallel=FALSE, cores=3, maxSteps=50, delta=10^(-6), numBootStraps=NULL, enableMC=NULL){
	if(missing(data1) || missing(data2))
		stop("Two valid data sets are required.")
	
	### Fix any old arguments
	if(!is.null(numBootStraps)){
		warning("'numBootStraps' is deprecated. It has been replaced with numPerms. View the help files for details.")
		numPerms <- numBootStraps
	}
	if(!is.null(enableMC)){
		warning("'enableMC' is deprecated. It has been replaced with parallel. View the help files for details.")
		parallel <- enableMC
	}
	
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
