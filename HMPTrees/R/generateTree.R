generateTree <-
function(data, numReadsPerSamp, theta=NULL, level="genus", split="."){
	if(missing(data) || missing(numReadsPerSamp))
		stop("data and/or numReadsPerSamp missing.")
	
	### Take a full tree and pull out a single level
	tempdata <- trimToTaxaLevel(data, level, FALSE, split=split)
	tempdata <- transformHMPTreetoHMP(tempdata, TRUE)
	
	### Get our starting shape
	dirfit <- dirmult::dirmult(tempdata)
	if(is.null(theta)){
		dirgamma <- dirfit$gamma
	}else{
		dirgamma <- dirfit$pi * ((1 - theta)/theta)
	}
	
	### Generate the data using HMP
	gendata <- HMP::Dirichlet.multinomial(numReadsPerSamp, dirgamma)
	colnames(gendata) <- colnames(tempdata)
	
	### Rotate back into HMPTree format
	gendata <- as.data.frame(t(gendata))
	
	### Build a full tree back out
	gendata <- buildTree(gendata, split)
	
	return(gendata)
}
