getBranchSizes <-
function(edgeLength, colors, divisions){
	edgeColor <- NULL
	edgeWidth <- NULL
	
	### Catch an all 0 tree
	if(max(edgeLength) == 0){ 
		retData <- list(edgecol=rep(0, length(edgeLength)), edgewid=rep(0, length(edgeLength)))
		return(retData)
	}
	
	if(is.null(divisions))
		divisions <- c(.1, 1, 10, 100, 1000, 10000)
	divisions <- sort(divisions)
	
	if(is.null(colors))
		colors <- c("red", "orange", "yellow", "green" , "cyan", "blue")
	
	### Check if we need more colors and add them
	if(length(divisions) > (length(colors)+1)) 
		colors <- c(colors, rep(colors[length(colors)], (length(divisions) - length(colors)-1)))
	palette(colors)
	
	for(i in 1:length(edgeLength)){ 	
		if(edgeLength[i] == 0){ #0 value so make it white
			edgeColor <- c(edgeColor, 0)
			edgeWidth <- c(edgeWidth, 0)
		}else{
			for(j in 1:length(divisions)){
				if(edgeLength[i] <= divisions[j]){
					edgeColor <- c(edgeColor, j)
					edgeWidth <- c(edgeWidth, floor(4*(j+1)/length(divisions)) )
					break
				}
			}
		}
	}
	
	retData <- list(edgecol=edgeColor, edgewid=edgeWidth)
	return(retData)
}
