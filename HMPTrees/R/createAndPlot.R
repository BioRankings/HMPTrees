createAndPlot <-
function(data, samples=NULL, level="genus", colors=NULL, divisions=NULL, main=NULL, sub="", 
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
