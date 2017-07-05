plotTree <-
function(treeList, colors=NULL, divisions=NULL, main=NULL, sub="", showTipLabel=TRUE, showNodeLabel=FALSE, displayLegend=TRUE, trees=NULL){
	### Fix any old arguments
	if(!is.null(trees)){
		warning("'trees' is deprecated. It has been replaced with treeList. View the help files for details.")
		treeList <- trees
	}
	
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
