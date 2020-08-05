
source("Tree Functions.R")


saliva <- read.csv("saliva.csv", row.names=1)
stool <- read.csv("stool.csv", row.names=1)
throat <- read.csv("throat.csv", row.names=1)


### ~~~~~~~~~~~~~~~~~~~~~
### pvalue functions
### ~~~~~~~~~~~~~~~~~~~~~
compareTwoDataSets.Test <- function(){
	data(saliva)
	data(stool)
	
	### We use 1 for the number of permutations for computation time
	### This value should be at least 1000 for an accurate result
	numPerms <- 1
	
	pval <- compareTwoDataSets(saliva, stool, numPerms)
	pval
}

getMLEandLoglike.Test <- function(){
	data(saliva)
	
	### We use 1 for the maximum number of steps for computation time
	### This value should be much higher to ensure an accurate result
	numSteps <- 1
	
	mle <- getMLEandLoglike(saliva, numSteps)$mleTree
}

pairedCompareTwoDataSets.Test <- function(){
	data(saliva)
	data(stool)
	
	### We use 1 for the number of permutations for computation time
	### This value should be at least 1000 for an accurate result
	numPerms <- 1
	
	pval <- pairedCompareTwoDataSets(saliva, stool, numPerms)
	pval
}



### ~~~~~~~~~~~~~~~~~~~~~
### transformation functions
### ~~~~~~~~~~~~~~~~~~~~~
trimToTaxaLevel.Test <- function(){
	data(saliva)
	
	### Trims saliva to only have the class level
	salivaClass <- trimToTaxaLevel(saliva, "class")
}

formatData.Test <- function(){
	data(saliva)
	
	saliva2 <- formatData(saliva, 1000, 10000)
}

mergeDataSets.Test <- function(){
	data(saliva)
	data(stool)
	
	dataComb <- mergeDataSets(list(saliva, stool))
}



### ~~~~~~~~~~~~~~~~~~~~~
### plotting functions
### ~~~~~~~~~~~~~~~~~~~~~
plotTree.Test <- function(){
	data(saliva)
	
	### Creates a tree for the 4th sample in 'Saliva' then plots it
	salivaTree <- createTrees(saliva[,4, drop=FALSE])
	plotTree(salivaTree, displayLegend=FALSE)
}

plotTreeDataMDS.Test <- function(){
	data(saliva)
	data(stool)
	
	plotTreeDataMDS(list(Saliva=saliva, Stool=stool))
}

createAndPlot.Test <- function(){
	data(saliva)
	
	### Plots the trees in column 2 and 3 in 'Saliva'
	createAndPlot(saliva[,2:3])
}

displayLegend.Test <- function(){
	displayLegend(c("red", "orange", "blue"), c(.1, 100, 10000))
}



### ~~~~~~~~~~~~~~~~~~~~~
### other functions
### ~~~~~~~~~~~~~~~~~~~~~
createTrees.Test <- function(){
	data(saliva)
	
	### Creates a tree for the 4th sample in 'Saliva'
	salivaTree <- createTrees(saliva[,4, drop=FALSE])
}

checkTreeValidity.Test <- function(){
	data(saliva) 
	
	validTree <- checkTreeValidity(saliva[,1, drop=FALSE])
	validTree
}

generateTree.Test <- function(){
	data(saliva)
	
	### Generate a the number of reads per sample
	### The first number is the number of reads and the second is the number of subjects
	nrs <- rep(10000, 2)
	
	gendata <- generateTree(saliva, nrs)
}






