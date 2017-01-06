
source("Tree Functions.R")


saliva <- read.csv("saliva.csv", row.names=1)

##########################################################
### generateTree
##########################################################
data(saliva)

gendata <- generateTree(saliva, 7000, 2)

##########################################################
### checkTreeValidity
##########################################################
data(saliva) 

validTree <- checkTreeValidity(saliva, 1)
validTree

##########################################################
### compareTwoDataSets
##########################################################
data(saliva)
data(stool)

### We use 1 for the number of boot straps for computation time
### This value should be at least 1000 for an accurate result
numBootStraps <- 1
pval <- compareTwoDataSets(saliva, stool, numBootStraps)
pval

##########################################################
### createAndPlot
##########################################################
data(saliva)

### Plots the trees in column 2 and 3 in 'Saliva'
createAndPlot(saliva, c(2:3))

##########################################################
### createTrees
##########################################################
data(saliva)

### Creates a object of type 'phylo' for the 4th tree in 'Saliva'
salivaTree <- createTrees(saliva, 4)

##########################################################
### displayLegend
##########################################################
displayLegend(c("red", "orange", "blue"), c(.1, 100, 10000))

##########################################################
### formatData
##########################################################
data(throat)

throat <- formatData(throat, 1000, 10000)

##########################################################
### getMLEandLoglike
##########################################################
data(saliva)

### We use 1 for the maximum number of steps for computation time
### This value should be much higher to ensure an accurate result
numSteps <- 1
mle <- getMLEandLoglike(saliva, numSteps)$mleTree

##########################################################
### mergeDataSets
##########################################################
data(saliva)
data(stool)

dataComb <- mergeDataSets(list(saliva, stool), FALSE, TRUE)

##########################################################
### plotTree
##########################################################
data(saliva)

### Creates a object of type 'phylo' for the 4th tree in 'Saliva'
### Then plots it
salivaTree <- createTrees(saliva, 4)
plotTree(salivaTree, displayLegend=FALSE)

##########################################################
### plotTreeDataMDS
##########################################################
data(saliva)
data(stool)

plotTreeDataMDS(list(saliva, stool), mleTitles=c("Saliva", "Stool"))

##########################################################
### transformHMPtoHMPTree
##########################################################
data(saliva)

### Trims saliva to only contain the class level
salivaClass <- trimToTaxaLevel(saliva, "class", TRUE)

### This transforms the saliva data set but retains
### any zero rows that may exist. 
transSaliva <- transformHMPTreetoHMP(salivaClass, FALSE, 0)

### saliva2 should be the same as salivaClass
saliva2 <- transformHMPtoHMPTree(transSaliva)

##########################################################
### transformHMPTreetoHMP
##########################################################
data(saliva)

### Trims saliva to only contain the class level
salivaClass <- trimToTaxaLevel(saliva, "class", TRUE)

### This transforms the saliva data set but retains
### any zero rows that may exist. 
transSaliva <- transformHMPTreetoHMP(salivaClass, FALSE, 0)

##########################################################
### trimToTaxaLevel
##########################################################
data(saliva)

### Trims saliva to only have the class level
salivaClass <- trimToTaxaLevel(saliva, "class", TRUE)

##########################################################
### pairedcompareTwoDataSets
##########################################################
data(saliva)
data(stool)

### We use 1 for the number of boot straps for computation time
### This value should be at least 1000 for an accurate result
numPerms <- 1
pval <- compareTwoDataSets(saliva, stool, numPerms)
pval



