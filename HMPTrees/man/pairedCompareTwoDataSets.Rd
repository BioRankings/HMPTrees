\name{pairedCompareTwoDataSets}
\alias{pairedCompareTwoDataSets}
\title{Likelihood-Ratio-Test Statistics to Compare the Distribution of 2 Paired Sets of RDP-Based Taxonomic Trees}
\description{
	This functions compares the distribution of two paired sets of RDP-based taxonomic trees using Likelihood-Ratio-Test statistics and a 
	p-value is computed using permutation.
}

\usage{
	pairedCompareTwoDataSets(data1, data2, numPerms = 1000, parallel = FALSE, 
		cores = 3, maxSteps=50, delta=10^(-6))
}

\arguments{
  \item{data1, data2}{Data frames in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{numPerms}{Number of permutations.  In practice this should be at least 1,000.}
  \item{parallel}{When this is 'TRUE' it allows for parallel calculation of the permutations. Requires the package \code{doParallel}.}
  \item{cores}{The number of parallel processes to run if parallel is 'TRUE'.}
  \item{maxSteps}{The maximum number of times to iterate though for the MLE.}
  \item{delta}{The minimum threshold of change in f to stop the search for the MLE.}
}

\details{
	Note: Both data sets should be standardized to the same number of reads.
}

\value{A p-value for the similarity of the two data sets based on the permutation test.}

\author{Patricio S. La Rosa, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	data(saliva)
	data(stool)
	
	### We use 1 for the number of permutations for computation time
	### This value should be at least 1000 for an accurate result
	numPerms <- 1
	
	pval <- pairedCompareTwoDataSets(saliva, stool, numPerms)
	pval
}