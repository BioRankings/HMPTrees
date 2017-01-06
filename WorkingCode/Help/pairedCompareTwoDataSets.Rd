\name{pairedCompareTwoDataSets}
\alias{pairedCompareTwoDataSets}
\title{Likelihood-Ratio-Test Statistics to Compare the Distribution of 2 Paired Sets of RDP-Based Taxonomic Trees}
\description{
	This functions compares the distribution of two paired sets of RDP-based taxonomic trees using Likelihood-Ratio-Test statistics and a 
	p-value is computed using permutation.
}

\usage{pairedCompareTwoDataSets(data1, data2, numPerms = 1000, parallel = FALSE, cores = 3)}

\arguments{
  \item{data1, data2}{Data frames in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{numPerms}{The number of times to run the permutation test, the default is 1000.}
  \item{parallel}{When this is 'TRUE' it allows for parallel calculation of the bootstraps. (See Note 1 in Details).}
  \item{cores}{The number of parallel processes to run if parallel is 'TRUE'.}
}

\details{
	Notes: 
	\enumerate{
		\item Enabling parallel calculation requires the package \code{doParallel}.
		\item Both data sets should be standardized to the same number of reads.
	}
}

\value{A p-value for the similarity of the two data sets based on the bootstrapping.}

\author{Patricio S. La Rosa, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	data(saliva)
	data(stool)
	
	### We use 1 for the number of boot straps for computation time
	### This value should be at least 1000 for an accurate result
	numPerms <- 1
	pval <- compareTwoDataSets(saliva, stool, numPerms)
	pval
}