\name{transformHMPtoHMPTree}
\alias{transformHMPtoHMPTree}
\title{Transform data from HMP to HMPTrees}
\description{This function will take data in the format needed in the package 'HMP' to the format used by 'HMPTrees'.}

\usage{transformHMPtoHMPTree(data)}

\arguments{
  \item{data}{A data frame of values in which each column is a separate taxa level and each row is a separate sample.}
}
\details{This function does the opposite of \code{transformHMPTreetoHMP}.}

\value{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	### Trims saliva to only contain the class level
	salivaClass <- trimToTaxaLevel(saliva, "class", TRUE)
	
	### This transforms the saliva data set but retains
	### any zero rows that may exist. 
	transSaliva <- transformHMPTreetoHMP(salivaClass, FALSE, 0)
	
	### saliva2 should be the same as salivaClass
	saliva2 <- transformHMPtoHMPTree(transSaliva)
}