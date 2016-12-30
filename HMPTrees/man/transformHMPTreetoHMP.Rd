\name{transformHMPTreetoHMP}
\alias{transformHMPTreetoHMP}
\title{Transform data from HMPTrees to HMP}
\description{This function will take data in the format needed in the package 'HMPTrees' to the format used by 'HMP'.}

\usage{transformHMPTreetoHMP(data, elimZero = FALSE, zeroValue = 1e-05)}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{elimZero}{If 'TRUE' this will eliminate taxa levels that contain only zero values and 'zeroValue' is ignored.}
  \item{zeroValue}{If 'elimZero' is 'FALSE' this parameter will fill in one zero with the number provided, at each all zero taxa level.  
  		The default is 0.00001. (See Details)}
}
\details{
	This function does the opposite of \code{transformHMPtoHMPTree}.\cr

	Some functions will not work properly in the HMP package if a taxa level exists with all zeros. 'elimZero' and 'zeroValue' 
	allow for an easy way to remove such levels or modify them with a small value in order to use them.  Generally the HMP package requires
	data sets that contain only 1 taxa level and therefore 'trimToTaxaLevel' should also be used with this function.
}

\value{A data set of values in which each column is a separate taxa level and each row is a separate sample.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	### Trims saliva to only contain the class level
	salivaClass <- trimToTaxaLevel(saliva, "class", TRUE)
	
	### This transforms the saliva data set but retains
	### any zero rows that may exist. 
	transSaliva <- transformHMPTreetoHMP(salivaClass, FALSE, 0)
}