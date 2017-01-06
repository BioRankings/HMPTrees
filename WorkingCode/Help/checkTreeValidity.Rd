\name{checkTreeValidity}
\alias{checkTreeValidity}
\title{Check Validity of an RDP-Based Taxonomic Tree}
\description{
	This function goes through every node in the tree and for each node it checks that the sum of that nodes children 
	are less than or equal to the value of that node.
}

\usage{checkTreeValidity(data, samples = 1, epsilon = 0.0001, split = ".")}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{samples}{A numeric vector that is used to select a single (or several) tree(s) out of a data set. A value of 0 will select all trees.}
  \item{epsilon}{This value allows for rounding problems or other such small errors in the data such that the (parent + epsilon > sum(children)).}
  \item{split}{This is the character that separates the taxa levels in the row names.}
}

\value{A boolean vector that indicates the validity of every tree tested.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva) 
	
	validTree <- checkTreeValidity(saliva, 1)
	validTree
}