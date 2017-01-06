\name{throat}
\alias{throat}
\docType{data}
\title{Throat Data Set}
\description{A data set containing all taxa from 22 subjects.}

\usage{data(throat)}

\format{
	The format is a data frame of 529 rows by 22 columns, with each column being a separate subject and each 
	row being a different taxa denoted by the row names.  The taxanomical levels are separated by a '.' in 
	their names (Bacteria.Phylum.Class....).  The values in each column are the sum of values that each taxa
	had in an RDP file. It should also be noted that the samples have not been normalized and should be 
	used with 'formatData'.  Also any level that ends with a U was unclassified in the RDP file.
}

\keyword{datasets}