\name{stool}
\alias{stool}
\docType{data}
\title{Stool Data Set}
\description{A data set containing all taxa from 24 subjects.}

\usage{data(stool)}

\format{
	The format is a data frame of 371 rows by 24 columns, with each column being a separate subject and each 
	row being a different taxa denoted by the row names.  The taxanomical levels are separated by a '.' in 
	their names (Bacteria.Phylum.Class....).  The values in each column are the sum of values that each taxa
	had in an RDP file.  It should also be noted that the samples are normalized to 7000 reads and any level 
	that ends with a U was unclassified in the RDP file.
}
	
\keyword{datasets}