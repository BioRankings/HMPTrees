\name{formatData}
\alias{formatData}
\title{Formats a Data Set}
\description{This function will take a data set and format it by removing low count trees, and/or normalizing counts.}

\usage{formatData(data, countThreshold = 1000, normalizeThreshold = 10000)}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{countThreshold}{A cut off threshold for reads - all trees with fewer than this number of reads will be removed.}
  \item{normalizeThreshold}{All the trees that are not removed will be normalized to this many reads.}
}

\details{
	When removing trees with too few reads, the cuts off is based on the value of the top level node, not the sum of all the reads in a sample.
}

\value{A new data set that is trimmed and standardized based on the specified parameters. The new data is also reordered alphabetically according to row labels.}

\author{Patricio S. La Rosa, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	data(throat)
	
	throat <- formatData(throat, 1000, 10000)
}