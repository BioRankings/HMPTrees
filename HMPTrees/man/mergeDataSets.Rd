\name{mergeDataSets}
\alias{mergeDataSets}
\title{Merge Multiple Data Sets}
\description{This function can take any number of data sets, calculate their individual and combined MLEs and then merge them.}

\usage{mergeDataSets(dataList, calcMLE = FALSE, uniqueNames = FALSE, data = NULL)}

\arguments{
  \item{dataList}{A list of data frames in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{calcMLE}{If 'FALSE' the MLEs for the data sets will not be calculated, otherwise they are added to the end.}
  \item{uniqueNames}{If 'TRUE' the column names in the combined data set will be appended to insure uniqueness, otherwise the column names
  		will follow the naming process from the \code{merge} function.}
  \item{data}{Deprecated. Replaced with dataList.}
}

\details{
	Although not required, all data sets should be standardized to the same number of reads before merging.
}

\value{A single data set containing all the data from the input data sets, in addition to their individual MLEs and a combined MLE if requested.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	data(stool)
	
	dataComb <- mergeDataSets(list(saliva, stool))
}