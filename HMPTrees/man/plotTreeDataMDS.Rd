\name{plotTreeDataMDS}
\alias{plotTreeDataMDS}
\title{Plot an MDS Plot of a Group of Trees}
\description{This function can take any number of data sets and plots them on an MDS plot to show relative closeness to one another.}
\usage{
	plotTreeDataMDS(dataList, main = "Tree MDS Comparisons", calcMLE = TRUE, 
		mleTitles = NULL, dotColors = NULL, dotSizes = NULL, showNames = FALSE, 
		returnCoords = FALSE, data = NULL)
}

\arguments{
  \item{dataList}{A list of a data frames in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{main}{A title for the MDS plot.}
  \item{calcMLE}{If 'FALSE' the MLEs for the data sets will not be calculated and plotted.}
  \item{mleTitles}{Deprecated. Replaced with the names in 'dataList'.}
  \item{dotColors}{The colors to be used when plotting the points and MLE points on the MDS plot.}
  \item{dotSizes}{A vector in which the first value is the data points CEX and the second value is the MLEs CEX.}
  \item{showNames}{When 'TRUE' the column name will be plotted above each corresponding point.}
  \item{returnCoords}{When 'TRUE' this function will return the x and y coordinates for every plotted point.}
  \item{data}{Deprecated. Replaced with dataList.}
}

\value{A MDS plot of the data.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	data(stool)
	
	plotTreeDataMDS(list(Saliva=saliva, Stool=stool))
}