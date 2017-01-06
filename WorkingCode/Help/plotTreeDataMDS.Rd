\name{plotTreeDataMDS}
\alias{plotTreeDataMDS}
\title{Plot an MDS Plot of a Group of Trees}
\description{This function can take any number of data sets and plots them on an MDS plot to show relative closeness to one another.}
\usage{
	plotTreeDataMDS(data, main = "Tree MDS Comparisons", calcMLE = TRUE, mleTitles, 
		dotColors = c("red", "orange", "blue", "green", "yellow", "purple"), 
		dotSizes = c(1, 2), showNames = FALSE, returnCoords = FALSE)
}

\arguments{
  \item{data}{A single data frame or a list of a data frames in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{main}{A title for the MDS plot.}
  \item{calcMLE}{If 'FALSE' the MLEs for the data sets will not be calculated and plotted.}
  \item{mleTitles}{A vector of labels for each MLE data point on the MDS plot.}
  \item{dotColors}{The colors to be used when plotting the points and/or MLE points on the MDS plot. (See Notes in Details)}
  \item{dotSizes}{A vector in which the first value is the data points CEX and the second value is the MLEs CEX.}
  \item{showNames}{When 'TRUE' the column name will be plotted above each corresponding point.}
  \item{returnCoords}{When 'TRUE' this function will return the x and y coordinates for every plotted point.}
}
\details{
	Notes: 
	\enumerate{
		\item If the length of dot colors is at least twice the length of data, the colors will correspond in the following manner:
			(data 1's points, data 1's MLE, data 2's points, data 2's MLE, etc)
		\item If the length of dot colors is less than twice the length of data, the colors will correspond in the following manner:
			(data 1's points & data 1's MLE, data 2's points & data 2's MLE, etc)
		\item If the length of dot colors is less than the length of data, dotColors will be set by the \code{rainbow} function and 
			the colors will correspond as above
	}
}

\value{A MDS plot of the data.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	data(stool)
	
	plotTreeDataMDS(list(saliva, stool), mleTitles=c("Saliva", "Stool"))
}