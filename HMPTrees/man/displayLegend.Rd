\name{displayLegend}
\alias{displayLegend}
\title{Displays Tree Plot Legend}
\description{This function displays a legend that shows the tree branch sizes/colors divisions.}

\usage{displayLegend(colors, divisions, title = "Confidence Value")}

\arguments{
  \item{colors}{A vector of colors to be used in the plot from lowest ranking to highest ranking.}
  \item{divisions}{A vector of numbers from lowest to highest to separate the tree branches into the color ranking.}
  \item{title}{The title for the legend.  The default is 'Confidence Value'.}
}

\details{
	The values for division should directly relate to the values of your data, i.e. if your data ranges from 0 
	to 50000 reads you should adjust the divisions to fit your data.
}

\value{A blank plot that contains a legend.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	displayLegend(c("red", "orange", "blue"), c(.1, 100, 10000))
}