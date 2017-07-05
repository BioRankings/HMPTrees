\name{createAndPlot}
\alias{createAndPlot}
\title{Create and Plot a Tree from a Data Set}
\description{This function combines the createTrees and plotTree functions to create and plot a set of trees.}

\usage{
	createAndPlot(data, samples = NULL, level = "genus", colors = NULL, 
		divisions = NULL, main = NULL, sub = "", showTipLabel = TRUE, 
		showNodeLabel = FALSE, displayLegend = TRUE, onePerPage = FALSE, 
		split = ".")
}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{samples}{Deprecated. Only send the columns in data to plot.}
  \item{level}{The depth the tree creation will go down to (kingdom, phylum, class, order, family, genus, species, subspecies).}
  \item{colors}{A vector of colors to be applied to the branches in the plot.}
  \item{divisions}{A vector of numbers to be used as break points to assign different colors.}
  \item{main}{A custom title(s) for the plot(s).}
  \item{sub}{A custom subtitle for the plot.}
  \item{showTipLabel}{Hides the tip labels if 'FALSE' otherwise it shows all non-zero tip labels.}
  \item{showNodeLabel}{Hides the interior node labels if 'FALSE' otherwise it shows all non-zero node labels.}
  \item{displayLegend}{Enables the display of a legend of the branch colors and divisions when 'TRUE'.}
  \item{onePerPage}{If 'TRUE' one tree will be plotted per page, if 'FALSE' four will be displayed per page.}
  \item{split}{This is the character that separates the taxa levels in the row names.}
}

\details{
	Notes:
	\enumerate{
		\item For 'level' k, p, c, o, f, g, s and ss can be used in place of kingdom, phylum, class, order, family, genus, species 
			and subspecies respectively.
		\item The values for division should directly relate to the values of your data, i.e. if your data ranges from 0 to 50000 reads you 
			should adjust the divisions to fit your data.
	}
}

\value{A plot of the tree(s).}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	### Plots the trees in column 2 and 3 in 'Saliva'
	createAndPlot(saliva[,2:3])
}