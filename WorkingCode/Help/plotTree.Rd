\name{plotTree}
\alias{plotTree}
\title{Plots a Tree Object}
\description{This function takes one or more 'phylo' objects and plots them.}

\usage{
	plotTree(trees, colors, divisions, main, sub, showTipLabel = TRUE, 
		showNodeLabel = FALSE, displayLegend = TRUE)
}

\arguments{
  \item{trees}{A list that contains at least one tree of type 'phylo'.}
  \item{colors}{A vector of colors to be applied to the branches in the plot.}
  \item{divisions}{A vector of numbers to be used as break points to assign different colors.}
  \item{main}{A custom title(s) for the plot(s).  The default value is the column names.}
  \item{sub}{A custom subtitle for the plot.}
  \item{showTipLabel}{Hides the tip labels if 'FALSE' otherwise it shows all non-zero tip labels.}
  \item{showNodeLabel}{Hides the interior node labels if 'FALSE' otherwise it shows all non-zero node labels.}
  \item{displayLegend}{Enables the display of a legend of the branch colors and divisions when 'TRUE'.}
}

\details{
	Notes: 
	\enumerate{
		\item The \code{phylo} type is a product of the \code{ape} package and the \code{createTrees} function in this package 
			produces a list of \code{phylo} type objects for use with this function.  
		\item The values for division should directly relate to the values of your data, i.e. if your data ranges from 
			0 to 50000 reads you should adjust the divisions to fit your data.
	}
}

\value{A plot of the tree(s).}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	### Creates a object of type 'phylo' for the 4th tree in 'Saliva'
	### Then plots it
	salivaTree <- createTrees(saliva, 4)
	plotTree(salivaTree, displayLegend=FALSE)
}