\name{createTrees}
\alias{createTrees}
\title{Create a Tree Object}
\description{This function creates a list tree objects of type 'phylo' for use in plotting the trees.}

\usage{createTrees(data, samples = 1, level = "genus")}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{samples}{A numeric vector that is used to select a single (or several) tree(s) out of a data set. A value of 0 will select all trees.}
  \item{level}{The depth the tree creation will go down to (kingdom, phylum, class, order, family, genus, species, subspecies). 'genus' is the default.}
}

\details{
	For 'level' k, p, c, o, f, g, s and ss can be used in place of kingdom, phylum, class, order, family, genus, 
	species and subspecies respectively.
}

\value{A list of 'phylo' objects that can be passed to plotTree to plot them.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	mySalivaTree <- createTrees(saliva, 4)
}