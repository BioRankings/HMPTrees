\name{createTrees}
\alias{createTrees}
\title{Create a Tree Object}
\description{This function creates a list tree objects of type 'phylo' for use in plotting the trees.}

\usage{createTrees(data, samples = NULL, level = "genus", split = ".")}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{samples}{Deprecated. Only send the columns in data to create.}
  \item{level}{The depth the tree creation will go down to (kingdom, phylum, class, order, family, genus, species, subspecies).}
  \item{split}{This is the character that separates the taxa levels in the row names.}
}

\details{
	For 'level' k, p, c, o, f, g, s and ss can be used in place of kingdom, phylum, class, order, family, genus, 
	species and subspecies respectively.
}

\value{A list of 'phylo' objects that can be passed to plotTree to plot them.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	### Creates a tree for the 4th sample in 'Saliva'
	salivaTree <- createTrees(saliva[,4, drop=FALSE])
}