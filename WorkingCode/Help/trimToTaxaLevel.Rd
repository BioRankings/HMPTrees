\name{trimToTaxaLevel}
\alias{trimToTaxaLevel}
\title{Trim a Tree to a Given Level}
\description{
	This function will take a tree and either remove all nodes lower than the given 
	level or will remove all nodes not of the given level.
}

\usage{
	trimToTaxaLevel(data, level = "genus", eliminateParentNodes = FALSE, 
		trimBelow = NULL, split = ".")
}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{level}{The depth the tree will go down to (kingdom, phylum, class, order, family, genus, species, subspecies). Genus is the default.}
  \item{eliminateParentNodes}{
  		If 'TRUE' the data set returned will only contain rows at the level specified by 'myTaxaLevel'.
  		If 'FALSE' the data set returned will contain all the nodes up to the level specified by 'myTaxaLevel'.}
  \item{trimBelow}{
  		If 'NULL' the function will pull out only the data at the level specified by 'myTaxaLevel'.
  		If 'TRUE' the function will remove all the levels below the specified level.
  		If 'FALSE' the function will remove all the levels above the specified level.
  }
  \item{split}{This is the character that separates the taxa levels in the row names.}
}

\details{
	Notes: 
	\enumerate{
		\item For 'level' k, p, c, o, f, g, s and ss can be used in place of kingdom, phylum, class, order, family, genus, species 
			and subspecies respectively.
		\item Numbers can also be used for 'level', with no maximum limit.
		\item The option to 'eliminateParentNodes' only works when 'trimBelow' is NULL.
	}
}
\value{A new data set that has been trimmed to the level selected.}

\author{Berkley Shands, Patricio S. La Rosa, Elena Deych, William D. Shannon}

\examples{
	data(saliva)
	
	### Trims saliva to only contain the class level 
	salivaClass <- trimToTaxaLevel(saliva, "class", TRUE)
}