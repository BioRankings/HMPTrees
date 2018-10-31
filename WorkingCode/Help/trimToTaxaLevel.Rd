\name{trimToTaxaLevel}
\alias{trimToTaxaLevel}
\title{Trim a Tree to a Given Level}
\description{
	This function will take a tree and remove all nodes above and/or below the given level
	and adjust the node names accordingly.
}

\usage{
	trimToTaxaLevel(data, level = "genus", trimNodes = "all", 
		trimNames = TRUE, split = ".")
}

\arguments{
  \item{data}{A data frame in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{level}{The depth the tree will go down to (kingdom, phylum, class, order, family, genus, species, subspecies).}
  \item{trimNodes}{
  		If 'all' only the rows at the level specified by 'myTaxaLevel' will be kept.
  		If 'below' all rows below the level specified by 'myTaxaLevel' will be removed.
  		If 'above' all rows above the level specified by 'myTaxaLevel' will be removed.
  }
  \item{trimNames}{
		This only applies when 'trimNodes' is 'all'.
  		If 'TRUE' the row names will be reduced to just the specified level.
  		If 'FALSE' the entire row name will be kept.
  }
  \item{split}{This is the character that separates the taxa levels in the row names.}
}

\details{
	Notes: 
	\enumerate{
		\item For 'level' k, p, c, o, f, g, s, ss, v, and fm can be used in place of kingdom, phylum, 
			class, order, family, genus, species, subspecies, variety, and form respectively.
		\item Numbers greater than 0 can also be used for 'level', with 1 being equivalent to kingdom 
			and no maximum limit.
	}
}
\value{A new data set that has been trimmed to the level selected.}

\examples{
	data(saliva)
	
	### Trims saliva to only contain the class level 
	salivaClass <- trimToTaxaLevel(saliva, "class")
}