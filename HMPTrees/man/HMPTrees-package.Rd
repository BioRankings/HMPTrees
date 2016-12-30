\name{HMPTrees-package}
\alias{HMPTrees-package}
\alias{HMPTrees}
\docType{package}
\title{Object Oriented Data Analysis of Taxonomic Trees}

\description{
	Object Oriented Data Analysis (OODA) methods to analyze Human Microbiome taxonomic trees directly. 
	We provide tools to model, compare, and visualize populations of taxonomic trees.
}

\details{
	HMP metagenomic sequences in a sample can be represented as a rooted taxonomic tree. Using supervised taxonomic methods a sequence is matched to a 
	hierarchical taxa or taxonomy bins defined in a bacterial-taxonomy library such as, for example, the Ribosomal Database Project (RDP) (Cole, 2005). The 
	supervised taxonomic analysis allows us to represent each sample (set of sequences) by a rooted taxonomic tree where the root corresponds to taxon at the 
	Kingdom level, i.e., bacteria, and the leaves correspond to the taxa at the Genus level, and the width of the edges (paths) between taxonomic levels correspond 
	to the 'abundances' of the descending taxon. \cr
	
	In particular, we combine RDP matches by adding RDP values of common taxon, which allows us to provide a measure 
	of taxa abundance weighting on the confidence of each taxa assignment. The resulting taxonomic trees satisfy the following conditions: i) branches closer to 
	the root have higher 'abundance' values than branches closer to leaves, and ii) the sum of the 'abundances' of all descending taxa under a common parent taxon 
	cannot be larger than the 'abundance' of the corresponding parent taxon. \cr
	
	It is important to note that due to how the \code{ape} package works the following naming conventions apply to taxa names: 
	\enumerate{
		\item Colons cannot be used in the taxa names at all.
		\item Each taxa name must be unique - you cannot have two seperate branches both have a child named 'unclassified' for example .
			(We took the parent name and added a 'U' to the end to signify an unclassified in our data sets)
		\item There can only be one top level node. (Bacteria and Archaea cannot both exist unless there is an additional single level above them for example)
	}
}

\references{
	\enumerate{
		\item Cole JR, Chai B, Farris RJ, Wang Q, Kulam SA, McGarrell DM, Garrity GM, Tiedje JM. The Ribosomal Database Project (RDP-II): 
			sequences and tools for high-throughput rRNA analysis. Nucleic Acids Research 2005; 33: D294-D296.
	
		\item P. S. La Rosa, Berkley Shands, Elena Deych, Yanjiao Zhou, Erica Sodergren, George Weinstock, and William D. Shannon, 
			"Object data analysis of taxonomic trees from human microbiome data,"PLoS ONE 7(11): e48996. doi:10.1371/journal.pone.0048996. Nov. 2012.
			
		\item Banks D, Constantine GM. Metric Models for Random Graphs. Journal of Classification 1998; 15: 199-223.
		
		\item Shannon WD, Banks D. Combining classification trees using MLE. Stat Med 1999; 18: 727-740.
	}
}

\author{
	Patricio S. La Rosa, Elena Deych, Berkley Shands, William D. Shannon
}

\keyword{ package }