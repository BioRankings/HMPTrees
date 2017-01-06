\name{compareTwoDataSets}
\alias{compareTwoDataSets}
\title{Likelihood-Ratio-Test Statistics to Compare the Distribution of 2 Sets of RDP-Based Taxonomic Trees}
\description{
	This functions compares the distribution of two sets of RDP-based taxonomic trees using Likelihood-Ratio-Test statistics and a 
	p-value is computed using bootstrapping.
}

\usage{compareTwoDataSets(data1, data2, numBootStraps = 1000, enableMC = FALSE, cores = 3)}

\arguments{
  \item{data1, data2}{Data frames in which each column contains the rdp read counts for every taxa given in the row names.}
  \item{numBootStraps}{The number of times to run the bootstrapping, the default is 1000.}
  \item{enableMC}{When this is 'TRUE' it allows for parallel calculation of the bootstraps. (See Note 1 in Details).}
  \item{cores}{The number of parallel processes to run if enableMC is 'TRUE'.}
}

\details{
	Notes: 
	\enumerate{
		\item Enabling parallel calculation requires the package \code{doParallel}.
		\item Both data sets should be standardized to the same number of reads.
	}
	
	We are interested in assessing whether the distributions from two metagenomic populations are the same or different, which is equivalent to evaluating 
	whether their respective parameters are the same or different. The corresponding hypothesis is given as follows:
	\deqn{H_{\mathrm{o}}: (g_{1}^{*},\tau_{1}) =  (g_{2}^{*},\tau_{2}) = (g_{0}^{*},\tau_{0}) vs H_{\mathrm{A}}: (g_{1}^{*},\tau_{1}) \neq (g_{2}^{*},\tau_{2}) ,}
	where \eqn{(g_{0}^{*},\tau_{0})} is the unknown common parameter vector. To evaluate this hypothesis we use the likelihood-ratio test (LRT) which is given by,
	\deqn{\lambda = -2 \log\left(\frac{L(g_{o}^{*},\tau_{o};{S_{1n},S_{2m}})}{L(g_{1}^{*},\tau_{1};{S_{1n}})+L(g_{2}^{*},\tau_{2};{S_{2m}})} \right),}
	where \eqn{S_{1n}} and \eqn{S_{2m}} are the sets containing \eqn{n} and \eqn{m} random samples of trees from each metagenomic population, respectively. 
	We assume that the model parameters are unknown under both the null and alternative hypothesis, therefore, we estimate these using the MLE procedure proposed 
	in La Rosa et al (see reference 2), and compute the corresponding p-value using non-parametric bootstrap.
}

\value{A p-value for the similarity of the two data sets based on the bootstrapping.}

\author{Patricio S. La Rosa, Elena Deych, Berkley Shands, William D. Shannon}

\examples{
	data(saliva)
	data(stool)
	
	### We use 1 for the number of boot straps for computation time
	### This value should be at least 1000 for an accurate result
	numBootStraps <- 1
	pval <- compareTwoDataSets(saliva, stool, numBootStraps)
	pval
}