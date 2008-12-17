\name{runBatchHansen}
\alias{runBatchHansen}
\title{Analyze Ornstein-Uhlenbeck models over a set of trees}
\description{
  \code{runBatchHansen} simulates the evolution of a single character on an \code{ape}-style user tree under arbitrarily
  complex Ornstein-Uhlenbeck models.
}
\usage{
runBatchHansen(ouchTrees, characterStates, cladeMembersList, filePrefix = NULL, di = NULL, nodeNames = NULL, maxNodes = length(cladeMembersList), regimeTitles = NULL, brown = F, ...)
}
\arguments{
  \item{ouchTrees}{
	A list of \code{ouchtree} style trees, or a single tree
  }
  \item{characterStates}{
	A vector of character states to analyze. As written, the function accepts a named vector, with names corresponding to 
	tip names in the trees, a vector of length equal to the number of terminal branches, or a vector of length equal
	to the number of nodes in the trees
  }
  \item{cladeMembersList}{
	A list of tip-name vectors, each defining a node of interest in the tree. Because these vectors are checked for monophyly on each tree, each
	node must by defined by all tips descendent from that node
    }
  \item{filePrefix}{
  	A file prefix for saving hansen objects. If left at default (\code{NULL}), no hansen objects will be saved. Typically 
  	this is fine, but if you want to extract parameters from individual hansen objects, or do any bootstrapping, 
  	it is helpful to save the individual analyses.
    }
  \item{di}{
	A directory in which to save hansen objects. Not required, but strongly recommended if you are saving hansen objects.  
    }
  \item{nodeNames}{
	A vector of names for the nodes being analyzed. Not implemented.
    }
  \item{maxNodes}{
  	The maximum number of nodes at which the regime is allowed to change. Defaults to the number of nodes specified.
    }
  \item{regimeTitles}{
	Names for regimes; only enter if you know how many regimes you will be getting and have reason to name them. Optional.
    }
  \item{brown}{
	Enter \code{T} or \code{F} according to whether you want to run the Brownian motion model.
    }
  \item{...}{
    Arguments passed along to \code{hansen}. As written now, \code{sigma} and \code{alpha} are set at 1
    for the optimizer to begin with. Other arguments are described in documentation for \code{hansen}.
    }
}
\details{
    This function is the primary function for estimating the probability of a change in with changes in regime at the nodes
    specified, with changes corresponding to permutations of the nodes fed to the function. See vignette for additional details.
    
    Note: summary estimates of the probability of change on a branch are only valid if 
    (1) the maximum number of regimes permitted approximates the actual maximum;
    (2) nodes at which changes actually occurred were included among the nodes being tested; and
    (3) any matrix you may have utilized to conduct analysis was balanced, 
        such that all nodes are present in the same number of models.
}
\value{
  A list of class \code{hansenBatch} with the following items:
  \item{hansens} {
  A matrix with rows equal to the number of models analyzed, and four columns: \code{loglik} = log-likelihood of the model;
  \code{dof} = number of free parameters in the model; \code{sigma.squared} = estimate of the variance parameter;
  \code{theta / alpha} = estimate of either the root character state (for the brownian motion model) or the \code{alpha}
  parameter (for Ornstein-Uhlenbeck models).
  }
  thetas = thetas
  regList = ar$regList
  regMatrix = ar$regMatrix
  nodeMatrix = ar$nodeMatrix
  brown = brown
  N = ouchTrees[[i]]@nterm
  analysisDate = date()
  call = match.call())
  }
\author{Andrew Hipp <ahipp@mortonarb.org>}
\seealso{
   \code{\link{batchHansen}} for running analysis on a single tree under an explicit regime,
   \code{\link{summary.hansenBatch}} for summarizing data.
   \code{\link{carex}} for an example
}
\examples{
## See \code{\link{carex}} for example
}
