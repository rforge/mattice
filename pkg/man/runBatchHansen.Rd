\name{runBatchHansen}
\alias{runBatchHansen}
\title{Analyze Ornstein-Uhlenbeck models over a set of trees, with changes in regime corresponding to permutations
       of the nodes fed to it.}
\description{
  \code{runBatchHansen} simulates the evolution of a single character on an \code{ape}-style user tree under arbitrarily
  complex Ornstein-Uhlenbeck models.
}
\usage{
runBatchHansen(ouchTrees, characterStates, cladeMembersList, nodeNames = NULL, maxNodes = NULL, regimeTitles = NULL, brown = F, rescale = 1, ...)}
\arguments{
  \item{ouchTrees}{
	A list of \code{ouchtree} style trees, or a single tree
  }
  \item{characterStates}{
	A vector of character states to analyze. As written, the function accepts a named vector, with names corresponding to 
	tip names in the trees, a vector of length equal to the number of terminal branches, or a vector of length equal
	to the number of nodes in the trees. 
  }
  \item{cladeMembersList}{
	A list of tip-name vectors, each defining a node of interest in the tree. Because these vectors are checked for monophyly on each tree, each
	node must by defined by all tips descendent from that node.
    }
  \item{nodeNames}{
	A vector of names for the nodes being analyzed. Not implemented. Optional.
    }
  \item{maxNodes}{
  	The maximum number of nodes at which the regime is allowed to change. Optional.
    }
  \item{regimeTitles}{
	Names for regimes; only enter if you know how many regimes you will be getting and have reason to name them. Optional.
    }
  \item{brown}{
	Enter \code{T} or \code{F} according to whether you want to run the Brownian motion model.
    }
  \item{rescale}
	Enter a value > 0 to rescale all trees to length = \code{rescale}.
    }
  \item{...}
    Arguments passed along to \code{hansen}. As written now, you have to include a beginning value for \code{sigma}
    and \code{alpha} for the optimizer to begin with. Other arguments are described in documentation for \code{hansen}.
    To begin, it is easiest to set \code{sigma} = 1, \code{alpha} = 1.
    }
}
\details{
    This function is the primary function for estimating the 
}
\value{
  XXX
  }
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{batchHansen} for running analysis on a single tree under an explicit regime,
   \code{summary.hansenBatch} for summarizing data.
}
\examples{
XXX
}
\keyword{models}
