\name{regimeMatrix}
\Rdversion{1.1}
\alias{regimeMatrix}
\title{
  Create a matrix that defines models for batch analysis
}
\description{
  \code{regimeMatrix} is a recursive function that returns a matrix of n columns by k rows, where n is a number
  of items and k is the number of all possible permutations of those items, up to a maximum number of items
  defined by \code{maxNodes}, and including a row with zero items. 
  }
\usage{
  regimeMatrix(n, maxNodes)
  }
\arguments{
  \item{n}{
    An integer between 1 and \code{maxNodes}
    }
  \item{maxNodes}{
    An integer between 1 and n; at \code{n} = \code{maxNodes}, \code{regimeMatrix} will return a matrix of 2^\code{n} rows.
    }
  }
\details{
  Typically, this function is called by \code{regimeVectors},
  but it may be useful to call this function if you are interested in defining an arbitrary set of models
  that doesn't include all possible models in the 
  }
\value{
  A matrix of n columns by k rows, where n is a number
  of items and k is the number of all possible permutations of those items, up to a maximum number of items
  defined by \code{maxNodes}, and including a row with zero items
  }
\author{
  Andrew Hipp \email{ahipp@mortonarb.org}
  }
\note{
  Use of \code{regimeMatrix} or \code{regimeMaker} to generate tailored sets of regimes is not currently integrated smoothly with \code{runBatchHansen},
  but for the time being they can be used to create regimes for manual analysis.
  }
\seealso{
  \code{\link{regimeMaker}}, \code{\link{paintBranches}}
  }
\examples{
  # allowing all possible permutations of changes at 15 nodes will be quite time-consuming:
  print(paste("Number of models =", dim(regimeMatrix(15,15))[1]))
  # paring that down to all possible permutations of changes at 15 nodes, 
  # up to a maximum of 5 changes, is much more manageable:
  print(paste("Number of models =", dim(regimeMatrix(15,5))[1]))
  }
\keyword{models}