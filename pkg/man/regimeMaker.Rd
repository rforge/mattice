\name{regimeMaker}
\Rdversion{1.1}
\alias{regimeMaker}
\title{
  Paint ouch-style regimes over a set of trees based on a user-defined matrix
}
\description{
  This function utilizes a matrix of the format returned by \code{regimeMatrix} to define the models it is going to create. Then, it paints these models onto the trees provided and returns a list suitable for further analysis.
}
\usage{
regimeMaker(ouchTrees, regMatrix, nodeMembers)
}
\arguments{
  \item{ouchTrees}{
  A list of \pkg{ouch}-style trees
  }
  \item{regMatrix}{
  A model by nodes matrix of the format returned by \code{regimeMatrix}
  }
  \item{nodeMembers}{
  A list of individuals defining the nodes specified in \code{regMatrix}
  }
}
\details{
  The easiest way to use this function will usually be to utilize \code{regimeMatrix} to create the matrix input as \code{regMatrix} and manually prune or add rows. See the example below.
}
\value{
  A list with three items:
  \item{regList}{
    A list of lists of regimes, with the number of first-order items corresponding to the number of trees in \code{ouchTrees}, 
    the number of second-order items corresponding to the number of models being tested. 
    In the example below, \code{regs$regList[[2]][[4]]} is the fourth model defined on the second tree.
    }
  \item{nodeMatrix}{
    A boolean matrix with rows corresponding to the trees in \code{ouchTrees} and columns corresponding to the nodes specified in \code{regMatrix}; 
    cells indicate by \code{TRUE} or \code{FALSE} whether a node is present in each of the trees being analyzed.
    }
  \item{regMatrix}{
    A list of \code{regimeMatrix}-format matrices that define the models applicable for each tree in \code{ouchTrees}, with missing nodes designated by 0s and 
    duplicate rows blanked out with NA.
    }
}
\author{
  Andrew L. Hipp <ahipp@mortonarb.org>
}
\note{
  Use of \code{regimeMatrix} or \code{regimeMaker} to generate tailored sets of regimes is not currently integrated smoothly with \code{runBatchHansen},
  but for the time being they can be used to create regimes for manual analysis.
}
\seealso{
  \code{\link{regimeMatrix}}, \code{\link{paintBranches}}
}
\examples{
  data(carex)
  # tree comes in in phylo format, but we need an ouchtree object
  trees <- lapply(carex$ovales.bayesTrees[1:10], ape2ouch) 
  # if we are interested in comparing the single-change models that allow changes only at nodes 1, 2, 4, 5, and 8, as well as the no-change model, we could do:
  mat <- regimeMatrix(n = 8, maxNodes = 1)[c(1,2,4,5,8,9), ] # the 9th row is the no-change model
  #... then we make regimes based on that regime matrix...
  regs <- regimeMaker(trees, mat, carex$ovales.nodes)
  #... and paint the regimes from tree 1 onto tree 1 either singly...
  plot(trees[[1]], regime = regs$regList[[1]][[1]])
  #... or in a layout:
  plot(trees[[1]], regime = regs$regList[[1]])
}
\keyword{ models }