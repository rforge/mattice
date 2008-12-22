\name{regimeVectors}
\alias{regimeVectors}
\title{Paint all possible regimes}
\description{
  Creates all regimes possible for the nodes specified on the trees provided, returning for each tree only those regimes 
  that are present on that tree, and 
}
\usage{
  regimeVectors(ouchTrees, cladeMembersList, maxNodes = NULL)
}
\arguments{
\item{ouchTrees}{
  A list of S4 \code{ouchtree} objects
  }
\item{cladeMembersList}{
  A list of taxa defining nodes. 
  }
\item{maxNodes}{
  The maximum number of nodes at which a regime may change.
  } 
}
\details{
  \code{regimeVectors} is the basic call to get all regime permutations for analysis in batchHansen.
}
\value{
A list with three items:
\item{regList} A list (by tree) of lists of regimes that can be input directly into \code{hansen}.
\item{nodeMatrix} A matrix indicating which nodes are presentin each tree.
\item{regMatrix} A list (by tree) of matrices indicating which nodes are present in each model.
}
\seealso{
\code{\link{paintBranches}}
}