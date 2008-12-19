\name{nodeDescendents}
\alias{nodeDescendents}
\title{Find all descendents of a node}
\description {
  Recursive function to find all the descendents of a node on an \code{ouchtree} object; returns only terminals, not
  internal nodes descendent from \code{startNode}.
}
\usage {
  nodeDescendents(tree, startNode) 
}
\arguments{
  \item{tree} {
    An \code{ouchtree} object.
    }
  \item{startNode} {
    The node of interest; nomenclature follows \code{tree@nodes}.
    }
}
\value{
  A vector of tip labels.
}
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{\link{isMonophyletic}} \code{\link{mrcaOUCH}}
}
/keyword{tree}