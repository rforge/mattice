\name{mrcaOUCH}
\alias{mrcaOUCH}
\title{Find the most recent common ancestor for a vector of taxa}
\description {
  Looks at the intersection of ancestor lines for all taxa in \code{cladeVector}, and returns the most recent as
  the most recent common ancestor.
}
\details{
  As written, a zero-length branch could cause the function to return something other than the m.r.c.a.; this will be
  fixed in subsequent versions.
}
\usage {
  mrcaOUCH(cladeVector, tree)
}
\arguments{
  \item{cladeVector} {
    A vector of tip labels.
    }
  \item{tree} {
    An \code{ouchtree} object.
    }
}
\value{
  A single node label; one of the items in \code{tree@nodes}.
}
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{\link{isMonophyletic}} \code{\link{nodeDescendents}}
}
