\name{isMonophyletic}
\alias{isMonophyletic}
\title{Test whether a set of taxa is monophyletic on an ouch-style tree.}
\usage{
  isMonophyletic(tree, taxa)
}
\arguments{
  \item{tree} {
    An \code{ouchtree} object.
    }
  \item{taxa} {
    A vector of taxon names.
    }
}
\value{
    Boolean, \code{TRUE} if taxa form a monophyletic group, \code{FALSE} if not.
    }
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{\link{nodeDescendents}} \code{\link{mrcaOUCH}}
}
/keyword{tree}