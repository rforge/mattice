\name{paintBranches}
\alias{paintBranches}
\title{Paint branches for Ornstein-Uhlenbeck analysis}
\description{
  Paints branches with regimes changing at nodes specified. Called by 
  \code{regimeVectors}.
}
\usage{
  paintBranches(tree, regimeShiftNodes, regimeTitles)
}
\arguments{
  \item{tree}{
   An S4 ouchtree object, generated using \code{ape2ouch} or \code{ouchtree}. 
   Old-style ouch trees (data frames) will not work but can be translated using \code{ouchtree}.
   }
  \item{regimeShiftNodes}{
   Either a vector of nodes (by number) at which the selective regime changes, 
   or a list of taxa defining nodes. This allows you to use paintBranches with the same
   nodes you use to call runBatchHansen. This parameter may include the root but does
   not need to.
   }
  \item{regimeTitles}{
   Titles for the regimes.
   } 
}
\details{
  \code{regimeTitles} is a vector of titles for the regimes that begin at the root and at the nodes indicated in 
  \code{regimeShiftNodes} in order of description in "regimeShiftNodes", except that the root is listed first in 
  \code{regimeTitles} but may not be included at all in \code{regimeShiftNodes}.
}
\value{
  A regime in which shifts in optimum are coded as occuring at all nodes specified in \code{regimeShiftNodes}.
  The regime is formated to input directly into \code{hansen}
}
\seealso{
  \code{\link{regimeVectors}}, \code{\link{isMonophyletic}}, \code{\link{nodeDescendents}}, \code{\link{mrcaOUCH}}
}
