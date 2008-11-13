\name{paintBranches}
\alias{paintBranches}
\title{Paint branches for Ornstein-Uhlenbeck analysis}
\description{
Paints branches with regimes changing at nodes specified. Called by 
\code{regimeVectors}.
  }
\usage{
          function(tree, regimeShiftNodes, regimeTitles)
          }
\arguments{
\item{tree}{An S4 ouchtree object, generated using 
\code{ape2ouch}
or
\code{ouchtree}. Old-style ouch trees (data frames) will not work but can be translated using 
\code{ouchtree}.
}
\item{regimeShiftNodes}{Nodes at which the selective regime changes.}
\item{regimeTitles}{Titles for the regimes.} 
}
\details{
\code{regimeTitles} is a vector of titles for the regimes that begin at the root and at the nodes indicated in "regimeShiftNodes",
#                   in order of description in "regimeShiftNodes", except that the root is listed first in "regimeTitles"
#                   but not at all in "regimeShiftNodes"
}
\examples{
}

\keyword{paint}