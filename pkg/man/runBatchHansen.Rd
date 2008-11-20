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

  }
  \item{characterStates}{

  }
  \item{cladeMembersList}{

    }
  \item{nodeNames}{

    }
  \item{maxNodes}{
    }
  \item{regimeTitles}{

    }
  \item{brown}{

    }
  \item{rescale}

    }
  \item{...}
    Arguments passed along to \code{hansen}. As written now, you have to include a beginning value for \code{sigma}
    and \code{alpha} for the optimizer to begin with. Other arguments are described in documentation for \code{hansen}.
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
   \code{batchHansen} running analysis on a single tree, \code{summary.hansenBatch} for summarizing data.
}
\examples{
XXX
}
\keyword{models}
