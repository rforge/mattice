\name{plot.ouSim}
\alias{plot.ouSim}
\title{Visualizing discrete-time simulation of Ornstein-Uhlenbeck models.}
\description{
  \code{plot.ouSim} is a \code{plot} method for an \code{ouSim} object.
}
\usage{
  plot.ouSim(x, nodeColor = "blue", nodeDotSize = 1.4, colors = NULL, ...)
}
\arguments{
  \item{x}{
    An \code{ouSim} object.
    }
  \item{nodeColor}{
    A single color for the dot indicating nodes, or a vector of colors corrsponding to node numbers.
    }
  \item{nodeDotSize}{
    Size of the dot at the nodes.
    }
  \item{colors}{
    A vector correponding to branches. If not provided, color is extracted from the regime in a \code{hansentree} object.
    This vector can be the same format as a regime, so you can color branches with a call to \code{\link{paintBranches}}.
    }
  \item{...}{
    Additional arguments passed along to \code{lines} (e.g., \code{lwd} for line width, \code{lty} for line type).
    }
}
\value{
  None. Called for its side-effect (plotting).
}
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{\link{ouSim}} for simulations; \code{carex} for examples
}
