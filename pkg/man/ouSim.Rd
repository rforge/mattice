\name{ouSim}
\alias{ouSim}
\title{Discrete-time simulation of Ornstein-Uhlenbeck models on a user tree.}
\description{
  \code{ouSim} simulates the evolution of a single character on an \code{ape}-style user tree under arbitrarily
  complex Ornstein-Uhlenbeck models.
}
\usage{
(phy, rootState = 0, shiftBranches = NULL, shiftStates = NULL, alpha = 0, variance = 1, theta = rootState, model = "OU", branchMeans = NULL, steps = 1000)
}
\arguments{
  \item{phy}{
    An object of class \code{phylo}.
  }
  \item{rootState}{
    The character state at the root of the tree
  }
  \item{shiftBranches}{
    An optional vector indicating any branches at which an OU model has a determined shift in ancestral state. 
    Same order as \code{shiftStates}. This argument and \code{shiftStates} are only needed if you want to specify a model with
    a disjunction in phenotype, similar to restarting evolution at a new point.
    }
  \item{shiftStates}{
    An optional vector of length = \code{length(shiftStates)} indicating the ancestrat states for the branches at which the state shifts.
    }
  \item{alpha}{
    The strength of attraction toward an equilibrium or optimum. This term is refered to as the rate of evolution by Hansen (1997)
    and the strength of selection by Butler and King (2004). It is a multiplier by the difference between the character state and the
    character state optimum. Alpha can be submitted as a single value applied to all branches or as a vector corresponding to branches
    in the \code{phylo} object. At \code{alpha = 0}, the simulation approximates a Brownian motion process.
    }
  \item{variance}{
    The variance on the stochastic portion of the Ornstein-Uhlenbeck model.
    }
  \item{theta}{
    The character state optimum. \code{theta} can be submitted as a single value or, like \code{alpha}, as a vector corresponding to
    branches of the tree.
    }
  \item{branchMeans}
    The mean for each branch, utilized only in the Platt model. 
    }
  \item{steps}
    The number of slices into which the tree is divided for simulation.
    }
}
\details{
    ouSim as currently written is a very basic simulation engine, but also flexible. As written, the user has to specify the 
    model using two vectors that correspond to the branches in an \code{ape} format tree (a \code{phylo} object). The next
    iteration of this function will take an ouchtree as input; interface with the \code{regimeVectors} function for painting
    regimes; accept a \code{hansentree} or \code{ouchtree} object as input, using the actual regimes and estimated paramters 
    to run the simulation; and return results that are scaled correctly for the tree (I believe at this point that the 
    sigma and alpha parameters don't return results that match the distributions you get using \code{simulate} on a \code{hansentree}
    or \code{browntree} object, which you should.
}
\value{
  A list of class 'ouSim' that describes the phenotype at the beginning and end of each branch segment, as well as the model.
  }
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{plot.ouSim} for visualizing simulation
}
\examples{
XXX
}
\keyword{models}
