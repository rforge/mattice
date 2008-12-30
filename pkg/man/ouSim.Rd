\name{ouSim}
\alias{ouSim}
\alias{ouSim.ouchtree}
\alias{ouSim.brownHansen}
\alias{ouSim.phylo}
\alias{ouSim.hansenBatch}
\title{Discrete-time simulation of Ornstein-Uhlenbeck models on a user tree.}
\description{
  \code{ouSim} simulates the evolution of a single character for visualization purposes; for parametric bootstrapping,
  utilize the \code{simulate} methods in \code{ouch}.
}
\usage{
  ouSim(object, ...) ## Generic function
  ouSim.phylo(phy, rootState = 0, shiftBranches = NULL, shiftStates = NULL, alpha = 0, variance = 1, theta = rootState, model = "OU", branchMeans = NULL, steps = 1000)
  ouSim.ouchtree(tree, rootState = 0, alpha = 0, variance = 1, theta = rootState, steps = 1000)
  ouSim.brownHansen(analysis, ...)
  ouSim.hansenBatch(analysis, tree, treeNum = 1, rootState = NULL, ...) 
}
\details{
  A call to \code{ouSim} detects the class of \code{object} and behaves as follows:
  \item{phylo}{
    With an \pkg{ape}-style tree, all parameters of the model are specified by branch. Most flexibility, least convenience.
    }
  \item{ouchtree}{
    With an \pkg{ouch}-style (S4) tree, again, all parameters of the model are specified by branch. 
    Somewhat less flexibility than \code{ouSim.phylo} but all the pain. Typically just called by the next two methods.
    }
  \item{browntree}{
    All Brownian motion model parameters and tree are extracted from the \code{browntree} object.
    }
  \item{hansentree}{
    All Ornstein-Uhlenbeck model parameters and tree are extracted from the \code{hansentree} object.
    Only the first regime is utilized.
    }
  \item{hansenBatch}{
    Model-averaged parameters from the \code{hansenBatch} object are used for analysis. 
    One of the trees used for analysis must be provided, and a corrsponding tree number must be provided so that 
    branches are indexed correctly.
    }
  \item{hansenSummary}{
    Behaves exactly like a \code{hansenBatch} object.
    }
  \code{ouSim.phylo} is a very basic simulation engine, but also the most flexible. As written, the user has to specify the 
  model using two vectors that correspond to the branches in an \code{ape}-format tree. It is important to note that
  this simulation method is really  heuristic device, not appropriate for estimating parameter distributions. For
  analysis purposes, you should utilize the \code{simulate} and \code{bootstrap} methods in \code{ouch}.
}
\arguments{
  \item{object}{
    In a call to the generic function, an object of class ouchtree, browntree, hansentree, hansenBatch, or hansenSummary.
    }
  \item{phy}{
    An \code{ape}-style tree, class \code{phylo}.
    }
  \item{rootState}{
    The character state at the root of the tree. In a {browntree} object, this value is provided. In a {hansentree} or
    {batchHansen} object, the value at the root is not provided, but it is taken to be the equilibrium or optimum 
    (\code{theta}) at the root of the tree.
    }
  \item{shiftBranches}{
    For a \code{phylo} tree only. An optional vector indicating any branches at which an OU model has a determined shift in ancestral state. 
    Same order as \code{shiftStates}. This argument and \code{shiftStates} are only needed if you want to specify a model with
    a disjunction in phenotype, similar to restarting evolution at a new point.
    }
  \item{shiftStates}{
    For a \code{phylo} tree only. An optional vector of length = \code{length(shiftStates)} indicating 
    the ancestral states for the branches at which the state shifts.
    }
  \item{alpha}{
    The rate of evolution toward an equilibrium or optimum. This term is refered to as the rate of evolution by Hansen (1997)
    and the strength of selection by Butler and King (2004). It is a multiplier by the difference between the character state and the
    character state optimum. Alpha can be submitted as a single value applied to all branches or as a vector corresponding to branches
    in the \code{phylo} object. At \code{alpha = 0}, the simulation approximates a Brownian motion process.
    This parameter is taken from the analysis results for \code{browntree} (\code{alpha} = 0), \code{hansentree} (point estimate), 
    \code{hansenBatch} or \code{hansenSummary} (model-averaged).
    }
  \item{variance}{
    The variance on the stochastic portion of the Ornstein-Uhlenbeck model.
    This parameter is taken from the analysis results for \code{browntree} (\code{alpha} = 0), \code{hansentree} (point estimate), 
    \code{hansenBatch} or \code{hansenSummary} (model-averaged).
    }
  \item{theta}{
    The character state optimum. \code{theta} can be submitted as a single value or, like \code{alpha}, 
    as a vector corresponding to branches of the tree.
    This parameter is taken from the analysis results for \code{browntree} (\code{alpha} = 0), \code{hansentree} (point estimate), 
    \code{hansenBatch} or \code{hansenSummary} (model-averaged).
    }
  \item{model}{
    For a \code{phylo} tree only. Specify "OU" for Ornstein-Uhlenbeck and Brownian motion models, "meanVar" for a model
    in which the variance is constant across the tree, but mean varies by branch, and distribution in each generation 
    depends on only these parameters (ancestry is not considered, only current mean and variance).
    }
  \item{branchMeans}{
    For a \code{phylo} tree only. The mean for each branch, utilized only in the "meanVar" model. 
    }
  \item{steps}{
    The number of slices into which the tree is divided for simulation.
    }
}
\value{
  A list of class 'ouSim' that describes the phenotype at the beginning and end of each branch segment, as well as the model.
  }
\author{Andrew Hipp ahipp@mortonarb.org}
\seealso{
   \code{plot.ouSim} for visualizing simulation; \code{carex} for examples
}
\references{
  Hansen, T.F. (1997)
  Stabilizing selection and the comparative analysis of adaptation.
  Evolution 51:1341-1351.

  Butler, M. and A.A. King. (2004) 
  Phylogenetic comparative analysis: a modeling approach for adaptive evolution. 
  American Naturalist 164:683-695.
}
