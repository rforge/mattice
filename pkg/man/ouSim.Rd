\name{ouSim}
\alias{ouSim}
\title{Discrete-time simulation and visualization of Ornstein-Uhlenbeck models on a user tree.}
\description{
  \code{ouSim} simulates the evolution of a single character on an {ape}-style user tree under arbitrarily
  complex Ornstein-Uhlenbeck models.
}

ouSim <- function(phy, rootState = 0, shiftBranches = NULL, shiftStates = NULL, alpha = 0, variance = 1, theta = rootState, model = "OU", branchMeans = NULL, steps = 1000) {
## function to plot a simulated dataset under brownian motion or Ornstein-Uhlenbeck (OU) model
## Arguments:
##   phy is an ape-style tree
##   alpha and theta are either single values or vectors of length (length(branchList))
##   shiftBranches is a vector indicating any branches at which an OU or brownian motion model has a determined shift in ancestral state
##   shiftStates is a vector of length = length(shiftBranches) indicaing the ancestral states for the determined break points
## Models:
##  "OU" is a brownian motion or OU model 
##  "platt" is the Platt-Molitor model, in which the only phylogenetic effect is the mean and variance for a given branch
## Andrew Hipp (ahipp@mortonarb.org), January 2008 
## July 2008: modified to accomodate a vector of alpha and theta corresponding to branches


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
  \item{branch}{
    An optional named vector specifying the end nodes of branches.
    The unique branch that terminates at the named node will be painted with the specified regime.
  }
  \item{which}{
    integer;
    if \code{tree} is a \code{hansentree}, start not with a blank canvas but with the regime specifications \code{tree} contains for the character indicated by \code{which}.
  }
}
\details{
  The names of \code{subtree} and \code{branch} must be the names of nodes of \code{tree}.
  The painting proceeds in a particular order: one can overpaint a branch.
  The subtrees indicated by the elements of \code{subtree} are painted first, in order.
  Then the branches indicated by \code{branch} are painted.
  If \code{tree} is a simple \code{ouchtree} object, then \code{paint} begins with a blank canvas, i.e., a tree painted with the single regime "nonspec".
  If \code{tree} inherits class \code{hansentree}, then \code{paint} begins with the regimes specified in the \code{regimes} slot of \code{tree}.
  Note that, if \code{tree} is a multivariate \code{hansentree}, then there are multiple regime specifications contained in \code{tree}.
  In this case, the argument \code{which} lets you pick which one you wish to begin with;
  by default, the first is used.
}
\value{
  A vector of class 'factor' with names corresponding to the nodes in \code{tree}, specifying selective regimes.
}
\author{Aaron A. King kingaa at umich dot edu}
\seealso{\code{ouchtree}, \code{hansen}}
\examples{
data(bimac)
x <- with(bimac,ouchtree(nodes=node,times=time/max(time),ancestors=ancestor,labels=species))
r <- paint(x,subtree=c("1"="medium","9"="large","2"="small"),branch=c("38"="large","2"="medium"))
plot(x,regimes=r,node.names=TRUE)
# compare to bimac['OU.LP']
h5 <- hansen(data=log(bimac['size']),tree=x,regimes=bimac['OU.LP'],alpha=1,sigma=1,reltol=1e-5)
r <- paint(h5,branch=c("18"="large"),subtree=c("9"="small"))
plot(x,regimes=r,node.names=TRUE)
}
\keyword{models}
