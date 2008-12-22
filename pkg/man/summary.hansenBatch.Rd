\name{summary.hansenBatch}
\alias{summary.hansenBatch}
\title{Summarize output from runBatchHansen}
\description{
  \code{summary.hansenBatch} is a summary method for a \code{hansenBatch} object.
}
\usage{
summary(hansenBatch)
}
\arguments{
  \item{hansenBatch}{
  Output from \code{runBatchHansen}.
  }
}
\details{
    This function averages parameters over models and over trees, and provides the cumulative information-criterion weight
    for models containing each node as an estimate of the support for the hypothesis that there was a character shift
    at that node.
}
\value{
  A list of class \code{hansenSummary} with the following items:
  \item{modelsMatrix}{
    A list of matrices, indexed by tree, providing the information criterion weights for each model.
  }
  \item{nodeWeightsMatrix}{
    A list of length two with two matrices: 
    \item{unnormalized}{
      The cumulative information criterion weight for each node, averaged only over trees that possess that node
      }
    \item{allNodes}{
      The cumulative information criterion weight for each node, averaged over all trees.
      }
    }
  \item{modelAvgAlpha}{
    The model-averaged estimate of \code{alpha}, the rate of evolution toward the equilibrium / optimum
    }
  \item{modelAvgSigmaSq}{
    The model-averaged estimate of \code{sigma.squared}, the variance or rate of stochastic evolution.
    }
  \item{thetaMatrix}{
    A matrix with trees as the rows, branch number as the columns, containing the model-averaged value of \code{theta} 
    (the character equilibrium or optimum value) for each branch on each tree
    }
  }
\author{Andrew Hipp <ahipp@mortonarb.org>}
\seealso{
   \code{\link{print.hansenSummary}} for a readable screen output summarizing most of the above; \code{\link{carex}} for an example
}
