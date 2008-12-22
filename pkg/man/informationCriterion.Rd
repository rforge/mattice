\name{informationCriterion}
\alias{informationCriterion}
\alias{informationCriterion.hansenBatch}
\alias{print.informationCriterion}
\title{Information criterion and weights for a set of models}
\description {
  Returns AIC, AICc, and BIC values and weights for a set of models.
}
\usage {
  informationCriterion(u = NULL, lnL = NULL, K, n = 0, names = NULL)
  informationCriterion.hansenBatch(hansenBatch)
  print.informationCriterion(ic)
}
\details {
  At the minimum, a vector of either the model log-likelihoods (\code{lnL}) or deviances (\code{u} = -2 * lnL) and a vector
  of number of free parameters for each model (\code{K}) must be provided for the function to work. If the sample size
  (\code{n}) is not provided, the function calculates AICc and BIC assuming \code{n} = 1. Information criterion statistics
  are calculated following Burnham and Anderson (2002).
}
\arguments{
  \item{u}{
     A vector of deviances, indexed by model.
     }
  \item{lnL}{
     A vector of log-likelihoods, indexed by model.
     }
  \item{K}{
     A vector of degrees-of-freedom / number of free parameters, indexed by model.
     }
  \item{n}{
     Sample size; for a phylogenetic comparative analysis, \code{n} = the number of tips.
     }
  \item{names}{
     Optional vector of model names, indexed by model.
     }
   \item{hansenBatch}{
     Output from \{runBatchHansen}.
   \item{ic}{
     Output from \{informationCriterion}.
}
\value{
  A list with the following vectors, all indexed by model number:
    \item{names}{
    Model names; if not provided, a vector of 1:length(u).
    }
    \item{u}{
    Deviance.
    }
    \item{K}{
    Degrees of freedom.
    }
    \item{AIC}{
    Akaike information criterion.
    }
    \item{AICc}{
    Small-sample AIC.
    }
    \item{BIC}{
    Bayes information criterion.
    }
    \item{AICwi}{
    AIC weight.
    }
    \item{AICcwi}{
    AICc weight.
    }
    \item{BICwi}{
    BIC weight.
    }
    
}
\author{Andrew Hipp ahipp@mortonarb.org}
\references{
  Burnham, K. P., and D. R. Anderson (2002) 
  Model selection and multimodel inference: a practical information-theoretic approach. 
  Springer, New York.
}