\name{isMonophyletic}
\alias{isMonophyletic}
\alias{mrcaOUCH}
\alias{nodeDescendents}
\title{Test whether a set of taxa is monophyletic on an ouch-style tree.}
\usage{
  isMonophyletic(tree, taxa)
  nodeDescendents(tree, startNode)
  mrcaOUCH(cladeVector, tree)
}
\arguments{
  \item{tree}{
    An \code{ouchtree} object.
    }
  \item{taxa}{
    A vector of taxon names.
    }
  \item{startNode}{
    The node of interest; nomenclature follows \code{tree@nodes}.
    }
  \item{cladeVector}{
    A vector of tip labels.
    }
}
\value{
    \item{isMonophyletic}{
    Boolean, TRUE if taxa form a monophyletic group, FALSE if not.
    }
    \item{mrcaOUCH}{
    A single node label; one of the items in \code{tree@nodes}.
    }
    \item{nodeDescendents}{
    A vector of tip labels.
    }
  }
\author{Andrew Hipp ahipp@mortonarb.org}
