\name{carex}
\alias{carex}
\docType{data}
\title{Carex section Ovales tree, chromosome data, and nodes}
\description{
  This is the Carex dataset used in Hipp (2007). It is a list with three items: a tree of 53 taxa for which chromosome
  counts are known; a vector of log-transformed mean chromosome counts for all taxa; and a list of eight taxon vectors.
}
\usage{data(carex)}
\format{
  A list with three items:
    \item{\code{ovales.tree}}{
     An ultrametric tree in \code{phylo} format with 53 tips
     }
    \item{\code{ovales.data}}{
     A named numeric vector of log-transformed mean chromosome counts for all 53 taxa
     }
    \item{\code{ovales.nodes}}{
     A list of eight taxon vectors defining the eight nodes studied.
     }
}
\details{
  Phylogeny (\code{ovales.tree}) was estimated for approximately 80 species (Hipp 2006), branch lengths optimized using
  penalized likelihood in \code{ape}, and then the tree was pruned down to only those taxa for which chromosome counts were
  available. Chromosome counts (\code{ovales.data}) were averaged by population, and the means log-transformed for
  analysis. Nodes were selected for analysis based on clade support and number of taxa in the clade.
}
\source{
  Hipp, A.L. (2007) 
  Non-uniform processes of chromosome evolution in sedges (Carex: Cyperaceae). 
  \emph{Evolution} 61:2175-2194.
}
\references{
  Hipp, A.L., P.E. Rothrock, A.A. Reznicek, and J.A. Weber (2006) 
  Phylogeny and classification of Carex section Ovales (Cyperaceae). 
  \emph{International Journal of Plant Sciences} 167:1029-1048.
  
  Hipp, A.L. (2007) 
  Non-uniform processes of chromosome evolution in sedges (Carex: Cyperaceae). 
  \emph{Evolution} 61:2175-2194.
}

\examples{
  library(maticce)
  data(carex)
  ovales.tree <- ape2ouch(carex$ovales.tree) # tree comes in in phylo format, but we need an ouchtree object
  trial <- runBatchHansen(carex$ovales.tree, carex$ovales.data, carex$ovales.nodes[1:4], maxNodes = 2) # for expedience, only tests for changes at up to 2 of the first 4 nodes
  summary(trial) # summarizes results using summary.hansenBatch and displays the results using print.hansenSummary
  multiModel(carex$ovales.tree, carex$ovales.data, carex$ovales.nodes[[2]]) # compares five different models of character change at node 2
  trialSim <- ouSim(trial, carex$ovales.tree) # simulates the evolution of the chromosome number under the model-averaged values
  plot(trialSim) # plots the character simulation, with all branches black
  plot(trialSim, colors = paintBranches(carex$ovales.nodes[1:4], carex$ovales.tree)) # plots the character simulation, with branch colors changing at all 8 nodes
  plot(trialSim, colors = paintBranches(list(carex$ovales.nodes[[2]]), carex$ovales.tree)) # plots the character simulation, with branch colors changing only at node 2
}
\author{Andrew L. Hipp <ahipp@mortonarb.org>}
\keyword{datasets}                                            