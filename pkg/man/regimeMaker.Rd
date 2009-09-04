\name{regimeMaker}
\Rdversion{1.1}
\alias{regimeMaker}
\title{
  Returns ouch-style regimes defined by a matrix, with each row specifying which nodes permit character transitions
}
\description{
  This function utilizes a matrix of the format returned by \code{regimeMatrix}
}
\usage{
regimeMaker(ouchTrees, regMatrix, nodeMembers)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ouchTrees}{
%%     ~~Describe \code{ouchTrees} here~~
}
  \item{regMatrix}{
%%     ~~Describe \code{regMatrix} here~~
}
  \item{nodeMembers}{
%%     ~~Describe \code{nodeMembers} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{regimeMatrix}}, \code{\link{paintBranches}}
}
\examples{
  data(carex)
  # tree comes in in phylo format, but we need an ouchtree object
  trees <- lapply(carex$ovales.bayesTrees[1:10], ape2ouch) 
  # if we are interested in comparing the single-change models that allow changes only at nodes 1, 2, 4, 5, and 8, as well as the no-change model, we could do:
  mat <- regimeMatrix(n = 8, maxNodes = 1)[c(1,2,4,5,8,9), ] # the 9th row is the no-change model
  regs <- regimeMaker(trees, mat, carex$ovales.nodes)
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
