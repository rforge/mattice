\name{regimeMatrix}
\Rdversion{1.1}
\alias{regimeMatrix}
\title{
%%  ~~function to do ... ~~
}
\description{
  \code{regimeMatrix} is a recursive function that returns a matrix of n columns by k rows, where n is a number
  of items and k is the number of all possible permutations of those items, up to a maximum number of items
  defined by \code{maxNodes}, and including a row with zero items. 
}
\usage{
regimeMatrix(n, maxNodes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{maxNodes}{
%%     ~~Describe \code{maxNodes} here~~
}
}
\details{
  Typically, this function is called by \code{regimeVectors},
  but it may be useful to 
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, maxNodes) 
{
    if (n == 1) 
        return(matrix(1:0, nrow = 2, ncol = 1))
    outmat <- matrix(NA, nrow = 0, ncol = n)
    for (i in 1:(n - 1)) {
        temp <- c(rep(0, (i - 1)), 1)
        remainder <- n - i
        if (maxNodes > 1 && remainder > 0) {
            nextMat <- regimeMatrix(remainder, maxNodes - 1)
            temp <- cbind(matrix(temp, dim(nextMat)[1], length(temp), 
                byrow = TRUE), nextMat)
        }
        else temp[(i + 1):n] <- rep(0, length((i + 1):n))
        outmat <- rbind(outmat, temp)
    }
    outmat <- rbind(outmat, c(rep(0, n - 1), 1))
    outmat <- rbind(outmat, rep(0, n))
    dimnames(outmat) = list(seq(dim(outmat)[1]), seq(dim(outmat)[2]))
    return(outmat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
