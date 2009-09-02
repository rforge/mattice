\name{regimeMaker}
\Rdversion{1.1}
\alias{regimeMaker}
\title{
  Returns ouch-style regimes defined by a matrix, with each row specifying which nodes permit character transitions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (ouchTrees, regMatrix, nodeMembers) 
{
    numTrees <- length(ouchTrees)
    numNodes <- length(nodeMembers)
    if (numNodes != dim(regMatrix)[2]) 
        stop("Number of nodes (columns) in regMatrix must equal number of items in nodeMembers list")
    nodeMatrix <- matrix(NA, nrow = numTrees, ncol = numNodes, 
        dimnames = list(seq(numTrees), dimnames(regMatrix)[[2]]))
    changeNodes <- list(numTrees)
    regList <- list(numTrees)
    regMatrices <- list(numTrees)
    for (i in seq(numNodes)) nodeMatrix[, i] <- unlist(lapply(ouchTrees, 
        isMonophyletic, taxa = nodeMembers[[i]]))
    for (i in seq(numTrees)) {
        tree <- ouchTrees[[i]]
        regMatrices[[i]] <- regMatrix * as.numeric(matrix(nodeMatrix[i, 
            ], dim(regMatrix)[1], dim(regMatrix)[2], byrow = TRUE))
        regMatrices[[i]][1:(dim(regMatrices[[i]])[1] - 1), ][which(apply(regMatrices[[i]][1:(dim(regMatrices[[i]])[1] - 
            1), ], 1, sum) == 0), ] <- rep(NA, numNodes)
        regMatrices[[i]][duplicated(apply(regMatrices[[i]], 1, 
            as.decimal)), ] <- rep(NA, numNodes)
        dimnames(regMatrices[[i]]) <- list(seq(dim(regMatrices[[i]])[1]), 
            dimnames(regMatrices[[i]])[[2]])
        numTreeRegs <- dim(regMatrices[[i]])[1]
        treeRegs <- list(numTreeRegs)
        nodesVector <- unlist(lapply(nodeMembers, mrcaOUCH, tree = ouchTrees[[i]]))
        for (j in seq(numTreeRegs)) {
            if (any(is.na(regMatrices[[i]][j, ]))) 
                treeRegs[[j]] <- NA
            else {
                treeRegs[[j]] <- as.factor(paintBranches(c("1", 
                  nodesVector[as.logical(regMatrices[[i]][j, 
                    ])]), tree))
                names(treeRegs[[j]]) <- tree@nodes
            }
        }
        regList[[i]] <- treeRegs
    }
    regMatrices$overall <- regMatrix
    outdata <- list(regList = regList, nodeMatrix = nodeMatrix, 
        regMatrix = regMatrices)
    return(outdata)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
