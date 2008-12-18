ouSim.hansenBatch <- function(analysis, tree, treeNum = 1, rootState = NULL, ...) {
## runs ouSim.ouchtree for a hansenBatch or hansenSummary object, using the model-averaged alpha, sigma.squared, and theta vector from one tree
  if(class(analysis) == "hansenBatch") analysis <- summary(analysis)
  if(identical(rootState, NULL)) rootState <- analysis$thetaMatrix[treeNum, ][tree@root] # rootstate taken to be the optimum at the root
  outdata <- ouSim(tree, rootState, alpha = analysis$modelAvgAlpha, variance = analysis$modelAvgSigmaSq, theta = analysis$thetaMatrix[treeNum, ], ...)
  class(outdata) <- "ouSim"
  return(outdata)
}

ouSim.brownHansen <- function(analysis, ...) {
  su <- summary(analysis)
  if(length(analysis@regimes) > 1) warning("Theta is based on analysis@regimes[[1]]")
  if(dim(su$alpha)[1] != 1) stop("This is a one-character simulation; analysis appears to be based on > 1 character")
  if(class(analysis) == "browntree") {
    alpha <- 0
    theta <- 0
    rootState <- su$theta[[1]]
    }
  if(class(analysis) == "hansentree") {
    alpha <- as.vector(su$alpha)
    theta <- su$optima[[1]][analysis@regimes[[1]]]
    rootState <- theta[analysis@root] # rootstate taken to be the optimum at the root
    }
  variance <- as.vector(su$sigma.squared)
  tree <- ouchtree(analysis@nodes, analysis@ancestors, analysis@times) 
  outdata <- ouSim.ouchtree(tree, rootState, alpha, variance, theta, ...)
  outdata$colors <- analysis@regimes[[1]]
  class(outdata) <- "ouSim"
  return(outdata)
}