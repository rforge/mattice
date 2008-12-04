ouSim.hansenBatch <- function(analysis, tree, treeNum = 1, rootState = NULL, ...) {
## runs ouSim.ouchtree for a hansenBatch or hansenSummary object, using the model-averaged alpha, sigma.squared, and theta vector from one tree
## tree, rootState = 0, alpha = 0, variance = 1, theta = rootState, steps = 1000
  if(class(analysis) == "hansenBatch") analysis <- summary(analysis)
  if(identical(rootState, NULL)) rootState <- mean(analysis$thetaMatrix[treeNum, ])
  outdata <- ouSim(tree, rootState, alpha = analysis$modelAvgAlpha, theta = analysis$thetaMatrix[treeNum, ], ...)
  return(outdata)
}