# ---------------------------------------------------------------------
# FUNCTIONS FOR SUMMARIZING ANALYSES
# ---------------------------------------------------------------------

summary.hansenBatch <- function(hansenBatch){
## items in output: hansens, regimeList, regimeMatrix
  icObject <- informationCriterion.hansenBatch(hansenBatch) # Get information criterion weights for all models
  nmodels <- dim(hansenBatch$hansens[[1]])[1] # number of models per tree (ignores the fact that models may not be present in all trees)
  ntrees <- length(hansenBatch$hansens) # number of trees
  nodeSums <- colSums(hansenBatch$nodeMatrix) # number of trees possessing each node
  nnodes <- length(nodeSums) # number of nodes being studied
  nodes <- dimnames(hansenBatch$regMatrix$overall)[[2]] # grab the overall regMatrix, which includes all possible nodes
  sigmaSqVector <- numeric(ntrees) # vector to capture model-averaged sigma^2 for each tree
  alphaVector <- numeric(ntrees) # vector to capture model-averaged alpha for each tree
  modelsMatrix <- vector('list', ntrees) # list of matrices, indexed by tree, holding the weight for each model
  matrixRows <- c('AIC.weight', 'AICc.weight', 'BIC.weight') # rows in the matrix
  nodeWeightsSummed <- matrix(0, nrow = length(matrixRows), ncol = nnodes, dimnames = list(matrixRows, nodes)) # holds node weights summed
  thetaMatrix <- matrix(NA, nrow = ntrees, 
                            ncol = dim(hansenBatch$thetas[[1]])[2], 
                            dimnames = list(1:ntrees, dimnames(hansenBatch$thetas[[1]])[[2]])
                       )
  for(tree in 1:ntrees) {
    modelsMatrix[[tree]] <- cbind(icObject[[tree]]$AICwi, icObject[[tree]]$AICcwi, icObject[[tree]]$BICwi)
    bic <- icObject[[tree]]$BICwi
    for(i in seq(nnodes)) {
      modelsMatrixSubset <- modelsMatrix[[tree]][hansenBatch$regMatrix$overall[, nodes[i]] == 1, ] # subset models that contain node i
      if(identical(dim(modelsMatrixSubset), NULL)) # is modelsMatrixSubset a 1-d vector? if so then:
        nodeWeightsSummed[, nodes[i]] <- nodeWeightsSummed[, nodes[i]] + replace.matrix(modelsMatrixSubset, NA, 0) # because extracting a single row yields a vector, and dim returns NULL for a vector
      else nodeWeightsSummed[, nodes[i]] <- nodeWeightsSummed[, nodes[i]] + colSums(modelsMatrixSubset, na.rm = T)
	  }
    sigmaSqVector[tree] <- weighted.mean(hansenBatch$hansens[[tree]][, 'sigma.squared'], bic, na.rm = T)
    if(hansenBatch$brown) bicOU <- bic[1: (length(bic) - 1)]
    alphaVector[tree] <- ifelse(hansenBatch$brown, 
                                weighted.mean(hansenBatch$hansens[[tree]][1:(nmodels - 1), 'theta / alpha'], bicOU, na.rm = T),
                                weighted.mean(hansenBatch$hansens[[tree]][ , 'theta / alpha'], bic, na.rm = T) 
                                )
    if(hansenBatch$brown) w <- bicOU else w <- bic
    thetaMatrix[tree, ] <- apply(hansenBatch$thetas[[tree]], 2, 
                                 weighted.mean, 
                                 w = w, 
                                 na.rm = T
                                 )
                                 
  }
  # in this matrix, the weight for each node is averaged only over trees that possess that node
  nodeWeightsMatrix.unnormalized <- nodeWeightsSummed / matrix(nodeSums, nrow = dim(nodeWeightsSummed)[1], ncol = nnodes, byrow = T)
  # in this matrix, the weight for each node is averaged over all trees
  nodeWeightsMatrix.allNodes <- nodeWeightsSummed / ntrees 
  
  # sum over number of parameters
  # create a vector of sums that tells us how many categories there are for each model: dof = sum(nodes) + 1 [because a node indicates a change in 
  #   regime, thus the total number of thetas = nodes + 1] + alpha + sigma = sum(nodes) + 3; for Brownian motion model, dof = 2
  
  #nodeSums <- apply(hansenBatch$regMatrix$overall, 1, sum) + 3
  #if(hansenBatch$brown) nodeSums['brown'] <- 2
  #kCats <- sort(unique(nodeSums)) # and just give us the unique degree-of-freedom categories, sorted
  #kMatrix <- matrix(NA, nrow = length(matrixRows), ncol = length(kCats), dimnames = list(matrixRows, as.character(kCats))) # make the empty kMatrix
  #for(i in as.character(kCats)) {
  #  modelsMatrixSubset <- weightsMatrix.allNodes[nodeSums == i, ] 
  #  if(identical(dim(modelsMatrixSubset), NULL)) kMatrix[, i] <- modelsMatrixSubset # is modelsMatrixSubset a 1-d vector?
  #  else kMatrix[, i] <- apply(modelsMatrixSubset, 2, sum) 
  #}
  modelAvgAlpha <- mean(alphaVector, na.rm = T)
  modelAvgSigmaSq <- mean(sigmaSqVector, na.rm = T)
  outdata <- list(modelsMatrix = modelsMatrix, nodeWeightsMatrix = list(unnormalized = nodeWeightsMatrix.unnormalized, allNodes = nodeWeightsMatrix.allNodes), modelAvgAlpha = modelAvgAlpha, modelAvgSigmaSq = modelAvgSigmaSq, thetaMatrix = thetaMatrix)
  class(outdata) <- 'hansenSummary'
  return(outdata)
}

replace.matrix <- function (x, oldValue, newValue) {
  if(is.na(oldValue)) x[is.na(x)] <- newValue
  else x[x == oldValue] <- newValue
  return(x)
}

print.hansenSummary <- function(hansenSummary) {
## This just formats a hansenSummary object so that it is readable on the screen; you can still store the summary object and extract elements as needed
  message(paste("\nSummarizing hansenBatch analyses over", length(hansenSummary$modelsMatrix), "trees and", dim(hansenSummary$modelsMatrix[[1]])[1], "models"))
  message("-----------------------------------------------------------")
  message("ESTIMATED SUPPORT FOR CHANGES OCCURRING AT DESIGNATED NODES")
  message("Averaged over all trees, thus multiplying the clade distribution by the mean support for the model:")
  print(hansenSummary$nodeWeightsMatrix$allNodes)
  message("\nSupport conditioned on trees that possess the node")
  print(hansenSummary$nodeWeightsMatrix$unnormalized)
  #message("\nESTIMATED SUPPORT FOR NUMBER OF PARAMETERS IN THE MODEL")
  #message("The properties of this support value have not been studied and are likely to be biased strongly toward the median value of\nK, as K is largest at the median values (they are distributed according to Stirling numbers of the first kind).")
  #print(hansenSummary$kMatrix)
  message("\nMODEL-AVERAGED PARAMETERS")
  message(paste("alpha =", hansenSummary$modelAvgAlpha))
  message(paste("sigma^2 =", hansenSummary$modelAvgSigmaSq))
  if(any(dim(hansenSummary$thetaMatrix) > 12)) message(paste("theta matrix is too long to display; access through the summary object"))
  else {
    message("theta matrix, with branches as the columns and trees as the rows:")
    print(hansenSummary$thetaMatrix)
    }
}