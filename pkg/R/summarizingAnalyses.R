# ---------------------------------------------------------------------
# FUNCTIONS FOR SUMMARIZING ANALYSES
# ---------------------------------------------------------------------

# March 2011 - summary.hansenBatch found to give incorrect answers with multiple trees; corrected. Also:
#  - changed the weighting so that model averaging can be by AICc, AIC, or BIC
#  - added support for nodeNames

summary.hansenBatch <- function(object, ic = 'AICc', ...){
## items in output: hansens, regimeList, regimeMatrix
## ic = choice of information criterion weight to use in model averaging
  hansenBatch <- object
  icObject <- informationCriterion.hansenBatch(object) # Get information criterion weights for all models
  nmodels <- dim(hansenBatch$hansens[[1]])[1] # number of models per tree (ignores the fact that models may not be present in all trees)
  ntrees <- length(hansenBatch$hansens) # number of trees
  nodeSums <- colSums(hansenBatch$nodeMatrix) # number of trees possessing each node
  nnodes <- length(nodeSums) # number of nodes being studied
  nodes <- dimnames(hansenBatch$regMatrix$overall)[[2]] # grab the overall regMatrix, which includes all possible nodes
  sigmaSqVector <- numeric(ntrees) # vector to capture model-averaged sigma^2 for each tree
  alphaVector <- numeric(ntrees) # vector to capture model-averaged alpha for each tree
  modelsMatrix <- vector('list', ntrees) # list of matrices, indexed by tree, holding the weight for each model
  matrixRows <- c('AICwi', 'AICcwi', 'BICwi') # rows in the matrix
  nodeWeightsSummed <- matrix(0, nrow = length(matrixRows), ncol = nnodes, dimnames = list(matrixRows, nodes)) # holds node weights summed; zero-filled b/c it is a sum?
  icMats <- vector('list', length(matrixRows))
  names(icMats) <- matrixRows
  for(i in matrixRows) icMats[[i]] <- matrix(NA, nrow = ntrees, ncol = nnodes, dimnames = list(NULL, nodes))
  thetaMatrix <- matrix(NA, nrow = ntrees, 
                            ncol = dim(hansenBatch$thetas[[1]])[2], 
                            dimnames = list(1:ntrees, dimnames(hansenBatch$thetas[[1]])[[2]])
                       )
  for(tree in 1:ntrees) {
    modelsMatrix[[tree]] <- cbind(icObject[[tree]]$AICwi, icObject[[tree]]$AICcwi, icObject[[tree]]$BICwi)
    dimnames(modelsMatrix[[tree]]) <- list( dimnames(hansenBatch$hansens[[1]])[[1]], c("AICwi", "AICcwi", "BICwi"))
    icWeight <- switch(ic,
    			BIC = icObject[[tree]]$BICwi,
			AIC = icObject[[tree]]$AICwi,
			AICc = icObject[[tree]]$AICcwi
			)
     for(i in matrixRows) {
	  icMats[[i]][tree, ] <- colSums(modelsMatrix[[tree]][, i] * hansenBatch$regMatrix$overall, na.rm = T)
      nodeWeightsSummed[i, ] <- icMats[[i]][tree, ] + nodeWeightsSummed[i, ]
	  } # close i
    sigmaSqVector[tree] <- weighted.mean(hansenBatch$hansens[[tree]][, 'sigma.squared'], icWeight, na.rm = TRUE)
    if(hansenBatch$brown) icOU <- icWeight[1: (length(icWeight) - 1)]
    alphaVector[tree] <- ifelse(hansenBatch$brown, 
                                weighted.mean(hansenBatch$hansens[[tree]][1:(nmodels - 1), 'theta / alpha'], icOU, na.rm = TRUE),
                                weighted.mean(hansenBatch$hansens[[tree]][ , 'theta / alpha'], icWeight, na.rm = TRUE)
                                )
    if(hansenBatch$brown) w <- icOU else w <- icWeight
    thetaMatrix[tree, ] <- apply(hansenBatch$thetas[[tree]], 2, 
                                 weighted.mean, 
                                 w = w, 
                                 na.rm = TRUE
                                 )
                                 
  }
  # in this matrix, the weight for each node is averaged only over trees that possess that node
  nodeWeightsMatrix.unnormalized <- nodeWeightsSummed / matrix(nodeSums, nrow = dim(nodeWeightsSummed)[1], ncol = nnodes, byrow = TRUE)
  # in this matrix, the weight for each node is averaged over all trees
  nodeWeightsMatrix.allNodes <- nodeWeightsSummed / ntrees 
  
  if(!identical(object$nodeNames, NULL)) dimnames(nodeWeightsMatrix.unnormalized)[[2]] <- dimnames(nodeWeightsMatrix.allNodes)[[2]] <- object$nodeNames
  modelAvgAlpha <- c(mean(alphaVector, na.rm = TRUE), quantile(alphaVector, c(0.025, 0.975)))
  modelAvgSigmaSq <- c(mean(sigmaSqVector, na.rm = TRUE),quantile(sigmaSqVector, c(0.025, 0.975)))
  names(modelAvgAlpha) <- names(modelAvgSigmaSq) <- c('mean', 'lower.CI', 'upper.CI')
  outdata <- list(modelsMatrix = modelsMatrix, nodeWeightsMatrix = list(unnormalized = nodeWeightsMatrix.unnormalized, allNodes = nodeWeightsMatrix.allNodes), modelAvgAlpha = modelAvgAlpha, modelAvgSigmaSq = modelAvgSigmaSq, thetaMatrix = thetaMatrix, icMats = icMats, ic = ic)
  class(outdata) <- 'hansenSummary'
  return(outdata)
}

replace.matrix <- function (x, oldValue, newValue) {
  if(is.na(oldValue)) x[is.na(x)] <- newValue
  else x[x == oldValue] <- newValue
  return(x)
}

print.hansenSummary <- function(x, ...) {
## This just formats a hansenSummary object so that it is readable on the screen; you can still store the summary object and extract elements as needed
  hansenSummary <- x
  cat(paste("\nSummarizing hansenBatch analyses over", length(hansenSummary$modelsMatrix), "trees and", dim(hansenSummary$modelsMatrix[[1]])[1], "models"))
  cat("\n-----------------------------------------------------------")
  cat("\nESTIMATED SUPPORT FOR CHANGES OCCURRING AT DESIGNATED NODES")
  cat("\nAveraged over all trees:\n")
  print(hansenSummary$nodeWeightsMatrix$allNodes)
  cat("\nSupport conditioned on trees that possess the node\n")
  print(hansenSummary$nodeWeightsMatrix$unnormalized)
  cat("\nMODEL-AVERAGED PARAMETERS BASED ON", x$ic, "WEIGHTS\n")
  cat("\nalpha:\n")
  print(hansenSummary$modelAvgAlpha)
  cat("\nsigma^2:\n")
  print(hansenSummary$modelAvgSigmaSq)
  if(any(dim(hansenSummary$thetaMatrix) > 12)) message(paste("\ntheta matrix is too long to display; access through the summary object"))
  else {
    cat("\ntheta matrix, with branches as the columns and trees as the rows:\n")
    print(hansenSummary$thetaMatrix)
    cat('\n')
    }
}