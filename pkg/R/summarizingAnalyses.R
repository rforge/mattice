# ---------------------------------------------------------------------
# FUNCTIONS FOR SUMMARIZING ANALYSES
# ---------------------------------------------------------------------
# functions included in this file:
# 1. summary.hansenBatch

summary.hansenBatch <- function(hansenBatch){
## items in output: hansens, regimeList, regimeMatrix
## the summary will eventually sum weights over all nodes over all trees
## for now, only doing first tree
  # 0. Get information criterion weights for all models
  icObject <- informationCriterion.hansenBatch(hansenBatch)
  matrixRows <- c('aic', 'aicc', 'bic')
  
  tree <- 1 # for now ignore all but the first tree
  
  # 1. sum over nodes
  aicList <- icObject[[tree]]$AICwi
  aiccList <- icObject[[tree]]$AICcwi
  bicList <- icObject[[tree]]$BICwi
  modelsMatrix <- cbind(aicList, aiccList, bicList) # value: modelsMatrix
  dimnames(modelsMatrix) = list(icObject[[tree]]$names, matrixRows)
  
  nonBrownWI <- 1 # defaults to no brownian motion model weights in case none are present
  if(hansenBatch$brown) {
    brownWeights <- modelsMatrix['brown', ] # value: brownWeights
    nonBrownWI <- 1-brownWeights # this is just to make it easy to normalize the non-Brownian OU models
    }

  nodes <- dimnames(hansenBatch$regimeMatrices[[tree]])[[2]]
  nodeWeightsMatrix <- matrix(NA, nrow = length(matrixRows), ncol = length(nodes), dimnames = list(matrixRows, nodes)) # value: nodeWeightsMatrix
  for(i in 1:length(nodes)) {
    modelsMatrixSubset <- modelsMatrix[hansenBatch$regimeMatrices[[tree]][, i] == 1, ]
    if(identical(dim(modelsMatrixSubset), NULL)) nodeWeightsMatrix[, i] <- modelsMatrixSubset / nonBrownWI # because extracting a single row yields a vector, and dim returns NULL for a vector
    else nodeWeightsMatrix[, i] <- apply(modelsMatrixSubset, 2, sum) / nonBrownWI
    }

  # 2. sum over number of parameters
  
  kCats <- sort(unique(hansenBatch$hansens[[tree]][, 'dof']))
  kMatrix <- matrix(NA, nrow = length(matrixRows), ncol = length(kCats), dimnames = list(matrixRows, as.character(kCats))) #value: kMatrix
  for(i in as.character(kCats)) {
    modelsMatrixSubset <- modelsMatrix[hansenBatch$hansens[[tree]][, 'dof'] == i, ]
    if(identical(dim(modelsMatrixSubset), NULL)) kMatrix[, i] <- modelsMatrixSubset
    else kMatrix[, i] <- apply(modelsMatrixSubset, 2, sum) 
    }
  
  outdata <- list(brownWeights = brownWeights, modelsMatrix = modelsMatrix, nodeWeightsMatrix = nodeWeightsMatrix, kMatrix = kMatrix)
  print(brownWeights)
  print(nodeWeightsMatrix)
  return(outdata)
  }