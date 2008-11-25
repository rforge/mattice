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
  nmodels <- dim(hansenBatch$hansens[[1]])[1]
  ntrees <- length(hansenBatch$hansens)
  modelsMatrix <- vector('list', ntrees)
  matrixRows <- c('AIC.weight', 'AICc.weight', 'BIC.weight')
  outMatrix <- matrix(0, nrow = length(icObject[[1]]$AICwi), ncol = length(matrixRows), dimnames = list(icObject[[1]]$names, matrixRows))
  nmodelsMatrix <- outMatrix # a matrix to divide outMatrix by so that average BICwi are only based on trees that have a node
  for(tree in 1:ntrees) {
    aicList <- icObject[[tree]]$AICwi
    aiccList <- icObject[[tree]]$AICcwi
    bicList <- icObject[[tree]]$BICwi
    modelsMatrix[[tree]] <- cbind(aicList, aiccList, bicList) # value: modelsMatrix
    temp <- modelsMatrix[[tree]]; temp[!is.na(temp)] <- 1
    nmodelsMatrix <- nmodelsMatrix + replace.matrix(temp, NA, 0) # nmodelsMatrix gets a 1 for each model present, a 0 for each model not present
    outMatrix <- outMatrix + replace.matrix(modelsMatrix[[tree]], NA, 0) # replace NAs with 0 so that sum works correctly at the end
  
    ## the lines below made the weights on branches ignore the fact that the Brownian motion model was part of the
    ##   model set; however, I've removed them b/c support for the Brownian motion model does (and should) contribute 
    ##   to reduced probability of change at any of the nodes. You can uncomment them if you feel differently.
    
    # nonBrownWI <- 1 # defaults to no brownian motion model weights in case none are present
    # if(hansenBatch$brown) {
    #   brownWeights <- modelsMatrix['brown', ] # value: brownWeights
    #   nonBrownWI <- 1-brownWeights # this is just to make it easy to normalize the non-Brownian OU models
    # }
  }
  weightsMatrix.unnormalized <- outMatrix / nmodelsMatrix # in this matrix, the weight for each model is averaged only over trees
                                                          # that possess that node; weights may not sum to 1.0
  weightsMatrix.allNodes <- outMatrix / ntrees # in this matrix, the weight for each model is averaged over all trees, setting weight
                                               # equal to zero in any trees that lack that model. Weights sum to 1.0, and they
                                               # factor in the posterior probability (or bootstrap proportion) for each model. 
                                               # Which to use? they are both informative, but this one has the desirable property of
                                               # acting being a proper probability distribution, albeit one that confounds clade support
                                               # with model support.
  
  # 1. sum over nodes
  nodes <- dimnames(hansenBatch$regMatrix$overall)[[2]] # grab the overall regMatrix, which includes all possible nodes, no matter which trees do or don't have them
  nodeWeightsMatrix.unnormalized <- matrix(NA, nrow = length(matrixRows), ncol = length(nodes), dimnames = list(matrixRows, nodes)) # value: nodeWeightsMatrix
  nodeWeightsMatrix.allNodes <- nodeWeightsMatrix.unnormalized # same dimensions, another empty matrix to fill up
  for(i in 1:length(nodes)) {
    modelsMatrixSubset <- weightsMatrix.allNodes[hansenBatch$regMatrix$overall[, nodes[i]] == 1, ]
    if(identical(dim(modelsMatrixSubset), NULL)) # is modelsMatrixSubset a 1-d vector? if so then:
      nodeWeightsMatrix.allNodes[, nodes[i]] <- modelsMatrixSubset # because extracting a single row yields a vector, and dim returns NULL for a vector
    else nodeWeightsMatrix.allNodes[, nodes[i]] <- apply(modelsMatrixSubset, 2, sum)
    
    modelsMatrixSubset <- weightsMatrix.unnormalized[hansenBatch$regMatrix$overall[, nodes[i]] == 1, ]
    if(identical(dim(modelsMatrixSubset), NULL)) # again, is modelsMatrixSubset a 1-d vector? 
      nodeWeightsMatrix.unnormalized[, nodes[i]] <- modelsMatrixSubset # because extracting a single row yields a vector, and dim returns NULL for a vector
    else nodeWeightsMatrix.unnormalized[, nodes[i]] <- apply(modelsMatrixSubset, 2, sum)
  }

  # 2. sum over number of parameters
  # create a vector of sums that tells us how many categories there are for each model: dof = sum(nodes) + 1 [because a node indicates a change in 
  #   regime, thus the total number of thetas = nodes + 1] + alpha + sigma = sum(nodes) + 3; for Brownian motion model, dof = 2
  nodeSums <- apply(hansenBatch$regMatrix$overall, 1, sum) + 3
  if(hansenBatch$brown) nodeSums['brown'] <- 2
  kCats <- sort(unique(nodeSums)) # and just give us the unique degree-of-freedom categories, sorted
  kMatrix <- matrix(NA, nrow = length(matrixRows), ncol = length(kCats), dimnames = list(matrixRows, as.character(kCats))) # make the empty kMatrix
  for(i in as.character(kCats)) {
    modelsMatrixSubset <- weightsMatrix.allNodes[nodeSums == i, ] # note this is right!
    if(identical(dim(modelsMatrixSubset), NULL)) kMatrix[, i] <- modelsMatrixSubset # is modelsMatrixSubset a 1-d vector?
    else kMatrix[, i] <- apply(modelsMatrixSubset, 2, sum) 
  }
  warning('the K-matrix as currently implemented is calculated over the all-nodes (normalized) weights.')
  outdata <- list(modelsMatrix = modelsMatrix, weightsMatrix = list(unnormalized = weightsMatrix.unnormalized, allNodes = weightsMatrix.allNodes), nodeWeightsMatrix = list(unnormalized = nodeWeightsMatrix.unnormalized, allNodes = nodeWeightsMatrix.allNodes), kMatrix = kMatrix)
  return(outdata)
}

replace.matrix <- function (x, oldValue, newValue) {
  if(is.na(oldValue)) x[is.na(x)] <- newValue
  else x[x == oldValue] <- newValue
  return(x)
}
