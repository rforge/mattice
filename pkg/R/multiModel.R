multiModel <- function(phy, dat, node, models = c('whole.brown', 'whole.ou1', 'whole.ou2','part.brown', 'part.ou')) {
# test the support for alternative models on simple and partitioned trees
# currently only works on one tree; eventually should be modified so it runs on a set of trees, conditioned on those trees 
#   that have the node of interest and returning percent of trees possessing that node as an additional value
  paramHeader <- c('loglik', 'dof', 'sigma.squared', 'alpha', 'theta', 'optimum', 'optimum.uptree', 'optimum.downtree')
  paramsAll <- c('loglik', 'dof', 'sigma.squared')
  paramSets <- list(brown = c(paramsAll, 'theta'), 
                  ou1 = c(paramsAll, 'alpha', 'optimum'), 
                  ou2 = c(paramsAll, 'alpha', 'optimum.uptree', 'optimum.downtree')
                  )
  modelsAll = c('whole.ou2', 'whole.ou1', 'whole.brown', 'part.ou.uptree', 'part.ou.downtree', 'part.ou.summed', 'part.brown.uptree', 'part.brown.downtree', 'part.brown.summed')
  pSum <- c('loglik', 'dof') # parameters to sum for evaluating partitioned trees

  aboveNodeTaxa <- phy$tip.label[-which(phy$tip.label %in% node)]
  wholeTree <- ape2ouch(phy)

  upTree <- ape2ouch(drop.tip(phy, node))
  downTree <- ape2ouch(drop.tip(phy, aboveNodeTaxa))
  partialTree <- list(upTree = upTree, downTree = downTree)

  outMatrix <- matrix(NA, nrow = length(modelsAll), ncol = length(paramHeader), dimnames = list(modelsAll, paramHeader))
  compareK <- rep(NA, length(models)); names(compareK) <- models
  comparelnL <- rep(NA, length(models)); names(comparelnL) <- models
  modelsSplit <- strsplit(models, ".", fixed = TRUE)

  for(i in modelsSplit) {
   model <- paste(i[1], ".", i[2], sep = "")
   if (i[1] == 'whole') {
     outMatrix[model,  ] <- wholeModel(wholeTree, dat, i[2], node, paramSets[[i[2]]], paramHeader)$params
     comparelnL[model] <- outMatrix[model, 'loglik']
     compareK[model] <- outMatrix[model, 'dof']
     }
   if (i[1] == 'part') {
     outMatrix[c(paste(model, '.uptree', sep = ""), paste(model, '.downtree', sep = ""), paste(model, '.summed', sep = "")), ] <- partialModel(partialTree, dat, i[2], c('uptree', 'downtree'), paramSets[[i[2]]], pSum, paramHeader)$params                                 
     compareK[model] <- outMatrix[paste(model, '.summed', sep = ""), 'dof'] 
     comparelnL[model] <- outMatrix[paste(model, '.summed', sep = ""), 'loglik'] 
     }
   }
  IC <- informationCriterion(lnL = comparelnL, K = compareK, n = length(phy$tip.label), names = models)
  outdata <- list(IC = IC, modelMatrix = outMatrix)
  return(outdata)
  }

wholeModel <- function(phy, dat, model, node, parameterVector, paramHeader) {
  dat <- dat[phy@nodelabels]; names(dat) <- phy@nodes
  if(model == "brown") analysis <- brown(dat, phy)
  if(model == "ou1") analysis <- hansen(dat, 
  					phy, 
                                        regimes = structure(rep(phy@root, phy@nnodes), names = phy@nodes, levels = 1, class = 'factor'),
                                        sqrt.alpha = 1, 
                                        sigma = 1
                                        )
  if(model == "ou2") {
    regime <- paintBranches(list(node), phy)
    uptreeNum <- as.character(phy@root)
    downtreeNum <- as.character(unique(regime))[unique(regime) != phy@root]
    analysis <- hansen(dat, phy, regime, sqrt.alpha = 1, sigma = 1)
    }
  params <- unlist(summary(analysis)[parameterVector])[paramHeader]
  names(params) <- paramHeader
  if(model == "ou2") {
    params['optimum.uptree'] <- summary(analysis)$optima[[1]][uptreeNum]
    params['optimum.downtree'] <- summary(analysis)$optima[[1]][downtreeNum]
    }
  if(model == "ou1") params['optimum'] <- summary(analysis)$optima[[1]]
  out <- list(analysis = analysis, params = params)
  return(out)
  }

partialModel <- function(phyList, dat, model, treeNames, parameterVector = NULL, pSum = NULL, paramHeader) {
# sums a subset of parameters (indicated by pSum) and leaves the others separate
  allParams <- pSum
  analysis <- vector('list', length(phyList))
  if(model == "brown") {
    for (i in seq(length(phyList))) analysis[[i]] <- brown(data = structure(dat[phyList[[i]]@nodelabels], .Names = phyList[[i]]@nodes), tree = phyList[[i]])
    }
  else {
    for (i in seq(length(phyList))) analysis[[i]] <- hansen(data = structure(dat[phyList[[i]]@nodelabels], .Names = phyList[[i]]@nodes), tree = phyList[[i]], 
                                                            regimes = structure(rep(1, phyList[[i]]@nnodes), 
                                                                                names = phyList[[i]]@nodes, 
                                                                                levels = 1, class = 'factor'),
                                                            sigma = 1, sqrt.alpha = 1)
    } 
  names(analysis) <- treeNames
  params <- matrix(NA, nrow = length(c(treeNames, 'summed')), ncol = length(paramHeader), dimnames = list(c(treeNames, 'summed'), paramHeader))
  for (i in treeNames) {
    temp <- summary(analysis[[i]])
    params[i, ] <- unlist(temp[paramHeader])[paramHeader]
    if(model != "brown") params[i, 'optimum'] <- temp$optima[[1]]
    }
  params['summed', pSum] <- colSums(params[treeNames, pSum])
  out <- list(analysis = analysis, params = params)
  return(out)
  }
  