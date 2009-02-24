multiModel <- function(phy, dat, node, models = c('whole.brown', 'whole.ou1', 'whole.ou2','part.brown', 'part.ou'), ic = "BIC") {
# test the support for alternative models on simple and partitioned trees
# currently only works on one tree; fix so it runs on a set of trees, conditioned on those trees that have the node of interest;
# return percent of trees possessing that node as an additional value
  paramHeader <- c('AICc', 'BIC', 'AICc.weight', 'BIC.weight', 'loglik', 'dof', 'sigma.squared', 'alpha', 'theta', 'optimum', 'optimum.uptree', 'optimum.downtree')
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
  modelsSplit <- strsplit(models, ".", fixed = T)

  for(i in modelsSplit) {
   model <- paste(i[1], ".", i[2], sep = "")
   if (i[1] == 'whole') outMatrix[model,  ] <- wholeModel(wholeTree, dat, i[2], node, paramSets[[i[2]]], paramHeader)$params
   if (i[1] == 'part') outMatrix[c(paste(model, '.uptree', sep = ""), paste(model, '.downtree', sep = ""), paste(model, '.summed', sep = "")), ] <- partialModel(partialTree, dat, i[2], c('uptree', 'downtree'), paramSets[[i[2]]], pSum, paramHeader)$params                                 
   }
  # MAKE A WEIGHTS COLUMN using informationCriterion
  return(outMatrix)
  }

wholeModel <- function(phy, dat, model, node, parameterVector, paramHeader) {
  if(model == "brown") analysis <- brown(dat, phy)
  if(model == "ou1") analysis <- hansen(dat, 
  					phy, 
                                        regimes = structure(rep(phy@root, phy@nnodes), names = phy@nodes, levels = 1, class = 'factor'),
                                        alpha = 1, 
                                        sigma = 1
                                        )
  if(model == "ou2") {
    regime <- paintBranches(list(node), phy)
    uptreeNum <- as.character(phy@root)
    downtreeNum <- as.character(unique(regime))[unique(regime) != phy@root]
    analysis <- hansen(dat, phy, regime, alpha = 1, sigma = 1)
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
#  for(i in parameterVector) {
#   if(i %in% pSum) next
#   else allParams <- c(allParams, paste(treeNames, i, sep = ""))
#   }
  if(model == "brown") analysis <- lapply(phyList, brown, data = dat)
  else {
    analysis <- vector('list', length(phyList))
    for (i in seq(length(phyList))) analysis[[i]] <- hansen(data = dat, tree = phyList[[i]], 
                                                            regimes = structure(rep(1, phyList[[i]]@nnodes), 
                                                                                names = phyList[[i]]@nodes, 
                                                                                levels = 1, class = 'factor'),
                                                            sigma = 1, alpha = 1)
    } # close else
  names(analysis) <- treeNames
  params <- matrix(NA, nrow = length(c(treeNames, 'summed')), ncol = length(paramHeader), dimnames = list(c(treeNames, 'summed'), paramHeader))
  for (i in treeNames) {
    temp <- summary(analysis)
    params[i, ] <- unlist(temp[paramHeader])[paramHeader]
    if(model != "brown") params[i, 'optimum'] <- summary(temp)$optima[[1]]
    }
  params['summed', pSum] <- colSums(params[treeNames, pSum])
  out <- list(analysis = analysis, params = params)
  return(out)
  }
  