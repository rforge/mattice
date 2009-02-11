multiModels <- function(phy, dat, node, models = c('whole.brown', 'whole.ou1', 'whole.ou2','part.brown', 'part.ou'), ic = "BIC") {
# test the support for alternative models on simple and partitioned trees
# currently only works on one tree; fix so it runs on a set of trees
  parameters <- c('loglik', 'dof', 'alpha', 'sigma.squared', 'theta', 'optima.nodeAncestor', 'optima.nodeDescendent')
  pSum <- c('loglik', 'dof') # parameters to sum for evaluating partitioned trees
  aboveNodeTaxa <- phy$tip.label[-which(phy$tip.label %in% node)]
  wholeTree <- ape2ouch(phy)
  upTree <- ape2ouch(drop.tip(phy, node))
  downTree <- ape2ouch(drop.tip(phy, aboveNodeTaxa))
  partialTree <- list(upTree = upTree, downTree = downTree)
  outMatrix <- matrix(NA, nrow = length(models), ncol = length(parameters), dimnames = list(models, parameters))
  modelsSplit <- strsplit(models, ".", fixed = T)
  for(i in modelsSplit) {
   model <= paste(i[1], ".", i[2], sep = "")
   # as written, the following won't correctly chunk parameters into the matrix... mismatch between partitioned- and whole-tree tests -- NOW MAY WORK... DOUBLE CHECK RELATIVE TO NAMES SPIT OUT IN WHOLE MODEL, AND CHEKC PARTIALMODEL NAMES
   outMatrix[model,  ] <- ifelse(i[1] == 'whole', 
                                 wholeModel(wholeTree, dat, i[2], node, parameters)$params,
                                 partialModel(partialTree, dat, i[2], c('uptree', 'downtree'), parameters)$params
                                 )
   }
  MAKE A WEIGHTS COLUMN using informationCriterion
  return(outMatrix)
  }

wholeModel <- function(phy, dat, model, node, parameterVector) {
  if(model == "brown") analysis <- brown(dat, phy)
  if(model == "ou1") analysis <- hansen(dat, phy, 
                                        regimes = structure(rep(1, phy@nnodes), 
                                        names = phy@nodes, levels = 1, 
                                        class = 'factor'),
                                        alpha = 1, sigma = 1)
  if(model = "ou2") {
    regime <- paintBranches(list(node), phy)
    ancNum <- as.character(tree@root)
    descNum <- as.character(unique(regime))[unique(regime) != tree@root]
    analysis <- hansen(dat, phy, regime, alpha = 1, sigma = 1)
    } # close if
  params <- unlist(summary(analysis)[parameterVector])[parameterVector]
  params['optima.nodeAncestor'] <- summary(analysis)$optima[[1]][ancNum]
  params['optima.nodeDescendent'] <- summary(analysis)$optima[[1]][descNum]
  out <- list(analysis = analysis, params = params)
  return(out)
  }

partialModel <- function(phyList, dat, model, treeNames, parameterVector = NULL, pSum = NULL) {
# sums a subset of parameters (indicated by pSum) and leaves the others separate
  allParams <- pSum
  for(i in parameterVector) {
   if(i %in% pSum) next
   else allParams <- c(allParams, paste(treeNames, i, sep = ""))
   }
  if(model == "brown") analysis <- lapply(phyList, brown, data = dat),
  else {
    analysis <- list(length(phyList))
    for (i in seq(length(phyList))) analysis[[i]] <- hansen(data = dat, tree = phyList[[i]], 
                                                            regimes = structure(rep(1, phyList[[i]]@nnodes), 
                                                                                names = phyList[[i]]@nodes, 
                                                                                levels = 1, class = 'factor'),
                                                            sigma = 1, alpha = 1)
    } # close else
  params <- lapply(analysis, function(x) {unlist(summary(x)[pV])[pV]}, pV = allParams)
  rawMat <- matrix(unlist(params), nrow = length(params), ncol = length(allParams), byrow = T, dimnames = list(treeNames, allParams))
  params <- colSums(rawMat, na.rm = T)
  out <- list(analysis = analysis, rawMat = rawMat, params = params)
  return(out)
  }
  