multiModels <- function(phy, dat, node, models = c('whole.brown', 'whole.ou1', 'whole.ou2','part.brown', 'part.ou'), ic = "BIC") {
# test the support for alternative models on simple and partitioned trees
# currently only works on one tree; fix so it runs on a set of trees
  parameters <- c('loglik', 'dof', 'alpha', 'sigma.squared', 'theta', 'optima.dat.1', 'optima.dat.n')
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
   # as written, the following won't correctly chunk parameters into the matrix... mismatch between partitioned- and whole-tree tests
   outMatrix[model,  ] <- ifelse(i[1] == 'whole', 
                                 wholeModel(wholeTree, dat, i[2], node, parameters)$params,
                                 partialModel(partialTree, dat, i[2], c('uptree', 'downtree'), parameters)$params
                                 )
   }
  MAKE A WEIGHTS COLUMN using informationCriterion
  return(outMatrix)
  }

wholeModel <- function(phy, dat, model, nodes, parameterVector) {
  if(model == "brown") analysis <- brown(dat, phy)
  if(model == "ou1") analysis <- hansen(dat, phy, 
                                        regimes = structure(rep(1, phy@nnodes), 
                                        names = phy@nodes, levels = 1, 
                                        class = 'factor'),
                                        alpha = 1, sigma = 1)
  if(model = "ou2") analysis <- hansen(dat, phy, regimes = paintBranches(nodes, phy), alpha = 1, sigma = 1)
  params <- unlist(summary(analysis)[parameterVector])[parameterVector]
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
  analysis <- ifelse(model == "brown",
                     lapply(phyList, brown, data = dat),
                     lapply(phyList, hansen, data = dat, 
                            regimes = structure(rep(1, phyList[[1]]@nnodes), 
                                                names = phyList[[1]]@nodes, 
                                                levels = 1, class = 'factor'),       ## won't work... multiple subtrees, different regimes!
                            sigma = 1, alpha = 1)
                     )
  params <- lapply(analysis, function(x) {unlist(summary(x)[pV])[pV]}, pV = allParams)
  rawMat <- matrix(unlist(params), nrow = length(params), ncol = length(allParams), dimnames = list(treeNames, allParams))
  params <- colSums(rawMat, na.rm = T)
  out <- list(analysis = analysis, rawMat = rawMat, params = params)
  return(out)
  }
  