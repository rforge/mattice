# ---------------------------------------------------------------------
# FUNCTIONS FOR PERFORMING A SERIES OF OU ANALYSES ON A BATCH OF TREES
# ---------------------------------------------------------------------
## Changes needed:
## - measurement error portions need to be fixed
## - In a better world, allow graphical selection of subtrees to test on a single tree, then extract defining taxa
##    based on those nodes, using locator() or something like it.

runBatchHansen <-
# 11 nov 08: renamed to runBatchHansen
# Runs batchHansenFit and brown over a list of ouchTrees
# Arguments:
#  "ouchTrees" = list of OUCH-style trees
#  "characterStates" = vector of character states, either extracted from an ouch-style tree data.frame or a named vector
#  REMOVED: "SEM"= standard error of the mean, vector extracted from an ouch-style tree data.frame
#  "rescale" = factor to multiply against (times / max(times)) -- choose based on trial analyses; set at <= 0 if you don't want to rescale trees
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
#  "brown" = whether to analyse the data under a Brownian motion model
#  "..." = additional arguments to pass along to hansen
function(ouchTrees, characterStates, cladeMembersList, filePrefix = NULL, di = NULL, nodeNames = NULL, maxNodes = NULL, regimeTitles = NULL, brown = F, rescale = 1, alpha = 1, sigma = 1, ...) {
  ## do all the objects in ouchTrees inherit ouchtree?
  if(is(ouchTrees,'ouchtree')) ouchTrees <- list(ouchTrees)
  treeCheck <- unlist(lapply(ouchTrees, function(x) is(x,'ouchtree')))
  if(F %in% treeCheck) 
        stop(paste('This function has been rewritten to use the new S4 ', sQuote('ouchtree'), ' class.',
	'\nYou can generate a tree of this class by calling ', sQuote('ouchtree()'), '.', sep = ""))
  
  ## Check character states to make sure that they are either named and match names in the trees, or are the same length as the tips
  for (i in 1:length(ouchTrees)) {
    dataFlag <- NULL
    stopFlag <- F
    tree <- ouchTrees[[i]]
    terminals <- tree@nodelabels[(tree@nnodes - tree@nterm + 1):tree@nnodes]
    if(any(F %in% (terminals %in% names(characterStates)))) {
      message(paste("Not every terminal branch in tree", i, "has a corresponding name in", sQuote("characterStates")))
      if(length(characterStates) == tree@nterm) {
        message("Data assumed to be in the same order as terminals")
        dataFlag <- 'sameOrderTerminals' 
        }
      if(length(characterStates) == tree@nnodes) {
        message("Data assumed to be in the same order as nodes;\nany data not associated with a terminal branch will be ignored")
        dataFlag <- 'sameOrderNodes'
        }
      if(identical(dataFlag, NULL)) stopFlag <- T
      message("-------------------\n")
      }
    else dataFlag <- 'named'
    if(stopFlag) stop("Correct discrepancies between trees and data and try again!")
    }
  if(!identical(di, NULL)) dir.create(di)
  ar = regimeVectors(ouchTrees, cladeMembersList, maxNodes)
  hansenBatch <- list(length(ouchTrees))
  thetas <- list(length(ouchTrees))
  for (i in 1:length(ouchTrees)) {
    fP <- NULL
    if(!identical(filePrefix, NULL)) fP <- paste(filePrefix, ".t", i, ".", sep = "")
    if(!identical(di, NULL)) fP <- paste(di, "/", fP, sep = "")
    tree <- ouchTrees[[i]]
    if(identical(regimeTitles, NULL)) {
      regimeTitles <- as.character(1:length(ar$regList[[i]]))
      if(brown) regimeTitles <- c(regimeTitles, 'brown')
      }
    
    ## rescale tree if requested
    if(rescale>0) tree@times <- rescale * tree@times / max(tree@times) 
    
    ## make sure data fits the tree
    dataIn <- NULL
    if(dataFlag == 'sameOrderTerminals') dataIn <- c(rep(NA, tree@nnodes - tree@nterm), characterStates)
    if(dataFlag == 'sameOrderNodes') dataIn <- characterStates
    if(dataFlag == 'named') dataIn <- characterStates[match(tree@nodelabels, names(characterStates))]
    if(identical(dataIn, NULL)) stop(paste("There is a problem with your data that I failed to catch at the outset of", sQuote('runBatchHansen()')))
    else names(dataIn) <- tree@nodes
    
    ## send it off to batchHansen and just stick the results in hansenBatch... this won't work as the number of regimes gets large, 
    ##   so there should be some option here to just hang onto the coefficients for each run (i.e., hang onto 'coef(hansen(...))' rather than 'hansen(...)')
    ##   there could also be an option to save the entire object as a series of files in addition to hanging onto 
    hb <- batchHansen(tree, dataIn, ar$regList[[i]], regimeTitles, brown, fP, alpha, sigma, ...)
    hansenBatch[[i]] <- hb$treeData
    thetas[[i]] <- hb$thetas
    message(paste("Tree",i,"of",length(ouchTrees),"complete", "\n-----------------------------"))
  }
    
    ## right now no summary is returned; perhaps one is needed, summarizing over trees what is summarized for each tree in batchHansen
  outdata <- list(hansens = hansenBatch, thetas = thetas, regList = ar$regList, regMatrix = ar$regMatrix, nodeMatrix = ar$nodeMatrix, brown = brown, N = ouchTrees[[i]]@nterm, analysisDate = date(), call = match.call())
  class(outdata) <- 'hansenBatch'
  return(outdata)}

batchHansen <-
# Runs hansen and brown on a tree over a batch of selective regimes
# Arguments:
#  "tree" = the standard OUCH-style (S4) tree
#  "regimesList" = list of regime-paintings as output from regimeVectors
#  "scalingFactor" = factor to multiply against (times / max(times)) -- choose based on trial analyses
# Value: a matrix with nrow = regimes (+ 1 if brownian model is included) and columns for u, d.f., all estimated parameters, LRvsBM, AIC, and AIC weight
function(tree, data, regimesList, regimeTitles, brown, filePrefix = NULL, alpha, sigma, ...) {
  n <- tree@nterm
  ## set up a matrix that returns lnL, K, sigmasq, theta0, and alpha for every model
  ## thetas go into a models-by-branch matrix
  hansenOptima <- list(length(regimeTitles))
  variables <- c("loglik", "dof", "sigma.squared", "theta / alpha") # only display variables... set the selecting variables in the next two lines
  brVars <- c("loglik", "dof", "sigma.squared", "theta")
  haVars <- c("loglik", "dof", "sigma.squared", "alpha")
  if(brown) thetaModels <- regimeTitles[1: (length(regimeTitles) - 1)]
  else thetaModels <- regimeTitles
  thetas <- matrix(NA, 
                   nrow = length(thetaModels), 
                   ncol = tree@nnodes, 
                   dimnames = list(thetaModels, tree@nodes))
  treeData <- matrix(data = NA, nrow = length(regimeTitles), ncol = length(variables), dimnames = list(regimeTitles,variables))
  if(brown) {
    br <- brown(data, tree)
    if(!identical(filePrefix, NULL)) save(br, file = paste(filePrefix, 'b.Rdata', sep = ""))
    treeData["brown", ] <- unlist(summary(br)[brVars])
    }
  for (i in seq(regimesList)) {
    if(any(is.na(regimesList[[i]]))) {
      message(paste("skipping regime", i))
      treeData[i, ] <- rep(NA, dim(treeData)[[2]])
      }
    else {
      message(paste("Running regime",i))
      ## at this point, the user has to give an initial alpha and sigma for hansen to search on... this should be relaxed
      ha = hansen(data = data, tree = tree, regimes = regimesList[[i]], alpha = alpha, sigma = sigma, ...)
      treeData[i, ] <- unlist(summary(ha)[haVars])
      thetas[i, ] <- ha@theta$data[ha@regimes[[1]]]
      if(!identical(filePrefix, NULL)) save(ha, file = paste(filePrefix, 'r', i, '.Rdata', sep = ""))
      }
  }
  outdata <- list(treeData = treeData, thetas = thetas)
  return(outdata) }