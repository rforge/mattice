# ---------------------------------------------------------------------
# FUNCTIONS FOR PERFORMING A SERIES OF OU ANALYSES ON A BATCH OF TREES
# ---------------------------------------------------------------------
# Copied from functions used in Hipp 2007 Evolution paper
# *** To run a hansen (OU) analysis, call runBatchHansenFit ***
# Utilizes ouch v 1.2-4
# This is the original set of functions utilized in Hipp 2007 (Evolution 61: 2175-2194) as modified
#  for Lumbsch, Hipp et al. 2008 (BMC Evolutionary Biology 8: 257). At the time of uploade to r-forge, 
#  they have not been checked for compatibility with subsequent versions of ouch.


## Project details
## maticce: Mapping Transitions In Continuous Character Evolution
## full project name: Continuous Character Shifts on Phylogenies
## unix name: mattice
## Project page requested from R-forge on 7 november 2008

## Changes needed:
## 1. calls should be to hansen rather than hansen.fit
## 2. measurement error portions need to be fixed
## 3. Analysis should be conducted over multiple trees, summarizing only over trees for which a given node is present;
##    node presence should be checked on each tree by looking to see whether the defining group is monophyletic,
##    and probably a matrix created for each multiple-tree analysis that makes summarizing quicker.
## 4. DONE -- Max number of simultaneous nodes should be set 
## 5. In a better world, allow graphical selection of subtrees to test on a single tree, then extract defining taxa
##    based on those nodes, using locator() or something like it.
## 6. IT statistics should use informationCriterion or something else to clean up the code

runBatchHansen <-
# 11 nov 08: renamed to runBatchHansen
# Runs batchHansenFit and brown.fit over a list of ouchTrees
# Arguments:
#  "ouchTrees" = list of OUCH-style trees; if a single tree, send as 'list(TREE)'
#  "characterStates" = vector of character states, either extracted from an ouch-style tree data.frame or a named vector
#  REMOVED: "SEM"= standard error of the mean, vector extracted from an ouch-style tree data.frame
#  "rescale" = factor to multiply against (times / max(times)) -- choose based on trial analyses; set at <= 0 if you don't want to rescale trees
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
#  "brown" = whether to analyse the data under a Brownian motion model
#  "..." = additional arguments to pass along to hansen
function(ouchTrees, characterStates, cladeMembersList, maxNodes = NULL, regimeTitles = NULL, logData = FALSE, brown = F, rescale = 1, ...) {
  ## do all the objects in ouchTrees inherit ouchtree?
  treeCheck <- unlist(lapply(ouchTrees, function(x) is(x,'ouchtree')))
  if(F %in% treeCheck) 
        stop(paste('This function has been rewritten to use the new S4 ', sQuote('ouchtree'), ' class.',
	'\nYou can generate a tree of this class by calling ', sQuote('ouchtree()'), '.', sep = ""))
  
  ## Check character states to make sure that they are either named and match names in the trees, or are the same length as the tips
  for (i in 1:length(ouchTrees)) {
    stopFlag <- F
    tree <- ouchTrees[[i]]
    terminals <- tree@nodelabels[(tree@nnodes - tree@nterm + 1):tree@nnodes]
    if(any(F %in% (terminals %in% names(characterStates)))) {
      message(paste("Not every terminal branch in tree", i, "has a corresponding name in", sQuote("characterStates")))
      if(length(characterStates) == tree@nterm) {
        message("Data assumed to be in the same order as terminals")
        dataFlag <- 'sameOrderTerminals' 
        }
      else if (length(characterStates) == tree@nnodes) {
        message("Data assumed to be in the same order as nodes;\nany data not associated with a terminal branch will be ignored")
        dataFlag <- 'sameOrderNodes'
        }
      else stopFlag <- T}
      message("-------------------\n")
      }
    else dataFlag <- 'named'
  if(stopFlag) stop("Correct discrepancies between trees and data and try again!")
  stopFlag <- F  

  if(logData) characterStates <- log(characterStates)
  if(identical(regimeTitles, NULL)) {
    regimeTitles <- as.character(1:length(regimesList))
    if(brown) regimeTitles <- c(regimeTitles, 'brown')
  }
  hansenBatch = vector("list", length(ouchTrees))
  treeCounter = 0
  for (i in 1:length(ouchTrees)) {
    tree <- ouchTrees[[i]]
    regimesList = regimeVectors(tree, cladeMembersList, maxNodes)
    
    ## rescale tree if requested
    if(rescale>0) tree@times <- rescale * tree@times / max(tree@times) 
    
    ## need to revise regimeVectors so that it only returns regimes for nodes that are supported in the tree
    ## for now, assume nodes of interest are present in all trees

    ## make sure data fits the tree
    dataIn <- NULL
    if(dataFlag) == 'sameOrderTerminals' dataIn <- c(rep(NA, tree@nnodes - tree@nterm), characterStates)
    if(dataFlag) == 'sameOrderNodes' dataIn <- characterStates
    if(dataFlag) == 'named' dataIn <- characterStates[match(tree@nodelabels), characterStates]
    if(identical(dataIn, NULL)) stop(paste("There is a problem with your data that I failed to catch at the outset of", sQuote('runBatchHansen()')))
    
    ## send it off to batchHansen and just stick the results in hansenBatch... this won't work as the number of regimes gets large, 
    ##   so there should be some option here to just hang onto the coefficients for each run (i.e., hang onto 'coef(hansen(...))' rather than 'hansen(...)')
    ##   there could also be an option to save the entire object as a series of files in addition to hanging onto 
    hansenBatch[[i]] = batchHansen(tree, dataIn, regimesList, regimeTitles, brown, ...)
    message(paste("Tree",i,"of",length(ouchTrees),"complete"))}
    
    ## right now no summary is returned; one is needed, summarizing over trees what is summarized for each tree in batchHansen
  
  return(list(hansens = hansenBatch, regimes = regimesList)) }

batchHansen <-
# Runs hansen.fit and brown.fit on a tree over a batch of selective regimes
# Arguments:
#  "node" "ancestor" "times" "data" = the standard tree specification vectors of the OUCH-style tree
#  "regimesList" = list of regime-paintings as output from regimeVectors
#  "scalingFactor" = factor to multiply against (times / max(times)) -- choose based on trial analyses
# Value: a matrix with nrow = regimes (+ 1 if brownian model is included) and columns for u, d.f., all estimated parameters, LRvsBM, AIC, and AIC weight
function(tree, data, regimesList, regimeTitles, brown, ...) {
  n <- tree@nterm
  ## set up a matrix that returns lnL, K, sigmasq, theta0, and alpha for every model; thetas will go along into a list that is indexed by model
  hansenOptima <- list(length(regimeTitles))
  variables <- c("loglik", "dof", "sigma.squared", "theta / alpha") # it's important that 'alpha' go last so that the matrix fills up right when the brownian motion model is used
  brVars <- c("loglik", "dof", "sigma.squared", "theta")
  haVars <- c("loglik", "dof", "sigma.squared", "alpha")
  treeData <- matrix(data = NA, nrow = length(regimeTitles), ncol = length(variables), dimnames = list(regimeTitles,variables))
  if(brown) br <- brown(data, tree)
  treeData["brown", ] <- unlist(summary(br)[brVars])
  for (i in seq(regimesList)) {
    message(paste("Running regime",i))
    ## at this point, the user has to give an initial alpha and sigma for hansen to search on... this should be relaxed
    ha = hansen(data, tree, regimesList[[i]], ...)
    treeData[i, ] <- unlist(summary(ha)[haVars])
    hansenOptima[[i]] <- summary(ha)$optima[[1]]
  }
  return(treeData) }