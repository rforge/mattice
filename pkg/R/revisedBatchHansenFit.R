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
## 4. Max number of simultaneous nodes should be set 
## 5. In a better world, allow graphical selection of subtrees to test on a single tree, then extract defining taxa
##    based on those nodes. This shouldn't be too tricky using locator() or something like it.
## 6. IT statistics should use informationCriterion or something else to clean up the code

runBatchHansen <-
# 11 nov 08: renamed to runBatchHansen
# Runs batchHansenFit and brown.fit over a list of ouchTrees
# Arguments:
#  "ouchTrees" = list of OUCH-style trees; if a single tree, send as 'list(TREE)'
#  "characterStates" = vector of character states, either extracted from an ouch-style tree data.frame or a named vector
#  REMOVED: "SEM"= standard error of the mean, vector extracted from an ouch-style tree data.frame
#  REMOVED: "scalingFactor" = factor to multiply against (times / max(times)) -- choose based on trial analyses
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
function(ouchTrees, characterStates, cladeMembersList, regimeNames = NULL, logData = FALSE) {
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
      if(length(characterStates) == tree@nterm) message("Data assumed to be in the same order as terminals")
      else if (length(characterStates) == tree@nnodes) message("Data assumed to be in the same order as nodes;\nany data not associated with a terminal branch will be ignored")
      else stopFlag <- T}
      message("-------------------\n") }
  if(stopFlag) stop("Correct discrepancies between trees and data and try again!")
  stopFlag <- F  

  if(logData) characterStates <- log(characterStates)
  if(identical(regimeNames, NULL)) regimeNames <- c(as.character(1:length(regimesList)), "brown")
  
  hansenBatch = vector("list", length(ouchTrees))
  treeCounter = 0
  for (i in 1:length(ouchTrees)) {
    tree <- ouchTrees[[i]]
    regimesList = regimeVectors(tree, cladeMembersList)
    
    ## need to revise regimeVectors so that it only returns regimes for nodes that are supported in the tree
    ## for now, assume nodes of interest are present in all trees

    hansenBatch[[i]] = batchHansen(tree, characterStates, regimesList, regimeNames, numberOfTermini = tree@nterm)
    message(paste("Tree",i,"of",length(ouchTrees),"complete"))}
    
    ## right now no summary is returned; it should be
  
  return(list(hansens = hansenBatch, regimes = regimesList)) }

batchHansen <-
# Runs hansen.fit and brown.fit on a tree over a batch of selective regimes
# Arguments:
#  "node" "ancestor" "times" "data" = the standard tree specification vectors of the OUCH-style tree
#  "regimesList" = list of regime-paintings as output from regimeVectors
#  "scalingFactor" = factor to multiply against (times / max(times)) -- choose based on trial analyses
# Value: a matrix with nrow = regimes + 1 and columns for u, d.f., all estimated parameters, LRvsBM, AIC, and AIC weight
# ADDED "error" on 2 march 07 to accomodate the "me" in Pienaar's ouch
function(tree, data, regimesList, regimeTitles, numberOfTermini) {
  variables = c("loglik", "df", "alpha", "sigma", "theta0", "theta1", "theta2", "theta3", "theta4", "theta5","theta6","theta7","theta8","theta9", "LRvsBM", "AIC", "AICweight", "deltaAIC", "exp(-0.5*deltaAIC)", "AICc", "AICcWeight", "deltaAICc", "exp(-0.5*deltaAICc)")
  treeData = matrix(data = NA, nrow = (length(regimesList) + 1), ncol = length(variables),
             dimnames = list(regimeTitles,variables))
  brownianResults = brown.fit(data, error, node, ancestor, (times/max(times))*scalingFactor)
  treeData[length(regimesList) + 1, "AIC"] = brownianResults$aic
  treeData[length(regimesList) + 1, "loglik"] = brownianResults$loglik
  treeData[length(regimesList) + 1, "sigma"] = brownianResults$sigma
  treeData[length(regimesList) + 1, "theta0"] = brownianResults$theta
  treeData[length(regimesList) + 1, "df"] = brownianResults$df
  treeData[length(regimesList) + 1, "AICc"] = brownianResults$aic + ((2 * brownianResults$df * (brownianResults$df + 1)) / (numberOfTermini - brownianResults$df - 1))
  for (i in 1:length(regimesList)) {
    message(paste("Running regime",i))
    tempHansen = hansen.fit(data, error, node, ancestor, (times/max(times))*scalingFactor, regimesList[[i]])
    treeData[i,"loglik"] = tempHansen$loglik
    treeData[i,"df"] = tempHansen$df
    treeData[i,"alpha"] = tempHansen$alpha
    treeData[i,"sigma"] = tempHansen$sigma
    treeData[i,"LRvsBM"] = 2 * ((tempHansen$loglik / 2) - (treeData[length(regimesList) + 1, "loglik"] / 2))
    treeData[i,"AIC"] = tempHansen$aic
    treeData[i,"AICc"] = tempHansen$aic + ((2 * tempHansen$df * (tempHansen$df + 1)) / (numberOfTermini - tempHansen$df - 1))
    for (j in 0:(length(names(tempHansen$theta))-1)) {
      treeData[i, paste("theta",j,sep = "")] = tempHansen$theta[[j+1]]}}
  for (i in 1:(length(regimesList) + 1)) {
    treeData[i, "deltaAIC"] = treeData[i, "AIC"] - min(treeData[1:(length(regimesList) + 1), "AIC"])
    treeData[i, "exp(-0.5*deltaAIC)"] = exp(-0.5 * treeData[i,"deltaAIC"])
    treeData[i, "deltaAICc"] = treeData[i, "AICc"] - min(treeData[1:(length(regimesList) + 1), "AICc"])
    treeData[i, "exp(-0.5*deltaAICc)"] = exp(-0.5 * treeData[i,"deltaAICc"]) }
  for (i in 1:(length(regimesList) + 1)) {
    treeData[i, "AICweight"] = treeData[i, "exp(-0.5*deltaAIC)"] / sum(treeData[1:(length(regimesList) + 1), "exp(-0.5*deltaAIC)"])
    treeData[i, "AICcWeight"] = treeData[i, "exp(-0.5*deltaAICc)"] / sum(treeData[1:(length(regimesList) + 1), "exp(-0.5*deltaAICc)"]) }
  return(treeData) }
