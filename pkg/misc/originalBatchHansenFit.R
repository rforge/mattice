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

runBatchHansenFit <-
# Runs batchHansenFit and brown.fit over a list of ouchTrees
# Arguments:
#  "ouchTrees" = list of OUCH-style trees; if a single tree, send as 'list(TREE)'
#  "characterStates" = vector of character states, extracted from an ouch-style tree data.frame
#  "SEM"= standard error of the mean, vector extracted from an ouch-style tree data.frame
#  "scalingFactor" = factor to multiply against (times / max(times)) -- choose based on trial analyses
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
function(ouchTrees, characterStates, cladeMembersList, SEM = rep(0, length(characterStates)), scalingFactor = 1, logData = FALSE) {
  hansenBatch = vector("list", length(ouchTrees))
  treeCounter = 0
  for (i in ouchTrees) {
    treeCounter = treeCounter + 1
    i = list(node = as.character(i$node), ancestor = as.character(i$ancestor), times = as.numeric(i$times), species = as.character(i$species))
    regimesList = regimeVectors(i$node, i$ancestor, i$times, i$species, cladeMembersList)
    write.table(regimesList, 'regimesList.txt')
    if (logData) hansenBatch[[treeCounter]] = batchHansenFit(i$node, i$ancestor, i$times, characterStates, error = SEM, scalingFactor, regimesList, c(as.character(1:length(regimesList)), "brown"), numberOfTermini = sum(!is.na(i$species)))
      else hansenBatch[[treeCounter]] = batchHansenFit(i$node, i$ancestor, i$times, characterStates, error = SEM, scalingFactor, regimesList, c(as.character(1:length(regimesList)), "brown"), numberOfTermini = sum(!is.na(i$species)))
    message(paste("Tree",treeCounter,"of",length(ouchTrees),"complete"))}
  return(list(hansens = hansenBatch, regimes = regimesList)) }

batchHansenFit <-
# Runs hansen.fit and brown.fit on a tree over a batch of selective regimes
# Arguments:
#  "node" "ancestor" "times" "data" = the standard tree specification vectors of the OUCH-style tree
#  "regimesList" = list of regime-paintings as output from regimeVectors
#  "scalingFactor" = factor to multiply against (times / max(times)) -- choose based on trial analyses
# Value: a matrix with nrow = regimes + 1 and columns for u, d.f., all estimated parameters, LRvsBM, AIC, and AIC weight
# ADDED "error" on 2 march 07 to accomodate the "me" in Pienaar's ouch
function(node, ancestor, times, data, error, scalingFactor = 1, regimesList, regimeTitles, numberOfTermini = 53) {
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
