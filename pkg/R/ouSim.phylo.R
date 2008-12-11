ouSim.phylo <- function(phy, rootState = 0, shiftBranches = NULL, shiftStates = NULL, alpha = 0, variance = 1, theta = rootState, model = "OU", branchMeans = NULL, steps = 1000) {
## function to plot a simulated dataset under brownian motion or Ornstein-Uhlenbeck (OU) model
## Arguments:
##   phy is an ape-style tree
##   alpha and theta are either single values or vectors of length (length(branchList))
##   shiftBranches is a vector indicating any branches at which an OU or brownian motion model has a determined shift in ancestral state
##   shiftStates is a vector of length = length(shiftBranches) indicaing the ancestral states for the determined break points
## Models:
##  "OU" is a brownian motion or OU model 
##  "meanVar" is a model in which the only phylogenetic effect is the mean and variance for a given branch
## Andrew Hipp (ahipp@mortonarb.org), January 2008 
## July 2008: modified to accomodate a vector of alpha and theta corresponding to branches
## Dec 2008: This function I'm leaving as is for the time being and just letting the phylo method be as raw as always.
##           New developments will be in the ouchtree, brown, hansen, and hansenBatch methods

preorderOU <- function(branchList, phy, startNode, startState, alpha, theta) {
## Recursive function to generate the data under a Brownian motion or OU model (not needed in the Platt model)
  startBranch = which(phy$edge[,2] == startNode)
  if(!identical(shiftStates, NULL)) {
    if(startBranch %in% shiftBranches) startState <- shiftStates[match(startBranch, shiftBranches)] }
  message(paste('Working on branch',startBranch,'with starting state',startState))
  branchList[[startBranch]][1] <- startState
  for (i in 2:length(branchList[[startBranch]])) {
    branchList[[startBranch]][i] <- branchList[[startBranch]][i - 1] + branchList[[startBranch]][i] + alpha[startBranch] / steps * (theta[startBranch] - branchList[[startBranch]][i - 1]) }
  endState = branchList[[startBranch]][length(branchList[[startBranch]])]
  daughterBranches <- phy$edge[which(phy$edge[, 1] == startNode), 2]
  if(!identical(as.integer(daughterBranches), integer(0))) {
    for(i in daughterBranches) branchList <- preorderOU(branchList, phy, i, endState, alpha, theta) }
  return(branchList) }  

## 1. initialize

if(length(alpha) == 1) alpha <- rep(alpha, length(phy$edge.length))
if(length(theta) == 1) theta <- rep(theta, length(phy$edge.length))
## The following creates a list of random draws from the normal distribution, with standard deviation scaled by total tree length and the number of draws for each branch equal to the number of steps in that branch. If there is a separate variance for each branch, I assume the variance is expressed in tree-length units, not branch-length units, so the scaling is the same for all branches (viz., sd = sqrt(variance / steps))
if(model == "OU") {
	if(length(variance) == 1) branchList <- lapply(round(phy$edge.length*steps), rnorm, sd = sqrt(variance/steps))
	  else {
	    branchList <- vector("list", length(phy$edge.length))
	    brLengths <- round(phy$edge.length*steps)
	    brSD <- sqrt(variance/steps)
	    for(i in seq(length(brLengths))) branchList[[i]] <- rnorm(n = brLengths[i], mean = 0, sd = brSD[i]) }}
if(model == "meanVar") {
	if(length(variance) == 1) variance <- rep(variance, length(phy$edge.length))
	branchList <- vector("list", length(phy$edge.length))
	if(length(branchMeans) == 1) branchMeans <- rep(branchMeans, length(phy$edge))
	brLengths <- round(phy$edge.length*steps)
	brSD <- sqrt(variance)
	for(i in seq(length(brLengths))) branchList[[i]] <- rnorm(n = brLengths[i], mean = branchMeans[i], sd = brSD[i])
	}
timesList <- lapply(branchList, seq)
rootNode = length(phy$tip.label) + 1

startTimes <- max(branching.times(phy)) - branching.times(phy) ## WORKS, BUT ASSUMES ULTRAMETRICITY
startTimes <- round(startTimes[match(phy$edge[,1], names(startTimes))] * steps)
for (i in 1:length(timesList)) timesList[[i]] <- timesList[[i]] + startTimes[i]

## 3. traverse
if(model == "OU") {
	for(i in which(phy$edge[, 1] == rootNode)) {
	  branchList <- preorderOU(branchList, phy, phy$edge[i,2], rootState, alpha, theta) }}
if(model == "meanVar") branchList <- branchList

value <- (list(branchList = branchList, timesList = timesList, steps = steps, parameters = list(rootState = rootState, alpha = alpha, variance = variance, theta = theta))) 
class(value) <- "ouSimPhylo"
return(value)
}