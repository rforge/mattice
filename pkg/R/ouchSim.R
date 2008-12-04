ouchSim.ouchtree <- function(tree, rootState = 0, shiftBranches = NULL, shiftStates = NULL, alpha = 0, variance = 1, theta = rootState, steps = 1000) {
## function to plot a simulated dataset under brownian motion or Ornstein-Uhlenbeck (OU) model
## Arguments:
##   tree is an ouch-style (S4) tree
##   alpha and theta are either single values or vectors of length (length(branchList))
##   shiftBranches is a vector indicating any branches at which an OU or brownian motion model has a determined shift in ancestral state
##   shiftStates is a vector of length = length(shiftBranches) indicaing the ancestral states for the determined break points

## 1. initialize

if(length(alpha) == 1) alpha <- rep(alpha, tree@nnodes)
if(length(theta) == 1) theta <- rep(theta, tree@nnodes)
## The following creates a list of random draws from the normal distribution, with standard deviation scaled by total tree length and the number of draws for each branch equal to the number of steps in that branch. If there is a separate variance for each branch, I assume the variance is expressed in tree-length units, not branch-length units, so the scaling is the same for all branches (viz., sd = sqrt(variance / steps))
if(length(variance) == 1) branchList <- lapply(round(phy$edge.length*steps), rnorm, sd = sqrt(variance/steps))
else {
  branchList <- vector("list", length(phy$edge.length))
  brLengths <- round(phy$edge.length*steps)
  brSD <- sqrt(variance/steps)
  for(i in seq(length(brLengths))) branchList[[i]] <- rnorm(n = brLengths[i], mean = 0, sd = brSD[i])
  }
timesList <- lapply(branchList, seq)
rootNode = length(phy$tip.label) + 1

startTimes <- max(branching.times(phy)) - branching.times(phy) ## WORKS, BUT ASSUMES ULTRAMETRICITY
startTimes <- round(startTimes[match(phy$edge[,1], names(startTimes))] * steps)
for (i in 1:length(timesList)) timesList[[i]] <- timesList[[i]] + startTimes[i]

## 3. traverse
for(i in which(phy$edge[, 1] == rootNode)) {
  branchList <- preorderOU(branchList, phy, phy$edge[i,2], rootState, alpha, theta) 
  }

value <- (list(branchList = branchList, timesList = timesList, steps = steps, parameters = list(rootState = rootState, alpha = alpha, variance = variance, theta = theta))) 
class(value) <- "ouSim"
return(value)}

preorderOU <- function(branchList, phy, startNode, startState, alpha, theta) {
## Recursive function to generate the data under a Brownian motion or OU model
## modified for ouchtree (s4) Dec 08
## branch times back from each tip are in tree@epochs, indexed by tip number
## nodes back from each node (including tips) are in tree@lineages, indexed by node number
## a branch length is the time of a node - the time of its ancestor
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

branchLength <- function(tree, endNode) {
  endNode <- as.character(endNode)
  ancestor <- tree@ancestors[which(tree@nodes == endNode)]
  endTime <- ot@times[which(ot@nodes == endNode)]
  startTime <- ot@times[which(ot@nodes == ancestor)]
  return(endTime - startTime)
}

plot.ouSim <- function(ouSim, nodeColor = "blue", nodeDotSize = 1.4, colors = rep("black", length(ouSim$branchList)), ...) {
## Plot the data by calling plot(ouSimObject)
## To plot different clades, set the colors vector according to the branches in the original 
## only passes the ... along to lines
  branches = length(ouSim$branchList)
  plot(1:ouSim$steps, ylim = range(unlist(ouSim$branchList)), type = "n", ylab = "Trait value", xlab = "Time")
  for(i in 1:branches) lines(ouSim$timesList[[i]], ouSim$branchList[[i]], col = colors[i], ...)
  for(i in 1:branches) points(ouSim$timesList[[i]][1], ouSim$branchList[[i]][1], pch = 19, col = nodeColor, cex = nodeDotSize) }
