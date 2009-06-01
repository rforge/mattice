ouSim.ouchtree <- function(object, rootState = 0, alpha = 0, variance = 1, theta = rootState, steps = 1000, ...) {
## function to plot a simulated dataset under brownian motion or Ornstein-Uhlenbeck (OU) model
## Arguments:
##   object is an ouch-style (S4) tree
##   alpha and theta are either single values or vectors of length (length(branchList))
tree <- object
message(paste("running sim with root =", rootState, ", alpha =", mean(alpha), ", var =", variance, "theta =", mean(theta)))

	## embedded function---------------------
	## could be released to the wild, but more arguments will need to be passed around
	preorderOU <- function(branchList, tree, startNode, startState, alpha, theta) {
	  ## Recursive function to generate the data under a Brownian motion or OU model
	  ## modified for ouchtree (s4) Dec 08
	  ## branch times back from each tip are in tree@epochs, indexed by tip number
	  ## nodes back from each node (including tips) are in tree@lineages, indexed by node number
	  ## a branch length is the time of a node - the time of its ancestor
	  startBranch <- startNode ## startNode really means start branch... it's the end node of the branch starting this process
	  workingBranch <- branchList[[startBranch]]
	  daughterBranches <- tree@nodes[which(tree@ancestors == startNode)]
	  if(length(workingBranch) == 0) endState <- startState
	  else {
	    for (brStep in 1:length(workingBranch)) {
	      workingBranch[brStep] <- 
	        startState + workingBranch[brStep] + alpha[startBranch] / steps * (theta[startBranch] - startState) # denom was mult'd by steps... should be? 
	      startState <- workingBranch[brStep] 
	      }
	    branchList[[startBranch]] <- workingBranch
	    endState <- branchList[[startBranch]][length(branchList[[startBranch]])]
	    }	  
	  if(!identical(as.integer(daughterBranches), integer(0))) {
	    for(i in daughterBranches) branchList <- preorderOU(branchList, tree, i, endState, alpha, theta) } 
	  return(branchList) 
	}  
	## --------------------------------------

 #debug(preorderOU)

  ## 1. initialize
  if(length(alpha) == 1) alpha <- rep(alpha, tree@nnodes)
  if(length(theta) == 1) theta <- rep(theta, tree@nnodes)
  brLengths <- c(0, unlist(lapply(2:tree@nnodes, branchLength, tree = tree))) # assumes first node is root; this should be relaxed
  names(brLengths) <- tree@nodes # branches are indexed by end node
  names(alpha) <- tree@nodes
  names(theta) <- tree@nodes

  ## 2. The following creates a list of random draws from the normal distribution, with standard deviation scaled by total 
  ##   tree length and the number of draws for each branch equal to the number of steps in that branch. 
  ##   If there is a separate variance for each branch, I assume the variance is expressed in terms of total tree-length, 
  ##   so the scaling is the same for all branches
  sdStep <- sqrt(variance / steps)
  if(length(variance) == 1) branchList <- lapply(round(brLengths*steps), rnorm, sd = sdStep)
  else {
    branchList <- vector("list", length(tree@nnodes)) 
    for(i in seq(tree@nnodes)) branchList[[i]] <- rnorm(n = brLengths[i], mean = 0, sd = sdStep[i])
    }
  names(branchList) <- tree@nodes # branches are indexed by their end nodes
  timesList <- lapply(lapply(branchList, length), seq) # creates a list of sequential numbers corresponding to the rnorm vectors
  for(i in which(sapply(branchList, length) == 0)) timesList[[i]] <- numeric(0) # added b/c seq(numeric(0)) == 1:0
  startTimes <- round(unlist(sapply(tree@nodes, ancestorTime, tree = tree)) * steps) ## ASSUMES ULTRAMETRICITY
  for (i in 1:length(timesList)) timesList[[i]] <- timesList[[i]] + startTimes[i]

  ## 3. traverse
  for(i in which(tree@ancestors == tree@root)) { ## calls preorderOU for each descendent from the root.
    branchList <- preorderOU(branchList, tree, tree@nodes[i], rootState, alpha, theta) 
    }
		
  value <- list(branchList = branchList, timesList = timesList, steps = steps, parameters = list(rootState = rootState, alpha = alpha, variance = variance, theta = theta))
  class(value) <- "ouSim"
  return(value)
}

nodeTime <- function(node, tree) {
  return(tree@times[which(tree@nodes == node)])
}

branchLength <- function(endNode, tree) {
  return(nodeTime(endNode, tree) - ancestorTime(endNode, tree))
}

ancestorTime <- function(endNode, tree) {
  endNode <- as.character(endNode)
  ancestor <- tree@ancestors[which(tree@nodes == endNode)]
  if(is.na(ancestor)) startTime <- NA
  else startTime <- tree@times[which(tree@nodes == ancestor)]
  return(startTime)
}
