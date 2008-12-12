# ---------------------------------------------------
# FUNCTIONS FOR PAINTING REGIMES ON AN S4 OUCH TREE #
# ---------------------------------------------------
# Adapted from functions used in Hipp 2007 Evolution paper
# Initially written for ouch v 1.2-4
# updated to ouch >= 2.4-2 Nov 2008
# updated to accommodate multiple trees Nov 2008

regimeVectors <-
# This is the basic call to get the full range of regimes over a set of trees
# Generates the list of painted branches representing all possible selective regimes for OU analyses, taking as argument
# species vectors that describe the clades at the bases of which regimes are specified to change.
# Arguments:
#  "tree" = the standard tree specification vectors of the OUCH-style tree
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
#  "maxNodes" = maximum number of nodes at which regime is permitted to change
# Value: 
#  "regList" = list of vectors that can each be plugged directly into OU analysis as the "regimes" argument
#  "nodeMatrix" = matrix of trees (rows) by nodes (columns) indicating what nodes are present on which trees
# 19 nov 08: changing to accept a list of trees and trimmed down greatly
function(ouchTrees, cladeMembersList, maxNodes = NULL) {
  nnode <- length(cladeMembersList)
  regMatrix <- regimeMatrix(n = nnode, maxNodes = maxNodes)
  apr = regimeMaker(ouchTrees, regMatrix, cladeMembersList)
  outdata <- list(regList = apr$regList, regMatrix = apr$regMatrix, nodeMatrix = apr$nodeMatrix)
  return(outdata) 
}

paintBranches <-
# Paints branches with regimes changing at nodes specified
# arguments
#  "tree" = OUCH-style (S4) tree
#  "regimeShiftNodes" = a vector of nodes at which selective regimes shift: root must be included, but tips are meaningless in this context
#  "regimeTitles" = a vector of titles for the regimes that begin at the root and at the nodes indicated in "regimeShiftNodes",
#                   in order of description in "regimeShiftNodes", except that the root is listed first in "regimeTitles"
#                   but not at all in "regimeShiftNodes"... defaults to "node[x]regime
# Value: a vector of regimes that can be plopped right into an OUCH-style tree data frame
function(regimeShiftNodes, tree, regimeTitles = NULL) {
  ## ------------------ begin ouchtree block -----------------
  ## check to see if tree inherits 'ouchtree'
  if (!is(tree,'ouchtree')) 
	stop(paste('This function has been rewritten to use the new S4 ', sQuote('ouchtree'), ' class.',
	'\nYou can generate a tree of this class by calling ', sQuote('ouchtree()'), '.', sep = ""))
  ## get the vectors we need:
  ancestor <- tree@ancestors # class = "character"
  node <- tree@nodes # class = "character"
  species <- tree@nodelabels # class = "character" -- note that nodelabels is more general than this indicates and the name should be changed throughout at some point
  times <- tree@times # class = "numeric"
  ## ------------------ end ouchtree block -------------------
  
  if(identical(regimeTitles, NULL)) regimeTitles <- as.character(regimeShiftNodes)
  names(regimeTitles) = as.character(regimeShiftNodes)
  colorsVector = character(length(node))
  for (i in 1:length(ancestor)) {
    # First three lines fill up the vector for nodes that are hit in order
    if (is.na(ancestor[i])) {
      colorsVector[i] = regimeTitles["1"]
      next }
    if (as.character(ancestor[i]) %in% as.character(regimeShiftNodes)) {
      colorsVector[i] = regimeTitles[as.character(ancestor[i])]
      next }
    if (colorsVector[as.integer(ancestor[i])] != "") {
      colorsVector[i] = colorsVector[as.integer(ancestor[i])]
      next }
    # These lines fill up the vector for nodes run reached before their immediate ancestor
    nodeQ = integer(length(node))
    ii = i
    repeat {
      nodeQ = c(ii, nodeQ)
      ii = as.numeric(ancestor[ii])
      if (as.character(ancestor[ii]) %in% as.character(regimeShiftNodes)) {
        colorsVector[ii] = colorsVector[as.integer(ancestor[ii])]
        break}
      if (colorsVector[as.integer(ancestor[ii])] != "") {
        colorsVector[ii] = colorsVector[as.integer(ancestor[ii])]
        break} }
    
    for(j in nodeQ) {
      colorsVector[j] = colorsVector[as.integer(ancestor[j])] } 
      
      } # closes for(i in 1:length(ancestor)) loop 

      # a little hack to fix a problem I don't understand... with the undesired side effect that it colors the stem of some subtrees rather than the crown as originally written
      for(i in 1:length(colorsVector)) if(colorsVector[i] == "") colorsVector[i] <- as.character(i) 
      
      # colors terminal branches if any terminal branches are in the regimeShiftNodes
      for(i in regimeShiftNodes) if(i %in% tree@term) colorsVector[as.numeric(i)] <- as.character(i)
  return(colorsVector) }

regimeMaker <- function(ouchTrees, regMatrix, nodeMembers) {
## supplants the old 'allPossibleRegimes'
## takes a list of ouchtree objects, a regimeMatrix ouput, and a list of nodeMembers (the taxa definining each node of interest)
## Value:
##  regList = a list of nodes defining the change points for each tree (i.e., a list of lists)
##  nodeMatrix = a matrix of trees (rows) by nodes (columns) indicating whether the node is present in each tree
  
  # set up variables
  numTrees <- length(ouchTrees)
  numNodes <- length(nodeMembers)
  if(numNodes != dim(regMatrix)[2]) stop('Number of nodes (columns) in regMatrix must equal number of items in nodeMembers list')
  nodeMatrix <- matrix(NA, nrow = numTrees, ncol = numNodes, dimnames = list(seq(numTrees), dimnames(regMatrix)[[2]]))
  changeNodes <- list(numTrees)
  regList <- list(numTrees)
  regMatrices <- list(numTrees)
  
  # fill outdata
  for(i in seq(numNodes)) nodeMatrix[, i] <- unlist(lapply(ouchTrees, isMonophyletic, taxa = nodeMembers[[i]]))
  for(i in seq(numTrees)) {
    tree <- ouchTrees[[i]]
    regMatrices[[i]] <- regMatrix * as.numeric(matrix(nodeMatrix[i, ], dim(regMatrix)[1], dim(regMatrix)[2], byrow = T)) # multiplies regMatrix by nodes present
    regMatrices[[i]][1:(dim(regMatrices[[i]])[1] - 1), ][which(apply(regMatrices[[i]][1:(dim(regMatrices[[i]])[1] - 1), ], 1, sum) == 0), ] <- rep(NA, numNodes) # set to NA regimes that have no nodes, except for OU1 model
    regMatrices[[i]][duplicated(apply(regMatrices[[i]], 1, as.decimal)), ] <- rep(NA, numNodes) ## set to NA non-unique regimes
    dimnames(regMatrices[[i]]) <- list(seq(dim(regMatrices[[i]])[1]), dimnames(regMatrices[[i]])[[2]])
    numTreeRegs <- dim(regMatrices[[i]])[1]
    treeRegs <- list(numTreeRegs) # this will be assigned to regList[[i]]
    nodesVector <- unlist(lapply(nodeMembers, mrcaOUCH, tree = ouchTrees[[i]])) # as written, gets the MRCA for even invalid nodes just so indexing stays right
    for(j in seq(numTreeRegs)) {
      if(any(is.na(regMatrices[[i]][j, ]))) treeRegs[[j]] <- NA
      else {
        treeRegs[[j]] <- as.factor(paintBranches(c("1", nodesVector[as.logical(regMatrices[[i]][j, ])]), tree))
        names(treeRegs[[j]]) <- tree@nodes
      }
    }
    regList[[i]] <- treeRegs
  }
  regMatrices$overall <- regMatrix # this is the matrix that includes all regimes without regard to any tree
  outdata <- list(regList = regList, nodeMatrix = nodeMatrix, regMatrix = regMatrices)
  return(outdata)
}

regimeMatrix <- function(n = NULL, nodeNames = NULL, regimeNames = NULL, maxNodes = NULL) {
## a brute-force approach, very inefficient as n and maxNodes diverge
## I think there's a recursive approach that would be efficient for small maxNodes, but it would probably be slower
##   than this approach as maxNodes -> n
  if(identical(n, NULL) && identical(nodeNames, NULL)) stop("You have to give regimeMatrix the number of nodes, a vector of node names, or both")
  if(identical(nodeNames, NULL)) nodeNames <- as.character(seq(n))
  else n <- length(nodeNames)
  if(identical(maxNodes, NULL)) maxNodes <- n
  outmatrix <- matrix(NA, nrow = 0, ncol = n, dimnames = list(NULL, nodeNames))
  maxNumberOfRegimes <- ifelse(n == 1, 2, 2^n)
  counter <- 1
  repeat {
    temp <- as.binary(counter, digits = n)
    if(sum(temp) <= maxNodes) outmatrix <- rbind(outmatrix, temp) 
    if(counter == maxNumberOfRegimes - 1) break
    counter <- counter + 1
    }
  outmatrix <- rbind(outmatrix, as.binary(0, digits = n))
  dimnames(outmatrix)[[1]] <- seq(dim(outmatrix)[1])
  return(outmatrix)
}

#regMatRec <- function(n, maxNodes, dat) {
#  for (i in 1:n) {
#    place a 1 in the ith position
#    temp = regMatRec on the zeros, thus with n = n-1, maxNodes = maxNodes - 1
#    make a matrix of the results, with the temp and the ith-position 1 concatenated
#  }
#  make an all-zero row
#  return results
#}
    

as.binary <- function(n, base = 2, r = FALSE, digits = NULL)
# Robin Hankin <initialDOTsurname at soc.soton.ac.uk (edit in obvious way; spam precaution)>
# submitted to R listserv Thu Apr 15 12:27:39 CEST 2004
# AH added 'digits' to make it work with regimeMatrix
# https://stat.ethz.ch/pipermail/r-help/2004-April/049419.html
{
   out <- NULL
   while(n > 0) {
     if(r) {
       out <- c(out , n%%base)
     } else {
       out <- c(n%%base , out)
     }   
     n <- n %/% base
   }
   if(!identical(digits, NULL) && !r) out <- c(rep(0, digits-length(out)), out)
   if(!identical(digits, NULL) && r) out <- c(out, rep(0, digits-length(out)))
   return(out)
}

as.decimal <- function(n) {
# takes a binary vector and makes it a decimal
  digits <- length(n)
  result <- 0
  for(i in digits:1) result <- result + n[i] * 2 ^ (digits - i)
  result
}