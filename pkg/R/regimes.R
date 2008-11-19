# ---------------------------------------------------
# FUNCTIONS FOR PAINTING REGIMES ON AN S4 OUCH TREE #
# ---------------------------------------------------

# Modified from functions used in Hipp 2007 Evolution paper
# Initially written for ouch v 1.2-4
# updated to ouch >= 2.4-2 Nov 2008

# FINISH REGIME MAKER


regimeVectors <-
# Generates the list of painted branches representing all possible selective regimes for OU analyses, taking as argument
# species vectors that describe the clades at the bases of which regimes are specified to change.
# Arguments:
#  "tree" = the standard tree specification vectors of the OUCH-style tree
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
# Value: list of vectors that can each be plugged directly into OU analysis as the "regimes" argument
# 19 nov 08: changing to accept a list of trees

function(ouchTrees, cladeMembersList, maxNodes = NULL) {
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
  
  nnode <- length(cladeMembersList)
  regMatrix <- regimeMatrix(n = nnode)
  apr = regimeMaker(ouchTrees, regMatrix, cladeMembersList) ## HOLD IT! NOW REGIME MAKER WORKS ON ALL TREES AT ONCE... RETHINK THIS
  # apr$regList = a list of vectors, each indicating changeNodes
  # apr$nodeMatrix = a matrix of trees (rows) by nodes (columns) indicating whether the node is present in each tree
  nodeMatri <- unlist(changeNodesList)
  #changeNodesVector = vector("character", length(changeNodesList))
  #for (i in 1:length(changeNodesList)) # Changing cladeMemberList to a 1-d vector
  #  {changeNodesVector[i] = changeNodesList[[i]]}
  allRegimes <- regimesList
  regimePaintings = vector("list", length(allRegimes))
  for (i in 1:length(allRegimes)) {
    allRegimes[[i]] <- c("1", allRegimes[[i]])
    regimePaintings[[i]] <- as.factor(paintBranches(tree, allRegimes[[i]], as.character(allRegimes[[i]])))
    names(regimePaintings[[i]]) <- tree@nodes
    message(paste('Created regime',i))}
  outdata <- list(regimeList = regimePaintings, regimeMatrix = regMatrix)
  return(outdata) }

paintBranches <-
# Paints branches with regimes changing at nodes specified
# arguments
#  "node" "ancestor" "times" = the standard tree specification vectors of the OUCH-style tree
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
  return(colorsVector) }

allPossibleRegimes <-
# Generates a list of vectors of all possible 2^n regimes for a given list of n ancestor nodes, assuming that each node
#  can either be a change node or not, and that all combinations are meaningful.
# NOTE: This routine returns all possible permutations, ignoring the fact that the root must be included as a changeNode
#  for paintBranches to work properly. Exclude the root when calling this routine.
# Arguments:
#  "changeNodes" = vector of all change nodes in all possible scenarios, or number of regimes assumed if nodeMatrix = T
#  "maxNodes" = single number indicating the maximum number of nodes at which a regime can change
#  "nodeMatrix" = indicates whether to return a binary table for interp or a list of changeNode vectors for analysis
# Value:
#    regimeList = a list of changeNode vectors (assumes type = "character"), one for each possible scenario
#    regimeMatrix = : a binary table indicating whether a regime node is present or absent based on allPossibleRegimes output; 
#                     presumes nodes are labeled 1:n
# 10 nov 08: this function now takes over regimeNodes
function(changeNodes, maxNodes = NULL, nodeMatrix = F) {
    if(!identical(maxNodes, NULL) && maxNodes > length(changeNodes)) warning(paste(sQuote('maxNodes'), 'cannot be larger than the number of nodes; maxNodes ignored'))
    numberOfRegimes = ifelse(length(changeNodes) == 1, 2, 2^length(changeNodes))
    regime = vector("list", numberOfRegimes)
    for (i in 1:(numberOfRegimes - 1)) {
      remainder = i
      n = NULL
      for (j in as.integer(log2(i)):0) {
        if (2^j > remainder) n[j+1] = NA
        else {
          n[j+1] = changeNodes[j+1]
          remainder = remainder %% 2^j 
        }
      }
      regime[[i]] = sort(n[!is.na(n)]) 
    }
    regime[[numberOfRegimes]] = rep("0", times = as.integer(log2(i)) + 1) 

    ## make node matrix
    regimesNameMatrix = matrix(
        data = NA, nrow = numberOfRegimes, ncol = length(changeNodes), dimnames = list(
          as.character(seq(numberOfRegimes)), as.character(seq(length(changeNodes)))
          )
      )
      for (i in seq(numberOfRegimes)) {
        for (j in seq(length(changeNodes))) {
          if (is.na(match(changeNodes[j],regime[[i]]))) regimesNameMatrix[i,j] = 0 # changed this so that j indexes changeNodes
          else regimesNameMatrix[i,j] = 1 
        }
      }
      outmatrix <- regimesNameMatrix
      if(!identical(maxNodes, NULL)) {
        outmatrix <- outmatrix[apply(outmatrix,1,sum) <= maxNodes, ]
        dimnames(outmatrix)[[1]] = as.character(seq(dim(outmatrix)[1]))
        }

    ## prune list
    outlist <- regime 
      if(!identical(maxNodes, NULL)) {
        outlist <- outlist[sapply(outlist, length) <= maxNodes]
        outlist[[length(outlist) + 1]] <- rep("0", length(changeNodes))
        }
  #    }
  outdata <- list(regimeList = outlist, regimeMatrix = outmatrix)
  return(outdata) }
  
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
  
  # fill outdata
  for(i in seq(numNodes)) nodeMatrix[, i] <- unlist(lapply(ouchTrees, isMonophyletic, taxa = nodeMembers[[i]]))
  for(i in seq(numTrees)) {
    treeRegMatrix <- regMatrix * matrix(nodeMatrix[i, ], dim(regMatrix)[1], dim(regMatrix[2], byrow = T) # multiplies regMatrix by nodes present
    treeRegMatrix <- treeRegMatrix[which(apply(treeRegMatrix, 1, sum) > 0), ] # subset for regimes that still have nodes
    numTreeRegs <- dim(treeRegMatrix)[1]
    treeRegs <- list(numTreeRegs)
    for(j in seq(numTreeRegs)) {
      changeNodes[[i]] <- c("1", unlist(lapply(nodeMembers[as.logical(nodeMatrix[i, ])], mrcaOUCH, tree = ouchTrees[[i]]))) # adds the root as a change so that paintBranches will work correctly
      FINISH! 
      }
    numNodesTemp <- sum(nodeMatrix[i, ])
    regList[[i]] <- lapply(changeNodes, paintBranches, tree = ouchTrees[[i]])
  }
  outdata <- list(regList = regList, nodeMatrix = nodeMatrix)
  return(outdata)
}

regimeMatrix <- function(n = NULL, nodeNames = NULL, regimeNames = NULL, maxNodes = NULL) {
  if(identical(n, NULL) && identical(nodeNames, NULL)) stop("You have to give regimeMatrix the number of nodes, a vector of node names, or both")
  if(identical(nodeNames, NULL)) nodeNames <- as.character(seq(n))
  else n <- length(nodeNames)
  numberOfRegimes <- ifelse(n == 1, 2, 2^n)
  outmatrix <- matrix(NA, nrow = numberOfRegimes, ncol = n, dimnames = list(regimeNames, nodeNames))
  for(i in 1:(numberOfRegimes - 1)) outmatrix[i, ] <- as.binary(i, digits = n)
  outmatrix[numberOfRegimes, ] <- as.binary(0, digits = n)
  if(!identical(maxNodes, NULL)) {
    outmatrix <- outmatrix[apply(outmatrix,1,sum) <= maxNodes, ]
    dimnames(outmatrix)[[1]] = as.character(seq(dim(outmatrix)[1]))
  }
  return(outmatrix)
}

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