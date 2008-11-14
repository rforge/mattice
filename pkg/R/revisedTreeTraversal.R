# ---------------------------------------------------------------
# FUNCTIONS FOR TRAVERSING AN S4 OUCH TREE AND PAINTING REGIMES #
# ---------------------------------------------------------------

# Modified from functions used in Hipp 2007 Evolution paper
# Initially written for ouch v 1.2-4
# 10 November 2008: changed everything to operate on an ouchtree object (ouch >= v2), otherwise functions the same;
#  checked on ouch v2.4-2
# functions included in this file:
# 1. paintBranches
# 2. mrcaOUCH
# 3. ancestorLine
# 4. allPossibleRegimes
# 5. regimeVectors


paintBranches <-
# Paints branches with regimes changing at nodes specified
# arguments
#  "node" "ancestor" "times" = the standard tree specification vectors of the OUCH-style tree
#  "regimeShiftNodes" = a vector of nodes at which selective regimes shift: root must be included, but tips are meaningless in this context
#  "regimeTitles" = a vector of titles for the regimes that begin at the root and at the nodes indicated in "regimeShiftNodes",
#                   in order of description in "regimeShiftNodes", except that the root is listed first in "regimeTitles"
#                   but not at all in "regimeShiftNodes"... defaults to "node[x]regime
# Value: a vector of regimes that can be plopped right into an OUCH-style tree data frame
function(tree, regimeShiftNodes, regimeTitles) {
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

mrcaOUCH <-
# Finds most recent common ancestor for a vector of tips by:
#  1. Creating a vector of ancestral nodes leading to each tip
#  2. Creating an intersection set of ancestral nodes for all taxa by intersecting taxa successively with the last intersection set
#  3. Returning the node of the final intersection set that has the highest time
# Arguments:
#  "node" "ancestor" "times" "species" = the standard tree specification vectors of the OUCH-style tree
#  "cladeVector" = vector of species for which you want to find the most recent common ancestor
# Value: the node number (as an integer) of the most recent common ancestor
# Works! 3-31-06
function(cladeVector, tree) {
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
  
  tips = match(cladeVector, species) 
  listOfAncestorLines = lapply(tips, ancestorLine, tree = tree) # 10 nov 08: this is identical to the appropriate subset of tree@lineages
  latestMatch = listOfAncestorLines[[1]]
  for (i in listOfAncestorLines) {
    latestMatch = i[match(latestMatch, i, nomatch = 0)] }
  timesVector = times[as.integer(latestMatch)]
  if(length(timesVector) == 1) {
    if (is.na(timesVector)) mrca = "1"
      else mrca = timesVector}
    else mrca = latestMatch[match(max(as.double(timesVector), na.rm = TRUE), timesVector)]
  return(mrca) }

ancestorLine <-
# Creates a vector of ancestral nodes for a tip
# Arguments:
#  CHANGED to "tree" 10 nov 08: "node" and "ancestor" = the standard tree specification vectors of the OUCH-style tree
#  "tip" = the tip node to trace back
# Value: a vector of nodes leading from a tip to the root
# 10 nov 08: changed to just grab the appropriate element from tree@lineages
function(tip, tree) {
  ## ------------------ begin ouchtree block -----------------
  ## check to see if tree inherits 'ouchtree'
  if (!is(tree,'ouchtree')) 
	stop(paste('This function has been rewritten to use the new S4 ', sQuote('ouchtree'), ' class.',
	'\nYou can generate a tree of this class by calling ', sQuote('ouchtree()'), '.', sep = ""))
  ## get the vectors we need:
  #ancestor <- tree@ancestors # class = "character"
  #node <- tree@nodes # class = "character"
  #species <- tree@nodelabels # class = "character" -- note that nodelabels is more general than this indicates and the name should be changed throughout at some point
  #times <- tree@times # class = "numeric"
  ## ------------------ end ouchtree block -------------------
  tip <- as.numeric(tip)
  #nodesVector = vector("character")
  #counter = 0
  #repeat {
  #  if (is.na(tip)) break
  #  counter = counter + 1
  #  nodesVector[counter] = ancestor[tip]
  #  tip = ancestor[tip] }
  nodesVector <- c(as.character(tree@lineages[[tip]][2:length(tree@lineages[[tip]])]), NA)
  return(nodesVector) 
  }

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
#    if nodeMatrix = F: a list of changeNode vectors (assumes type = "character"), one for each possible scenario
#    if nodeMatrix = T: a binary table indicating whether a regime node is present or absent based on allPossibleRegimes output; 
#                       presumes nodes are labelled 1:n
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
    if(nodeMatrix == T) {
      #n <- ifelse(length(changeNodes) == 1, as.numeric(changeNodes), length(changeNodes))
      regimesNameMatrix = matrix(
        data = NA, nrow = numberOfRegimes, ncol = length(changeNodes), dimnames = list(
          as.character(seq(numberOfRegimes)), as.character(seq(length(changeNodes)))
          )
      )
      for (i in seq(numberOfRegimes)) {
        for (j in seq(length(changeNodes))) {
          if (is.na(match(j,regime[[i]]))) regimesNameMatrix[i,j] = 0
          else regimesNameMatrix[i,j] = 1 
        }
      }
      outdata <- regimesNameMatrix
      if(!identical(maxNodes, NULL)) {
        outdata <- outdata[apply(outdata,1,sum) <= maxNodes, ]
        dimnames(outdata)[[1]] = as.character(seq(dim(outdata)[1]))
        }
    }
    else {
      outdata <- regime 
      if(!identical(maxNodes, NULL)) {
        outdata <- outdata[sapply(outdata, length) <= maxNodes]
        outdata[[length(outdata) + 1]] <- rep("0", length(changeNodes))
        }
      }
  return(outdata) }

regimeVectors <-
# Generates the list of painted branches representing all possible selective regimes for OU analyses, taking as argument
# species vectors that describe the clades at the bases of which regimes are specified to change.
# Arguments:
#  "node" "ancestor" "times" "species" = the standard tree specification vectors of the OUCH-style tree
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
# Value: list of vectors that can each be plugged directly into OU analysis as the "regimes" argument
function(tree, cladeMembersList, maxNodes = NULL) {
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
      
  changeNodesList <- lapply(cladeMembersList, mrcaOUCH, tree = tree) #Returns a list of length-1 character vectors, each containing a single changeNode -- the fact that this is a list causes problems in paintBranches if not changed to a 1-d vector
  changeNodesVector <- unlist(changeNodesList)
  #changeNodesVector = vector("character", length(changeNodesList))
  #for (i in 1:length(changeNodesList)) # Changing cladeMemberList to a 1-d vector
  #  {changeNodesVector[i] = changeNodesList[[i]]}
  allRegimes = allPossibleRegimes(changeNodesVector, maxNodes)
  regimePaintings = vector("list", length(allRegimes))
  for (i in 1:length(allRegimes)) {
    allRegimes[[i]] = c("1", allRegimes[[i]])
    regimePaintings[[i]] = as.factor(paintBranches(tree, allRegimes[[i]], as.character(allRegimes[[i]])))
    names(regimePaintings[[i]]) <- tree@nodes
    message(paste('Created regime',i))}
  return(regimePaintings) }
