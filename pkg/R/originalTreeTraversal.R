# ---------------------------------------------------------------------
# FUNCTIONS FOR TRAVERSING AN OUCH TREE AND PAINTING REGIMES
# ---------------------------------------------------------------------
# Copied from functions used in Hipp 2007 Evolution paper
# Last checked on ouch v 1.2-4
# functions included in this file:
# 1. paintBranches
# 2. mrcaOUCH
# 3. ancestorLine
# 4. allPossibleRegimes
# 5. regimeVectors
# 6. regimeNodes

paintBranches <-
# Reads down an ouch data.frame tree row by row and paints regimes according, with regimes changing at nodes specified
# arguments
#  "node" "ancestor" "times" = the standard tree specification vectors of the OUCH-style tree
#  "regimeShiftNodes" = a vector of nodes at which selective regimes shift: root must be included, but tips are meaningless in this context
#  "regimeTitles" = a vector of titles for the regimes that begin at the root and at the nodes indicated in "regimeShiftNodes",
#                   in order of description in "regimeShiftNodes", except that the root is listed first in "regimeTitles"
#                   but not at all in "regimeShiftNodes"... defaults to "node[x]regime
# Value: a vector of regimes that can be plopped right into an OUCH-style tree data frame
function(node, ancestor, times, regimeShiftNodes, regimeTitles) {
  names(regimeTitles) = regimeShiftNodes
  colorsVector = vector("character", length(node))
  for (i in 1:length(ancestor)) {
    # First three lines fill up the vector for nodes that are hit in order
    # uncomment the following line if you are using a newer version of OUCH:
    if (is.na(ancestor[as.integer(i)])) {
    # and comment the following line:
    # if (ancestor[as.integer(i)] == 0) {
      colorsVector[i] = regimeTitles["1"]
      next }
    if (any(ancestor[i] == regimeShiftNodes)) {
      colorsVector[i] = regimeTitles[as.character(ancestor[as.integer(i)])]
      next }
    if (colorsVector[as.integer(ancestor[i])] != "") {
      colorsVector[i] = colorsVector[as.integer(ancestor[i])]
      next }

    # These lines fill up the vector for nodes run reached before their immediate ancestor
    nodeQ = vector("integer", length(node))
    ii = i
    repeat {
      nodeQ = c(ii, nodeQ)
      ii = ancestor[ii]
      if (any(ancestor[ii] == regimeShiftNodes)) {
        colorsVector[ii] = colorsVector[as.integer(ancestor[ii])]
        break}
      if (colorsVector[ancestor[ii]] != "") {
        colorsVector[ii] = colorsVector[as.integer(ancestor[ii])]
        break} }

    for(j in nodeQ) {
      colorsVector[j] = colorsVector[as.integer(ancestor[j])] } }

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
function(node, ancestor, times, species, cladeVector) {
  tips = match(cladeVector, species)
  listOfAncestorLines = lapply(tips, ancestorLine, node = node, ancestor = ancestor)
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
#  "node" and "ancestor" = the standard tree specification vectors of the OUCH-style tree
#  "tip" = the tip node to trace back
# Value: a vector of nodes leading from a tip to the root
function(node, ancestor, tip) {
  nodesVector = vector("character")
  counter = 0
  repeat {
    if (is.na(tip)) break
    counter = counter + 1
    nodesVector[counter] = ancestor[as.integer(tip)]
    tip = ancestor[as.integer(tip)] }
  return(nodesVector) }

allPossibleRegimes <-
# Generates a list of vectors of all possible 2^n regimes for a given list of n ancestor nodes, assuming that each node
#  can either be a change node or not, and that all combinations are meaningful.
# NOTE: This routine returns all possible permutations, ignoring the fact that the root must be included as a changeNode
#  for paintBranches to work properly. Exclude the root when calling this routine.
# Arguments:
#  "changeNodes" = vector of all change nodes in all possible scenarios
# Value: a list of changeNode vectors (assumes type = "character"), one for each possible scenario
function(changeNodes) {
  numberOfRegimes = ifelse(length(changeNodes) == 1, 2, 2^length(changeNodes))
  regime = vector("list", numberOfRegimes)
  for (i in 1:(numberOfRegimes - 1)) {
    remainder = i
    n = NULL
    for (j in as.integer(log2(i)):0) {
      if (2^j > remainder) n[j+1] = NA
      else {n[j+1] = changeNodes[j+1]
            remainder = remainder %% 2^j }}
    regime[[i]] = sort(n[!is.na(n)]) }
  regime[[numberOfRegimes]] = rep("0", times = as.integer(log2(i)) + 1)
  return(regime) }

regimeVectors <-
# Generates the list of painted branches representing all possible selective regimes for OU analyses, taking as argument
# species vectors that describe the clades at the bases of which regimes are specified to change.
# Arguments:
#  "node" "ancestor" "times" "species" = the standard tree specification vectors of the OUCH-style tree
#  "cladeMembersList" = list of vectors containing names of the members of each clade (except for the root of the tree)
# Value: list of vectors that can each be plugged directly into OU analysis as the "regimes" argument
function(node, ancestor, times, species, cladeMembersList) {
  changeNodesList = lapply(cladeMembersList, mrcaOUCH, node = node, ancestor = ancestor, times = times, species = species) #Returns a list of length-1 character vectors, each containing a single changeNode -- the fact that this is a list causes problems in paintBranches if not changed to a 1-d vector
  changeNodesVector = vector("character", length(changeNodesList))
  for (i in 1:length(changeNodesList)) # Changing cladeMemberList to a 1-d vector
    {changeNodesVector[i] = changeNodesList[[i]]}
  allRegimes = allPossibleRegimes(changeNodesVector)
  regimePaintings = vector("list", length(allRegimes))
  for (i in 1:length(allRegimes)) {
    allRegimes[[i]] = c("1", allRegimes[[i]])
    regimePaintings[[i]] = paintBranches(node, ancestor, times, allRegimes[[i]], as.character(allRegimes[[i]])) }
  return(regimePaintings) }

regimeNodes <-
# Returns a binary table indicating whether a regime node is present or absent based on allPossibleRegimes output;
#  presumes nodes are labelled 1:n)
# ONLY NEEDED FOR LABELLING AND INTERPRETATION; not actually used in the analysis
function(n) {
  regimesList = allPossibleRegimes(1:n)
  regimesNameMatrix = matrix(data = NA, ncol = n, nrow = length(regimesList), dimnames = list(as.character(1:length(regimesList)), as.character(1:n)))
  for (i in 1:length(regimesList)) {
    for (j in 1:n) {
      if (is.na(match(j,regimesList[[i]]))) regimesNameMatrix[i,j] = 0
       else regimesNameMatrix[i,j] = 1 }}
  return(regimesNameMatrix) }