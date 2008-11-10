apeToOuch <-
## Functions for converting from APE tree format to OUCH tree format
## Andrew Hipp, February 2005
## Generalized to non-ultrametric trees Sept 2006
# Converts a "phylo" object from APE to a 3-column matrix ordered for OUCH format
function(tree, scalingFactor = 1) {
  numberOfNodes <- length(tree$edge.length)
  treeMatrix <- c(0*1:(5*(numberOfNodes + 1)))
  dim(treeMatrix) <- c(numberOfNodes + 1, 5)
  species = 1
  ancestor = 2
  times = 3
  apeAncestor = 4
  apeChild = 5
  # totalLengthUnscaled <- distFromRoot(tree, "1")
  # 'totalLengthUnscaled' was used in the original function to determine endtime of each terminal branch; assumes ultrametricity
  nodeCounter <- 1
  treeMatrix[1, 1:3] <- c(NA, 0, 0)
  # This version sets the ancestor of the root as 0, following the older version of ouch; change this by uncommenting the following line:
  # treeMatrix[1, 1:3] <- c(NA, NA, 0)
  for (lineNumber in 1:numberOfNodes) {
    parentNode = tree$edge[lineNumber, 1]
    childNode = tree$edge[lineNumber, 2]
    if (as.integer(childNode) < 0) {
      nodeCounter <- nodeCounter + 1
      treeMatrix[nodeCounter, species] <- NA
      treeMatrix[nodeCounter, times] <- distFromRoot(tree, childNode) * scalingFactor
      treeMatrix[nodeCounter, apeChild] <- childNode
      treeMatrix[nodeCounter, apeAncestor] <- parentNode }}
  numberOfInternalNodes <- nodeCounter

  for (lineNumber in 1:numberOfNodes) {
    parentNode = tree$edge[lineNumber, 1]
    childNode = tree$edge[lineNumber, 2]
    if (as.integer(childNode) > 0) {
      nodeCounter <- nodeCounter + 1
      treeMatrix[nodeCounter, species] <- tree$tip.label[as.integer(childNode)]
      treeMatrix[nodeCounter, times] <- distFromRoot(tree, childNode) * scalingFactor
         # This is revised from original function, which just dropped in the maximum time
      treeMatrix[nodeCounter, apeChild] <- childNode
      treeMatrix[nodeCounter, apeAncestor] <- parentNode }}

  for (lineNumber in (1:numberOfNodes + 1)) {
    if (treeMatrix[lineNumber, apeAncestor] == "-1") treeMatrix[lineNumber, ancestor] <- 1
      else {treeMatrix[lineNumber, ancestor] <- match(treeMatrix[lineNumber, apeAncestor], treeMatrix[1:numberOfInternalNodes, apeChild])}
    }

  matrixLength <- length(treeMatrix) / 5
  treeDataFrame = data.frame(node = I(as.integer(1:matrixLength)),
                             species = I(treeMatrix[1:matrixLength,1]),
                             ancestor = I(as.integer(treeMatrix[1:matrixLength,2])),
                             times = I(as.double(treeMatrix[1:matrixLength,3])))
  return(treeDataFrame) }

distFromRoot <-
# Finds the distance of a particular node from the root of the tree.
# Takes an ape tree (class "phylo") and child node label (class "character") as its arguments.
# GENERALLY ONLY CALLED BY apeToOuch
function(tree, childNode) {
  if (childNode == "-1") return(0)
  numberOfNodes <- length(tree$edge.length)
  childLine <- match(childNode, tree$edge[1:numberOfNodes, 2])
  ancestor <- tree$edge[childLine, 1]
  distance <- tree$edge.length[childLine] + distFromRoot(tree, ancestor)
  return(distance) }

mergeData <-
# Merges an OUCH-style tree as output from apeToOuch
function(ouchTree, dataMatrix) {
  newTree = merge(ouchTree, data.frame(dataMatrix), all.x = T)
  newTree = newTree[order(newTree$node), ]
  row.names(newTree) = 1:length(newTree$times)
  return(newTree) }