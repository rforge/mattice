plot.ouSim <- function(x, nodeColor = "blue", nodeDotSize = 1.4, colors = NULL, ...) {
## To plot different clades, set the colors vector according to the branches in the original 
## only passes the ... along to lines
  branches = length(x$branchList)
  if(length(nodeColor) == 1) nodeColor <- rep(nodeColor, branches)
  if(identical(colors, NULL)) {
    if("colors" %in% names(x)) colors <- x$colors
    else colors <- rep("black", branches)
    }
  plot(1:x$steps, ylim = range(unlist(x$branchList)), type = "n", ylab = "Trait value", xlab = "Time")
  for(i in 1:branches) lines(x$timesList[[i]], x$branchList[[i]], col = colors[i], ...)
  for(i in 1:branches) points(x$timesList[[i]][1], x$branchList[[i]][1], pch = 19, col = nodeColor[i], cex = nodeDotSize) 
}