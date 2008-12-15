plot.ouSimPhylo <- function(ouSim, nodeColor = "blue", nodeDotSize = 1.4, colors = NULL, ...) {
## To plot different clades, set the colors vector according to the branches in the original 
## only passes the ... along to lines
  if(identical(colors, NULL)) {
    if("colors" %in% names(ouSim)) colors <- ouSim$colors
    else colors <- rep("black", length(ouSim$branchList))
    }
  branches = length(ouSim$branchList)
  plot(1:ouSim$steps, ylim = range(unlist(ouSim$branchList)), type = "n", ylab = "Trait value", xlab = "Time")
  for(i in 1:branches) lines(ouSim$timesList[[i]], ouSim$branchList[[i]], col = colors[i], ...)
  for(i in 1:branches) points(ouSim$timesList[[i]][1], ouSim$branchList[[i]][1], pch = 19, col = nodeColor, cex = nodeDotSize) 
}