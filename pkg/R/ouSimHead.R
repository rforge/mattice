ouSim <- function(tree, ...) {
# right now this is just a switcher, but eventually these should be turned into proper methods of a generic ouSim
  switch(class(tree), 
         phylo = ouSim.phylo(tree, ...),
         ouchtree = ouSim.ouchtree(tree, ...),
         brown = ouSim.brownHansen(tree, ...),
         hansen = ouSim.brownHansen(tree, ...),
         batchHansen = ouSim.batchHansen(tree, ...),
         "Unrecognized tree class"
         )
}