ouSim <- function(tree, ...) {
# right now this is just a switcher, but eventually these should be turned into proper methods of a generic ouSim
  switch(class(tree), 
         phylo = ouSim.phylo(tree, ...), # original function
         ouchtree = ouSim.ouchtree(tree, ...), # completed
         brown = ouSim.brownHansen(tree, ...),
         hansen = ouSim.brownHansen(tree, ...),
         hansenBatch = ouSim.hansenBatch(tree, ...),
         hansenSummary = ouSim.hansenBatch(tree, ...),
         stop("Unrecognized tree class")
         )
}