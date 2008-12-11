ouSim <- function(tree, ...) {
# right now this is just a switcher, but eventually these should be turned into proper methods of a generic ouSim
  switch(class(tree), 
         phylo = ouSim.phylo(tree, ...), # original function
         ouchtree = ouSim.ouchtree(tree, ...), # completed
         browntree = ouSim.brownHansen(tree, ...), # completed
         hansentree = ouSim.brownHansen(tree, ...),
         hansenBatch = ouSim.hansenBatch(tree, ...), # completed
         hansenSummary = ouSim.hansenBatch(tree, ...), # completed
         stop("Unrecognized tree class")
         )
}