ouSim <- function(tree, ...) {
  switch(class(tree), 
         phylo = ouSim.phylo(tree, ...),
         ouchtree = ouSim.ouchtree(tree, ...),
         brown = ouSim.brown(tree, ...),
         hansen = ouSim.hansen(tree, ...),
         batchHansen = ouSim.batchHansen(tree, ...),
         "Unrecognized tree class"
         )
}