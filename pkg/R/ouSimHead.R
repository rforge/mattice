ouSim <- function(object, ...) {
# right now this is just a switcher, but eventually these should be turned into proper methods of a generic ouSim
  switch(class(object), 
         phylo = ouSim.phylo(object, ...), 
         ouchtree = ouSim.ouchtree(object, ...), 
         browntree = ouSim.brownHansen(object, ...), 
         hansentree = ouSim.brownHansen(object, ...),
         hansenBatch = ouSim.hansenBatch(object, ...), 
         hansenSummary = ouSim.hansenBatch(object, ...), 
         stop("Unrecognized tree class")
         )
}