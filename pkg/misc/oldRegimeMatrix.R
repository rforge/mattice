oldRegimeMatrix <- function(n = NULL, nodeNames = NULL, regimeNames = NULL, maxNodes = NULL) {
## a brute-force approach, very inefficient as n and maxNodes diverge
## I'm leaving this in the project only b/c i don't know whether the recursive function will fail at large maxNodes
  if(identical(n, NULL) && identical(nodeNames, NULL)) stop("You have to give regimeMatrix the number of nodes, a vector of node names, or both")
  if(identical(nodeNames, NULL)) nodeNames <- as.character(seq(n))
  else n <- length(nodeNames)
  if(identical(maxNodes, NULL)) maxNodes <- n
  outmatrix <- matrix(NA, nrow = 0, ncol = n, dimnames = list(NULL, nodeNames))
  maxNumberOfRegimes <- ifelse(n == 1, 2, 2^n)
  counter <- 1
  repeat {
    temp <- as.binary(counter, digits = n)
    if(sum(temp) <= maxNodes) outmatrix <- rbind(outmatrix, temp) 
    if(counter == maxNumberOfRegimes - 1) break
    counter <- counter + 1
    }
  outmatrix <- rbind(outmatrix, as.binary(0, digits = n))
  dimnames(outmatrix)[[1]] <- seq(dim(outmatrix)[1])
  return(outmatrix)
}


as.binary <- function(n, base = 2, r = FALSE, digits = NULL)
# Robin Hankin <initialDOTsurname at soc.soton.ac.uk (edit in obvious way; spam precaution)>
# submitted to R listserv Thu Apr 15 12:27:39 CEST 2004
# AH added 'digits' to make it work with regimeMatrix
# https://stat.ethz.ch/pipermail/r-help/2004-April/049419.html
{
   out <- NULL
   while(n > 0) {
     if(r) {
       out <- c(out , n%%base)
     } else {
       out <- c(n%%base , out)
     }   
     n <- n %/% base
   }
   if(!identical(digits, NULL) && !r) out <- c(rep(0, digits-length(out)), out)
   if(!identical(digits, NULL) && r) out <- c(out, rep(0, digits-length(out)))
   return(out)
}

