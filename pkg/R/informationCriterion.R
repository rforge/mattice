informationCriterion <- function(u = NULL, lnL = NULL, K, n = 1, names = NULL) {
  if(identical(n, 1)) warning("Information criterion values calculated assuming sample size = 1; if this is accurate, consider additional sampling for future projects.")
  if(identical(u,NULL)) u <- -2 * lnL # deviance (u) needed; take from lnL if not provided, ignore lnL if provided
  if(identical(names, NULL)) names <- seq(length(u))
  AIC <- vector("numeric", length(u))
  BIC <- vector("numeric", length(u))
  AICc <- vector("numeric", length(u))
  for(i in 1:length(u)) {
    AICc[i] <- u[i] + (2 * K[i] * (n / (n - K[i] - 1)))
    AIC[i] <- u[i] + (2 * K[i])
    BIC[i] <- u[i] + (log(n) * K[i]) }
  deltaAIC <- as.vector(lapply(AIC, function(x, allX) {x - min(allX, na.rm = T)}, allX = AIC), mode = "numeric")
  deltaAICc <- as.vector(lapply(AICc, function(x, allX) {x - min(allX, na.rm = T)}, allX = AICc), mode = "numeric")
  deltaBIC <- as.vector(lapply(BIC, function(x, allX) {x - min(allX, na.rm = T)}, allX = BIC), mode = "numeric")
  AICwi <- as.vector(lapply(deltaAIC, function(x, allDelta) {exp(-0.5 * x) / sum(exp(-0.5 * allDelta), na.rm = T)}, allDelta = deltaAIC), mode = "numeric")
  AICcwi <- as.vector(lapply(deltaAICc, function(x, allDelta) {exp(-0.5 * x) / sum(exp(-0.5 * allDelta), na.rm = T)}, allDelta = deltaAICc), mode = "numeric")
  BICwi <- as.vector(lapply(deltaBIC, function(x, allDelta) {exp(-0.5 * x) / sum(exp(-0.5 * allDelta), na.rm = T)}, allDelta = deltaBIC), mode = "numeric")
  outdata <- list(names = names, u = u, K = K, AIC = AIC, AICc = AICc, BIC = BIC, AICwi = AICwi, AICcwi = AICcwi, BICwi = BICwi)
  class(outdata) <- 'informationCriterion'
  return(outdata)
}

informationCriterion.hansenBatch <- function(hansenBatch) {
## call informationCriterion for a 'hansen.batch' object
## Just returns AIC, AICc, and BIC weights for each of the trees analyzed in a hansenBatch object
  outdata <- vector("list", length(hansenBatch$hansens))
  N = hansenBatch$N
  for(i in seq(outdata)) {
    temp <- hansenBatch$hansens[[i]]
    outdata[[i]] <- informationCriterion(lnL = temp[, 'loglik'], K = temp[, 'dof'], n = N, names = row.names(temp))
    }
  return(outdata)
}

print.informationCriterion <- function(ic, ...) {
  items <- c('u', 'K', 'AIC', 'AICc', 'BIC', 'AICwi', 'AICcwi', 'BICwi')
  out <- matrix(NA, nrow = length(ic$names), ncol = length(items), dimnames = list(ic$names, items))
  for(i in items) out[, i] <- ic[[i]]
  print(out) 
  return(NULL)
}