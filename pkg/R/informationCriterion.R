informationCriterion <- function(u = NULL, lnL = NULL, K, n = 0, names = NULL) {
## Returns information criterion values + weights for a vector of u or lnL, a vector of K (= df), and a single n (sample size); names for analyses are optional
  if(identical(u,NULL)) u <- -2 * lnL # deviance (u) needed; take from lnL if not provided, ignore lnL if provided
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
  return(list(names = names, u = u, K = K, AIC = AIC, AICc = AICc, BIC = BIC, AICwi = AICwi, AICcwi = AICcwi, BICwi = BICwi)) }

informationCriterion.hansenBatch <- function(hansenBatch) {
## call informationCriterion for a 'hansen.batch' object
## Just returns AIC, AICc, and BIC weights for each of the trees analyzed in a hansenBatch object
  outdata <- vector("list", length(hansenBatch$hansens))
  N = hansenBatch$N
  for(i in 1:length(outdata)) {
    temp <- hansenBatch$hansens[[i]]
    outdata[[i]] <- informationCriterion(lnL = temp[, 'loglik'], K = temp[, 'dof'], n = N, names = row.names(temp))
    }
  outdata
}