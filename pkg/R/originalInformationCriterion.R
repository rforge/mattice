informationCriterion <- function(u = NULL, lnL = NULL, K, n = 0, names = NULL) {
## Returns information criterion values + weights for a vector of u or lnL, K (= df), and n (sample size); names for analyses are optional
  if(identical(u,NULL)) u <- -2 * lnL
  AIC = vector("numeric", length(u))
  BIC = vector("numeric", length(u))
  AICc = vector("numeric", length(u))
  for(i in 1:length(u)) {
    AICc[i] = u[i] + (2 * K[i] * (n / (n - K[i] - 1)))
    AIC[i] = u[i] + (2 * K[i])
    BIC[i] = u[i] + (log(n) * K[i]) }
  deltaAIC = as.vector(lapply(AIC, function(x, allX) {x - min(allX)}, allX = AIC), mode = "numeric")
  deltaAICc = as.vector(lapply(AICc, function(x, allX) {x - min(allX)}, allX = AICc), mode = "numeric")
  deltaBIC = as.vector(lapply(BIC, function(x, allX) {x - min(allX)}, allX = BIC), mode = "numeric")
  AICwi = as.vector(lapply(deltaAIC, function(x, allDelta) {exp(-0.5 * x) / sum(exp(-0.5 * allDelta))}, allDelta = deltaAIC), mode = "numeric")
  AICcwi = as.vector(lapply(deltaAICc, function(x, allDelta) {exp(-0.5 * x) / sum(exp(-0.5 * allDelta))}, allDelta = deltaAICc), mode = "numeric")
  BICwi = as.vector(lapply(deltaBIC, function(x, allDelta) {exp(-0.5 * x) / sum(exp(-0.5 * allDelta))}, allDelta = deltaBIC), mode = "numeric")
  return(list(names = names, u = u, K = K, AIC = AIC, AICc = AICc, BIC = BIC, AICwi = AICwi, AICcwi = AICcwi, BICwi = BICwi)) }

