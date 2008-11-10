# ---------------------------------------------------------------------
# FUNCTIONS FOR PERFORMING SUMMARIZING ANALYSES
# ---------------------------------------------------------------------
# Copied from functions used in Hipp 2007 Evolution paper
# Last checked in ouch v 1.2-4
# functions included in this file:
# 1. hansenStats

hansenStats <-
# Returns a list of statistics for a runBatchHansenFit list
function(hansenRunList) {
  columnNames = c("maxAlpha","bestFitModel", "maxAICweight")
  rowNames = as.character(1:length(hansenRunList))
  hansenStatsTemp = matrix(data = NA, ncol = length(columnNames), nrow = length(hansenRunList), dimnames = list(rowNames, columnNames))
  counter = 0
  for (i in hansenRunList) {
    counter = counter + 1
    hansenStatsTemp[counter, "maxAlpha"] = max(i[0:16, "alpha"])
    hansenStatsTemp[counter, "bestFitModel"] = match(max(i[0:16, "AICweight"]), i[0:16, "AICweight"])
    hansenStatsTemp[counter, "maxAICweight"] = max(i[0:16, "AICweight"])}
  return(hansenStatsTemp)}

