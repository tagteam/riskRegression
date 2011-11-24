.packageName <- "riskRegression"

.First.lib <- function(lib, pkg) {
  library.dynam("riskRegression", pkg, lib)
}

.Last.lib <- function(lib){
  library.dynam.unload("riskRegression",lib)
}

