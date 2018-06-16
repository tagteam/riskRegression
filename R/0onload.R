### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 14 2018 (17:36) 
## Version: 
## Last-Updated: jun 14 2018 (17:37) 
##           By: Brice Ozenne
##     Update #: 2
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

.onAttach <- function(lib, pkg="riskRegression") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

survival_model.matrix <- get("model.matrix.coxph", envir = asNamespace("survival"), inherits = FALSE)


######################################################################
### 0onload.R ends here
