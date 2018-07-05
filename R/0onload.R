### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 14 2018 (17:36) 
## Version: 
## Last-Updated: jul  5 2018 (15:57) 
##           By: Brice Ozenne
##     Update #: 10
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

## copy from survival package
## get("model.matrix.coxph", envir = asNamespace("survival"), inherits = FALSE) ## (old version)
survival_model.matrix <- function (object, data = NULL, contrast.arg = object$contrasts, ...){
    if (is.null(data) && !is.null(object[["x"]])) 
        return(object[["x"]])
    Terms <- delete.response(object$terms)
    if (is.null(data)) 
        mf <- stats::model.frame(object)
    else {
        if (is.null(attr(data, "terms"))) 
            mf <- stats::model.frame(Terms, data, xlev = object$xlevels)
        else mf <- data
    }
    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        temp <- survival::untangle.specials(Terms, "cluster")
        dropterms <- temp$terms
    }
    else dropterms <- NULL
    attr(Terms, "intercept") <- 1
    adrop <- 0
    stemp <- survival::untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) > 0) {
        hasinteractions <- FALSE
        for (i in stemp$vars) {
            if (any(attr(Terms, "order")[attr(Terms, "factors")[i, 
                ] > 0] > 1)) 
                hasinteractions <- TRUE
        }
        if (!hasinteractions) 
            dropterms <- c(dropterms, stemp$terms)
        else adrop <- c(0, match(stemp$var, colnames(attr(Terms, 
            "factors"))))
    }
    if (length(dropterms)) {
        temppred <- attr(terms, "predvars")
        Terms2 <- Terms[-dropterms]
        if (!is.null(temppred)) {
            attr(Terms2, "predvars") <- temppred[-(1 + dropterms)]
        }
        X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
        renumber <- match(colnames(attr(Terms2, "factors")), 
            colnames(attr(Terms, "factors")))
        attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
    }
    else X <- model.matrix(Terms, mf, contrasts = contrast.arg)
    Xatt <- attributes(X)
    xdrop <- Xatt$assign %in% adrop
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    X
}


######################################################################
### 0onload.R ends here
