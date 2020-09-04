### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 14 2018 (17:36) 
## Version: 
## Last-Updated: sep  4 2020 (09:37) 
##           By: Brice Ozenne
##     Update #: 37
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * .onAttach (i.e. message when loading the package)
.onAttach <- function(lib, pkg="riskRegression") {
    desc <- utils::packageDescription(pkg)
    packageStartupMessage(desc$Package, " version ",desc$Version)
}

## * riskRegression.options
## adapted from the lava package https://github.com/kkholst/lava
riskRegression.env <- new.env()
assign("options",
       list(method.predictRisk = paste0("predictRisk.",c("ARR", "BinaryTree", "CauseSpecificCox", "Cforest", "coxph", "coxph.penal", "coxphTD", "cph", "CSCTD", "Ctree", "default", "double", "factor", "FGR", "flexsurvreg", "formula", "gbm", "glm", "integer", "lrm", "matrix", "numeric", "penfitS3", "prodlim", "psm", "randomForest", "ranger", "rfsrc", "riskRegression", "rpart", "selectCox", "SmcFcs", "SuperPredictor", "survfit","wglm")),
            method.predictRiskIID = paste0("predictRiskIID.",c("CauseSpecificCox", "coxph", "cph", "default", "glm", "phreg", "wglm"))),
       envir = riskRegression.env)

## cat(paste("c(\"",paste(gsub("predictRiskIID.","",as.character(utils::methods("predictRiskIID")), fixed=TRUE),collapse = "\", \""),"\")\n",sep=""))
## cat(paste("c(\"",paste(gsub("predictRisk.","",as.character(utils::methods("predictRisk")), fixed=TRUE),collapse = "\", \""),"\")\n",sep=""))

##' @title Global options for \code{riskRegression}
##'
##' @description Output and set global options for the \code{riskRegression} package.
##'
##' @param ... for now limited to \code{method.predictRisk} and \code{mehtod.predictRiskIID}.
##'
##' @details only used by the \code{ate} function.
##'
##' @examples
##' options <- riskRegression.options()
##'
##' ## add new method.predictRiskIID
##' riskRegression.options(method.predictRiskIID = c(options$method.predictRiskIID,"xx"))
##'
##' riskRegression.options()
##' @export
riskRegression.options <- function(...) {
    dots <- list(...)
    out <- get("options", envir = riskRegression.env)
    if (length(dots)==0){
        return(out)
    }else{
        if(any(names(dots) == "")){
            stop("All element must be named \n")
        }
        dots.name <- names(dots)
        dots.existing <- intersect(dots.name,names(out))
        dots.new <- setdiff(dots.name,names(out))
        if(length(dots.existing)>0){
            out[dots.existing] <- dots[dots.existing] ## assign
        }
        if(length(dots.new)){
            out <- c(out,dots[dots.new])
        }
        out <- out[order(names(out))] ## reorder names
        assign("options", out, envir = riskRegression.env)
        invisible(out)
    }
}

## * survival_model.matrix
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
