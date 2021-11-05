### 0onload.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 14 2018 (17:36) 
## Version: 
## Last-Updated: nov  5 2021 (16:11) 
##           By: Brice Ozenne
##     Update #: 44
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
       list(method.predictRisk = paste0("predictRisk.",c("ARR", "BinaryTree", "CauseSpecificCox", "Cforest", "cox.aalen", "coxph", "coxph.penal", "coxphTD", "cph", "CSCTD", "Ctree", "default", "double", "factor", "FGR", "flexsurvreg", "formula", "gbm", "glm", "hal9001", "integer", "lrm", "matrix", "multinom", "numeric", "penfitS3", "prodlim", "psm", "randomForest", "ranger", "rfsrc", "riskRegression", "rpart", "selectCox", "singleEventCB", "SmcFcs", "SuperPredictor", "survfit", "wglm", "aalen")),
            method.predictRiskIID = paste0("predictRiskIID.",c("CauseSpecificCox", "coxph", "cph", "default", "glm", "phreg", "wglm","multinom"))),
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


######################################################################
### 0onload.R ends here
