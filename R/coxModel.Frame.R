### coxModel.Frame.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:33) 
## Version: 
## Last-Updated: feb 16 2026 (09:48) 
##           By: Thomas Alexander Gerds
##     Update #: 8
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * coxModelFrame
#' @title Extract the design matrix used to train a Cox model
#' @description Extract the design matrix used to train a Cox model. Should contain the time of event, the type of event, 
#' the variable for the linear predictor, the strata variables and the date of entry (in case of delayed entry).
#' @name coxModelFrame 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package), \code{cph}
#'     (rms package), or \code{phreg} (mets package).
#' @param center [logical] Should the variables of the linear predictor be added ?
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxModelFrame
#' @export
coxModelFrame <- function(object, center){
  UseMethod("coxModelFrame") 
} 

## ** coxModelFrame.coxph
#' @rdname coxModelFrame
#' @method coxModelFrame coxph
#' @export
coxModelFrame.coxph <- function(object, center = FALSE){

    default.start <- 0
    
    if("y" %in% names(object) == FALSE){
        stop("invalid object \n",
             "set y=TRUE in the call to ",class(object)[1]," \n")
    }
    
    if("x" %in% names(object) == FALSE){
        rhs.formula <- update(formula(object),"0~.")
        if(length(all.vars(rhs.formula))==0 || (!is.null(object$nevent) && object$nevent==0)){
            ## name.var <- setdiff(all.vars(object$terms), unlist(SurvResponseVar(coxFormula(object))[-1]))
            name.var <- names(object$coef)
            object$x <- matrix(0, nrow = NROW(object$y), ncol = length(name.var),
                               dimnames = list(NULL, name.var))
        }else{
            stop("invalid object \n",
                 "set x=TRUE in the call to ",class(object)[1]," \n")
        }
    }
 
    ## ** add x
    if(NCOL(object[["x"]])!=0){
        if(center){
            dt <- as.data.table(rowCenter_cpp(object[["x"]], center = coxCenter(object)))
        }else{
            dt <- as.data.table(object[["x"]])
        }        
        if(any(names(dt) != names(coef(object))) && NCOL(dt) == length(coef(object))){
            colnames(dt) <- names(coef(object))
        }
    }else{
        dt <- NULL
    }
    if("strata" %in% names(dt)){
        stop("The variables in the linear predictor should be named \"strata\" \n")
    }

    ## ** add y
    if(is.null(dt)){
        dt <- data.table(status = object[["y"]][,"status"])
    }else{
        dt[,c("status") := object[["y"]][,"status"]]
    }
    if("start" %in% colnames(object$y) == FALSE){
        dt[,c("start","stop") := list(default.start, object[["y"]][,"time"])]
    }else{        
        dt[,c("start","stop") := list(object[["y"]][,"start"], object[["y"]][,"stop"])]
    }
  
    ## ** add strata
    if("strata" %in% names(object)){
        dt[,c("strata") := object[["strata"]]]
    }else{
        dt[,c("strata") := factor(1)]
    }
    
    ## ** export
    first.col <- c("start","stop","status")
    data.table::setcolorder(dt, c(first.col,setdiff(names(dt),first.col)))
    dt[]
}


## ** coxModelFrame.cph
#' @rdname coxModelFrame
#' @method coxModelFrame cph
#' @export
coxModelFrame.cph <- coxModelFrame.coxph

## ** coxModelFrame.phreg
#' @rdname coxModelFrame
#' @method coxModelFrame phreg
#' @export
coxModelFrame.phreg <- function(object, center = FALSE){

    default.start <- 0

    ## ** add y
    dt <- as.data.table(unclass(object$model.frame[,1]))
    if("start" %in% names(dt) == FALSE){
        dt[, c("start") := default.start]
    }
  
    ## normalize names
    name.old <- names(dt)
    name.new <- gsub("entry","start",gsub("time","stop",name.old))
    setnames(dt, old = name.old, new = name.new)

    ## ** add x
    if(center){
        M.X <- rowCenter_cpp(object$X, center = coxCenter(object))
        colnames(M.X) <- colnames(object$X)
        dt <- cbind(dt, as.data.table(M.X))
    }else{
        dt <- cbind(dt, as.data.table(object$X))
    }
    if("strata" %in% names(dt)){
        stop("The variables in the linear predictor should be named \"strata\" \n")
    }

    ## ** add strata
    if(!is.null(object$strata.name)){
        dt[,c("strata") := as.factor(object$model.frame[[object$strata.name]])]
    }else{
        dt[,c("strata") := factor(1)]
    }

    ## ** export
    first.col <- c("start","stop","status")
    data.table::setcolorder(dt, c(first.col,setdiff(names(dt),first.col)))
    return(dt)
    
}
# }}}

## ** coxModelFrame.prodlim
#' @rdname coxModelFrame
#' @method coxModelFrame prodlim
#' @export
coxModelFrame.prodlim <- function(object, center = FALSE){
    ## argument center not used
    default.start <- 0
    originalDataOrder <- object$originalDataOrder

    ## ** add y
    dt <- object$model.response[originalDataOrder,,drop=FALSE]
    if(!is.data.table(dt)){
        if(inherits(dt,"Surv")){
            dt <- as.matrix(dt)
        }
        dt <- as.data.table(dt)
    }
    if("entry" %in% names(dt) == FALSE){
        dt[, c("entry") := default.start]
    }
    if(object$reverse){
        dt[,c("status") := 1-.SD$status]
    }
    ## normalize names
    data.table::setnames(dt, old = c("entry","time"), new = c("start","stop"))
    data.table::setcolorder(x = dt, neworder = c("start","stop","status"))
    ## ** add x
    if(length(object$continuous.predictors)>0){
        stop("Method coxModelFrame for prodlim object only implemented for discrete predictor and not for continuous predictors. \n")
    }

    strata.var <- object$discrete.predictors
    if(length(strata.var)==0){
        dt[,c("strata") := factor(1)]
    }else{
        if("strata" %in% strata.var){
            stop("The variables in the linear predictor should not be named \"strata\" \n")
        }
        dt.X <- as.data.table(object$model.matrix[originalDataOrder,,drop=FALSE])
        dt[,c("strata") := interaction(dt.X, sep = ", ")]
    }   
    return(dt)
}


## ** coxModelFrame.GLMnet
#' @rdname coxModelFrame
#' @method coxModelFrame GLMnet
#' @export
coxModelFrame.GLMnet <- function(object,center = FALSE){
    # FIXME: argument center not used? 
    default.start <- 0
    dt <- data.table(start = default.start,
                     stop = object$y[, "time"],
                     status = object$y[, "status"])
    dt <- cbind(dt, object$X)
    if("strata" %in% names(attributes(object$y))){
        dt$strata <- attributes(object$y)$strata
    } else{
        dt[, strata := as.factor(1)]
    }
    dt
}



######################################################################
### coxModel.Frame.R ends here
