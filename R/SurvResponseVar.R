### SurvResponseVar.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:37) 
## Version: 
## Last-Updated: Apr 27 2025 (07:37) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## SurvResponseVar
#' @title Extract the time and event variable from a Cox model
#' @description Extract the time and event variable from a Cox model
#' @name SurvResponseVar
#' 
#' @param formula a formula
#'
#' @author Brice Ozenne broz@@sund.ku.dk
#'
#' @examples
#' \dontrun{
#' SurvResponseVar(Surv(time,event)~X1+X2)
#' SurvResponseVar(Hist(time,event==0)~X1+X2)
#' SurvResponseVar(Surv(start,time, status,type="counting") ~ X3+X5)
#' SurvResponseVar(Surv(start,event=status, time2=time,type="counting") ~ X3+X5)
#'
#' SurvResponseVar(survival::Surv(start,event=status, time2=time,type="counting") ~ X3+X5)
#' SurvResponseVar(status ~ X3+X5)
#' SurvResponseVar(I(status == 1) ~ X3+X5)
#' SurvResponseVar(list(Hist(time, event) ~ X1+X6,Hist(time, event) ~ X6))
#' }
SurvResponseVar <- function(formula){

    if(inherits(x=formula,what="call")){
        formula <- eval(formula)
    }
    
    if(inherits(x=formula,what="list")){
        if(any(unlist(lapply(formula,class))!="formula")){
            stop("When a list, the elements of argument \'formula\' must all be a formula \n")
        }
        
        ls.formula <- lapply(formula, SurvResponseVar)
        if(length(ls.formula)>1){
            test.identical <- sapply(ls.formula[-1], function(x){identical(ls.formula[[1]],x)})
            if(all(test.identical)){
                return(ls.formula[[1]])
            }else{
                stop("Different left hand side in the formulae \n")
            }
        }else{
            return(ls.formula[[1]])
        }
        
    }else if(!inherits(formula,"formula")){
        stop("Argument \'formula\' can either be a formula or a list of formula \n")
    }

    ls.formula <- as.list(formula[[2]])
    if(length(ls.formula)==1){ ## deal with the case of simple formula, e.g. Y~X
        return(list(status = all.vars(formula)[1]))
    }
    
    ## extra names of the inputs
    n.input <- length(ls.formula) - 1
    all.input <- names(ls.formula)[-1]
    if(is.null(all.input)){
        all.input <- rep("",n.input)
    }
    

    ## extract arguments
    operator <- deparse(ls.formula[[1]])
    if(operator %in% c("survival::Surv","Surv")){
        ## names(as.list(args(survival::Surv)))
        if( ("time2" %in% all.input) || ( sum(all.input=="") >= (1 + ("time" %in% all.input == FALSE) + ("event" %in% all.input == FALSE)) )){
            all.args <- c(entry = "time", time = "time2", status = "event", type = "type", origin = "origin")
        }else{
            all.args <- c(time = "time", status = "event", type = "type", origin = "origin")
        }
        
    }else if(operator %in% c("riskRegression::Hist","Hist")){
        ## names(as.list(args(Hist)))
        all.args <- c(time = "time", status = "event", entry = "entry", id = "id", cens.code = "cens.code", addInitialState = "addInitialState")
    }else if(operator %in% "I"){
        xx <- all.vars(formula)[1]
        attr(xx, "call") <- deparse(formula[[2]])
        return(list(status = xx))
    }else{
        stop("The left side of the formula must contain Hist() or Surv() \n") 
    }

    ## associate named inputs
    remaining.input <- all.input
    remaining.args <- all.args
    index.input <- 1:n.input
    if(any(all.input!="")){
        for(iNames in all.input[all.input!=""]){
            remaining.args <- remaining.args[remaining.args != iNames]
            index.input <- index.input[remaining.input != iNames]
            remaining.input <- remaining.input[remaining.input != iNames]
        }
    }
    ## set names to unnamed inputs
    all.input[index.input] <- remaining.args[1:length(index.input)]

    all.args <- all.args[all.args %in% all.input]
    for(iNewname in 1:length(all.args)){ ## iNewname <- 1
        all.input[all.input == all.args[iNewname]] <- names(all.args)[iNewname]
    }
    
    ## output
    names(ls.formula) <- c("operator",all.input)
    
    for(iArg in 1:length(ls.formula)){ ## iArg <- 3
        if(inherits(x=ls.formula[[iArg]],what="name")){
            ls.formula[[iArg]] <- deparse(ls.formula[[iArg]])
        }else if(inherits(x=ls.formula[[iArg]],what="call")){
            xx <- as.character(ls.formula[[iArg]])[2]
            attr(xx,"call") <- deparse(ls.formula[[iArg]])
            ls.formula[[iArg]] <- xx
        }
    }

  ## export
  return(ls.formula)
  
}

######################################################################
### SurvResponseVar.R ends here
