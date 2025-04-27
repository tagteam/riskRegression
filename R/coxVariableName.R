### coxVariableName.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Apr 27 2025 (07:29) 
## Version: 
## Last-Updated: Apr 27 2025 (07:30) 
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
# {{{ coxVariableName
## * coxVariableName
#' @title Extract variable names from a model
#' @description Extract the name of the variables belonging to the linear predictor or used to form the strata
#' @name coxVariableName
#' 
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param model.frame [data.frame] dataset containing all the relevant variables (entry, time to event, type of event, variables in the linear predictor, strata).
#' Output from \code{coxModelFrame}.
#' 
#' @author Brice Ozenne broz@@sund.ku.dk

#' @rdname coxVariableName
#' @export
coxVariableName <- function(object, model.frame){

    ## ** get formula
    ff <- coxFormula(object)

    ## ** get special operators
    special.object <- coxSpecial(object)

    ## ** get names of the start/stop/status variables
    out <- SurvResponseVar(ff)

    ## ** identify special terms
    xterms <- delete.response(terms(ff,
                                    special = unlist(special.object),
                                    data = model.frame))    
    ls.specials <- attr(xterms,"specials")
    if(any(attr(xterms,"order")>1) && length(ls.specials$strata)>0){
        ## the presence of an interaction makes attr(xterms,"specials")$strata consider the 'main' effects which may not automatically be there
        ## e.g. coxph(Surv(time,event)~X1*X3+strata(X5),data=train, y=TRUE, x = TRUE) lead to design matrix [X1,X3,strata(X5),X1:X3] and attr(xterms,"specials")$strata being 3 (as expected)
        ## e.g. coxph(Surv(time,event)~X1:X3+strata(X5),data=train, y=TRUE, x = TRUE) lead to design matrix [strata(X5),X1:X3] and attr(xterms,"specials")$strata being 3 (incorrect)
        ## This is a fix
        term.strata <- names(which(colSums(attr(xterms,"factors")[attr(xterms,"specials")$strata,,drop=FALSE])>0))
        ls.specials$strata <- intersect(which(attr(xterms,"order")==1), ## strata terms are always main terms (handle X2*strata(X1)--> strata(X1), X2:strata(X1))
                                        which(attr(xterms,"term.labels") %in% term.strata))
    }
    n.specials <- length(unlist(ls.specials))
    
    n.xterms <- length(attr(xterms,"term.labels"))

    ## ** linear predictor
    if(n.xterms>n.specials && !inherits(object,"prodlim")){
        out$lpvars <- names(coef(object)) ##attr(xterms.lp, "term.labels") not ok for categorical variables with coxph
        if(n.specials>0){
            out$lp.sterms <- stats::drop.terms(xterms,unlist(ls.specials))
        }else{
            out$lp.sterms <- xterms
        }
        out$lpvars.original <- all.vars(out$lp.sterms)
        
    }else{
        out["lpvars"] <- list(NULL)
        out["lpvars.original"] <- list(NULL)
        out["lp.sterms"] <- list(NULL)
    }

    ## ** strata variables
    out$strataspecials <- special.object$strata
    if(inherits(object,"prodlim")){
        out$is.strata <- n.xterms>0
    }else{
        out$is.strata <- length(ls.specials[[special.object$strata]])>0
    }

    if(out$is.strata){
        if(!inherits(object,"prodlim") && length(ls.specials[[special.object$strata]])!=n.xterms){
            out$strata.sterms <- stats::drop.terms(xterms,setdiff(1:n.xterms,ls.specials[[special.object$strata]]))
        }else{
            out$strata.sterms <- xterms
        }
        out$stratavars <- attr(out$strata.sterms, "term.labels")
        out$stratavars.original <- all.vars(out$strata.sterms)
        out$strata.levels <- coxStrataLevel(object)

        if(length(out$strata.levels)==1){
            out$is.strata <- FALSE
        }
    }else{
        out["stratavars"] <- list(NULL)
        out["stratavars.original"] <- list(NULL)
        out["strata.sterms"] <- list(NULL)
        out$strata.levels <- factor("1") ## not NULL - coherent with coxStrata
    }

    ## ** export
    return(out)
} 
# }}}

######################################################################
### coxVariableName.R ends here
