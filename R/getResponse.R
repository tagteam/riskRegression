getResponse <- function(formula,cause,data,vars){
    ## case 1: continuous 
    ## case 2: binary
    ## case 3: ordinal
    ## case 4: survival (Surv, Hist) 
    ## case 5: competing risks (Hist)
    ## case 6: updating (?)
    if (length(vars)==1){
        ## case 1,2,3
        m <- stats::model.frame(formula=formula,data=data,na.action=na.fail)
        response <- stats::model.response(m)
        rlevs <- unique(response)
        if (is.factor(response) || length(rlevs)==2){
            if (is.factor(response)) rlevs <- levels(response)
            if (length(rlevs)==2) {
                if (missing(cause)||is.null(cause)){
                    if (is.factor(response)){
                        cause <- levels(response)[2]
                    } else{
                        cause <- NULL
                        if (all(rlevs %in% c("0","1"))) cause=1
                        if (all(rlevs %in% c("1","2"))) cause=2
                        if (is.null(cause)){
                            warning("Outcome cause of interest not specified. Using first occuring value for now.")
                            cause <- rlevs[1]
                        }
                    }
                }
                ## coercing to 0/1 variable
                response <- data.table(ReSpOnSe=as.numeric(response==cause))
                data.table::setattr(response,"states",c("0","1"))
                data.table::setattr(response,"event","1")
                data.table::setattr(response,"model","binary")
            }
            else{
                data.table::setattr(response,"model","multi.level")
                stop("Methods for factors with more than two classes are not\n (not yet) implemented.")
            }
        }
        else{
            data.table::setattr(response,"model","continuous")
            attr(response,"event") <- NULL
        }
    }
    else{
        responseNames <- all.names(formula)
        if (responseNames[2] %in% c("Surv","Hist")){
            ## case 4,5
            ## replace Surv by Hist to use the attributes of
            ## the Hist object
            requireNamespace("prodlim")
            if (responseNames[2]=="Surv"){
                formula[[2]][[1]] <- as.name("Hist")
            }
            response <- unclass(eval(formula[[2]],data))
            ## m <- stats::model.frame(formula=formula,
            ## data=data,
            ## na.action=na.fail)
            ## response <- unclass(stats::model.response(m))
        }
        else{
            stop("Cannot assign response type.")
        }
    }
    response
}
