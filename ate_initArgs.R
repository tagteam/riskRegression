### ate_initArgs.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 14 2025 (08:49) 
## Version: 
## Last-Updated: May 14 2025 (08:55) 
##           By: Thomas Alexander Gerds
##     Update #: 2
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * ate_initArgs
ate_initArgs <- function(object.event,
                         object.treatment,
                         object.censor,
                         landmark,
                         formula = NULL,
                         estimator,
                         mydata, ## instead of data to avoid confusion between the function data and the dataset when running the bootstrap
                         data.index,
                         cause,
                         times,
                         handler,
                         product.limit,
                         store){


    ## ** user-defined arguments
    ## handler
    handler <- match.arg(handler, c("foreach","mclapply","snow","multicore"))
    ## data.index
    if(is.null(data.index)){data.index <- 1:NROW(mydata)}

    ## ** fit a regression model when the user specifies a formula
    ## Event
    if(missing(object.event) || is.null(object.event)){
        formula.event <- NULL
        model.event <- NULL
    }else if(inherits(object.event,"formula")){ ## formula
        formula.event <- object.event
        if(any(grepl("Hist(",object.event, fixed = TRUE))){
            model.event <- do.call(CSC, args = list(formula = object.event, data = mydata))
        }else if(any(grepl("Surv(",object.event, fixed = TRUE))){
            model.event <- do.call(survival::coxph, args = list(formula = object.event, data = mydata, x = TRUE, y = TRUE))
        }else{
            model.event <- do.call(stats::glm, args = list(formula = object.event, data = mydata, family = stats::binomial(link = "logit")))
        }
    }else if(is.list(object.event) && all(sapply(object.event, function(iE){inherits(iE,"formula")}))){   ## list of formula
        formula.event <- object.event
        model.event <- CSC(object.event, data = mydata)
    }else if(inherits(object.event,"glm") || inherits(object.event,"wglm") || inherits(object.event,"CauseSpecificCox")){ ## glm / CSC
        formula.event <- stats::formula(object.event)
        model.event <- object.event
    }else if(inherits(object.event,"coxph") || inherits(object.event,"cph") || inherits(object.event,"phreg")  || inherits(object.event,"prodlim")){ ## Cox
        formula.event <- coxFormula(object.event)
        model.event <- object.event
    }else if(is.character(object.event)){ ## character (time,event) i.e. no model
        formula.event <- NULL
        model.event <- NULL
    }else{
        message("Unknown event model. \n",
                "Recognized event models: stats::glm, riskRegression::wglm, survival::coxph, rms::cph, mets::phreg, prodlim::prodlim, riskRegression::wglm. \n")
        formula.event <- stats::formula(object.event)
        model.event <- object.event
    }

    ## Treatment
    if(missing(object.treatment) || is.null(object.treatment)){
        formula.treatment <- NULL
        model.treatment <- NULL
    }else if(inherits(object.treatment,"formula")){
        model.treatment <- do.call(stats::glm, args = list(formula = object.treatment, data = mydata, family = stats::binomial(link = "logit")))        
        formula.treatment <- object.treatment
    }else if(inherits(object.treatment,"glm") || inherits(object.treatment,"nnet")){
        formula.treatment <- stats::formula(object.treatment)
        model.treatment <- object.treatment
    }else if(is.character(object.treatment)){
        formula.treatment <- NULL
        model.treatment <- NULL
    }else{
        message("Unknown treatment model: \"",paste(class(object.treatment), collapse="\", \""),"\". \n",
                "Recognized treatment models: stats::glm, nnet::multinom. \n")
        formula.treatment <- stats::formula(object.treatment)
        model.treatment <- object.treatment
    }

    ## Censoring
    if(missing(object.censor) || is.null(object.censor)){
        if(inherits(object.event,"wglm") && any(object.event$n.censor>0)){
            formula.censor <- coxFormula(object.event$model.censor)
            model.censor <- object.event$model.censor
        }else{
            formula.censor <- NULL
            model.censor <- NULL
        }
    }else if(inherits(object.censor,"formula")){
        formula.censor <- object.censor
        model.censor <- do.call(survival::coxph, args = list(formula = object.censor, data = mydata, x = TRUE, y = TRUE))
    }else if(inherits(object.censor,"coxph") || inherits(object.censor,"cph") || inherits(object.censor,"phreg") || inherits(object.censor,"prodlim")){
        formula.censor <- coxFormula(object.censor)
        model.censor <- object.censor
    }else if(is.character(object.censor)){
        formula.censor <- NULL
        model.censor <- NULL
    }else{
        message("Unknown censoring model. \n",
                "Recognized censoring models: survival::coxph, rms::cph, mets::phreg, prodlim::prodlim. \n")
        formula.censor <- stats::formula(object.censor)
        model.censor <- object.censor
    }

    ## ** times
    if(missing(times) && ((inherits(model.event,"glm") || inherits(model.event,"multinom")) || (inherits(model.treatment,"glm") && is.null(model.event) && is.null(model.censor)))){
        times <- NA
    }
    
    ## ** identify event, treatment and censoring variables
    if(inherits(model.censor,"coxph") || inherits(model.censor,"cph") || inherits(model.censor,"phreg") || inherits(model.censor,"prodlim")){
        censoringMF <- coxModelFrame(model.censor)
        test.censor <- censoringMF$status == 1
        n.censor <- sapply(times, function(t){sum(test.censor * (censoringMF$stop < t))})
        info.censor <- SurvResponseVar(coxFormula(model.censor))
        censorVar.status <- info.censor$status
        censorVar.time <- info.censor$time
        if(inherits(model.censor,"prodlim")){
            level.censoring <- attr(model.censor$model.response,"cens.code")
            if(is.numeric(mydata[[censorVar.status]])){
                level.censoring <- as.numeric(level.censoring)
            }
            ## Not reliable due to re-ordering and reverse arguments
            ## unique(mydata[[censorVar.status]][model.censor$model.response[,"status"]==1])
        }else{
            mydata.censor <- try(eval(model.censor$call$data), silent = TRUE)
            if(!inherits(mydata.censor, "try-error") && (inherits(mydata.censor,"list") || inherits(mydata.censor,"data.frame"))){
                level.censoring <- unique(mydata.censor[[censorVar.status]][model.censor$y[,2]==1])
            }else if(is.numeric(mydata[[censorVar.status]])){
                level.censoring <- 0
            }else if(is.character(mydata[[censorVar.status]])){
                level.censoring <- sort(unique(mydata[[censorVar.status]]))[1]
            }else if(is.factor(mydata[[censorVar.status]])){
                level.censoring <- levels(mydata[[censorVar.status]])[1]
            }            
        }
    }else{ ## G-formula or IPTW (no censoring)
        censorVar.status <- NA
        censorVar.time <- NA
                
        if(inherits(model.event,"CauseSpecificCox")){
            test.censor <- model.event$response[,"status"] == 0        
            n.censor <- sapply(times, function(t){sum(test.censor * (model.event$response[,"time"] < t))})
            level.censoring <- attr(model.event$response,"cens.code")
        }else if(inherits(model.event,"coxph") || inherits(model.event,"cph") || inherits(model.event,"phreg") || inherits(model.event,"prodlim")){ 
            censoringMF <- coxModelFrame(model.event)
            test.censor <- censoringMF$status == 0
            n.censor <- sapply(times, function(t){sum(test.censor * (censoringMF$stop < t))})
            level.censoring <- 0
        }else if(inherits(model.event,"wglm")){
            n.censor <- object.event$n.censor
            level.censoring <- setdiff(unique(object.event$data[[object.event$var.outcome]]), object.event$causes)
        }else{
            n.censor <- rep(0,length(times))
        }
        if(all(n.censor==0)){level.censoring <- NA}
    }

    ## treatment
    if(missing(object.treatment) || is.null(object.treatment)){
        treatment <- NULL
    }else if(!is.null(formula.treatment)){
        treatment <- all.vars(formula.treatment)[1]
    }else if(is.character(object.treatment)){
        treatment <- object.treatment
    }

    ## event
    if(inherits(model.event,"CauseSpecificCox")){
        responseVar <- SurvResponseVar(formula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status
        if(is.na(cause)){
            cause <- model.event$theCause
        }
        if(is.null(product.limit)){product.limit <- TRUE}
        level.states <- model.event$cause
    }else if(inherits(model.event,"prodlim")){
        responseVar <- SurvResponseVar(formula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status
        level.states <- attr(model.event$model.response,"state")
        if(is.numeric(mydata[[eventVar.status]])){
            level.states <- as.numeric(level.states)
        }
        if(is.na(cause)){ ## handle Hist(time,event > 0) ~ ...
            cause <- level.states[1]
        }
        if(is.null(product.limit)){product.limit <- TRUE}
    }else if(inherits(model.event,"coxph") || inherits(model.event,"cph") || inherits(model.event,"phreg")){
        responseVar <- SurvResponseVar(formula.event)
        eventVar.time <- responseVar$time
        eventVar.status <- responseVar$status

        if(is.na(cause)){ ## handle Surv(time,event > 0) ~ ...
            if(any(model.event$y[,NCOL(model.event$y)]==1)){
                modeldata <- try(eval(model.event$call$data), silent = TRUE)
                if(inherits(x=modeldata,what="try-error")){
                    if(NROW(mydata)==NROW(model.event$y)){
                        cause <- unique(mydata[[eventVar.status]][model.event$y[,NCOL(model.event$y)]==1])
                    }else{
                        stop("Cannot guess which value of the variable \"",eventVar.status,"\" corresponds to the outcome of interest \n",
                             "Number of rows differs between the input data and the model frame (object$y). \n")
                    }                    
                }else{
                    if(NROW(modeldata)==NROW(model.event$y)){
                        cause <- unique(modeldata[[eventVar.status]][model.event$y[,NCOL(model.event$y)]==1])
                    }else{
                        stop("Cannot guess which value of the variable \"",eventVar.status,"\" corresponds to the outcome of interest \n",
                             "Number of rows differs between the training data (",deparse(model.event$call$data),") and the model frame (object$y), possibly due to missing values. \n")
                    }
                }
            }
        }
        if(is.null(product.limit)){product.limit <- FALSE}
        level.states <- 1
    }else if(inherits(model.event,"glm") || inherits(model.event,"multinom")){
        eventVar.time <- as.character(NA)
        eventVar.status <- all.vars(formula.event)[1]
        if(is.na(cause)){ ## handle I(Y > 0) ~ ...
            if(inherits(model.event,"glm")){
                cause <- unique(stats::model.frame(model.event)[[eventVar.status]][model.event$y==1])
            }else if(inherits(model.event,"multinom")){
                stop("Argument \'cause\' should be specified for \"multinom\" objects. \n")
            }
        }
        if(is.null(product.limit)){product.limit <- FALSE}
        level.states <- unique(mydata[[eventVar.status]])
    }else if(inherits(object.event,"wglm")){
        eventVar.time <- object.event$var.time
        eventVar.status <- object.event$var.outcome
        if(is.na(cause)){ ## handle I(Y > 0) ~ ...
            cause <- object.event$theCause
        }
        level.states <- object.event$causes
        if(is.null(product.limit)){product.limit <- (length(level.states)>1)}
    }else{ ## no outcome model
        if(identical(names(object.event),c("time","status"))){
            eventVar.time <- object.event[1]
            eventVar.status <- object.event[2]
        }else if(identical(names(object.event),c("status","time"))){
            eventVar.status <- object.event[1]
            eventVar.time <- object.event[2]
        }else if(length(object.event)==2){            
            eventVar.time <- object.event[1]
            eventVar.status <- object.event[2]
        }else if(length(object.event)==1){
            eventVar.time <- NA
            eventVar.status <- object.event[1]
            if(is.na(cause)){
                if(any(mydata[[eventVar.status]] %in% 0:2 == FALSE)){
                    stop("The event variable must be an integer variable taking value 0, 1, or 2. \n")
                }
                cause <- 1
            }
        }
        level.states <- setdiff(unique(mydata[[eventVar.status]]), level.censoring)
        if(is.na(cause)){
            cause <- sort(level.states)[1]
        }
        if(is.null(product.limit)){product.limit <- FALSE}
    }

    ## presence of time dependent covariates
    TD <- switch(class(object.event)[[1]],
                 "coxph"=(attr(object.event$y,"type")=="counting"),
                 "CauseSpecificCox"=(attr(object.event$models[[1]]$y,"type")[1]=="counting"),
                 FALSE)
    browser(skipCalls=1L)
    if(TD){
        fct.pointEstimate <- ATE_TD
    }else{
        landmark <- NULL
        formula <- NULL
        fct.pointEstimate <- ATE_TI
    }

    ## estimator
    if(is.null(estimator)){
        if(!is.null(model.event) && is.null(model.treatment)){
            if(TD){
                estimator <- "GFORMULATD"
            }else{
                estimator <- "GFORMULA"
            }
        }else if(is.null(model.event) && !is.null(model.treatment)){
            if(any(n.censor>0)){
                estimator <- "IPTW,IPCW"
            }else{
                estimator <- "IPTW"
            }
        }else if(!is.null(model.event) && !is.null(model.treatment)){
            if(any(n.censor>0)){
                estimator <- "AIPTW,AIPCW"
            }else{
                estimator <- "AIPTW"
            }
        }else{
            estimator <- NA
        }
        test.monotone <- FALSE
    }else{
        estimator <- toupper(estimator)

        ## should montonicity constraint be enforced on estimators based on IPW
        mestimator <- estimator
        estimator <- gsub("MONOTONE","",estimator)
        
        index.westimator <- which(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW"))
        if(length(index.westimator)>0){
            test.monotone <- unique(grepl(pattern = "MONOTONE",mestimator[index.westimator]))
        }else{
            test.monotone <- FALSE
        }
        if(any(estimator == "IPTW") && any(n.censor>0)){
            estimator[estimator == "IPTW"] <- "IPTW,IPCW"
        }
        if(any(estimator == "AIPTW") && any(n.censor>0)){
            estimator[estimator == "AIPTW"] <- "AIPTW,AIPCW"
        }
        if(any(estimator == "G-FORMULA")){
            estimator[estimator == "G-FORMULA"] <- "GFORMULA"
        }
        if(any(estimator == "G-FORMULATD")){
            estimator[estimator == "G-FORMULATD"] <- "GFORMULA"
        }
    }

    ## which terms are used by the estimator
    estimator.output <- unname(sapply(estimator, switch,
                                      "GFORMULA" = "GFORMULA",
                                      "GFORMULATD" = "GFORMULA",
                                      "IPTW" = "IPTW",
                                      "IPTW,IPCW" = "IPTW",
                                      "AIPTW" = "AIPTW",
                                      "AIPTW,AIPCW" = "AIPTW"
                                      ))

    if(TD){
        attr(estimator.output,"TD") <- TRUE
    }else{
        attr(estimator.output,"TD") <- FALSE
    }

    ## term 1: F_1(\tau|A=a,W) or F_1(\tau|A=a,W) (1-1(A=a)/Prob[A=a|W])
    if(any(estimator %in% c("GFORMULA","GFORMULATD","AIPTW","AIPTW,AIPCW"))){
        attr(estimator.output,"GFORMULA") <- TRUE
    }else{
        attr(estimator.output,"GFORMULA") <- FALSE
    }

    ## term 2 (treatment): 1(A=a)Y(tau)/Prob[A=a|W] or 1(A=a)Y(tau)/Prob[A=a|W] 1(\Delta!=0)/Prob[\Delta!=0]
    if(any(estimator %in% c("IPTW","IPTW,IPCW","AIPTW","AIPTW,AIPCW"))){
        attr(estimator.output,"IPTW") <- TRUE
    }else{
        attr(estimator.output,"IPTW") <- FALSE
    }

    ## term 2 (censoring): 1(A=a)Y(tau)/Prob[A=a|W] 1(\Delta!=0)/Prob[\Delta!=0]
    if(any(estimator %in% c("IPTW,IPCW","AIPTW,AIPCW"))){
        attr(estimator.output,"IPCW") <- TRUE
    }else{
        attr(estimator.output,"IPCW") <- FALSE
    }

    ## term 3: 1(A=a)Y(tau) \int () dM
    if(any(estimator %in% c("AIPTW,AIPCW"))){
        attr(estimator.output,"integral") <- !inherits(object.event,"wglm")
    }else{
        attr(estimator.output,"integral") <- FALSE
    }

    if(any(estimator %in% c("AIPTW","AIPTW,AIPCW"))){
        attr(estimator.output,"augmented") <- TRUE
    }else{
        attr(estimator.output,"augmented") <- FALSE
    }

    ## which estimates should be output?
    if(any(estimator %in% c("GFORMULA","GFORMULATD"))){
        attr(estimator.output,"export.GFORMULA") <- TRUE
    }else{
        attr(estimator.output,"export.GFORMULA") <- FALSE
    }
    if(any(estimator %in% c("IPTW","IPTW,IPCW"))){
        attr(estimator.output,"export.IPTW") <- TRUE
    }else{
        attr(estimator.output,"export.IPTW") <- FALSE
    }
    if(any(estimator %in% c("AIPTW","AIPTW,AIPCW"))){
        attr(estimator.output,"export.AIPTW") <- TRUE
    }else{
        attr(estimator.output,"export.AIPTW") <- FALSE
    }
    attr(estimator.output,"full") <- estimator
    if(!attr(estimator.output,"integral") && any(estimator == "AIPTW,AIPCW")){
        attr(estimator.output,"full")[attr(estimator.output,"full")=="AIPTW,AIPCW"] <- "AIPTW,IPCW"
    }else{
        
    }
    attr(estimator.output,"monotone") <- test.monotone

    ## ** sample size
    n.obs <- c(data = NROW(mydata),
               model.event = if(!is.null(model.event)){coxN(model.event)}else{NA},
               model.treatment = if(!is.null(model.treatment)){stats::nobs(model.treatment)}else{NA},
               model.censor = if(!is.null(model.censor)){coxN(model.censor)}else{NA}
               )

    ## ** store
    store.data <- NULL
    store.iid <- "full"
    store.size.split <- NULL
    if(!is.null(store)){
        if(length(store) > 2){
            stop("Argument \'store\' should contain at most two elements. \n",
                 "For instance store = c(data = \"full\", iid = \"full\") or store = c(data = \"minimal\", iid = \"minimal\").\n")
        }
        if(is.null(names(store)) || any(names(store) %in% c("data","iid", "size.split") == FALSE)){
            stop("Incorrect names for argument \'store\': should be \"data\", \"iid\", \"size.split\". \n",
                 "For instance store = c(data = \"full\", iid = \"full\") or store = c(data = \"minimal\", iid = \"minimal\").\n")
        }
        if("data" %in% names(store) && !is.null(store[["data"]])){
            if(store[["data"]] %in% c("minimal","full") == FALSE){
                stop("Element \"data\" in argument \'store\' should take value \'minimal\' or \'full\'.\n",
                     "For instance store = c(data = \"full\") or store = c(data = \"minimal\").\n")
            }
            store.data <- store[["data"]]
        }
        if("iid" %in% names(store) && !is.null(store[["iid"]])){
            if(store[["iid"]] %in% c("minimal","full") == FALSE){
                stop("Element \"iid\" in argument \'store\' should take value \'minimal\' or \'full\'.\n",
                     "For instance store = c(iid = \"full\") or store = c(iid = \"minimal\").\n")
            }
            store.iid <- store[["iid"]]
        }
        if("size.split" %in% names(store)){
            if(!is.numeric(store[["size.split"]]) || store[["size.split"]] < 0 || (store[["size.split"]] %% 1>0)){
                stop("Element \"size.split\" in argument \'store\' should be a positive integer.\n")
            }
            store.size.split <- store[["size.split"]]
        }
        
    }
    store <- list(data = store.data, iid = store.iid, size.split = store.size.split)

    ## ** output
    return(list(object.event = model.event,
                object.treatment = model.treatment,
                object.censor = model.censor,
                times = times,
                handler = handler,
                estimator = estimator.output,
                formula = formula,
                landmark = landmark,
                fct.pointEstimate = fct.pointEstimate,
                n.obs = n.obs,
                n.censor = n.censor,
                treatment = treatment,
                cause = cause,
                level.censoring = level.censoring,
                level.states = level.states,
                eventVar.time = eventVar.time,
                eventVar.status = eventVar.status,
                censorVar.time = censorVar.time,
                censorVar.status = censorVar.status,
                product.limit = product.limit,
                data.index = data.index,
                store = store
                ))
}

######################################################################
### ate_initArgs.R ends here
