### ate_checkArgs.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: May 14 2025 (08:49) 
## Version: 
## Last-Updated: May 14 2025 (15:26) 
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

## * ate_checkArgs
ate_checkArgs <- function(call,
                          object.event,
                          object.treatment,
                          object.censor,
                          mydata,
                          formula,
                          treatment,
                          strata,
                          contrasts,
                          allContrasts,
                          times,
                          cause,
                          se,
                          iid,
                          data.index,
                          return.iid.nuisance,
                          band,
                          B,
                          confint,
                          seed,
                          handler,
                          mc.cores,
                          cl,
                          verbose,
                          n.censor,
                          level.censoring,
                          level.states,
                          estimator,
                          eventVar.time,
                          eventVar.status,
                          censorVar.time,
                          censorVar.status,
                          n.obs){

    options <- riskRegression.options()

    ## ** times
    if(inherits(object.event,"glm") && length(times)!=1){
        stop("Argument \'times\' has no effect when using a glm object \n",
             "It should be set to NA \n")        
    }

    ## ** status
    if(eventVar.status %in% names(mydata) == FALSE){
        stop("The data set does not seem to have a variable called \'",eventVar.status,"\' (argument: event[2]). \n")
    }
    if(inherits(object.event,"glm") && is.na(cause)){
        stop("Argument \'cause\' should not be specified when using a glm model in argument \'event\'")
    }
    if(length(cause)>1 || is.na(cause)){
        if(!is.null(object.event$na.action)){
            stop("Cannot guess which value of the variable \"",eventVar.status,"\" corresponds to the outcome of interest \n",
                 "Likely due to missing data in the dataset used to fit the survival model. \n")
        }else{
            stop("Cannot guess which value of the variable \"",eventVar.status,"\" corresponds to the outcome of interest \n",
                 "Please specify the argument \'cause\'. \n")
        }
    }
    if(cause %in% mydata[[eventVar.status]] == FALSE){
        stop("Value of the argument \'cause\' not found in column \"",eventVar.status,"\" \n")
    }

    ## ** estimator
    if(any(is.na(estimator))){
        stop("No model for the outcome/treatment/censoring has been specified. Cannot estimate the average treatment effect. \n")
    }
    
    valid.estimator <- c("GFORMULA","IPTW,IPCW","IPTW","AIPTW,AIPCW","AIPTW")
    if(any(estimator %in% valid.estimator == FALSE)){
        stop("Incorrect value(s) for argument \'estimator\': \"",paste(estimator[estimator %in% valid.estimator == FALSE], collapse = "\" \""),"\"\n",
             "Valid values: \"G-formula\", \"IPTW\", \"AIPTW\" \n")
    }
    if(attr(estimator,"GFORMULA")){
        if(is.null(object.event)){
            stop("Using a G-formula/AIPTW estimator requires to specify a model for the outcome (argument \'event\') \n")
        }
    }
    if(attr(estimator,"IPTW")){
        if(is.null(object.treatment)){
            stop("Using a ITPW/AIPTW estimator requires to specify a model for the treatment allocation (argument \'treatment\') \n")
        }
    }
    if(attr(estimator,"IPCW")){
        if(is.null(object.censor)){
            stop("Using a ITPW/AIPTW estimator requires to specify a model for the censoring mechanism (argument \'censor\') in presence of right-censoring \n")
        }
    }
    if(length(attr(estimator,"monotone"))>1){
        stop("monotonicity constrain should be request for all or none of the IPW estimators \n")
    }
    
    ## ** object.event
    if(!is.null(object.event)){
        ## look for all function in riskRegression
        ## unclass(lsf.str(envir = asNamespace("riskRegression"), all = T))
        vec.candidateMethods <- paste("predictRisk",class(object.event),sep=".")
        for(iCandidate in vec.candidateMethods){
            test <- try(is.function(eval(parse(text=iCandidate))), silent = TRUE)            
            if(!inherits(test, "try-error") && test==TRUE){candidateMethods <- iCandidate; break}
        }
        if (inherits(test, "try-error") || test==FALSE){
            stop(paste("Could not find predictRisk S3-method for ",paste(class(object.event),collapse=", "),sep=""))
        }
        if((inherits(object.event,"wglm") || inherits(object.event,"CauseSpecificCox")) && (cause %in% object.event$causes == FALSE)){
            stop("Argument \'cause\' does not match one of the available causes: ",paste(object.event$causes,collapse=" "),"\n")
        }
        if(inherits(object.event,"wglm") && (cause != object.event$theCause)){
            stop("Argument \'cause\' does not match the one of \'object.event\' \n",
                 "Consider re-fitting the model with appropriate cause.")
        }
        if(any(is.na(level.censoring)) && any(na.omit(level.censoring) == cause)){
            stop("The cause of interest must differ from the level indicating censoring (in the outcome model) \n")
        } 
        if(any(is.na(coef(object.event)))){
            stop("Cannot handle missing values in the model coefficients (event) \n")
        }
        if(!is.null(treatment) && identical(eventVar.status, treatment)){
            stop("The treatment variable has the same name as the event variable. \n")
        }
    }
    
    ## ** object.treatment    
    if(!is.null(object.treatment)){
        if(!inherits(object.treatment,"multinom") && (!inherits(object.treatment,"glm") || object.treatment$family$family!="binomial")){
            stop("Argument \'treatment\' must be a logistic regression (glm object) or a multinomial regression (nnet::multinom)\n",
                 " or a character variable giving the name of the treatment variable. \n")
        }

        if(any(levels(mydata[[treatment]]) %in% mydata[[treatment]] == FALSE)){
            mistreat <- levels(mydata[[treatment]])[levels(mydata[[treatment]]) %in% mydata[[treatment]] == FALSE]
            stop("Argument \'mydata\' needs to take all possible treatment values when using IPTW/AIPTW \n",
                 "missing treatment values: \"",paste(mistreat,collapse="\" \""),"\"\n")
        }

        if(any(is.na(coef(object.treatment)))){
            stop("Cannot handle missing values in the model coefficients (object.treatment) \n")
        }

        if(inherits(object.treatment,"glm") && object.treatment$family$family == "binomial" && length(unique(mydata[[treatment]]))>2){
            stop("Cannot use a logistic regression in argument \"treatment\" when there are more than two treatment categories \n",
                 "Consider subsetting the dataset to have only two treatment categories or using a multinomial regression (nnet::multinom).")
        }
    }

    ## ** object.censor
    if(attr(estimator,"IPCW")){
        ## require censoring model
        if(!inherits(object.censor,"coxph") && !inherits(object.censor,"cph") && !inherits(object.censor,"phreg") && !inherits(object.censor,"prodlim")){
            stop("Argument \'object.censor\' must be a Cox model \n")
        }
        
        if(inherits(object.event,"glm") && any(n.censor > 0) && all(object.event$prior.weights==1)){
            stop("Argument \'object.event\' should not be a standard logistic regression in presence of censoring \n")
        }

        if(!identical(censorVar.time,eventVar.time)){
            stop("The time variable should be the same in \'object.event\' and \'object.censor\' \n")
        }
        if(any(cause == level.censoring)){
            stop("The level indicating censoring should differ between the outcome model and the censoring model \n",
                 "maybe you have forgotten to set the event type == 0 in the censoring model \n")
        }
        if(!identical(censorVar.status,eventVar.status)){
            if(length(censorVar.status)==1 && length(eventVar.status) == 1 && censorVar.status == eventVar.status){
                ## ok only difference is in attributes
            }else if(sum(diag(table(mydata[[censorVar.status]] == level.censoring, mydata[[eventVar.status]] %in% level.states == FALSE)))!=NROW(mydata)){
                stop("The status variables in object.event and object.censor are inconsistent \n")
            }
        }
        if(any(is.na(coef(object.censor)))){
            stop("Cannot handle missing values in the model coefficients (object.treatment) \n")
        }
    }

    ## ** bootstrap
    if(B>0){
        if(se==FALSE){
            warning("Argument 'se=0' means 'no standard errors' so number of bootstrap repetitions is forced to B=0.")
        }
        if(iid==TRUE){
            stop("Influence function cannot be computed when using the bootstrap approach \n",
                 "Either set argument \'iid\' to FALSE to not compute the influence function \n",
                 "or set argument \'B\' to 0 \n")
        }
        if(band==TRUE){
            stop("Confidence bands cannot be computed when using the bootstrap approach \n",
                 "Either set argument \'band\' to FALSE to not compute the confidence bands \n",
                 "or set argument \'B\' to 0 to use the estimate of the asymptotic distribution instead of the bootstrap\n")
        }
        if(!is.null(object.event) && is.null(object.event$call)){
            stop("Argument \'event\' does not contain its own call, which is needed to refit the model in the bootstrap loop.")
        }
        if(!is.null(object.treatment) && is.null(object.treatment$call)){
            stop("Argument \'treatment\' does not contain its own call, which is needed to refit the model in the bootstrap loop.")
        }
        if(!is.null(object.censor) && is.null(object.censor$call)){
            stop("Argument \'censor\' does not contain its own call, which is needed to refit the model in the bootstrap loop.")
        }
        ## if(any(data.index != 1:na.omit(n.obs)[1])){
        ##     stop("Argument \'data.index\' cannot be used when using bootstrap resampling \n")
        ## }
        ## if((Sys.info()["sysname"] == "Windows") && (handler == "mclapply") && (mc.cores>1) ){
        ## stop("mclapply cannot perform parallel computations on Windows \n",
        ## "consider setting argument handler to \"foreach\" \n")
        ## }
        max.cores <- parallel::detectCores()
        if(mc.cores > max.cores){
            stop("Not enough available cores \n","available: ",max.cores," | requested: ",mc.cores,"\n")
        }
    }

    ## ** functional delta method
    if(B==0 && (se|band|iid)){
        if(length(unique(stats::na.omit(n.obs[-1])))>1){
            stop("Arguments \'",paste(names(stats::na.omit(n.obs[-1])), collapse ="\' "),"\' must be fitted using the same number of observations \n")
        }
        test.data.index <- any(data.index %in% 1:max(n.obs[-1],na.rm = TRUE) == FALSE)
        if(is.null(data.index) || (is.null(call$data.index) && test.data.index)){
            stop("Incompatible number of observations between argument \'data\' and the dataset from argument(s) \'",paste(names(stats::na.omit(n.obs[-1])), collapse ="\' "),"\' \n",
                 "Consider specifying argument \'data.index\' \n \n")
        }
        if(!is.numeric(data.index) || any(duplicated(data.index)) || any(is.na(data.index)) || test.data.index){
            stop("Incorrect specification of argument \'data.index\' \n",
                 "Must be a vector of integers between 0 and ",max(n.obs[-1],na.rm = TRUE)," \n")
        }
        
        if(!is.null(object.event)){
            if ("average.iid" %in% names(formals(candidateMethods)) == FALSE){
                stop(paste("The method ",candidateMethods," is missing an argument average.iid to be used for the functional delta method. \n",
                           "Set argument \'B\' to a positive integer to use a bootstrap instead \n",sep=""))
            }
        }
        if(!is.null(object.treatment) && stats::nobs(object.treatment)!=NROW(mydata)){ ## note: only check that the datasets have the same size
            stop("Argument \'treatment\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                 "(discrepancy found in number of rows) \n")
        }

        if(!is.null(object.censor) && coxN(object.censor)!=NROW(mydata)){ ## note: only check that the datasets have the same size
                stop("Argument \'censor\' must be have been fitted using argument \'data\' for the functional delta method to work\n",
                     "(discrepancy found in number of rows) \n")
        }

    }

    
    ## ** treatment
    if(!is.null(treatment)){
        if(treatment %in% names(mydata) == FALSE){
            stop("The data set does not seem to have a variable ",treatment," (argument: object.treatment). \n")
        }
        if(is.numeric(mydata[[treatment]])){
            stop("The treatment variable must be a factor variable. \n",
                 "Convert treatment to factor, re-fit the object using this new variable and then call ate. \n")
        }
    }else if(is.null(strata)){
        stop("The treatment variable must be specified using the argument \'treatment\' \n")
    }

    ## ** contrasts/allContrasts
    if(!is.null(contrasts)){
        if(any(contrasts %in% unique(mydata[[treatment]]) == FALSE)){
            stop("Incorrect values for the argument \'contrasts\' \n",
                 "Possible values: \"",paste(unique(mydata[[treatment]]),collapse="\" \""),"\" \n")
        }
    
    }
    if(!is.null(allContrasts)){
        if(any(allContrasts %in% unique(mydata[[treatment]]) == FALSE)){
            stop("Incorrect values for the argument \'allContrasts\' \n",
                 "Possible values: \"",paste(unique(mydata[[treatment]]),collapse="\" \""),"\" \n")
        }
        if(!is.matrix(allContrasts) || NROW(allContrasts)!=2){
            stop("Argument \'allContrasts\' must be a matrix with 2 rows \n")
        }
    }

    ## ** strata
    if(!is.null(strata)){
        if(attr(estimator, "IPTW")){
            stop("Argument \'strata\' only compatible with the G-formula estimator \n")
        }
        if(any(strata %in% names(mydata) == FALSE)){
            stop("The data set does not seem to have a variable \"",paste0(strata, collapse = "\" \""),"\" (argument: strata). \n")
        }
        if(length(strata) != 1){
            stop("Argument strata should have length 1. \n")
        }
    }
    
    ## ** event time
    if(!is.na(eventVar.time)){
        if(eventVar.time %in% names(mydata) == FALSE){
            stop("The data set does not seem to have a variable ",eventVar.time," (argument: object.event[1]). \n")
        }
        if(is.na(cause)){
            stop("Argument \'cause\' not specified\n")
        }
        data.status <- mydata[[eventVar.status]]
        if(!is.null(treatment)){
            data.strata <- mydata[[treatment]]
        }else{
            data.strata <- mydata[[strata]]
        }
        if(is.factor(data.strata)){
            data.strata <- droplevels(data.strata)
        }
        freq.event <- tapply(data.status, data.strata, function(x){mean(x==cause)})
        count.event <- tapply(data.status, data.strata, function(x){sum(x==cause)})
        
        if(any(count.event < 5)  ){
            warning("Rare event \n")
        }
        if(any(mydata[[eventVar.time]]<0)){
            stop("The time to event variable should only take positive values \n")
        }
    }

    ## ** iid.nuisance    
    if((return.iid.nuisance == FALSE) & (iid == TRUE) & (attr(estimator, "augmented") == FALSE)){
        warning("Ignoring the uncertainty associated with the estimation of the nuisance parameters may lead to inconsistent standard errors. \n",
                "Consider using a double robust estimator or setting the argument \'known.nuisance\' to TRUE \n")
    }
    

    ## ** output
    return(TRUE)
}

######################################################################
### ate_checkArgs.R ends here
