### ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2018 (17:47) 
## Version: 
## Last-Updated: jul  4 2018 (16:48) 
##           By: Brice Ozenne
##     Update #: 316
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * ateRobust (documentation)
#' @title Compute the Average Treatment Effects using Double Robust Approach
#' @description Compute the average treatment effect (ATE) for a Cox regression model using five different approaches:
#' G-formula, inverse probability weighting (IPW), augmented inverse probability weighting (AIPW).
#' For IPW and AIPW both the non-efficient and efficient estimates are output. Experimental!
#' @name ateRobust
#'
#' @param data [data.frame or data.table] Data set in which to evaluate the ATE.
#' @param formula.event [formula] The Cox model for the event of interest. Typically \code{Surv(time,event)~treatment}.
#' @param formula.censor [formula] The Cox model for the censoring. Typically \code{Surv(time,event==0)~treatment}.
#' @param formula.treatment [formula] The logistic model for the treatment. Typically \code{treatment~1}.
#' @param times [numeric] Time point at which to evaluate average treatment effects.
#' @param treatment [character] Name of the treatment variable. Must be binary.
#' @param fitter [character] Routine to fit the Cox regression models.
#' If \code{coxph} use \code{survival::coxph} else use \code{rms::cph}.
#' @param efficient [logical] should the efficient estimates of IPW and AIPW be computed.
#' Can be time consuming for large dataset.
#' @param product.limit [logical] If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param iid [logical] If \code{TRUE} add the influence function to the output.
#' For now the influence function is computed treating the survival probability, censoring probability,
#' and probability of being treated as known.
#' 

## * ateRobust (example)
#' @rdname ateRobust
#' @examples
#'
#' library(survival)
#' library(lava)
#' library(data.table)
#' 
#' ## generative model
#' mSimSurv <- lava::lvm()
#' lava::distribution(mSimSurv, ~a) <- lava::binomial.lvm(p = c(0.5))
#' lava::distribution(mSimSurv, ~w) <- lava::binomial.lvm(p = c(0.5))
#' lava::distribution(mSimSurv, "eventtime") <- lava::coxExponential.lvm(scale = 1)
#' lava::distribution(mSimSurv, "censtime") <- lava::coxExponential.lvm(scale = 1)
#' mSimSurv <- lava::eventTime(mSimSurv, time ~ min(eventtime = 1, censtime = 0), "event")
#' lava::regression(mSimSurv) <- eventtime ~ alpha*a + beta*w
#'
#' ## settings
#' alpha <- 0.3 
#' beta <- 0.4 
#'
#' ## check bias 
#' set.seed(10)
#' dt <- as.data.table(sim(mSimSurv, n = 1e4, p = c(alpha = alpha, beta = beta)))
#' setkeyv(dt, c("a","w"))
#'
#' ## True value
#' psi.TRUE <- c("surv.0" = exp(-1*exp(alpha*0+beta*0))*0.5 + exp(-1*exp(alpha*0+beta*1))*0.5,
#'               "surv.1" = exp(-1*exp(alpha*1+beta*0))*0.5 + exp(-1*exp(alpha*1+beta*1))*0.5)
#' psi.TRUE
#' 
#' ## Approximate true value
#' dt[,.(surv = mean(eventtime>1)),by = c("a")]
#'
#' ## Estimated value
#' res1 <- ateRobust(data = dt, type = "survival",
#'             formula.event = Surv(time, event) ~ strata(a,w),
#'             formula.censor = Surv(time, event==0) ~ strata(a,w),
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             treatment = "a",
#'             efficient = FALSE,
#'             product.limit = FALSE)
#' res1
#'
#' if(FALSE){
#'   e.cox <- coxph(Surv(time, event) ~ strata(a,w), data = dt, x = TRUE, y = TRUE)
#'   ate(e.cox, data = dt, treatment = "a", times = 1)
#' }
#'
#' ## check efficient estimator
#' set.seed(10)
#' dtRed <- dt[sample.int(n = .N, size = 1e3, replace = FALSE)]
#' setkeyv(dtRed, c("a","w","time"))
#' dtRed[,max(time), by = c("a","w")]
#' 
#' res2 <- ateRobust(data = dtRed, type = "survival",
#'             formula.event = Surv(time, event) ~ strata(a,w), ## ~ a + w,##
#'             formula.censor = Surv(time, event==0) ~ strata(a,w), ## ~ a + w,##
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             treatment = "a",
#'             efficient = TRUE,
#'             product.limit = TRUE)
#' res2
#' 
#' res3 <- ateRobust(data = dtRed, type = "survival",
#'             formula.event = Surv(time, event) ~ a + w,##
#'             formula.censor = Surv(time, event==0) ~ a + w,##
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             treatment = "a",
#'             efficient = TRUE,
#'             product.limit = TRUE)
#' res3


## * ateRobust (code)
#' @rdname ateRobust
#' @export
ateRobust <- function(data, times, treatment, cause, type,
                      formula.event, formula.censor, formula.treatment, 
                      fitter = "coxph", efficient = TRUE, product.limit = FALSE, iid = FALSE){

    
    ## ** normalize arguments
    if(!is.data.table(data)){
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }
    n.obs <- NROW(data)
    if(length(times)!=1){
        stop("Argument \'times\' must have length 1 \n")
    }
    if(treatment %in% names(data) == FALSE){
        stop("Variable \"",treatment,"\" not found in data \n")
    }
    if(!is.factor(data[[treatment]])){
        data[, c(treatment) := as.factor(.SD[[treatment]])]
    }
    level.treatment <- levels(data[[treatment]])
    if(length(level.treatment)!=2){
        stop("only implemented for binary treatment variables \n")
    }
    type <- match.arg(type, c("survival","competing.risks"))
    if(product.limit){
        predictor.cox <- "predictCoxPL"
    }else{
        predictor.cox <- "predictCox"        
    }
    
    ## times since it is an argument that will be used in the data table
    reserved.name <- c("times","time.tau","status.tau","weights","prob.event","prob.event0","prob.event1","prob.censoring","prob.treatment","treatment.bin",
                       "Lterm0","Lterm1")
    if(any(reserved.name %in% names(data))){
        txt <- reserved.name[reserved.name %in% names(data)]
        stop("Argument \'data\' should not contain column(s) named \"",paste0(txt, collapse = "\" \""),"\"\n")
    }

    ## ** prepare dataset
    ## convert to binary
    data[, c("treatment.bin") := as.numeric(as.factor(.SD[[treatment]]))-1]

    ## counterfactual
    data0 <- copy(data)
    data0[,c(treatment) := factor(level.treatment[1], level.treatment)]

    data1 <- copy(data)
    data1[,c(treatment) := factor(level.treatment[2], level.treatment)]

    ## ** outcome model: conditional expectation
    if(type == "survival"){
        model.event <- do.call(fitter, args = list(formula = formula.event, data = data, x = TRUE, y = TRUE))
        coxMF <- coxModelFrame(model.event)

        ## [:change outcome:]
        data[, c("prob.event") := do.call(predictor.cox, args = list(model.event, newdata = data, times = times))$survival[,1]]
        data[, c("prob.event0") := do.call(predictor.cox, args = list(model.event, newdata = data0, times = times))$survival[,1]]
        data[, c("prob.event1") := do.call(predictor.cox, args = list(model.event, newdata = data1, times = times))$survival[,1]]
        ## data[, c("prob.event") := 1-do.call(predictor.cox, args = list(model.event, newdata = data, times = times))$survival[,1]]
        ## data[, c("prob.event0") := 1-do.call(predictor.cox, args = list(model.event, newdata = data0, times = times))$survival[,1]]
        ## data[, c("prob.event1") := 1-do.call(predictor.cox, args = list(model.event, newdata = data1, times = times))$survival[,1]]

    }else if(type=="competing.risks"){
        model.event <- CSC(formula.event, data = data, fitter = fitter, cause = cause)
        coxMF <- coxModelFrame(model.event$models[[paste0("Cause ",model.event$theCause)]])

        ## [:change outcome:]
        data[, c("prob.event") := 1-predict(model.event, newdata = data, times = times, cause = cause, product.limit = product.limit)$absRisk[,1]]
        data[, c("prob.event0") := 1-predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit)$absRisk[,1]]
        data[, c("prob.event1") := 1-predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit)$absRisk[,1]]
        ## data[, c("prob.event") := predict(model.event, newdata = data, times = times, cause = cause, product.limit = product.limit)$absRisk[,1]]
        ## data[, c("prob.event0") := predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit)$absRisk[,1]]
        ## data[, c("prob.event1") := predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit)$absRisk[,1]]
    }

    ## truncate at times
    data[,c("time.tau") := pmin(coxMF$stop,times)]
    data[,c("status.tau") := (coxMF$status==1)*(coxMF$stop<=times) - (coxMF$status==0)*(coxMF$stop<times)]    
    ## -1 censored, 0 at risk, 1 event

    ## ** Propensity score model: weights
    model.treatment <- do.call(glm, args = list(formula = formula.treatment, data = data, family = stats::binomial(link = "logit")))
    data[, c("prob.treatment") := as.double(predict(model.treatment, newdata = data, type = "response"))]

    ## ** Censoring model: weights
    if(all(coxMF$status==1)){
        message("no censoring")
        model.censor <- NULL
        data[,c("prob.censoring") := 1]
    }else{
        ## fit model
        model.censor <- do.call(fitter, args = list(formula = formula.censor, data = data, x = TRUE, y = TRUE)) ## , no.opt = TRUE
        ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok
        data[,c("prob.censoring") := do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), diag = TRUE))$survival[,1]]
    }
    ## (status.tau>=0) : not censored
    data[,c("weights") := (.SD$status.tau>=0)/.SD$prob.censoring]
    sumW <- data[,sum(.SD$weights)]
    
    ## ** Compute parameter of interest
    IF <- list()

    ## *** Gformula
    IF$Gformula <- data[,cbind(
        prob.event0,
        prob.event1
    )] / n.obs
    ## 1-colSums(IF$Gformula)
    
    ## *** naive IPW
    ## (status.tau==0) : still at risk
    ## data[treatment.bin==0,.(mean(status.tau==0),mean(prob.event0))]
    ## data[treatment.bin==1,.(mean(status.tau==0),mean(prob.event1))]
    IF$IPWnaive <- data[,cbind(
        ## [:change outcome:]
        .SD$weights * ((.SD$status.tau == 0) * (1-.SD$treatment.bin) / (1-.SD$prob.treatment)),
        .SD$weights * ((.SD$status.tau == 0) * (.SD$treatment.bin) / (.SD$prob.treatment))
        ## .SD$weights * ((.SD$status.tau == 1) * (1-.SD$treatment.bin) / (1-.SD$prob.treatment)),
        ## .SD$weights * ((.SD$status.tau == 1) * (.SD$treatment.bin) / (.SD$prob.treatment))
    )]/ sumW

    ## *** naive AIPW
    IF$AIPWnaive <- data[,cbind(
        ## [:change outcome:]
        .SD$weights * (((.SD$status.tau == 0) - .SD$prob.event0) * (1-.SD$treatment.bin) / (1-.SD$prob.treatment) + .SD$prob.event0),
        .SD$weights * (((.SD$status.tau == 0) - .SD$prob.event1) * (.SD$treatment.bin) / (.SD$prob.treatment) + .SD$prob.event1)
        ## .SD$weights * (((.SD$status.tau == 1) - .SD$prob.event0) * (1-.SD$treatment.bin) / (1-.SD$prob.treatment) + .SD$prob.event0),
        ## .SD$weights * (((.SD$status.tau == 1) - .SD$prob.event1) * (.SD$treatment.bin) / (.SD$prob.treatment) + .SD$prob.event1)
    )] / sumW
    ## lapply(IF, colSums)

    ## *** correction for efficiency
    if(efficient){
        ## data is updated within the function .calcLterm
        .calcLterm(data = data, data0 = data0, data1 = data1,
                   n.obs = n.obs, times = times, type = type, cause = cause,
                   model.censor = model.censor, model.event = model.event,
                   predictor.cox = predictor.cox, product.limit = product.limit)

        IF$IPWefficient <- IF$IPWnaive + data[, cbind(
        (1-treatment.bin)/(1-prob.treatment) * prob.event * Lterm,
        (treatment.bin/prob.treatment) * prob.event * Lterm
        )] / sumW
        
        IF$AIPWefficient <- IF$IPWefficient + data[, cbind(
                                                  prob.event0 * (1-(1-treatment.bin)/(1-prob.treatment)),
                                                  prob.event1 * (1-treatment.bin/prob.treatment)
                                              )] / sumW

        ## IF.tempo$correctionAIPW <- cbind(data[,(1-(1-.SD$treatment.bin)/(1-.SD$prob.treatment))*.SD$prob.event0*(1-.SD$weights)],
        ## data[,(1-(.SD$treatment.bin)/(.SD$prob.treatment))*.SD$prob.event1*(1-.SD$weights)])
        ## range((IF$AIPWnaive + (IF.tempo$correctionIPW + IF.tempo$correctionAIPW) / sumW) - IF$AIPWefficient)
    }
  
    ## ** export
    out <- list()
    
    ## value
    n.estimator <- length(IF)
    name.surv <- paste0("surv.",level.treatment)
    out$ate.value <- matrix(NA, nrow = 3, ncol = n.estimator,
                            dimnames = list(c(name.surv,"ate.diff"),
                                              names(IF)))

    for(iL in 1:n.estimator){ ## iL <- 1
        out$ate.value[name.surv,iL] <- colSums(IF[[iL]])
        IF[[iL]] <- rowCenter_cpp(IF[[iL]], center = out$ate.value[paste0("surv.",level.treatment),iL]/n.obs)
    }
    out$ate.value["ate.diff",] <- out$ate.value[name.surv[2],] - out$ate.value[name.surv[1],]
    ##    out$ate.value["ate.ratio",] <- out$ate.value[name.surv[2],] / out$ate.value[name.surv[1],]


    ## standard error
    out$ate.se <- do.call(cbind,lapply(IF, function(iIF){ ## iIF <- IF[[1]]
        sqrt(c(colSums(iIF^2), sum((iIF[,2]-iIF[,1])^2)))
    }))
    rownames(out$ate.se) <- rownames(out$ate.value)

    ## center influence function
    if(iid){
        out$iid <- do.call(cbind, IF)
        colnames(out$iid) <- unlist(lapply(paste0(names(IF),"_"),paste0,paste0("surv.",level.treatment)))
    }
    
    return(out)
    

}

## * .calcLterm
.calcLterm <- function(data, data0, data1,
                       n.obs, times, type,
                       model.censor, model.event,
                       predictor.cox, product.limit, cause){

    status <- stop <- NULL ## [:CRANtest:] data.table
    
    eXb.censor <- exp(coxLP(model.censor, data = NULL, center = FALSE))
    MF.censor <- coxModelFrame(model.censor)

    jump.times <- sort(MF.censor[status==1&stop<=times,stop]) ## only look at jump times before the prediction time    
    resBaseline.censor <- predictCox(model.censor,
                                     times = jump.times,
                                     type = "hazard",
                                     keep.infoVar = TRUE)
    is.strata <- resBaseline.censor$infoVar$is.strata

    calcLterm <- function(status, time, jump.time, hazard0, eXb, surv.event, surv.censor){
        ## jump times are before the prediction time
        index.beforeJump <- which(jump.time<=time)
        if(length(index.beforeJump)==0){
            return(0)
        }else{
            dN <- (jump.time[index.beforeJump]==time)*status
            dLambda <- hazard0[index.beforeJump]*eXb
            return(sum((dN-dLambda)/(surv.event[index.beforeJump]*surv.censor[index.beforeJump])))
        }
    }

    data[,c("Lterm") := as.numeric(NA)]
    
    if(!is.strata){
        hazard0 <- resBaseline.censor$hazard
        all.times <- resBaseline.censor$times

        ## survival = P[C>t] = P[Delta(t)==1] - ok
        pred.censor <- do.call(predictor.cox,
                               args = list(model.censor, newdata = data, times = all.times-(1e-10), type = "survival"))$survival

        if(type=="survival"){
            ## [:change outcome:]
            pred.event <- do.call(predictor.cox,
                                  args = list(model.event, newdata = data, times = all.times, type = "survival"))$survival
            ## pred.event <- 1-do.call(predictor.cox,
            ## args = list(model.event, newdata = data, times = all.times, type = "survival"))$survival
            pred.eventC <- colScale_cpp(pred.event, scale = data[["prob.event"]])
            
        }else if(type == "competing.risks"){
            pred.event <- predict(model.event, newdata = data, times = all.times,
                                  product.limit = product.limit, cause = cause)$absRisk
            pred.surv <- do.call(predictor.cox,
                                 args = list(model.event$models[[paste0("Cause ",cause)]], newdata = data, times = all.times, type = "survival"))$survival

            pred.eventC <- -colCenter_cpp(pred.event, center = data[["prob.event"]])/pred.surv

        }
        
        data[,"Lterm" := sapply(1:n.obs, function(iObs){ ## iObs <- 1

                calcLterm(status = MF.censor$status[iObs],
                          time = MF.censor$stop[iObs],
                          jump.time = all.times,
                          hazard0 = hazard0,
                          eXb = eXb.censor[iObs],
                          surv.event = pred.eventC[iObs,],
                          surv.censor = pred.censor[iObs,])
            
        })]

    }else{

        strata.levels.censor <- resBaseline.censor$infoVar$strata.levels
        n.strata.censor <- length(strata.levels.censor)
            
        for(iStrata in 1:n.strata.censor){ ## iStrata <- 2
            iIndex.strata <- which(MF.censor$strata == strata.levels.censor[iStrata])
            iN.strata <- length(iIndex.strata)

            iIndex.baseline <- intersect(which(resBaseline.censor$strata==strata.levels.censor[iStrata]),
                                         which(resBaseline.censor$times<=times))
            iHazard0 <- resBaseline.censor$hazard[iIndex.baseline]
            iAll.times <- resBaseline.censor$times[iIndex.baseline]

            ## survival = P[C>t] = P[Delta(t)==1] - ok
            iPred.censor <- do.call(predictor.cox,
                                    args = list(model.censor, newdata = data[iIndex.strata], times = iAll.times-(1e-10), type = "survival"))$survival

            if(type == "survival"){
                ## [:change outcome:]
                iPred.event <- do.call(predictor.cox,
                                       args = list(model.event, newdata = data[iIndex.strata], times = iAll.times, type = "survival"))$survival
                ## iPred.event <- 1-do.call(predictor.cox,
                ## args = list(model.event, newdata = data[iIndex.strata], times = iAll.times, type = "survival"))$survival
                iPred.eventC <- colScale_cpp(iPred.event, scale = data[iIndex.strata,.SD$prob.event])

            }else if(type == "competing.risks"){
                iPred.event <- predict(model.event, newdata = data[iIndex.strata], times = iAll.times,
                                       product.limit = product.limit, cause = cause)$absRisk
                iPred.surv <- do.call(predictor.cox,
                                      args = list(model.event$models[[paste0("Cause ",cause)]], newdata = data[iIndex.strata], times = iAll.times, type = "survival"))$survival

                iPred.eventC <- -colCenter_cpp(iPred.event, center = data[iIndex.strata,.SD$prob.event])/iPred.surv

            }

            data[iIndex.strata, Lterm := sapply(1:iN.strata, function(iObs){ ## iObs <- 10

                calcLterm(status = MF.censor$status[iIndex.strata[iObs]],
                          time = MF.censor$stop[iIndex.strata[iObs]],
                          jump.time = iAll.times,
                          hazard0 = iHazard0,
                          eXb = eXb.censor[iIndex.strata[iObs]],
                          surv.event = iPred.eventC[iObs,],
                          surv.censor = iPred.censor[iObs,])
            
            })]

        }
        
            
    }

    return(NULL)
}

######################################################################
### ateRobust.R ends here
