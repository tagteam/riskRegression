### ateTMLE.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2018 (17:47) 
## Version: 
## Last-Updated: jun 28 2018 (16:05) 
##           By: Brice Ozenne
##     Update #: 165
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
#' @param formula.time [formula] The Cox model for the event of interest. Typically \code{Surv(time,event)~treatment}.
#' @param formula.censor [formula] The Cox model for the censoring. Typically \code{Surv(time,event==0)~treatment}.
#' @param formula.treatment [formula] The logistic model for the treatment. Typically \code{treatment~1}.
#' @param times [numeric] Time point at which to evaluate average treatment effects.
#' @param treatment [character] Name of the treatment variable. Must be binary.
#' @param fitter.cox [character] Routine to fit the Cox regression models.
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
#' dt[,.(surv = mean(eventtime>=1)),by = c("a")]
#'
#' ## Estimated value
#' res <- ateRobust(data = dt,
#'             formula.time = Surv(time, event) ~ strata(a,w),
#'             formula.censor = Surv(time, event==0) ~ strata(a,w),
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             treatment = "a",
#'             efficient = FALSE,
#'             product.limit = FALSE)
#' res
#'
#' ## check efficient estimator
#' set.seed(10)
#' dtRed <- dt[sample.int(n = .N, size = 1e3, replace = FALSE)]
#' setkeyv(dtRed, c("a","w","time"))
#' dtRed[,max(time), by = c("a","w")]
#' 
#' res <- ateRobust(data = dtRed,
#'             formula.time = Surv(time, event) ~ strata(a,w), ## ~ a + w,##
#'             formula.censor = Surv(time, event==0) ~ strata(a,w), ## ~ a + w,##
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             treatment = "a",
#'             efficient = TRUE,
#'             product.limit = TRUE)
#' res
#' 


## * ateRobust (code)
#' @rdname ateRobust
#' @export
ateRobust <- function(data, formula.time, formula.censor, formula.treatment, times, treatment,
                    fitter.cox = "coxph", efficient = TRUE, product.limit = FALSE,
                    iid = FALSE){

    
    ## ** normalize arguments
    if(!is.data.table(data)){
        data <- data.table::as.data.table(data)
    }else{
        data <- data.table::copy(data)
    }
    n.obs <- NROW(data)
    if(!is.factor(data[[treatment]])){
        data[, c(treatment) := as.factor(.SD[[treatment]])]
    }
    level.treatment <- levels(data[[treatment]])
    if(length(level.treatment)!=2){
        stop("only implemented for binary treatment variables \n")
    }
    if(product.limit){
        predictor.cox <- "predictCoxPL"
    }else{
        predictor.cox <- "predictCox"        
    }
    
    ## times since it is an argument that will be used in the data table
    reserved.name <- c("times","time.tau","status.tau","weights","prob.survival","prob.survival0","prob.survival1","prob.survivalCensoring","prob.treatment","treatment.bin")
    if(any(reserved.name %in% names(data))){
        txt <- reserved.name[reserved.name %in% names(data)]
        stop("Argument \'data\' should not contain column(s) named \"",paste0(txt, collapse = "\" \""),"\"\n")
    }

    infoVar <- SurvResponseVar(formula.time)
    var.time <- infoVar$time
    var.status <- infoVar$status

    ## ** prepare
    ## *** dataset
    ## convert to binary
    data[, c("treatment.bin") := as.numeric(as.factor(.SD[[treatment]]))-1]

    ## truncate at times
    data[,c("time.tau") := pmin(.SD[[var.time]],times)]
    data[,c("status.tau") := (.SD[[var.status]]==1)*(.SD[[var.time]]<=times) - (.SD[[var.status]]==0)*(.SD[[var.time]]<times)]
    ## -1 censored, 0 at risk, 1 event
    ## table(data$status.tau)
    
    ## *** counterfactual
    data0 <- copy(data)
    data0[,c(treatment) := factor(level.treatment[1],level.treatment)]

    data1 <- copy(data)
    data1[,c(treatment) := factor(level.treatment[2],level.treatment)]

    ## *** Outcome model: conditional expectation
    model.event <- do.call(fitter.cox, args = list(formula = formula.time, data = data, x = TRUE, y = TRUE))

    data[, c("prob.survival") := as.data.table(do.call(predictor.cox, args = list(model.event, newdata = data, times = times)))$survival]
    data[, c("prob.survival0") := as.data.table(do.call(predictor.cox, args = list(model.event, newdata = data0, times = times)))$survival]
    data[, c("prob.survival1") := as.data.table(do.call(predictor.cox, args = list(model.event, newdata = data1, times = times)))$survival]

    ## *** Propensity score model: weights
    model.treatment <- do.call(glm, args = list(formula = formula.treatment, data = data, family = stats::binomial(link = "logit")))
    data[, c("prob.treatment") := as.double(predict(model.treatment, newdata = data, type = "response"))]

    ## *** Censoring model: weights
    if((model.event$n-model.event$nevent)==0){
        message("no censoring")
        model.censor <- NULL
        data[,c("prob.survivalCensoring") := 1]
    }else{
        ## fit model
        model.censor <- do.call(fitter.cox, args = list(formula = formula.censor, data = data, x = TRUE, y = TRUE)) ## , no.opt = TRUE
        ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok
        data[,c("prob.survivalCensoring") := do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), diag = TRUE))$survival[,1]]
    }
    ## (status.tau>=0) : not censored
    data[,c("weights") := (.SD$status.tau>=0)/.SD$prob.survivalCensoring]
    sumW <- data[,sum(.SD$weights)]
    
    ## ** Compute parameter of interest
    IF <- list()

    ## *** Gformula
    IF$Gformula <- cbind(data[,.SD$prob.survival0/.N],                                               
                         data[,.SD$prob.survival1/.N]
                         )
        
    ## *** naive IPW
    ## (status.tau==0) : still at risk
    ## data[treatment.bin==0,.(mean(status.tau==0),mean(prob.survival0))]
    ## data[treatment.bin==1,.(mean(status.tau==0),mean(prob.survival1))]
    
    IF$IPWnaive <- cbind(data[,.SD$weights * ((.SD$status.tau==0)*(1-.SD$treatment.bin)/(1-.SD$prob.treatment))] / sumW,
                         data[,.SD$weights * ((.SD$status.tau==0)*(.SD$treatment.bin)/(.SD$prob.treatment))] / sumW
                         )
    ## c(data[status.tau>=0,weighted.mean((status.tau==0)*(1-treatment.bin)/(1-prob.treatment), w = weights)],
    ##   data[status.tau>=0,weighted.mean((status.tau==0)*(treatment.bin)/(prob.treatment), w = weights)]
    ##   )
    
    ## *** naive AIPW
    IF$AIPWnaive <- cbind(data[,.SD$weights * (((.SD$status.tau==0)-.SD$prob.survival0)*(1-.SD$treatment.bin)/(1-.SD$prob.treatment)+.SD$prob.survival0)] / sumW,
                          data[,.SD$weights * (((.SD$status.tau==0)-.SD$prob.survival1)*(.SD$treatment.bin)/(.SD$prob.treatment)+.SD$prob.survival1)] / sumW
                          )
    ## c(data[status.tau>=0,weighted.mean(((status.tau==0)-prob.survival0)*(1-treatment.bin)/(1-prob.treatment)+prob.survival0, w = weights)],
      ## data[status.tau>=0,weighted.mean(((status.tau==0)-prob.survival1)*(treatment.bin)/(prob.treatment)+prob.survival1, w = weights)]
      ## )
    

    ## *** correction for efficiency
    if(efficient){
        IF.tempo <- .correctionIPW(data = data, n.obs = n.obs, times = times,
                                   model.censor = model.censor, model.event = model.event,
                                   fitter.cox = fitter.cox, predictor.cox = predictor.cox)

        IF$IPWefficient <- IF$IPWnaive + (IF.tempo$correctionIPW) / sumW
        IF$AIPWefficient <- IF$AIPWnaive + (IF.tempo$correctionIPW + IF.tempo$correctionAIPW) / sumW
    }
  

        ## out$AIPWefficient <- out$AIPWnaive + IPWcorrrection + AIPWcorrrection
        

    ## ** export
    out <- list()
    
    ## value
    out$ate.value <- do.call(cbind,lapply(IF, colSums))
    rownames(out$ate.value) <- paste0("surv.",level.treatment)
    out$ate.value <- rbind(out$ate.value,
                           "ate.diff" = out$ate.value[2,]-out$ate.value[1,],
                           "ate.ratio" = out$ate.value[2,]/out$ate.value[1,])

    ## standard error
    out$ate.se <- do.call(cbind,lapply(IF, function(iIF){ ## iIF <- IF[[1]]
        sqrt(c(colSums(iIF^2), sum((iIF[,2]-iIF[,1])^2), NA))
    }))
    rownames(out$ate.se) <- rownames(out$ate.value)

    if(iid){
        out$iid <- do.call(cbind, IF)
        colnames(out$iid) <- unlist(lapply(paste0(names(IF),"_"),paste0,paste0("surv.",level.treatment)))
    }
    
    return(out)
    

}

## * .correctionIPW
.correctionIPW <- function(data, n.obs, times,
                           model.censor, model.event,
                           fitter.cox, predictor.cox){

    status <- stop <- NULL ## [:CRANtest:] data.table
    
    eXb.censor <- exp(coxLP(model.censor, data = NULL, center = FALSE))
    MF.censor <- coxModelFrame(model.censor)

    jump.times <- sort(MF.censor[status==1&stop<=times,stop]) ## only look at jump times before the prediction time    
    resBaseline.censor <- predictCox(model.censor,
                                     times = jump.times,
                                     type = "hazard", keep.infoVar = TRUE)
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

    if(!is.strata){
        hazard0 <- resBaseline.censor$hazard
        all.times <- resBaseline.censor$times
        
        pred.event <- do.call(predictor.cox,
                              args = list(model.event, newdata = data, times = all.times, type = "survival"))$survival
        ## survival = P[C>t] = P[Delta(t)==1] - ok
        pred.censor <- do.call(predictor.cox,
                               args = list(model.censor, newdata = data, times = all.times-(1e-10), type = "survival"))$survival

        data[, c("Lterm") := sapply(1:n.obs, function(iObs){ ## iObs <- 1

            calcLterm(status = MF.censor$status[iObs],
                      time = MF.censor$stop[iObs],
                      jump.time = all.times,
                      hazard0 = hazard0,
                      eXb = eXb.censor[iObs],
                      surv.event = pred.event[iObs,],
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

            iPred.event <- do.call(predictor.cox,
                                   args = list(model.event, newdata = data[iIndex.strata], times = iAll.times, type = "survival"))$survival
            ## survival = P[C>t] = P[Delta(t)==1] - ok
            iPred.censor <- do.call(predictor.cox,
                                    args = list(model.censor, newdata = data[iIndex.strata], times = iAll.times-(1e-10), type = "survival"))$survival


            data[iIndex.strata, c("Lterm") := sapply(1:iN.strata, function(iObs){ ## iObs <- 10

                calcLterm(status = MF.censor$status[iIndex.strata[iObs]],
                          time = MF.censor$stop[iIndex.strata[iObs]],
                          jump.time = iAll.times,
                          hazard0 = iHazard0,
                          eXb = eXb.censor[iIndex.strata[iObs]],
                          surv.event = iPred.event[iObs,],
                          surv.censor = iPred.censor[iObs,])
            
            })]
        }
        
            
    }

    ## ** export
    out <- list()
    
    out$correctionIPW <- cbind(data[,(1-.SD$treatment.bin)/(1-.SD$prob.treatment)*.SD$prob.survival0*.SD$Lterm],
                               data[,(.SD$treatment.bin)/(.SD$prob.treatment)*.SD$prob.survival1*.SD$Lterm])
    ## colSums(out$correctionIPW)
    out$correctionAIPW <- cbind(data[,(1-(1-.SD$treatment.bin)/(1-.SD$prob.treatment))*.SD$prob.survival0*(1-.SD$weights)],
                                data[,(1-(.SD$treatment.bin)/(.SD$prob.treatment))*.SD$prob.survival1*(1-.SD$weights)])

    return(out)
}

######################################################################
### ateRobust.R ends here
