## ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2018 (17:47) 
## Version: 
## Last-Updated: Sep  5 2018 (11:04) 
##           By: Thomas Alexander Gerds
##     Update #: 1025
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
#' @param formula.event [formula] Cox model for the event of interest (outcome model).
#' Typically \code{Surv(time,event)~treatment}.
#' @param formula.censor [formula] Cox model for the censoring (censoring model).
#' Typically \code{Surv(time,event==0)~treatment}.
#' @param formula.treatment [formula] Logistic regression for the treatment (propensity score model).
#' Typically \code{treatment~1}.
#' @param times [numeric] Time point at which to evaluate average treatment effects.
#' @param fitter [character] Routine to fit the Cox regression models.
#' If \code{coxph} use \code{survival::coxph} else use \code{rms::cph}.
#' @param product.limit [logical] If \code{TRUE} the survival is computed using the product limit estimator.
#' Otherwise the exponential approximation is used (i.e. exp(-cumulative hazard)).
#' @param cause [numeric/character] The cause of interest. Defaults to the first cause.
#' @param type [character] When set to \code{"survival"} uses a cox model for modeling the survival,
#' otherwise when set to \code{"competing.risks"} uses a Cause Specific Cox model for modeling the absolute risk of the event.
#' @param nuisance.iid [logical] If \code{TRUE} take into accound the uncertainty related to
#' the estimation of outcome model and the propensity score model.
#' The uncertainty related to the estimation of the censoring model is never accounted for.
#' @param efficient [logical] Should the efficient IPW and AIPW estimator be used?
#' @param na.rm [logical] If \code{TRUE} ignore observations whose influence function is NA.
#'
#' @details Argument \code{nuisance.iid}: the asymptotic distribution of the ATE
#' when we estimate the nuisance parameters
#' (parameters from the outcome/propensity score/censoring model)
#' equals the one if we would know the nuisance parameters. Therefore in large sample size,
#' the value of the argument \code{nuisance.iid} should not matter. Setting it to \code{FALSE} will save some computation time.
#' 

## * ateRobust (example)
#' @rdname ateRobust
#' @examples
#'
#' library(survival)
#' library(lava)
#' library(data.table)
#'
#' #### Survival ####
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
#' dt <- as.data.table(lava::sim(mSimSurv, n = 1e3, p = c(alpha = alpha, beta = beta)))
#' setkeyv(dt, c("a","w"))
#'
#' ## True value
#' psi.TRUE <- c("risk.0" = (1-exp(-1*exp(alpha*0+beta*0)))*0.5 + (1-exp(-1*exp(alpha*0+beta*1)))*0.5,
#'               "risk.1" = (1-exp(-1*exp(alpha*1+beta*0)))*0.5 + (1-exp(-1*exp(alpha*1+beta*1)))*0.5)
#' psi.TRUE
#' 
#' ## Approximate true value
#' dt[,.(risk = mean(eventtime<1)),by = c("a")]
#'
#' ## Estimated using stratified Cox model
#' res1 <- ateRobust(data = dt, type = "survival",
#'             formula.event = Surv(time, event) ~ strata(a,w),
#'             formula.censor = Surv(time, event==0) ~ strata(a,w),
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             product.limit = FALSE)
#' print(res1, efficient = TRUE, nuisance.iid = TRUE)
#' print(res1, efficient = TRUE, nuisance.iid = FALSE)
#' print(res1, efficient = FALSE, nuisance.iid = TRUE)
#' print(res1, efficient = FALSE, nuisance.iid = FALSE)
#'
#' if(FALSE){
#'   e.cox <- coxph(Surv(time, event) ~ strata(a,w), data = dt, x = TRUE, y = TRUE)
#'   ate(e.cox, data = dt, treatment = "a", times = 1)
#' rbind(estimate = res1$ate.value[,"Gformula2"],
#'       se = res1$ate.se[,"Gformula2"])
#' }
#'
#' ## Estimate using Cox model
#' res2 <- ateRobust(data = dt, type = "survival",
#'             formula.event = Surv(time, event) ~ a + w,
#'             formula.censor = Surv(time, event==0) ~ a + w,
#'             formula.treatment = a ~ w,
#'             times = 1,
#'             product.limit = FALSE,
#'             nuisance.iid = TRUE)
#' 
#' print(res2, efficient = TRUE, nuisance.iid = TRUE)
#' print(res2, efficient = TRUE, nuisance.iid = FALSE)
#' print(res2, efficient = FALSE, nuisance.iid = TRUE)
#' print(res2, efficient = FALSE, nuisance.iid = FALSE)
#'
#' #### Competing risks ####
#' set.seed(10)
#' n <- 1e3
#' 
#' ## simulate data
#' alphaE.X <- 2
#' alphaCR.X <- 1
#' alphaE.Y <- 3
#' alphaCR.Y <- 2
#'
#' set.seed(10)
#' df <- rbind(data.frame(time1 = rexp(n, rate = alphaE.X), time2 = rexp(n, rate = alphaCR.X), group = "1"),
#'             data.frame(time1 = rexp(n, rate = alphaE.Y), time2 = rexp(n, rate = alphaCR.Y), group = "2"))
#' df$time <- pmin(df$time1,df$time2) ## first event
#' df$event <- (df$time2<df$time1)+1 ## type of event
#' df$eventC <- df$event
#' df$eventC[rbinom(n, size = 1, prob = 0.2)==1] <- 0
#' 
#' ## true value
#' tau <- 1
#' c(CIF.X = alphaE.X/(alphaE.X+alphaCR.X)*(1-exp(-(alphaE.X+alphaCR.X)*(tau))),
#'   CIF.Y = alphaE.Y/(alphaE.Y+alphaCR.Y)*(1-exp(-(alphaE.Y+alphaCR.Y)*(tau))))
#'
#' ## estimating using a CSC (no censoring)
#' tapply(df$time,df$group,max)
#' 
#' resCR <- ateRobust(data = df, type = "competing.risks",
#'             formula.event = Hist(time, event) ~ group, ## strata(group),
#'             formula.censor = Surv(time, event==0) ~ group,## strata(group),
#'             formula.treatment = group ~ 1,
#'             times = tau,
#'             nuisance.iid = FALSE,
#'             product.limit = FALSE,
#'             cause = 1)
#' resCR
#'
#'
#' ## estimating using a CSC (censoring)
#' resCRC <- ateRobust(data = df, type = "competing.risks",
#'             formula.event = Hist(time, eventC) ~  strata(group), ## group, ##
#'             formula.censor = Surv(time, eventC==0) ~  strata(group), ## group,##
#'             formula.treatment = group ~ 1,
#'             times = tau,
#'             nuisance.iid = FALSE,
#'             product.limit = FALSE,
#'             cause = 1)
#' print(resCRC, efficient = TRUE)
#' print(resCRC, efficient = FALSE)
#'
#'



## * ateRobust (code)
#' @rdname ateRobust
#' @export
ateRobust <- function(data, times, cause, type,
                      formula.event, formula.censor, formula.treatment, 
                      fitter = "coxph", product.limit = FALSE, efficient = TRUE,
                      nuisance.iid = TRUE, na.rm = FALSE){

    modeSurvival <- "risk"
    
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
    treatment <- all.vars(formula.treatment)[1]
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
    reserved.name <- c("times",
                       "time.tau",
                       "status.tau",
                       "censoring.tau",
                       "treatment.bin",
                       "prob.event",
                       "prob.event0",
                       "prob.event1",
                       "prob.treatment",
                       "weights",
                       "prob.censoring",
                       "prob.indiv.censoring",
                       "Lterm")
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

    ## random variables stopped at times
    if(type == "survival"){
        ## coxMF contains three columns: start, stop, status (0 censored, 1 event), + strata?
        data[,c("time.tau") := pmin(coxMF$stop,times)]
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$status==1)]
        ## 0 survival or censored, 1 event
        data[,c("Ncensoring.tau") := as.numeric(((coxMF$stop>=times) + (coxMF$status!=0)) > 0)]
        ## 0 censored, 1 survival or event
    }else if(type=="competing.risks"){
        ## coxMF contains three columns: stop, status (0 censored, 1 any event), event (1 event of interest, 2 competing event, 3 censoring)
        data[,c("time.tau") := pmin(coxMF$stop,times)]
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$event==cause)]
        ## 0 survival competing event or censored, 1 event of interest
        data[,c("Ncensoring.tau") := as.numeric(((coxMF$stop>=times) + (coxMF$status!=0)) > 0)]
        ## 0 censored, 1 survival or event (of interest or comepting)
    }
    
    ## ** outcome model: conditional expectation
    if(type == "survival"){
        model.event <- do.call(fitter, args = list(formula = formula.event, data = data, x = TRUE, y = TRUE))
        coxMF <- coxModelFrame(model.event)

        prediction.event <- do.call(predictor.cox, args = list(model.event, newdata = data, times = times, iid = nuisance.iid))
        prediction.event0 <- do.call(predictor.cox, args = list(model.event, newdata = data0, times = times, iid = nuisance.iid))
        prediction.event1 <- do.call(predictor.cox, args = list(model.event, newdata = data1, times = times, iid = nuisance.iid))

        if(modeSurvival == "survival"){
            data[, c("prob.event") := prediction.event$survival[,1]]
            data[, c("prob.event0") := prediction.event0$survival[,1]]
            data[, c("prob.event1") := prediction.event1$survival[,1]]

            if(nuisance.iid){
                iid.event0 <- prediction.event0$survival.iid[,1,]
                iid.event1 <- prediction.event1$survival.iid[,1,]
            }
        }else{
            data[, c("prob.event") := 1 - prediction.event$survival[,1]]
            data[, c("prob.event0") := 1 - prediction.event0$survival[,1]]
            data[, c("prob.event1") := 1 - prediction.event1$survival[,1]]

            if(nuisance.iid){
                iid.event0 <- - prediction.event0$survival.iid[,1,]
                iid.event1 <- - prediction.event1$survival.iid[,1,]
            }
        }

    }else if(type=="competing.risks"){
        model.event <- CSC(formula.event, data = data, fitter = fitter, cause = cause, surv.type = "hazard")
        coxMF <- unclass(model.event$response)

        prediction.event <- predict(model.event, newdata = data, times = times, cause = cause, product.limit = product.limit, iid = FALSE)
        prediction.event0 <- predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit, iid = nuisance.iid)
        prediction.event1 <- predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit, iid = nuisance.iid)

        data[, c("prob.event") := prediction.event$absRisk[,1]]
        data[, c("prob.event0") := prediction.event0$absRisk[,1]]
        data[, c("prob.event1") := prediction.event1$absRisk[,1]]

        if(nuisance.iid){
            iid.event0 <- prediction.event0$absRisk.iid[,1,]
            iid.event1 <- prediction.event1$absRisk.iid[,1,]
        }

    }
    ## ** Censoring model: weights
    ## fit model
    if(n.censor==0){
        
        data[,c("prob.censoring") := 1]
        data[,c("prob.indiv.censoring") := 1]
        data[,c("weights") := 1]
        
        if(nuisance.iid){
            iidIPW.censoring0 <- rep(0,n.obs)
            iidIPW.censoring1 <- rep(0,n.obs)
        }
        
    }else{
        
        ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok

        ## stopped at tau
        predTau.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = times, type = "survival"))
        ## predTau.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = times, type = "survival", iid = nuisance.iid))
        data[,c("prob.censoring") := predTau.censor$survival[,1]]
        
        ## at each time
        predIndiv.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), type = "survival", diag = TRUE))
        ## predIndiv.censor <- do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), type = "survival", diag = TRUE, iid = nuisance.iid))
        data[,c("prob.indiv.censoring") := predIndiv.censor$survival[,1]]
        
        ## store
        data[,c("weights") := as.numeric(NA)]
        data[coxMF$stop<=times, c("weights") := .SD$Ncensoring.tau / .SD$prob.indiv.censoring]
        data[coxMF$stop>times, c("weights") := .SD$Ncensoring.tau / .SD$prob.censoring]

        if(nuisance.iid){
            ## iid.tempo <- matrix(NA, nrow = n.obs, ncol = n.obs)
            ## iid.tempo[coxMF$stop<=times,] <- predIndiv.censor$survival.iid[coxMF$stop<=times,1,]
            ## iid.tempo[coxMF$stop>times,] <- predTau.censor$survival.iid[coxMF$stop>times,1,]

            ## weight0 <- data[, .SD$weights^2 * .SD$status.tau * (1-.SD$treatment.bin) / (1-.SD$prob.treatment)]
            ## weight1 <- data[, .SD$weights^2 * .SD$status.tau * (.SD$treatment.bin) / (.SD$prob.treatment)]
            
            ## iidIPW.censoring0 <- -apply(iid.tempo,2,function(iCol){mean(weight0*iCol)})
            ## iidIPW.censoring1 <- -apply(iid.tempo,2,function(iCol){mean(weight1*iCol)})

            iidIPW.censoring0 <- rep(0,n.obs)
            iidIPW.censoring1 <- rep(0,n.obs)
        }

    }else if(type=="competing.risks"){
        data[,c("time.tau") := pmin(coxMF[,"time"],times)]
        data[,c("status.tau") := (coxMF[,"event"]==1)*(coxMF[,"time"]<=times) - (coxMF[,"status"]==0)*(coxMF[,"time"]<times)]
        ## -1 censored, 0 survival or death (i.e. no event), 1 event
    }
    n.censor <- sum(data$status.tau==-1)
    ## table(data$status.tau)

    ## ** Propensity score model: weights
    model.treatment <- do.call(glm, args = list(formula = formula.treatment, data = data, family = stats::binomial(link = "logit")))
    
    if(nuisance.iid){
        prediction.treatment <- .predictGLM(model.treatment, newdata = data)
        iid.treatment <- attr(prediction.treatment, "iid")
        attr(prediction.treatment, "iid") <- NULL
    }else{
        prediction.treatment <- cbind(predict(model.treatment, newdata = data, type = "response", se = FALSE))
    }

    ## ** Censoring model: weights
    ## fit model
    model.censor <- suppressWarnings(do.call(fitter, args = list(formula = formula.censor, data = data, x = TRUE, y = TRUE))) ## , no.opt = TRUE
    if(model.censor$nevent==0){
        data[,c("prob.censoring") := 1]
    }else{
        ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok
        data[,c("prob.censoring") := do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), diag = TRUE))$survival[,1]]
    }
    ## (status.tau>=0) : not censored
    data[,c("weights") := (.SD$status.tau>=0)/.SD$prob.censoring]
    sumW <- data[,sum(.SD$weights)]

    ## ** correction for efficiency
    ## data is updated within the function .calcLterm
    if(efficient){
        if(n.censor==0){
            data[,c("Lterm") := 0]
        }else{
            .calcLterm(data = data, data0 = data0, data1 = data1, coxMF = coxMF,
                       n.obs = n.obs, times = times, type = type, cause = cause,
                       model.censor = model.censor, model.event = model.event,
                       predictor.cox = predictor.cox, product.limit = product.limit,
                       modeSurvival = modeSurvival)
        }
    }
    
    ## ** Compute parameter of interest
    IF <- list()

    ## *** Gformula
    ## eq:IF-Gformula (blue term)
    IF$Gformula <- data[,cbind(
        .SD$prob.event0,
        .SD$prob.event1
    )] / n.obs

    if(nuisance.iid){
        ## additional term: eq:IF-Gformula-full (green term)
        IF$Gformula2 <- IF$Gformula + cbind(colMeans(iid.event0),
                                            colMeans(iid.event1))
    }
    ## *** IPW
    ## eq:IF-IPWc (blue term)
    IF$IPWnaive <- data[,cbind(
        .SD$weights * .SD$status.tau * (1-.SD$treatment.bin) / (1-.SD$prob.treatment),
        .SD$weights * .SD$status.tau * (.SD$treatment.bin) / (.SD$prob.treatment)
    )]/ n.obs
    ## colSums(IF$IPWnaive)
    ## colSums(IF$Gformula)
    
    if(efficient){
        ## additional term: eq:IF-IPWeff (red term)
        IF$IPWefficient <- IF$IPWnaive + data[,cbind(
                                               .SD$Lterm * (1-.SD$treatment.bin) / (1-.SD$prob.treatment), ## .SD$prob.event
                                               .SD$Lterm * (.SD$treatment.bin) / (.SD$prob.treatment) ## .SD$prob.event
                                         )]/ sumW ## n.obs ##
    }
    if(nuisance.iid){
        ## additional term: eq:IF-IPWfull (green term)
        weight.tempo0 <- data[,  (1-.SD$treatment.bin) * (.SD$status.tau==1) / (1-.SD$prob.treatment)^2]
        weight.tempo1 <- data[, .SD$treatment.bin * (.SD$status.tau==1) / (.SD$prob.treatment)^2]
        IPWadd <- - t(apply(iid.treatment,2,function(iCol){
            c(-mean(weight.tempo0*iCol),mean(weight.tempo1*iCol))
        }))
        ## colSums(IPWadd * data$weights * n.obs / sumW)
        ## colSums(IPWadd)

        ## IPW
        IF$IPWnaive2 <- IF$IPWnaive + IPWadd * data$weights * n.obs / sumW
        if(efficient){
            IF$IPWefficient2 <- IF$IPWefficient + IPWadd * data$weights * n.obs / sumW
        }
    }

    ## *** AIPW
    ## additional term:
    AIPWadd <- data[,cbind(
        .SD$weights * .SD$prob.event0 * (1-(1-.SD$treatment.bin)/(1-.SD$prob.treatment)),
        .SD$weights * .SD$prob.event1 * (1-.SD$treatment.bin/(.SD$prob.treatment))
    )]/ sumW ## n.obs ##
        
    IF$AIPWnaive <- IF$IPWnaive + AIPWadd

    if(efficient){
        IF$AIPWefficient <- IF$IPWefficient + AIPWadd
    }
    
    if(nuisance.iid){
        ## additional term: eq:IF-AIPWfull (green terms)
        weightY.tempo0 <- data[, 1 - (1-.SD$treatment.bin) / (1-.SD$prob.treatment)]
        weightY.tempo1 <- data[, 1 - .SD$treatment.bin / .SD$prob.treatment]

        AIPWaddY <- cbind(
            apply(iid.event0,2,function(iCol){mean(weightY.tempo0*iCol)}),
            apply(iid.event1,2,function(iCol){mean(weightY.tempo1*iCol)})
        )

        weightE.tempo0 <- data[,(1-.SD$treatment.bin) * .SD$prob.event0 / (1-.SD$prob.treatment)^2]
        weightE.tempo1 <- data[,.SD$treatment.bin * .SD$prob.event1 / .SD$prob.treatment^2]
        
        AIPWaddE <- t(apply(iid.treatment,2,function(iCol){
            c(-mean(weightE.tempo0*iCol),mean(weightE.tempo1*iCol))
        }))

        IF$AIPWnaive2 <- IF$AIPWnaive + (IPWadd + AIPWaddY + AIPWaddE) * data$weights * n.obs / sumW
        if(efficient){
            IF$AIPWefficient2 <- IF$AIPWefficient + (IPWadd + AIPWaddY + AIPWaddE) * data$weights * n.obs / sumW
        }
    }
    
    ## ** export
    out <- list()

    ## value
    n.estimator <- length(IF)
    name.risk <- paste0("risk.",level.treatment)
    out$ate.value <- matrix(NA, nrow = 3, ncol = n.estimator,
                            dimnames = list(c(name.risk,"ate.diff"),
                                              names(IF)))
    for(iL in 1:n.estimator){ ## iL <- 1
        if(na.rm){
            IF[[iL]] <- IF[[iL]][which(rowSums(is.na(IF[[iL]]))==0),,drop=FALSE]
        }
        
        out$ate.value[name.risk,iL] <- colSums(IF[[iL]])
        IF[[iL]] <- rowCenter_cpp(IF[[iL]], center = out$ate.value[paste0("risk.",level.treatment),iL]/n.obs)
    }
    out$ate.value["ate.diff",] <- out$ate.value[name.risk[2],] - out$ate.value[name.risk[1],]
    ##    out$ate.value["ate.ratio",] <- out$ate.value[name.risk[2],] / out$ate.value[name.risk[1],]

    ## standard error
    out$ate.se <- do.call(cbind,lapply(IF, function(iIF){ ## iIF <- IF[[1]]
        sqrt(c(colSums(iIF^2), sum((iIF[,2]-iIF[,1])^2)))
    }))
    rownames(out$ate.se) <- rownames(out$ate.value)

    class(out) <- "ateRobust"
    out$nuisance.iid <- nuisance.iid
    out$efficient <- efficient
    return(out)
    

}

## * .calcLterm
.calcLterm <- function(data, data0, data1, coxMF,
                       n.obs, times, type,
                       model.censor, model.event,
                       predictor.cox, product.limit, cause,
                       modeSurvival){

    status <- stop <- NULL ## [:CRANtest:] data.table

    eXb.censor <- exp(coxLP(model.censor, data = NULL, center = FALSE))
    MF.censor <- coxModelFrame(model.censor)

    jump.times <- sort(MF.censor[status==1&stop<=times,stop]) ## only look at jump times before the prediction time
    resBaseline.censor <- predictCox(model.censor,
                                     times = jump.times,
                                     type = "hazard",
                                     keep.infoVar = TRUE)
    
    is.strata <- resBaseline.censor$infoVar$is.strata

    calcLterm <- function(status, time, jump.time, hazard0, eXb, Esp.event, surv.censor){
        ## jump times are before the prediction time
        index.beforeJump <- which(jump.time<=time)
        if(length(index.beforeJump)==0){
            return(0)
        }else{
            dN <- (jump.time[index.beforeJump]==time)*status
            dLambda <- hazard0[index.beforeJump]*eXb
            return(sum(Esp.event[index.beforeJump]*(dN-dLambda)/surv.censor[index.beforeJump]))
        }
    }

    data[,c("Lterm") := as.numeric(NA)]

    ## note
    if(!is.strata){
        hazard0 <- resBaseline.censor$hazard
        all.times <- resBaseline.censor$times ## same as jump.times

        ## survival = P[C>t] = P[Delta(t)==1] - ok
        pred.censor <- do.call(predictor.cox,
                               args = list(model.censor, newdata = data, times = all.times-(1e-10), type = "survival"))$survival

        if(type=="survival"){
            
            pred.surv <- do.call(predictor.cox,
                                 args = list(model.event, newdata = data, times = all.times, type = "survival"))$survival

            if(modeSurvival == "survival"){
                pred.eventC <- colMultiply_cpp(1/pred.surv, scale = data[["prob.event"]])
            }else{
                pred.event <- 1 - pred.surv
                pred.eventC <- -colCenter_cpp(pred.event, center = data[["prob.event"]]) / pred.surv
            }

        }else if(type == "competing.risks"){
            pred.event <- predict(model.event, newdata = data, times = all.times,
                                  product.limit = product.limit, cause = cause)$absRisk

            ## all cause survival
            if(predictor.cox=="predictCox"){
                pred.cumhazard <- matrix(0, nrow = n.obs, ncol = length(all.times))
                for(iC in 1:length(model.event$models)){
                    pred.cumhazard <- pred.cumhazard + predictCox(model.event$models[[iC]], newdata = data, times = all.times-(1e-10), type = "cumhazard")$cumhazard
                }
                pred.surv <- exp(-pred.cumhazard)

                ## check computation of the survival
                ## data[["status"]] <- unclass(model.event$response)[,"status"]                
                ## eval(parse(text = paste0("coxAll <- coxph(",deparse(formula(model.event$models[[1]])),", data = data, x = TRUE)")))               
                ## pred.surv - predictCox(coxAll, newdata = data, times = all.times, type = "survival")$survival

            }else if(predictor.cox=="predictCoxPL"){
                all.jump.times <- model.event$eventTimes[model.event$eventTimes<=times]
                pred.allHazard <- matrix(0, nrow = n.obs, ncol = length(all.jump.times))
                for(iC in 1:length(model.event$models)){
                    pred.allHazard <- pred.allHazard + predictCox(model.event$models[[iC]],
                                                                  newdata = data,
                                                                  times = all.jump.times,
                                                                  type = "hazard")$hazard
                }
                ## check computation of the survival
                ## data[["status"]] <- unclass(model.event$response)[,"status"]                
                ## eval(parse(text = paste0("coxAll <- coxph(",deparse(formula(model.event$models[[1]])),", data = data, x = TRUE)")))               
                ## hazAll <- predictCox(coxAll, newdata = data, times = all.jump.times, type = "hazard")$hazard
                ## hazAll - pred.hazard
                index.jump <- prodlim::sindex(eval.times = all.times,
                                              jump.times = c(0,all.jump.times))
                pred.allSurv <- t(apply(1-pred.allHazard,1,cumprod))
                pred.surv <- cbind(0,pred.allSurv)[,index.jump,drop=FALSE]
                
                ## pred.surv - predictCoxPL(coxAll, newdata = data, times = all.times, type = "survival")$survival
            }
            ## conditional risk
            pred.eventC <- -colCenter_cpp(pred.event, center = data[["prob.event"]]) / pred.surv
            ## same as using option landmark but can be done in one go
            ## iIndex <- 15
            ## pred.eventC2 <- predict(model.event, newdata = data, times = times, landmark = all.times[iIndex],
            ## product.limit = product.limit, cause = cause)$absRisk
            ## pred.eventC[,iIndex]-pred.eventC2[,1]

        }
        
        data[,"Lterm" := sapply(1:n.obs, function(iObs){ ## iObs <- 1
            
            calcLterm(status = MF.censor$status[iObs],
                      time = MF.censor$stop[iObs],
                      jump.time = all.times,
                      hazard0 = hazard0,
                      eXb = eXb.censor[iObs],
                      Esp.event = pred.eventC[iObs,],
                      surv.censor = pred.censor[iObs,])
            
        })]
        ## data[,mean(Lterm)]
        
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

                iPred.surv <- do.call(predictor.cox,
                                      args = list(model.event, newdata = data[iIndex.strata], times = iAll.times, type = "survival"))$survival

                if(modeSurvival=="survival"){
                    iPred.eventC <- colMultiply_cpp(1/iPred.surv, scale = data[iIndex.strata,.SD$prob.event])
                }else{
                    iPred.event <- 1 - iPred.surv                
                    iPred.eventC <- -colCenter_cpp(iPred.event, center = data[iIndex.strata,.SD$prob.event]) / iPred.surv 
                }
                
            }else if(type == "competing.risks"){
                iPred.event <- predict(model.event, newdata = data[iIndex.strata], times = iAll.times,
                                       product.limit = product.limit, cause = cause)$absRisk

                ## all cause survival
                if(predictor.cox=="predictCox"){
                    iPred.cumhazard <- matrix(0, nrow = length(iIndex.strata), ncol = length(iAll.times))
                    for(iC in 1:length(model.event$models)){
                        iPred.cumhazard <- iPred.cumhazard + do.call(predictor.cox,
                                                                     args = list(model.event$models[[iC]], newdata = data[iIndex.strata],
                                                                                 times = iAll.times-(1e-10), type = "cumhazard"))$cumhazard
                    }
                    iPred.surv <- exp(-iPred.cumhazard)
                }else if(predictor.cox=="predictCoxPL"){
                     iAll.jump.times <- model.event$eventTimes[model.event$eventTimes<=times]
                     iPred.allHazard <- matrix(0, nrow = length(iIndex.strata), ncol = length(all.jump.times))
                     for(iC in 1:length(model.event$models)){
                         iPred.allHazard <- iPred.allHazard + predictCox(model.event$models[[iC]],
                                                                         newdata = data[iIndex.strata],
                                                                         times = iAll.jump.times,
                                                                         type = "hazard")$hazard
                     }
                     iIndex.jump <- prodlim::sindex(eval.times = iAll.times,
                                                    jump.times = c(0,iAll.jump.times))
                     iPred.allSurv <- t(apply(1-iPred.allHazard,1,cumprod))
                     iPred.surv <- cbind(0,iPred.allSurv)[,iIndex.jump,drop=FALSE]
                
                     ## pred.surv - predictCoxPL(coxAll, newdata = data, times = all.times, type = "survival")$survival
                }

                iPred.eventC <- -colCenter_cpp(iPred.event, center = data[iIndex.strata,.SD$prob.event]) / iPred.surv
                ## same as using option landmark but can be done in one go
                ## iIndex <- 15
                ## iPred.eventC2 <- predict(model.event, newdata = data[iIndex.strata], times = times, landmark = iAll.times[iIndex],
                ## product.limit = product.limit, cause = cause)$absRisk
                ## range(iPred.eventC[,iIndex]-iPred.eventC2[,1])

            }

            ## table(is.na(iPred.eventC))
            ## table(is.na(iPred.censor))

            data[iIndex.strata, c("Lterm") := sapply(1:iN.strata, function(iObs){ ## iObs <- 10

                iOut <- calcLterm(status = MF.censor$status[iIndex.strata[iObs]],
                                  time = MF.censor$stop[iIndex.strata[iObs]],
                                  jump.time = iAll.times,
                                  hazard0 = iHazard0,
                                  eXb = eXb.censor[iIndex.strata[iObs]],
                                  Esp.event = iPred.eventC[iObs,],
                                  surv.censor = iPred.censor[iObs,])
                return(iOut)
            
            })]

        ##        data[iIndex.strata, mean(Lterm)]
        }
        
            
    }

    return(NULL)
}

######################################################################
### ateRobust.R ends here
