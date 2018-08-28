## ateRobust.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: jun 27 2018 (17:47) 
## Version: 
## Last-Updated: aug 28 2018 (11:07) 
##           By: Brice Ozenne
##     Update #: 917
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
#' n <- 1e3
#'
#' ## check bias 
#' set.seed(10)
#' dt <- as.data.table(lava::sim(mSimSurv, n = n, p = c(alpha = alpha, beta = beta)))
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
    reserved.name <- c("times","time.tau","status.tau","treatment.bin",
                       "prob.event","prob.event0","prob.event1",
                       "prob.treatment",
                       "weights","prob.censoring",
                       "Lterm0","Lterm1")
    if(any(reserved.name %in% names(data))){
        txt <- reserved.name[reserved.name %in% names(data)]
        stop("Argument \'data\' should not contain column(s) named \"",paste0(txt, collapse = "\" \""),"\"\n")
    }

    ## ** fit models
    ## event
    if(type == "survival"){
        model.event <- do.call(fitter, args = list(formula = formula.event, data = data, x = TRUE, y = TRUE))
        coxMF <- coxModelFrame(model.event)
        ##  start       stop status strata ...
    }else{
        model.event <- CSC(formula.event, data = data, fitter = fitter, cause = cause, surv.type = "hazard")
        coxMF <- as.data.frame(unclass(model.event$response))
        names(coxMF)[names(coxMF)=="time"] <- "stop"
        ## stop status event ...
    }
    n.censor <- sum(coxMF$status==0)

    
    ## treatment
    model.treatment <- do.call(glm, args = list(formula = formula.treatment, data = data, family = stats::binomial(link = "logit")))       
    data[, c("prob.treatment") := predict(model.treatment, newdata = data, type = "response", se = FALSE)]

    ## censoring
    if(n.censor>0){
        model.censor <- do.call(fitter, args = list(formula = formula.censor, data = data, x = TRUE, y = TRUE)) ## , no.opt = TRUE
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
        data[,c("time.tau") := pmin(coxMF$stop,times)]
        data[,c("status.tau") := (coxMF$stop<=times)*(coxMF$status==1) - (coxMF$stop<times)*(coxMF$status==0)]
        ## -1 censored, 0 survival, 1 event, 
    }else if(type=="competing.risks"){
        data[,c("time.tau") := pmin(coxMF$stop,times)]
        data[,c("status.tau") := (coxMF$event==cause)*(coxMF$stop<=times) - (coxMF$status==0)*(coxMF$stop<times)]
        ## -1 censored, 0 survival or death (i.e. no event), 1 event
    }
    ## table(data$status.tau)

    ## ** outcome model: conditional expectation
    if(nuisance.iid){
        ## Computation of the influence function (Gformula, AIPW)
        ## additional term: eq:IF-AIPWfull (green terms)
        ## this is sent to predictCox to multiply each individual IF before averaging

        nuisance.iid0 <- TRUE
        attr(nuisance.iid0, "factor") <- list("Gformula" = matrix(1, nrow = n.obs, ncol = 1),
                                              "AIPW" = cbind(data[, 1 - (1-.SD$treatment.bin) / (1-.SD$prob.treatment)])
                                              )
        nuisance.iid1 <- TRUE
        attr(nuisance.iid1, "factor") <- list("Gformula" = matrix(1, nrow = n.obs, ncol = 1),
                                              "AIPW" = cbind(data[, 1 - .SD$treatment.bin / .SD$prob.treatment])
                                              )
    }else{
        nuisance.iid0 <- FALSE
        nuisance.iid1 <- FALSE
    }
    
    if(type == "survival"){
        ## Estimation of the survival + IF
        prediction.event <- do.call(predictor.cox, args = list(model.event, newdata = data, times = times, type = "survival"))
        prediction.event0 <- do.call(predictor.cox, args = list(model.event, newdata = data0, times = times, type = "survival", average.iid = nuisance.iid0))
        prediction.event1 <- do.call(predictor.cox, args = list(model.event, newdata = data1, times = times, type = "survival", average.iid = nuisance.iid1))

        ## store results
        data[, c("prob.event") := 1 - prediction.event$survival[,1]]
        data[, c("prob.event0") := 1 - prediction.event0$survival[,1]]
        data[, c("prob.event1") := 1 - prediction.event1$survival[,1]]

        if(nuisance.iid){
            data[,iidG.event0 := -prediction.event0$survival.average.iid[[1]][,1]]
            data[,iidAIPW.event0 := -prediction.event0$survival.average.iid[[2]][,1]]
            data[,iidG.event1 := -prediction.event1$survival.average.iid[[1]][,1]]
            data[,iidAIPW.event1 := -prediction.event1$survival.average.iid[[2]][,1]]
        }

    }else if(type=="competing.risks"){
        
        ## Estimation of the survival + IF
        prediction.event <- predict(model.event, newdata = data, times = times, cause = cause, product.limit = product.limit)
        prediction.event0 <- predict(model.event, newdata = data0, times = times, cause = cause, product.limit = product.limit, iid = nuisance.iid)
        prediction.event1 <- predict(model.event, newdata = data1, times = times, cause = cause, product.limit = product.limit, iid = nuisance.iid)

        data[, c("prob.event") := prediction.event$absRisk[,1]]
        data[, c("prob.event0") := prediction.event0$absRisk[,1]]
        data[, c("prob.event1") := prediction.event1$absRisk[,1]]

        ## store results
        if(nuisance.iid){
            weight0 <- attr(nuisance.iid0, "factor")[["AIPW"]][,1]
            data[,iidG.event0 := colMeans(prediction.event0$absRisk.iid[,1,])]
            data[,iidAIPW.event0 := apply(prediction.event0$absRisk.iid[,1,],2,function(iCol){mean(weight0*iCol)})]

            weight1 <- attr(nuisance.iid1, "factor")[["AIPW"]][,1]
            data[,iidG.event1 := colMeans(prediction.event1$absRisk.iid[,1,])]
            data[,iidAIPW.event1 := apply(prediction.event1$absRisk.iid[,1,],2,function(iCol){mean(weight1*iCol)})]
        }

    }

    ## ** Propensity score model: weights
    if(nuisance.iid){

        factor <- cbind("IPW0" = data[,  (1-.SD$treatment.bin) * (.SD$status.tau==1) / (1-.SD$prob.treatment)^2],
                        "IPW1" = data[, .SD$treatment.bin * (.SD$status.tau==1) / (.SD$prob.treatment)^2],
                        "AIPW0" = data[, (1-.SD$treatment.bin) * .SD$prob.event0 / (1-.SD$prob.treatment)^2],
                        "AIPW1" = data[, .SD$treatment.bin * .SD$prob.event1 / .SD$prob.treatment^2])

        average.iid <- TRUE
        attr(average.iid, "factor") <- factor
        
        prediction.treatment.iid <- attr(.predictGLM(model.treatment, newdata = data, average.iid = average.iid), "iid")
        data[,iidIPW.treatment0 := prediction.treatment.iid[,1]]
        data[,iidIPW.treatment1 := -prediction.treatment.iid[,2]]
        data[,iidAIPW.treatment0 := -prediction.treatment.iid[,3]]
        data[,iidAIPW.treatment1 := prediction.treatment.iid[,4]]
        
    }
    
    ## ** Censoring model: weights
    ## fit model
    if(n.censor==0){
        data[,c("prob.censoring") := 1]
        data[,c("prob.indiv.censoring") := 1]
    }else{
        ## survival = P[C>min(T,tau)] = P[Delta(min(T,tau))==1] - ok

        ## stopped at tau
        data[,c("prob.censoring") := do.call(predictor.cox, args = list(model.censor, newdata = data, times = times))$survival[,1]]

        ## at each time
        data[,c("prob.indiv.censoring") := do.call(predictor.cox, args = list(model.censor, newdata = data, times = data$time.tau-(1e-10), diag = TRUE))$survival[,1]]
    }
    data[coxMF$stop<=times,c("weights") := (.SD$status.tau >= 0) / .SD$prob.indiv.censoring]
    data[coxMF$stop>times,c("weights") := (.SD$status.tau >= 0) / .SD$prob.censoring]

    ## ** correction for efficiency
    ## data is updated within the function .calcLterm
    if(efficient){
        if(n.censor==0){
            data[,c("Lterm") := 0]
        }else{
            data[Lterm := .calcLterm(data = data,
                                     n.obs = n.obs,
                                     times = times,
                                     model.censor = model.censor,
                                     model.event = model.event,
                                     type = type,
                                     predictor.cox = predictor.cox,
                                     product.limit = product.limit,
                                     cause = cause)]
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
        nuisanceEvent.Gformula <- cbind(data$iidG.event0, data$iidG.event1)
        IF$Gformula2 <- IF$Gformula + nuisanceEvent.Gformula
    }

    ## *** IPW
    ## eq:IF-IPWc (blue term)
    IF$IPWnaive <- data[,cbind(
        .SD$weights * (.SD$status.tau == 1) * (1-.SD$treatment.bin) / (1-.SD$prob.treatment),
        .SD$weights * (.SD$status.tau == 1) * (.SD$treatment.bin) / (.SD$prob.treatment)
    )]/ n.obs
    ## colSums(IF$IPWnaive)
    ## colSums(IF$Gformula)
    
    if(efficient){
        ## additional term: eq:IF-IPWeff (red term)
        efficient.IPW <- data[,cbind(
            .SD$Lterm * (1-.SD$treatment.bin) / (1-.SD$prob.treatment), ## .SD$prob.event
            .SD$Lterm * (.SD$treatment.bin) / (.SD$prob.treatment) ## .SD$prob.event
        )]/ n.obs
        IF$IPWefficient <- IF$IPWnaive + efficient.IPW
        ## colSums(efficient.IPW)
        ## sum(data$Lterm)/n.obs
        ## data$prob.treatment
    }
    
    if(nuisance.iid){
        ## additional term: eq:IF-IPWfull (green term)
        nuisanceTreatment.IPW <- colMultiply_cpp(cbind(data$iidIPW.treatment0, data$iidIPW.treatment1),
                                                 scale = data$weights)
        
        ## colMeans(IPWadd)
        IF$IPWnaive2 <- IF$IPWnaive + nuisanceTreatment.IPW
        if(efficient){
            IF$IPWefficient2 <- IF$IPWefficient + nuisanceTreatment.IPW
        }
    }

    ## *** AIPW
    ## additional term:
    AIPWadd <- data[,cbind(
        .SD$prob.event0 * (1-(1-.SD$treatment.bin)/(1-.SD$prob.treatment)),
        .SD$prob.event1 * (1-.SD$treatment.bin/(.SD$prob.treatment))
    )]/ n.obs
        
    IF$AIPWnaive <- IF$IPWnaive + AIPWadd

    if(efficient){
        IF$AIPWefficient <- IF$IPWefficient + AIPWadd
    }
    
    if(nuisance.iid){
        nuisanceEvent.AIPW <- cbind(data$iidAIPW.event0, data$iidAIPW.event1) ## no outcome here so no censoring
        nuisanceTreatment.AIPW <- cbind(data$iidAIPW.treatment0, data$iidAIPW.treatment1) ## no outcome here so no censoring
        
        IF$AIPWnaive2 <- IF$AIPWnaive + (nuisanceTreatment.IPW + nuisanceEvent.AIPW + nuisanceTreatment.AIPW)
        if(efficient){
            IF$AIPWefficient2 <- IF$AIPWefficient + (nuisanceTreatment.IPW + nuisanceEvent.AIPW + nuisanceTreatment.AIPW)
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
.calcLterm2 <- function(data, n.obs, times,
                        model.censor,
                        model.event, type, predictor.cox, product.limit, cause){

    info.censor <- coxVariableName(model.censor, data)
    timeVar.censor <- info.censor$time
    statusVar.censor <- info.censor$status
    
    X.censor <- coxModelFrame(model.censor)
    new.time <- X.censor$stop
    new.status <- X.censor$status
    
    jump.time <- sort(X.censor$stop[X.censor$status == 1]) ##  only select jumps
    jump.time <- jump.time[jump.time <= times] ## before time horizon
    njump <- length(jump.time)
  
    ## ** compute conditional risk
    riskTau <- matrix(data$prob.event, nrow = n.obs, ncol = njump, byrow = FALSE)
    if(type == "survival"){
        riskTime <- 1 - predictCox(model.event, newdata = data, times = jump.time, type = "survival")$survival
  
        riskConditional <- (riskTau - riskTime)/(1-riskTime)
    }else if(type == "competing.risk"){
        riskTime <- predict(model.event, newdata = data, times = jump.time,
                            cause = cause, product.limit = product.limit)$absRisk

        survTime <- predictSurv(model.event, newdata = data, times = jump.time, product.limit = product.limit)
        
        riskConditional <- (riskTau - riskTime)/(survTime)

        ## check        
        ## index.test <- 5
        ## GS <- predict(model.event, newdata = data, times = times, landmark = jump.time[index.test],
        ## cause = cause, product.limit = product.limit)$absRisk
        ## range(GS-riskConditional[,index.test])
        
    }
    
    ## ** at risk indicator
    atRisk <- do.call(rbind,
                      lapply(1:n.obs, function(iTime){
                          as.numeric(jump.time <= new.time[iTime])
                      }))
  
    ## ** counting process for the censoring 
    dN <- do.call(rbind,
                  lapply(1:n.obs, function(iTime){
                      (jump.time == new.time[iTime]) * new.status[iTime]
                  }))
  
    ## ** compensator for the censoring
    dLambda <- predictCox(model.censor, newdata = data, times = jump.time, type = "hazard")$hazard

    ## ** survival for censoring at t-
    if(njump>1){
        survCensoring <- do.call(predictor.cox,
                                 args = list(model.censor, newdata = data, times = jump.time, type = "survival")
                                 )$survival
        survCensoring <- cbind(1,survCensoring[,1:(njump-1)])
    }else{
        survCensoring <- matrix(1, nrow = n.obs, ncol = 1)
    }

    ## ** integral
    out <- rowSums(atRisk * riskConditional/survCensoring * (dN-dLambda))

    ## ** export
    return(out)
}


