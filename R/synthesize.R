### synthesize.R ---
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff Thomas Alexander Gerds
## Created: Apr 28 2021 (09:04)
## Version:
## Last-Updated: May 14 2025 (06:49) 
##           By: Thomas Alexander Gerds
##     Update #: 153
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##' Synthesizes survival data (also works for linear models and generalized linear models).
##' The idea is to be able to simulate new data sets that mimic the original data.
##' See the vignette \code{vignette("synthesize",package = "riskRegression")} for more details.
##'
##' The simulation engine is: lava.
##' @title Cooking and synthesizing survival data
##' @description Fit parametric regression models to the outcome distribution and optionally
##' also parametric regression models for the joint distribution of the predictors
##' structural equation models. 
##' Then the function \code{sim.synth} can be called on the  resulting object to
##' to simulate from the parametric model based on the machinery of the \code{lava} package
##' @aliases synthesize.formula synthesize.lvm 
##' @param object Specification of the synthesizing model structures. Either a \code{formula} or a \code{lvm} object. See examples.
##' @param data Data to be synthesized.
##' @param recursive Let covariates recursively depend on each other.
##' @param max.levels Integer used to guess which variables are categorical. When set to \code{10}, the default,
##'                   variables with less than 10 unique values in data are treated as categorical.
##' @param verbose Logical. If \code{TRUE} then more messages and warnings are provided.
##' @param logtrans Vector of covariate names that should be log-transformed. This is primarily for internal use.
##' @param fix.names Fix possible problematic covariate names. 
##' @param return_code Logical. If \code{TRUE} return the R-code instead of the lava object.
##' @param ... Not used yet.
##' @return lava object
##' @seealso lvm
##' @examples
##' # pbc data
##' library(survival)
##' library(lava)
##' data(pbc)
##' pbc <- na.omit(pbc[,c("time","status","sex","age","bili")])
##' pbc$logbili <- log(pbc$bili)
##' v_synt <- synthesize(object=Hist(time,status)~sex+age+logbili,data=pbc)
##' set.seed(8)
##' d <- simsynth(v_synt,38)
##' fit_sim <- coxph(Surv(time,status==1)~age+sex+logbili,data=d)
##' fit_real <- coxph(Surv(time,status==1)~age+sex+logbili,data=pbc)
##' # compare estimated log-hazard ratios between simulated and real data
##' cbind(coef(fit_sim),coef(fit_real))
##'
##' # return the simulation code instead of the object 
##' synthesize(object=Hist(time,status)~logbili+age+sex,data=pbc,return_code=TRUE)
##' # recursive
##' synthesize(object=Hist(time,status)~logbili+age+sex,
##'              data=pbc,
##'       return_code=TRUE,
##'         recursive=TRUE)
##'
##' synthesize(object=logbili~age+sex,
##'              data=pbc,
##'       return_code=TRUE,
##'         recursive=TRUE)
##'
##' u <- lvm()
##' distribution(u,~sex) <- binomial.lvm()
##' distribution(u,~age) <- normal.lvm()
##' distribution(u,~trt) <- binomial.lvm()
##' distribution(u,~logbili) <- normal.lvm()
##' u <-eventTime(u,time~min(time.cens=0,time.transplant=1,time.death=2), "status")
##' lava::regression(u,logbili~age+sex) <- 1
##' lava::regression(u,time.transplant~sex+age+logbili) <- 1
##' lava::regression(u,time.death~sex+age+logbili) <- 1
##' lava::regression(u,time.cens~1) <- 1
##' transform(u,logbili~bili) <- function(x){log(x)}
##' u_synt <- synthesize(object=u, data=na.omit(pbc))
##' saveSynth(u_synt)
##' set.seed(8)
##' d <- simsynth(u_synt,n=1000)
##' # note: synthesize may relabel status variable
##' fit_sim <- coxph(Surv(time,status==1)~age+sex+logbili,data=d)
##' fit_real <- coxph(Surv(time,status==1)~age+sex+log(bili),data=pbc)
##' # compare estimated log-hazard ratios between simulated and real data
##' cbind(coef(fit_sim),coef(fit_real))
##'
##' #
##' # Cancer data
##' #
##' data(cancer)
##' b <- lvm()
##' distribution(b,~rx) <- binomial.lvm()
##' distribution(b,~age) <- normal.lvm()
##' distribution(b,~resid.ds) <- binomial.lvm()
##' distribution(b,~ecog.ps) <- binomial.lvm()
##' lava::regression(b,time.death~age+rx+resid.ds) <- 1
##' b<-eventTime(b,futime~min(time.cens=0,time.death=1), "fustat")
##' b_synt <- synthesize(object = b, data = ovarian)
##' D <- simsynth(b_synt,1000)
##' fit_real <- coxph(Surv(futime,fustat)~age+rx+resid.ds, data=ovarian)
##' fit_sim <- coxph(Surv(futime,fustat)~age+rx+resid.ds, data=D)
##' cbind(coef(fit_sim),coef(fit_real))
##' w_synt <- synthesize(object=Surv(futime,fustat)~age+rx+resid.ds, data=ovarian)
##' D <- simsynth(w_synt,1000)
##' fit_sim <- coxph(Surv(futime,fustat==1)~age+rx+resid.ds,data=D)
##' fit_real <- coxph(Surv(futime,fustat==1)~age+rx+resid.ds,data=ovarian)
##' # compare estimated log-hazard ratios between simulated and real data
##' cbind(coef(fit_sim),coef(fit_real))
##'
##'
##' @export
##' @author Johan Sebastian Ohlendorff <johan.ohlendorff@@sund.ku.dk>  and Thomas A. Gerds <tag@@biostat.ku.dk> 
synthesize <- function(object, data,...){
  requireNamespace("lava")
  UseMethod(generic = "synthesize",object=object)
}

##' @export synthesize.formula
##' @export
#' @rdname synthesize
#' @method synthesize formula
synthesize.formula <- function(object, # a formula object Surv(time,event) or Hist(time,event)
                               data,
                               recursive=FALSE,
                               max.levels=10,
                               verbose=FALSE,
                               return_code = FALSE,
                               ...){
  requireNamespace("lava",quietly=TRUE)
  if (missing(return_code)) return_code <- FALSE
  # time to event outcome
  tt <- all.vars(update(object,".~1"))
  # covariates
  vv <- all.vars(formula(delete.response(terms(object))))
  fml <- paste(vv,collapse = "+")
  vv.withtrans <- attr(terms(formula(delete.response(terms(object)))),"term.labels")
  logtrans <- vv[vv!=vv.withtrans]
  # could be done more easily if the order of the covariates is not important,
  # but it might be due to recursive
  # potentially being set to TRUE
  # this changes the formula to be logtransformed
  hasLog <- length(logtrans)>0
  for (v in logtrans){
    if (is.null(data[[v]])){
      if (verbose) warning(paste0("Covariate ",v," is not in data set. Removing it from the formula."))
      fml <- gsub(paste0("[+]",v),"",fml)
      vv <- logtrans[logtrans!=v]
    }
    else {
      data[[paste0("log",v)]] <- log(data[[v]])
      s<-paste0("log",v)
      fml <- gsub(v,s,fml)
    }
  }
  # check if covariates are categorical
  if (hasLog){
    #include transformed variables in model
    vv1 <- vv
    vv <-gsub("[(]|[)]","",vv.withtrans)
  }
  object <- update(object,as.formula(paste0("~",fml)))
  object <- lava::lvm(object)
  # specify distributions of covariates in the lava object
  for (v in vv){
    if (is.null(data[[v]])){
      if(verbose) warning(paste0("Covariate ",v," is not in data set. Removing it from the formula"))
      fml <- gsub(paste0("[+]",v),"",fml)
      vv <- vv[vv!=v]
      object <- lava::rmvar(object,v)
    }
    else if (categorize(v,max.levels,data) == 2){
      data[[v]] <- factor(data[[v]])
      lava::distribution(object,v) <- lava::binomial.lvm()
    }
    else if (categorize(v,max.levels,data) == 1){
      data[[v]] <- factor(data[[v]])
      object <- lava::categorical(object, v, K=length(unique(data[[v]])))
    }
    else {
      lava::distribution(object,v) <- lava::normal.lvm()
    }
  }
  # include recursive structure of covariates if requested
  if (recursive){
    cv <- vv
    while (length(cv) > 1) {
      cv2 <- cv[2:length(cv)]
      #lacks implementation in case it is categorical
      #needs multinomial logistic regression to be implemented
      if (categorize(cv[1],max.levels,data) == 1){
        warning("Categorical variables not supported with recursive=TRUE. Skipping regression")
      }
      else {
        if (categorize(cv[1],max.levels,data) == 2){
          lava::distribution(object,cv[1]) <- lava::binomial.lvm()
        }
        regression(object) <- as.formula(paste0(cv[1],"~",paste(cv2,collapse = "+")))
      }
      cv <- cv2
    }
  }
  # outcome
  # binary/continuous outcome
  if (length(tt)==1){
    if (is.null(data[[tt[1]]])){stop("Response variable not in data set. Add it to the data, then try again.")}
    len.tt <- length(unique(data[[tt[[1]]]]))
    is.binary <- length(unique(data[[tt[[1]]]]))==2
    if (len.tt == 2){
      lava::distribution(object,tt[[1]]) <- lava::binomial.lvm()
    }
    else if (len.tt > 2 && len.tt < max.levels){
      warning("Categorical responses are not supported (yet).")
      object <- lava::categorical(object, tt[1], K=len.tt)
    }
    else {
      lava::distribution(object,tt[[1]]) <- lava::normal.lvm()
    }
  }else{
    # event time outcome
    if(is.null(data[[tt[1]]]) || is.null(data[[tt[2]]])) {
      stop("Response variable not in data set. Add it to the data, then try again.")
    }
    events <- sort(unique(data[[tt[[2]]]]))
    et.formula <- formula(paste0(tt[[1]]," ~ min(",paste(paste0("time.event.",events,"=",events),collapse=", "),")"))
    object <- lava::eventTime(object,et.formula, tt[[2]])
  }
  # call the next method, now the class of object is 'lvm' :)
  attr(x=object,which="from.formula") <- TRUE
  out <- synthesize(object=object,data=data,verbose=verbose,logtrans=logtrans,return_code = return_code,...)
  out
}


#' @export synthesize.lvm
#' @export
#' @rdname synthesize
#' @method synthesize lvm
synthesize.lvm <- function(object,
                           data,
                           max.levels = 10,
                           logtrans = NULL,
                           verbose=FALSE,
                           fix.names = FALSE,
                           return_code = FALSE,
                           ...){
  from.formula <- length(attr(object,"from.formula"))>0
  # check whether variables in model are in data set
  if (!from.formula && !all(object$attributes$eventHistory$time$names %in% names(data))) {
    stop("Time or status variable could not be found in data set.")
  }
  var.model <- colnames(object$M)
  # should deal with transformations check if non transformed are in data and if not add them
  # use grep to find transformed variables
  has.eventTime <- length(object$attributes$eventHistory)>0
  if (has.eventTime){
    timename <- names(object$attributes$eventHistory)
    timeplusevent <- object$attributes$eventHistory[[timename]]$names
    if (!(timeplusevent[1] %in% names(data)) || !(timeplusevent[2] %in% names(data))) {
      if (verbose) warning("time/status not found. Removing it")
      object <- lava::rmvar(object,timeplusevent[1])
      object <- lava::rmvar(object,timeplusevent[2])
    }
  }
  else {
    timename <- NULL
    timeplusevent <- NULL
  }
  # should check group transform
  # if (any(grepl("grp|group",var.model))) {
  #   stop("")
  # }
  # should check logtransform in data and add if necessary
  if (!from.formula && is.null(logtrans) && any(grepl("log",var.model))){
    #these have log in front
    trans <- var.model[grepl("log",var.model)]
    #log trans dont have log in front
    logtrans <- sub("log","",trans)
    for (v in logtrans){
      if (is.null(data[[v]]) && is.null(data[[paste0("log",v)]]) ){
        if(verbose) warning(paste0("Could not find the variable: ", v, "in the data"))
        object <- lava::rmvar(object,v)
        object <- lava::rmvar(object,paste0("log",v))
        trans <- trans[trans!=paste0("log",v)]
        logtrans <- logtrans[logtrans!=v]
      }
      else {
        data[[paste0("log",v)]] <- log(data[[v]])
      }
    }
  }
  # Check if all variables in data also occur in object
  # here we need to remove them if they are not
  if(!all(others <- names(data) %in% dimnames(object$M)[[1]])){
    if (verbose)
      warning("Some variables in dataset are not in object (or the names don't match).\n These variables are not synthesized:\n",
              paste0(names(data)[!others],collapse="\n"))
    # if(data.table::is.data.table(object))
    #     data <- data[,dimnames(object$M)[[1]],with=FALSE]
    # elsesub("log","",logtrans)
    #     data <- data[,dimnames(object$M)[[1]],drop=FALSE]
  }
  #check if variables in model are in data
  
  if (!all(others <- dimnames(object$M)[[1]] %in% names(data))){
    if (verbose) warning("Some variables in object are not in dataset (or the names don't match).\n These variables are not synthesized:\n",paste0(dimnames(object$M)[[1]][!others],collapse="\n"))
    #we have checked the logtransformed
    miss <- dimnames(object$M)[[1]][!others]
    if (length(logtrans) > 0) miss <- intersect(trans,intersect(logtrans,miss))
    if (has.eventTime) miss <- setdiff(miss, object$attributes$eventHistory[[timename]]$latentTimes)
    for (v in miss){
      object <- lava::rmvar(object,v)
    }
  }
  
  # note: will be a problem if there are NAs in data, should check at beginning of function
  
  # find intersection between variables in model with variables in data (these are the actual variables of the lava object)
  var.model <- colnames(object$M)
  var.model <- intersect(var.model,names(data))
  ismissingvar <- sapply(var.model, function(x) {anyNA(data[[x]])})
  if (any(ismissingvar )) {
      missing.var <- names(ismissingvar[ismissingvar])
      printvar <- paste(missing.var[1:min(length(missing.var),20)],collapse = ", ")
      stop(paste0("There should not be NAs for the variables in the model. The following variables have NAs: \n ",printvar))
  }
  if (return_code){
      simulation_code <- "sim_model <- lava::lvm()"
  }else{
      sim_model <- lava::lvm()
  }
  latent_vars <- endogenous(object)
  # should check factors in object are factors in data. If not, transform them
  if (!from.formula){
      for (v in var.model){
          #ignore transformed variables do not take the variables from the event history model
          if (!(v %in% object$attributes$eventHistory[[timename]]$names) && categorize(v,max.levels,data) != 0){
              data[[v]] <- factor(data[[v]])
          }
      }
  }
  
  # Set distributions:
  # 1. normal
  # 2. binomial
  # 3. categorical
  # deal with exogenous variables (which are without covariates)
  labels <- list()
  categorical.vars <- c()
  for(var in exogenous(object)){
    var_formula <- as.formula(paste0("~", var))
    #case: gaussian
    if("gaussian" %in% attributes(object$attributes$distribution[[var]])$family){
        if (return_code){
            simulation_code <- c(simulation_code,list(paste0("lava::distribution(sim_model,",paste0("~", var),") <- lava::normal.lvm(mean=",mean(data[[var]]),",","sd=",sd(data[[var]]),")")))
        }else{
            lava::distribution(sim_model,var_formula) <- lava::normal.lvm(mean=mean(data[[var]]),sd=sd(data[[var]]))
        }
    }
      #case: binary
    else if("binomial" %in% attributes(object$attributes$distribution[[var]])$family){
        labels[[var]] <- levels(data[[var]])
        if (return_code){
            simulation_code <- c(simulation_code,list(paste0("lava::distribution(sim_model,","~", var,") <- lava::binomial.lvm(p=",mean(factor(data[[var]])==levels(factor(data[[var]]))[2]),")")))
        }else{
            lava::distribution(sim_model, var_formula) <- lava::binomial.lvm(p=mean(factor(data[[var]])==levels(factor(data[[var]]))[2]))
        }
    }
      #case: categorical
    else if(var %in% names(object$attributes$nordinal)){
        categorical.vars <- c(categorical.vars, var)
        num_cat <- object$attributes$nordinal[[var]] #what if levels is NULL?
        cat_probs <- vector(mode="numeric", length=num_cat-1)
        for(i in 1:(num_cat-1)){
            cat_probs[i] <- mean(data[[var]] == levels(data[[var]])[i])
        }
        if (return_code){
            simulation_code <- c(simulation_code,list(paste0("sim_model <- lava::categorical(sim_model,",
                                                             paste0("~", var),
                                                             ",",
                                                             "labels=c('",
                                                             paste0(levels(data[[var]]),collapse = "','"),
                                                             "'),",
                                                             "K=",
                                                             num_cat,
                                                             ",",
                                                             "p=c(",
                                                             paste0(cat_probs,collapse = ","),
                                                             "))")))
        }else{
            sim_model <- lava::categorical(sim_model,var_formula,labels=levels(data[[var]]),K=num_cat,p=cat_probs)
        }
    }
    else {
        stop("Error. Distribution not supported!")
    }
  }
  # dichotomize categorical variables
  # should not do this if the cat. variable is not in lvm.object
  # note: name of dich. variables must correspond with names used in lava::regression formulas outside synthesize-func.
  dichotomized_variables <- c()
  for(var in colnames(data)[grepl('factor|character', sapply(data, class))]){
      #fix factors with weird characters including binary
      #note: we do not fix formulas and weirdly named variables
      #  can give problems if synthesized directly from lvm
      # don't change 0/1 variables even if they are factors
      if (from.formula || fix.names){
          # remove warning, i.e. don't compare if the number of levels is greater than 2
          if (length(levels(data[[var]])) ==2 && all(levels(data[[var]]) == c("0","1"))){
          }
          else {
              lvls <- make.names(levels(data[[var]]))
              levels(data[[var]]) <- lvls
          }
      }
      if(!any(grepl(var,  dimnames(object$M)[[1]]))){
          #variable not in lvm object
      } else if (length(levels(data[[var]]))>2){
          dichotomized_variables <- c(dichotomized_variables, var)
          for(lvl in levels(data[[var]])[-1]){
              data[[paste0(var, lvl)]] <- 1*(data[[var]] == lvl)
              formula <- as.formula(paste0(var, lvl,"~",var))
              # updates the formula in a local environment; otherwise the lvl will be tied to the last value of lvl in the for loop
              if (return_code){
                  simulation_code <- c(simulation_code,list(paste0("lava::transform(sim_model, ",
                                                                   var,
                                                                   lvl,
                                                                   "~",
                                                                   var,
                                                                   ") <- function(x){1*(c(x)=='",
                                                                   lvl,
                                                                   "')}")))
              }else{
                  sim_model <- local({
                      l <- lvl
                      lava::transform(sim_model, formula) <- function(x){1*(x==l)}
                      return(sim_model)})
              }
          }
      }
  }
  # get covariates as string and include dichotomized variables instead of the original ones
  get_covariates <- function(obj, v, dichotomized_variables) {
      if (!(v %in% dimnames(obj$M)[[2]])){
          return("1")
      }
      covariates <- dimnames(obj$M)[[2]][obj$M[,v] == 1]
      if (length(covariates) == 0){
          return("1")
      }
      else {
          fml <- paste(covariates,collapse = "+")
          for (var in dichotomized_variables){
              lvl <- paste0(var,levels(data[[var]])[-1])
              cv <- paste0(lvl,collapse = "+")
              fml <- gsub(var,cv,fml)
          }
          return(fml)
      }
  }
  
  # define latent event time variables
  # note: will be a problem if there are several eventTime variables. Should loop through all these
  has.eventTime <- length(object$attributes$eventHistory)>0
  if (has.eventTime){
    covariates <- get_covariates(object,timename,dichotomized_variables)
    response1 <- paste(object$attributes$eventHistory[[timename]]$names, collapse = ",")
    
    #include event time variables and fit the regression model
    for (latvar in object$attributes$eventHistory[[timename]]$latentTime){
        if (!from.formula){
            covariates <- get_covariates(object,latvar,dichotomized_variables)
        }
        latvar_formula <- as.formula(paste0("~", latvar))
        status_ind <- object$attributes$eventHistory[[timename]]$events[which(object$attributes$eventHistory[[timename]]$latentTimes %in% latvar)]
        response <- paste0(response1, "==", status_ind)
        surv_formula <- as.formula(paste0("Surv(", response, ")~", covariates))
        if (!(inherits(status_ind,"numeric"))) {stop("event or status variable has to be numeric")}
        G <- survreg(surv_formula, data = data)
        reg_formula <- as.formula(paste0(latvar, "~", covariates))
        if (return_code){
            simulation_code <- c(simulation_code,list(paste0("lava::distribution(sim_model,~",latvar,")"," <- lava::coxWeibull.lvm(scale=",exp(-coef(G)["(Intercept)"]/G$scale),",shape=",1/G$scale,")")))
            simulation_code <- c(simulation_code,list(paste0("lava::regression(sim_model,",latvar,"~",covariates,")"," <- c(",paste0(-coef(G)[-1]/G$scale,collapse = ","),")")))
        }else{
            lava::distribution(sim_model, latvar_formula) <- lava::coxWeibull.lvm(scale=exp(-coef(G)["(Intercept)"]/G$scale),shape=1/G$scale)
            lava::regression(sim_model, reg_formula) <- -coef(G)[-1]/G$scale
        }
    }
    #calculate time and status afterwards
    events <- object$attributes$eventHistory[[timename]]$events
    if (!from.formula){
        cens <- object$attributes$eventHistory[[timename]]$latentTime
        et_formula_string <- paste0(timename," ~ min(",paste0(cens,"=",events,collapse=", "),")")
        et.formula <- formula(paste0(timename," ~ min(",paste0(cens,"=",events,collapse=", "),")"))
    }
    else {
        et_formula_string <- paste0(timename," ~ min(",paste(paste0("time.event.",events,"=",events),collapse=", "),")")
        et.formula <- formula(paste0(timename," ~ min(",paste(paste0("time.event.",events,"=",events),collapse=", "),")"))
    }
    if (return_code){
        simulation_code <- c(simulation_code,list(paste0("sim_model <- lava::eventTime(sim_model,",et_formula_string,",'",object$attributes$eventHistory[[timename]]$names[2],"')")))
    }else{
        sim_model <- lava::eventTime(sim_model,et.formula, object$attributes$eventHistory[[timename]]$names[2])
    }
  }
  # Estimate regression coefficients in real data
  # and add them to the lvm object using lava::regression
  for(var in latent_vars[!(latent_vars %in% object$attributes$eventHistory$time$names)]){
      covariates <- get_covariates(object,var,dichotomized_variables)
      reg_formula <- as.formula(paste0(var, "~", covariates))
      # we have three types of regression to deal with now. Either
      # 1. linear regression
      # 2. logistic regression
      # 3. multinomial logistic regression
      if ("binomial" %in% attributes(object$attributes$distribution[[var]])$family){
          # might also lack levels here
          labels[[var]] <- levels(data[[var]])
          fit <- glm(reg_formula,data=data,family="binomial")
          p0<-exp(coef(fit)[1])/(1+exp(coef(fit)[1]))
          if (return_code){
              simulation_code <- c(simulation_code,list(paste0("lava::distribution(sim_model,~",var,")"," <- lava::binomial.lvm(p=",p0,")")))
              simulation_code <- c(simulation_code,list(paste0("lava::regression(sim_model,",var,"~",covariates,")"," <- ","c(",paste0(coef(fit)[-1],collapse=","),")")))
          }else{
              lava::distribution(sim_model,as.formula(paste0("~", var))) <- lava::binomial.lvm(p=p0)
              lava::regression(sim_model,reg_formula)<-coef(fit)[-1]
          }
      }
      #case: categorical
      else if (var %in% dichotomized_variables){
          categorical.vars <- c(categorical.vars, var)
          if (verbose) warning("Synthesize untested for categorical variables")
          for(lvl in levels(data[[var]])[-1]){
              reg_formula <- as.formula(paste0(var, lvl, "~", covariates))
              fit <- glm(reg_formula,data=data,family="binomial")
              p0<-exp(coef(fit)[1])/(1+exp(coef(fit)[1]))
              if (return_code){
                  simulation_code <- c(simulation_code,list(paste0("lava::distribution(sim_model,~",var,")"," <- lava::binomial.lvm(p=",p0,")")))
                  simulation_code <- c(simulation_code,list(paste0("lava::regression(sim_model,",
                                                                   var,
                                                                   lvl,
                                                                   "~",
                                                                   covariates,
                                                                   ")",
                                                                   " <- ",
                                                                   "c(",
                                                                   paste0(coef(fit)[-1],collapse=","),
                                                                   ")")))
              }else{
                  lava::distribution(sim_model,as.formula(paste0("~", var))) <- lava::binomial.lvm(p=p0)
                  lava::regression(sim_model,reg_formula)<-coef(fit)[-1]
              }
          }
      }
      # ignore them
      else if (has.eventTime && var %in% object$attributes$eventHistory[[timename]]$latentTime){}
      # case: gaussian
      else {
          fit <- lm(reg_formula,data=data)
          if (return_code){
              simulation_code <- c(simulation_code,list(paste0("lava::distribution(sim_model,~",var,")"," <- lava::normal.lvm(mean=",coef(fit)[1],",sd=",summary(fit)$sigma,")")))
              simulation_code <- c(simulation_code,list(paste0("lava::regression(sim_model,",var,"~",covariates,")"," <- ","c(",paste0(coef(fit)[-1],collapse=","),")")))
          }else{
              lava::distribution(sim_model,as.formula(paste0("~", var))) <- lava::normal.lvm(mean = coef(fit)[1],sd = summary(fit)$sigma)
              lava::regression(sim_model,reg_formula)<-coef(fit)[-1]
          }
      }
      #else {if ("gaussian"%in% attributes(object$attributes$distribution[[var]])$family){
      #    stop("Distribution is not supported.")
      #}
  }
  
  #transform logtransformed covariates back
  for (v in logtrans){
      if (return_code){
          simulation_code <- c(simulation_code,list(paste0("lava::transform(sim_model,~log",v,") <- ","function(x){exp(x)}")))
      }else{
          transform(sim_model,as.formula(paste0(v,"~log",v))) <- function(x){exp(x)}
      }
  }
  if (return_code){
      ## cat(paste(simulation_code,collapse = "\n"),"\n")
      res <- paste0(paste(simulation_code,collapse = "\n"),"\n")
      class(res) <- c("synth_code",class(res))
      res
  }else{
      res <- list(lava.object = sim_model,labels = labels,categories = categorical.vars)
      class(res) <- c("synth",class(res))
      return(res)
  }
}

#returns categorical values
# 0=numeric, 1=categorical with more levels, 2=binary
categorize <- function(v,max.levels,data){
  is.cat <- is.factor(data[[v]])
  nu <- length(unique(data[[v]]))
  if (!is.cat) {is.cat <- (nu < max.levels)}
  is.bin <- nu==2
  return(is.bin+is.cat)
}

#----------------------------------------------------------------------
### synthesize.R ends here
