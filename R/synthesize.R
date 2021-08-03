### synthesize.R ---
#----------------------------------------------------------------------
## Author: Johan Sebastian Ohlendorff & Vilde Hansteen Ung & Thomas Alexander Gerds
## Created: Apr 28 2021 (09:04)
## Version:
## Last-Updated: Jul 27 2021 (09:07)
##           By: Thomas Alexander Gerds
##     Update #: 52
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:
##' Synthesizing survival data
##'
##' The simulation engine is: lava.
##' @title Cooking and synthesizing survival data
##' @param object
##' @param data
##' @param ...
##' @return lava object
##' @seealso lvm
##' @examples
##' # pbc data
##' library(survival)
##' library(lava)
##' u <- lvm()
##' distribution(u,~sex) <- binomial.lvm()
##' distribution(u,~age) <- normal.lvm()
##' distribution(u,~trt) <- binomial.lvm()
##' distribution(u,~logbili) <- normal.lvm()
##' distribution(u,~logprotime) <- normal.lvm()
##' u <-eventTime(u,time~min(time.cens=0,time.transplant=1,time.death=2), "status")
##' lava::regression(u,logbili~age+sex) <- 1
##' lava::regression(u,protimegrp1~age+sex+logbili) <- 1
##' lava::regression(u,protimegrp2~age+sex+logbili) <- 1
##' lava::regression(u,stage3~age+sex+protimegrp1+protimegrp2+logbili) <- 1
##' lava::regression(u,stage4~age+sex+protimegrp1+protimegrp2+logbili) <- 1
##' lava::regression(u,time.transplant~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4) <- 1
##' lava::regression(u,time.death~sex+age+logbili+protimegrp1+protimegrp2+stage3+stage4) <- 1
##' lava::regression(u,time.cens~1) <- 1
##' u <- categorical(u,~stage,labels=c("1/2","3","4"), K=3)
##' transform(u,protime~logprotime) <- function(x){exp(x)}
##' transform(u,protimegrp~protime) <- function(x){cut(c(unlist(x)), c(-Inf,10,11,Inf), labels=c("0","1","2"))}
##' u_synt <- synthesize(object=u, data=na.omit(pbc))
##' set.seed(8)
##' d <- sim(u_synt,n=1000)
##' fit_sim <- coxph(Surv(time,status==1)~age+sex+logbili,data=d)
##' fit_real <- coxph(Surv(time,status==1)~age+sex+logbili,data=pbc)
##' # compare estimated log-hazard ratios between simulated and real data
##' cbind(coef(fit_sim),coef(fit_real))
##'
##' data(pbc)
##' pbc <- na.omit(pbc[,c("time","status","sex","age","bili")])
##' pbc$logbili <- log(pbc$bili)
##' v_synt <- synthesize(object=Surv(time,status)~sex+age+logbili,data=pbc)
##' d <- sim(v_synt,1000)
##' fit_sim <- coxph(Surv(time,status==1)~age+sex+logbili,data=d)
##' fit_real <- coxph(Surv(time,status==1)~age+sex+logbili,data=pbc)
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
##' D <- sim(b_synt,1000)
##' fit_real <- coxph(Surv(futime,fustat)~age+rx+resid.ds, data=ovarian)
##' fit_sim <- coxph(Surv(futime,fustat)~age+rx+resid.ds, data=D)
##' cbind(coef(fit_sim),coef(fit_real))
##' w_synt <- synthesize(object=Surv(futime,fustat)~age+rx+resid.ds, data=ovarian)
##' D <- sim(w_synt,1000)
##' fit_sim <- coxph(Surv(futime,fustat==1)~age+rx+resid.ds,data=D)
##' fit_real <- coxph(Surv(futime,fustat==1)~age+rx+resid.ds,data=ovarian)
##' # compare estimated log-hazard ratios between simulated and real data
##' cbind(coef(fit_sim),coef(fit_real))
##'
##'
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
synthesize <- function(object, data,...){
    requireNamespace("lava")
    UseMethod("synthesize",object=object)
}

##' @export synthesize.formula
##' @export
synthesize.formula <- function(object, # a formula object Surv(time,event) or Hist(time,event)
                               data,
                               covariates, # either 'lava' or 'empirical'
                               max.levels=10,
                               verbose=FALSE,
                               recursive=FALSE,
                               ...){
    requireNamespace("lava",quietly=TRUE)
    # time to event outcome
    tt <- all.vars(update(object,".~1"))
    # covariates
    vv <- all.vars(formula(delete.response(terms(object))))
    vv.withtrans <- attr(terms(formula(delete.response(terms(object)))),"term.labels")

    # could be done more easily if the order of the covariates is not important, but it might be due to recursive
    # potentially being set to TRUE
    # this changes the formula to be logtransformed
    form.string <- ""
    hasLog <- FALSE
    for (i in 1:length(vv)){
      if (vv[i]==vv.withtrans[i]){
        if (form.string == ""){
          form.string <- vv[i]
        }
        else {
          form.string <- paste(form.string, "+",vv[i])
        }
      }
      else {
        if (form.string == ""){
          data[[paste0("log",vv[i])]] <- log(data[[vv[i]]])
          form.string <-paste0("log",vv[i])
          hasLog <- TRUE
        }
        else {
          data[[paste0("log",vv[i])]] <- log(data[[vv[i]]])
          form.string <- paste(form.string,"+",paste0("log",vv[i]))
          hasLog <- TRUE
        }
      }
    }
    object <- update(object,as.formula(paste0("~",form.string)))
    # add all variables
    object <- lava::lvm(object)
    # check if covariates are categorical
    if (hasLog){
      #include transformed variables in model
      for (v in vv[vv!=vv.withtrans]){
        transform(object,as.formula(paste0(v,"~log",v))) <- function(x){exp(x)}
      }
      vv1 <- vv
      vv <-c(vv[vv==vv.withtrans],paste0("log",vv[vv!=vv.withtrans]))
    }

    #returns categorical values
    # 0=numeric, 1=categorical with more levels, 2=binary
    categorize <- function(v){
        is.cat <- is.factor(data[[v]])
        nu <- length(unique(data[[v]]))
        if (!is.cat) {is.cat <- (nu < max.levels)}
        is.bin <- nu==2
        return(is.bin+is.cat)
    }

    # specify distributions of covariates in the lava object
    for (v in vv){
      if (categorize(v) == 2){
        data[[v]] <- factor(data[[v]])
        lava::distribution(object,v) <- lava::binomial.lvm()
      }
      else if (categorize(v) == 1){
        data[[v]] <- factor(data[[v]])
        object <- lava::categorical(object, v, K=length(unique(data[[v]])))
      }
      else {
        lava::distribution(object,v) <- lava::normal.lvm()
      }
    }

    # include recursive structure of covariates if true
    if (recursive){
      cv <- vv
      while (length(cv) > 1) {
        cv2 <- cv[2:length(cv)]
        #lacks implementation in case it is categorical
        #needs multinomial logistic regression to be implemented
        if (categorize(cv[1]) == 1){
          warning("not implemented")
          object <- lava::categorical(object, cv[1], K=length(unique(data[[cv[1]]])))
        }
        else if (categorize(cv[1]) == 2){
          lava::distribution(object,cv[1]) <- lava::binomial.lvm()
        }
        regression(object) <- as.formula(paste0(cv[1],"~",paste(cv2,collapse = "+")))
        cv <- cv2
      }
    }

    # outcome
    # binary/continuous outcome
    if (length(tt)==1){
        is.binary <- length(unique(data[[tt[[1]]]]))==2
        if (is.binary){
            lava::distribution(object,tt[[1]]) <- lava::binomial.lvm()
        }
    }else{
        # event time outcome
        events <- sort(unique(data[[tt[[2]]]]))
        et.formula <- formula(paste0(tt[[1]]," ~ min(",paste(paste0("time.event.",events,"=",events),collapse=", "),")"))
        object <- lava::eventTime(object,et.formula, tt[[2]])
    }
    # call the next method, now the class of object is 'lvm' :)
    synthesize(object=object,data=data,verbose=verbose,logtrans=vv1[vv1!=vv],...)
}


#' @export synthesize.lvm
#' @export
synthesize.lvm <- function(object, data, verbose=FALSE,logtrans = c(),...){

    # note: will be a problem if there are NAs in data, should check at beginning of function
    if(anyNA(data)){stop("There should not be NAs in data.")}

    # Check if all variables in data also occur in object
    if(!all(others <- (names(data) %in% dimnames(object$M)[[1]]))){
        if (verbose)
        warning("Some variables in dataset are not in object (or the names don't match).\n These variables are not synthesized:\n",
                paste0(names(data)[!others],collapse="\n"))
        # if(data.table::is.data.table(object))
        #     data <- data[,dimnames(object$M)[[1]],with=FALSE]
        # else
        #     data <- data[,dimnames(object$M)[[1]],drop=FALSE]
    }

    # dichotomize categorical variables
    # should not do this if the cat. variable is not in lvm.object
    # note: name of dich. variables must correspond with names used in lava::regression formulas outside synthesize-func.

  dichotomized_variables <- c()

  sim_model <- lava::lvm()


  for(var in colnames(data)[grepl('factor|character', sapply(data, class))]){
    if(!any(grepl(var,  dimnames(object$M)[[1]]))){
      warning("should not dichotimize if the cat. variable is not in lvm.object")
    } else  if (length(levels(data[[var]]))==2){

    } else {
      dichotomized_variables <- c(dichotomized_variables, var)
      for(lvl in levels(data[[var]])[-1]){
        data[[paste0(var, lvl)]] <- 1*(data[[var]] == lvl)
        formula <- as.formula(paste0(var, lvl,"~",var))
        lava::transform(sim_model, formula) <- function(x){1*(x==lvl)}
      }
    }
  }



    latent_vars <- endogenous(object)

    # Set distributions:
    # 1. normal
    # 2. binomial
    # 3. categorical
    # deal with exogenous variables (which are without covariates)
    for(var in exogenous(object)){
      var_formula <- as.formula(paste0("~", var))
      #continous
      if("gaussian" %in% attributes(object$attributes$distribution[[var]])$family){
        lava::distribution(sim_model,var_formula) <- lava::normal.lvm(mean=mean(data[[var]]),
                                                                   sd=sd(data[[var]]))
      }
      #binary
      else if("binomial" %in% attributes(object$attributes$distribution[[var]])$family){
        lava::distribution(sim_model, var_formula) <- lava::binomial.lvm(p=mean(factor(data[[var]])==levels(factor(data[[var]]))[2]))
      }
      #categorical
      else if(var %in% names(object$attributes$nordinal)){
        num_cat <- object$attributes$nordinal[[var]] #what if levels is NULL?
        cat_probs <- vector(mode="numeric", length=num_cat-1)

        for(i in 1:(num_cat-1)){
          cat_probs[i] <- mean(data[[var]] == levels(data[[var]])[i+1])
        }
        sim_model <- lava::categorical(sim_model,
                                    var_formula,
                                    labels=levels(data[[var]]),
                                    K=num_cat,
                                    p=cat_probs
        )
      }
      else {
        stop("Error. Distribution not supported!")
      }
    }

    make_dichotomized_string <- function(var){
      s<- ""
      for(lvl in levels(data[[var]])[-1]){
        if (s == ""){
          s <- paste0(var,lvl)
        }
        else {
          s<- paste0(s, " + ", var, lvl)
        }
      }
      return(s)
    }

    # get covariates as string and include dichotomized variables instead of the original ones
    get_covariates <- function(obj, v, dichotomized_variables) {
      covariates <- dimnames(obj$M)[[2]][obj$M[,v] == 1]
      if (length(covariates) == 0){
        return("1")
      }
      else {
        s <- ""
        for (var in covariates){
          if (var %in% dichotomized_variables){
            if (s == ""){
              s <- make_dichotomized_string(var)
            }
            else {
              s <- paste0(s,"+",make_dichotomized_string(var))
            }
          }
          else {
            if (s == ""){
              s <- var
            }
            else {
              s <- paste0(s,"+",var)
            }
          }
        }
        return(s)
      }
    }


    # define latent event time variables
    # note: will be a problem if there are several eventTime variables. Should loop through all these
    has.eventTime <- length(object$attributes$eventHistory)>0
    if (has.eventTime){
        timename <- names(object$attributes$eventHistory)
        covariates <- get_covariates(object,timename,dichotomized_variables)
        response1 <- paste(object$attributes$eventHistory[[timename]]$names, collapse = ",")

        #include event time variables and fit the regression model
        for (latvar in object$attributes$eventHistory[[timename]]$latentTime){
            latvar_formula <- as.formula(paste0("~", latvar))
            status_ind <- object$attributes$eventHistory[[timename]]$events[which(object$attributes$eventHistory[[timename]]$latentTimes %in% latvar)]
            response <- paste0(response1, "==", status_ind)
            surv_formula <- as.formula(paste0("Surv(", response, ")~", covariates))
            G <- survreg(surv_formula, data = data)
            reg_formula <- as.formula(paste0(latvar, "~", covariates))
            lava::distribution(sim_model, latvar_formula) <- lava::coxWeibull.lvm(scale=exp(-G$coef["(Intercept)"]/G$scale),shape=1/G$scale)
            lava::regression(sim_model, reg_formula) <- -coef(G)[-1]/G$scale
        }

        #calculate time and status afterwards
        events <- object$attributes$eventHistory$time$events
        et.formula <- formula(paste0(timename," ~ min(",paste(paste0("time.event.",events,"=",events),collapse=", "),")"))
        sim_model <- lava::eventTime(sim_model,et.formula, object$attributes$eventHistory$time$names[2])
    }

    # Estimate regression coefficients in real data
    # and add them to the lvm object using lava::regression
    '%!in%' <- function(x,y)!('%in%'(x,y))
    for(var in latent_vars[latent_vars %!in% object$attributes$eventHistory$time$names]){
        covariates <- get_covariates(object,var,dichotomized_variables)
        reg_formula <- as.formula(paste0(var, "~", covariates))
        # we have three types of regression to deal with now. Either
        # 1. linear regression
        # 2. logistic regression
        # 3. multinomial logistic regression
        if("gaussian" %in% attributes(object$attributes$distribution[[var]])$family) {
          fit <- lm(reg_formula,data=data)
          lava::distribution(sim_model,as.formula(paste0("~", var))) <- lava::normal.lvm(mean = coef(fit[1]), sd = summary(fit)$sigma)
          lava::regression(sim_model,reg_formula)<-coef(fit)[-1]
        }
        else if ("binomial" %in% attributes(object$attributes$distribution[[var]])$family){
          # correct link function?
          fit <- glm(reg_formula,data=data,family="binomial")
          lava::distribution(sim_model,as.formula(paste0("~", var))) <- lava::binomial.lvm()
          lava::regression(sim_model,reg_formula)<-coef(fit)[-1]
          # does this do the correct thing?
          lava::intercept(sim_model, all.vars(reg_formula)[1])<-coef(fit)[1]
        }
        else if ("categorical" %in% m$attributes$type[[var]]){
            stop("not implemented")
        }
        else if (var %in% logtrans){
          #we don't do anything on the original scale
        }
        else {
          stop("unkown type of regression")
        }
    }
    return(sim_model)
  }

#----------------------------------------------------------------------
### synthesize.R ends here
