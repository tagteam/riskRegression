### synthesize.R ---
#----------------------------------------------------------------------
<<<<<<< HEAD
## Author: Vilde Hansteen Ung & Thomas Alexander Gerds
## Created: Apr 28 2021 (09:04)
## Version:
## Last-Updated: Jul 21 2021 (10:10)
=======
## Author: Johan Sebastian Ohlendorff & Vilde Hansteen Ung & Thomas Alexander Gerds
## Created: Apr 28 2021 (09:04)
## Version:
## Last-Updated: Jul 27 2021 (09:07) 
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7
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
<<<<<<< HEAD
                               verbose=TRUE,
                               recursive=FALSE,
                               ...){
    requireNamespace("lava",quietly=FALSE)
    requireNamespace("Publish",quietly=FALSE)
    #requireNamespace("formula.tools",quietly=FALSE)
=======
                               verbose=FALSE,
                               recursive=FALSE,
                               ...){
    requireNamespace("lava",quietly=TRUE)
    # requireNamespace("Publish",quietly=TRUE)
    # requireNamespace("formula.tools",quietly=FALSE)
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7
    # time to event outcome
    tt <- all.vars(update(object,".~1"))
    # covariates
    #vv <- get.vars(formula(delete.response(terms(object))))
    #vv <- all.vars(formula(delete.response(terms(object))))
    vv <- attr(terms(formula(delete.response(terms(object)))),"term.labels")
    vv1 <- all.vars(formula(delete.response(terms(object))))
    s <- ""
    # can be done more nicely
    hasLog <- FALSE
    for (i in 1:length(vv)){
      if (vv[i]==vv1[i]){
        if (s == ""){
          s <- vv[i]
        }
        else {
          s <- paste(s, "+",vv[i])
        }
      }
      else {
        if (s == ""){
          data[[paste0("log",vv1[i])]] <- log(data[[vv1[i]]])
          s <-paste0("log",vv1[i])
          hasLog <- TRUE
        }
        else {
          data[[paste0("log",vv1[i])]] <- log(data[[vv1[i]]])
          s <- paste(s,"+",paste0("log",vv1[i]))
          hasLog <- TRUE
        }
      }
    }
    s <- as.formula(paste0("~",s))
    object <- update(object,s)
    # add all variables
    object <- lava::lvm(object)
    # check if covariates are categorical
    if (hasLog){
      vv <-c(vv1[vv==vv1],paste0("log",vv1[vv!=vv1]))
    }

    is.categorical <- sapply(vv,function(v){
        is.cat <- is.factor(data[[v]])
        if (is.cat==FALSE)
            if (nu <- length(unique(data[[v]]))<max.levels){
                is.cat <- TRUE
                nu <- length(unique(data[[v]]))
            }
            else
                is.cat <- FALSE
        else
            nu <- length(levels(data[[v]]))
        is.bin <- nu==2
        is.bin+is.cat # 0=numeric, 1=categorical with more levels, 2=binary
    })
    for (v in vv[is.categorical==2]){
        # make into factor
        data[[v]] <- factor(data[[v]])
        lava::distribution(object,v) <- lava::binomial.lvm()
    }
    for (v in vv[is.categorical==1]){
        # make into factor
        data[[v]] <- factor(data[[v]])
        #distribution(object,v) <- categorical.lvm() #categorical.lvm not found
        object <- lava::categorical(object, v, K=length(unique(data[[v]])))
    }
    for (v in vv[is.categorical==0]){
        lava::distribution(object,v) <- lava::normal.lvm()
    }

    if (recursive){
      #don't change original covariates
      cv <- vv
      while (length(cv) > 1) {
        cv2 <- cv[2:length(cv)]

        if (!is.numeric(pbc[[cv[1]]])){
          #lacks implementation in case it is categorical
          stop("not implemented")
        }
        regression(object) <- as.formula(paste0(paste0(cv[1],"~"),paste(cv2,collapse = "+")))
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
<<<<<<< HEAD
=======
    
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7
    # call the next method, now the class of object is 'lvm' :)
    synthesize(object=object,data=data,verbose=verbose,...)
}

<<<<<<< HEAD
synthesize.lvm <- function(object, data, verbose,...){
=======
#' @export synthesize.lvm
#' @export
synthesize.lvm <- function(object, data, verbose=FALSE,...){
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7

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
    for(var in colnames(data)[grepl('factor|character', sapply(data, class))]){
        if(!any(grepl(var,  dimnames(object$M)[[1]]))){

        } else  if (length(levels(data[[var]]))==2){

        } else {
            for(lvl in levels(data[[var]])[-1]){
                data[[paste0(var, lvl)]] <- 1*(data[[var]] == lvl)
                formula <- as.formula(paste0(var, lvl,"~",var))
                lava::transform(object, formula) <- function(x){1*(x==lvl)}
            }
        }
    }
<<<<<<< HEAD

=======
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7
    # define latent event time variables
    # note: will be a problem if there are several eventTime variables. Should loop through all these
    has.eventTime <- length(object$attributes$eventHistory)>0
    if (has.eventTime){
        timename <- names(object$attributes$eventHistory)
        for (latvar in object$attributes$eventHistory[[timename]]$latentTime){
            latvar_formula <- as.formula(paste0("~", latvar))
            lava::distribution(object, latvar_formula) <- lava::coxWeibull.lvm()
        }
    }

<<<<<<< HEAD
  # Estimate coefficients
  for(var in dimnames(object$M)[[1]]){
    var_formula <- as.formula(paste0("~", var))
=======
    # Set distributions:
    # 1. normal
    # 2. binomial
    # 3. weibull
    # 4. 
    for(var in dimnames(object$M)[[1]]){
        var_formula <- as.formula(paste0("~", var))
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7

    if("gaussian" %in% attributes(object$attributes$distribution[[var]])$family){
        lava::distribution(object,var_formula) <- lava::normal.lvm(mean=mean(data[[var]]),
                                                                   sd=sd(data[[var]]))
    }
    else if("binomial" %in% attributes(object$attributes$distribution[[var]])$family){
        lava::distribution(object, var_formula) <- lava::binomial.lvm(p=mean(factor(data[[var]])==levels(factor(data[[var]]))[2]))
    }
    else if(has.eventTime && "weibull" %in% attributes(object$attributes$distribution[[var]])$family){
        covariates <- dimnames(object$M)[[2]][object$M[,var] == 1]
        response <- paste(object$attributes$eventHistory[[timename]]$names, collapse = ",")
        status_ind <- object$attributes$eventHistory[[timename]]$events[which(object$attributes$eventHistory[[timename]]$latentTimes %in% var)]
        response <- paste0(response, "==", status_ind)
        #timename <- names(object$attributes$eventHistory)
        if(length(covariates) == 0){
            surv_formula <- as.formula(paste0("Surv(", response,")", "~1"))
        } else{
            covariates <- paste(covariates, collapse = "+")
            surv_formula <- as.formula(paste0("Surv(", response,")", "~", covariates))
        }
        G <- survreg(surv_formula, data = data)
        lava::distribution(object, var_formula) <- lava::coxWeibull.lvm(scale=exp(-G$coef["(Intercept)"]/G$scale),shape=1/G$scale)
    }
    else if(var %in% names(object$attributes$nordinal)){
<<<<<<< HEAD
        num_cat <- object$attributes$nordinal[[var]] #what is levels is NULL?
        cat_probs <- vector(mode="numeric", length=num_cat-1)

        for(i in 1:(num_cat-1)){
          cat_probs[i] <- mean(data[[var]] == levels(data[[var]])[i+1])
=======
        num_cat <- object$attributes$nordinal[[var]] #what if levels is NULL?
        cat_probs <- vector(mode="numeric", length=num_cat-1)

        for(i in 1:(num_cat-1)){
            cat_probs[i] <- mean(data[[var]] == levels(data[[var]])[i+1])
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7
        }
        object <- lava::categorical(object,
                                    var_formula,
                                    labels=levels(data[[var]]),
                                    K=num_cat,
                                    p=cat_probs
<<<<<<< HEAD
        )
    }
  }
  #find lava::regression relationships
  for(var in dimnames(object$M)[[1]]){
    covariates <- dimnames(object$M)[[2]][object$M[,var] == 1]

    if (has.eventTime && length(dimnames(object$M)[[2]][object$M[,timename] == 1]) != 0 && length(covariates)==0 && "weibull" %in% attributes(object$attributes$distribution[[var]])$family){
        covariates <- dimnames(object$M)[[2]][object$M[,timename] == 1]
        covariates <- paste(covariates, collapse = "+")
        response <- paste(object$attributes$eventHistory[[timename]]$names, collapse = ",")
        status_ind <- object$attributes$eventHistory[[timename]]$events[which(object$attributes$eventHistory[[timename]]$latentTimes %in% var)]
        response <- paste0(response, "==", status_ind)
        surv_formula <- as.formula(paste0("Surv(", response, ")~", covariates))
        G <- survreg(surv_formula, data = data)
        reg_formula <- as.formula(paste0(var, "~", covariates))
        lava::regression(object, reg_formula) <- -coef(G)[-1]/G$scale
    }
    else if (length(covariates) > 0) {
      covariates <- paste(covariates, collapse = "+")
=======
                                    )
    }
    }
    # Estimate regression coefficients in real data 
    # and add them to the lvm object using lava::regression
    for(var in dimnames(object$M)[[1]]){
        covariates <- dimnames(object$M)[[2]][object$M[,var] == 1]
        if (has.eventTime && length(dimnames(object$M)[[2]][object$M[,timename] == 1]) != 0
            && length(covariates)==0
            && "weibull" %in% attributes(object$attributes$distribution[[var]])$family){
            covariates <- dimnames(object$M)[[2]][object$M[,timename] == 1]
            covariates <- paste(covariates, collapse = "+")
            response <- paste(object$attributes$eventHistory[[timename]]$names, collapse = ",")
            status_ind <- object$attributes$eventHistory[[timename]]$events[which(object$attributes$eventHistory[[timename]]$latentTimes %in% var)]
            response <- paste0(response, "==", status_ind)
            surv_formula <- as.formula(paste0("Surv(", response, ")~", covariates))
            G <- survreg(surv_formula, data = data)
            reg_formula <- as.formula(paste0(var, "~", covariates))
            lava::regression(object, reg_formula) <- -coef(G)[-1]/G$scale
        }
        else if (length(covariates) > 0) {
            covariates <- paste(covariates, collapse = "+")
>>>>>>> ce0adfebaddb41edfd567da8562d05cdea96f1d7

      reg_formula <- as.formula(paste0(var, "~", covariates))

      # lava::regression depends on response variable type
      # if variable has a specified distribution, do:
      if(var %in% names(object$attributes$distribution)){
          if("gaussian" %in% attributes(object$attributes$distribution[[var]])$family) {
              lava::regression(object,reg_formula)<-coef(lm(reg_formula,data=data))[-1]
          } else if ("binomial" %in% attributes(object$attributes$distribution[[var]])$family){
              #browser()
              lava::regression(object,reg_formula)<-coef(glm(reg_formula,data=data,family="binomial"))[-1]
              lava::intercept(object, all.vars(reg_formula)[1])<-coef(glm(reg_formula,data=data,family="binomial"))[1]
          } else if ("weibull" %in% attributes(object$attributes$distribution[[var]])$family)
          {
              if(length(covariates) == 0){
                  #reg_formula <- as.formula(paste0(var,"~", "1"))
                  #is anything needed here?
              }
              else{
                  if (has.eventTime){
                      response <- paste(object$attributes$eventHistory[[timename]]$names, collapse = ",")
                      status_ind <- object$attributes$eventHistory[[timename]]$events[which(object$attributes$eventHistory[[timename]]$latentTimes %in% var)]
                      response <- paste0(response, "==", status_ind)
                      surv_formula <- as.formula(paste0("Surv(", response, ")~", covariates))
                      G <- survreg(surv_formula, data = data)
                      reg_formula <- as.formula(paste0(var, "~", covariates))
                      lava::regression(object, reg_formula) <- -coef(G)[-1]/G$scale
                  }
              }
          }
          # if variable does not have a specified distribution, then:
      } else if (!(class(data[[var]]) == "numeric") & length(levels(factor(data[[var]]))) == 2) #variable is binary
      {
          lava::regression(object,reg_formula) <- coef(glm(reg_formula,data=data,family="binomial"))[-1]
      } else if (class(data[[var]]) == "numeric") #variable is continous
      {
          lava::regression(object,reg_formula) <- coef(lm(reg_formula,data=data))[-1]
      }
      #what if variable is categorical?
    }
  }
    return(object)
  }


#----------------------------------------------------------------------
### synthesize.R ends here
