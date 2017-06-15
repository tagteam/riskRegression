#' @title Extract i.i.d. decomposition from a Cox model
#' @description Compute the influence function for each observation used to estimate the model
#' @name iid
#' 
#' @param object object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param newdata Optional new data at which to do i.i.d. decomposition 
#' @param tauHazard the vector of times at which the i.i.d decomposition of the baseline hazard will be computed
#' @param keep.times Logical. If \code{TRUE} add the evaluation times to the output.
#' @param store.iid the method used to compute the influence function and the standard error.
#' Can be \code{"full"}, \code{"approx"} or \code{"minimal"}. See the details section.

#' @details If there is no event in a strata, the influence function for the baseline hazard is set to 0.
#'
#' \code{store.iid} equal to \code{"full"} exports the influence function for the coefficients
#' and the baseline hazard at each event time.
#' \code{store.iid} equal to \code{"approx"} does the same except that the terms that do not contributes
#' to the variance are not ignored (i.e. set to 0)
#' \code{store.iid} equal to \code{"minimal"} exports the influence function for the coefficients. For the
#' baseline hazard it only computes the quantities necessary to compute the influence function in order to save memory.
#' 
#' @return A list containing:
#'  \itemize{
#'  \item{IFbeta}{Influence function for the regression coefficient.}
#'  \item{IFhazard}{Time differential of the influence function of the hazard.}
#'  \item{IFcumhazard}{Influence function of the cumulative hazard.}
#'  \item{calcIFhazard}{Elements used to compute the influence function at a given time.}
#'  \item{time}{Times at which the influence function has been evaluated.}
#'  \item{etime1.min}{Time of first event (i.e. jump) in each strata.}
#'  \item{etime.max}{Last observation time (i.e. jump or censoring) in each strata.}
#'  \item{indexObs}{Index of the observation in the original dataset.}
#' }
#'                 
#' @examples
#' library(survival)
#' library(data.table)
#' set.seed(10)
#' d <- sampleData(100, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' 
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
#' system.time(IF.cox <- iidCox(m.cox))
#' system.time(IF.cox_approx <- iidCox(m.cox, store.iid = "approx"))
#'
#' 
#' IF.cox <- iidCox(m.cox, tauHazard = sort(unique(c(7,d$eventtime))))
#'  
#' 

#' @rdname iid
#' @export
iidCox <- function(object, newdata = NULL, tauHazard = NULL, 
                   keep.times = TRUE, store.iid = "full"){

    #### extract elements from object ####
    infoVar <- CoxVariableName(object)
    iInfo <- CoxVarCov(object)
    object.design <- CoxDesign(object)
  
    object.status <- object.design[,"status"]
    object.time <- object.design[,"stop"]
    object.strata <- CoxStrata(object, data = NULL, stratavars = infoVar$stratavars)
    object.levelStrata <- levels(object.strata)
    object.eXb <- exp(CoxLP(object, data = NULL, center = FALSE))
    object.LPdata <- as.matrix(object.design[,infoVar$lpvars,drop = FALSE])
    nStrata <- length(levels(object.strata))
  
    #### Extract new observations ####
    if(!is.null(newdata)){
    
        if("data.frame" %in% class(newdata) == FALSE){
            stop("class of \'newdata\' must inherit from data.frame \n")
        }
    
        # if(infoVar$status %in% names(newdata)){ # call Cox model with with event==1
        tempo <- with(newdata, eval(CoxFormula(object)[[2]]))
        new.status <- tempo[,2]
        new.time <- tempo[,1]
        # }else{ # Cox model from CSC 
        #    new.status <- newdata[[infoVar$status]]
        #    new.time <- newdata[[infoVar$time]]
        # }
    
        new.strata <- CoxStrata(object, data = newdata, 
                                sterms = infoVar$sterms, stratavars = infoVar$stratavars, levels = object.levelStrata, stratalevels = infoVar$stratalevels)
        new.eXb <- exp(CoxLP(object, data = newdata, center = FALSE))
        new.LPdata <- model.matrix(object, newdata)
    
    }else{
    
      new.status <-  object.status
      new.time <-  object.time
      new.strata <-  object.strata
      new.eXb <- object.eXb
      new.LPdata <- object.LPdata
    
    }
  
    #### tests ####
    ## time at which the influence function is evaluated
    if(is.list(tauHazard) && length(tauHazard)!=nStrata){
        stop("argument \"tauHazard\" must be a list with ",nStrata," elements \n",
             "each element being the vector of times for each strata \n")
    }

    if(store.iid %in% c("full","approx","minimal") == FALSE){
        stop("store.iid can only be \"full\", or \"approx\" or \"minimal\"\n")
    }
    

    #### Compute quantities of interest ####
    p <- NCOL(object.LPdata)
  
    ## baseline hazard
    lambda0 <- predictCox(object,
                          type = "hazard",
                          centered = FALSE,
                          keep.strata = TRUE)
    etime1.min <- rep(NA, nStrata)
        
    ## S0, E, jump times
    object.index_strata <- list() 
    object.order_strata <- list()
  
    object.eXb_strata <- list()
    object.LPdata_strata <- list()
    object.status_strata <- list()
    object.time_strata <- list()
  
    new.index_strata <- list()
    new.order_strata <- list()
  
    new.eXb_strata <- list()
    new.LPdata_strata <- list()
    new.status_strata <- list()
    new.time_strata <- list()
  
    Ecpp <- list()
    new.indexJump <- list()
    new.order <- NULL
  
    for(iStrata in 1:nStrata){
    
        ## reorder object data
        object.index_strata[[iStrata]] <- which(object.strata == object.levelStrata[iStrata])
        object.order_strata[[iStrata]] <- order(object.time[object.index_strata[[iStrata]]])
    
        object.eXb_strata[[iStrata]] <- object.eXb[object.index_strata[[iStrata]][object.order_strata[[iStrata]]]]
        object.LPdata_strata[[iStrata]] <- object.LPdata[object.index_strata[[iStrata]][object.order_strata[[iStrata]]],,drop = FALSE]
        object.status_strata[[iStrata]] <- object.status[object.index_strata[[iStrata]][object.order_strata[[iStrata]]]]
        object.time_strata[[iStrata]] <- object.time[object.index_strata[[iStrata]][object.order_strata[[iStrata]]]]
    
        ## reorder new data
        if(!is.null(newdata)){
            new.index_strata[[iStrata]] <- which(new.strata == object.levelStrata[iStrata])
            new.order_strata[[iStrata]] <- order(new.time[new.index_strata[[iStrata]]])
      
            new.eXb_strata[[iStrata]] <- new.eXb[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
            new.LPdata_strata[[iStrata]] <- new.LPdata[new.index_strata[[iStrata]][new.order_strata[[iStrata]]],,drop = FALSE]
            new.status_strata[[iStrata]] <- new.status[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
            new.time_strata[[iStrata]] <- new.time[new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
        }else{
            new.index_strata[[iStrata]] <- object.index_strata[[iStrata]]
            new.order_strata[[iStrata]] <- object.order_strata[[iStrata]]
      
            new.eXb_strata[[iStrata]] <- object.eXb_strata[[iStrata]]
            new.LPdata_strata[[iStrata]] <- object.LPdata_strata[[iStrata]]
            new.status_strata[[iStrata]] <- object.status_strata[[iStrata]]
            new.time_strata[[iStrata]] <- object.time_strata[[iStrata]]
        }
    
        ## E
        Ecpp[[iStrata]] <-  calcE_cpp(status = object.status_strata[[iStrata]], 
                                      eventtime = object.time_strata[[iStrata]],
                                      eXb = object.eXb_strata[[iStrata]],
                                      X = object.LPdata_strata[[iStrata]],
                                      p = p, add0 = TRUE)
    
        new.indexJump[[iStrata]] <- prodlim::sindex(Ecpp[[iStrata]]$Utime1, new.time) - 1
        # if event/censoring is before the first event in the training dataset 
        # then sindex return 0 thus indexJump is -1
        # the following 3 lines convert -1 to 0
        if(any(new.indexJump[[iStrata]]<0)){
            new.indexJump[[iStrata]][new.indexJump[[iStrata]]<0] <- 0
        }
    
        ## store order
        if(length(new.order>0)){
            new.order <- c(new.order, new.index_strata[[iStrata]][new.order_strata[[iStrata]]])
        }else{
            new.order <- new.index_strata[[iStrata]][new.order_strata[[iStrata]]]
        }
    
    }
  
    #### Computation of the influence function ####
    IFbeta <- NULL
    IFcumhazard <- NULL
    IFhazard <- NULL
    calcIFhazard <- list(delta_iS0 = NULL,
                         Elambda0 = NULL,
                         cumElambda0 = NULL,
                         lambda0_iS0= NULL,
                         cumLambda0_iS0= NULL,
                         time1 = NULL)
    ls.Utime1 <- NULL
  
    #### beta
    for(iStrata in 1:nStrata){
    
    new.indexJump_strata <- new.indexJump[[iStrata]][new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
    
    ## IF
    if(p>0){
      if(store.iid != "approx"){
        IFbeta_tempo <- IFbeta_cpp(newT = new.time_strata[[iStrata]],
                                   neweXb = new.eXb_strata[[iStrata]],
                                   newX = new.LPdata_strata[[iStrata]],
                                   newStatus = new.status_strata[[iStrata]], 
                                   newIndexJump = new.indexJump_strata, 
                                   S01 = Ecpp[[iStrata]]$S0,
                                   E1 = Ecpp[[iStrata]]$E,
                                   time1 = Ecpp[[iStrata]]$Utime1,
                                   iInfo = iInfo,
                                   p = p)
      }else{
        IFbeta_tempo <- IFbetaApprox_cpp(newX = new.LPdata_strata[[iStrata]],
                                         newStatus = new.status_strata[[iStrata]],
                                         newIndexJump = new.indexJump_strata,  
                                         E1 = Ecpp[[iStrata]]$E,
                                         iInfo = iInfo,
                                         p = p)
      }
    }else{
      IFbeta_tempo <- matrix(NA, ncol = 1, nrow = length(new.index_strata[[iStrata]]))
    }
    
    ## output
    IFbeta <- rbind(IFbeta, IFbeta_tempo)
    
  }
  
  ## set original order
  IFbeta <- IFbeta[order(new.order),,drop=FALSE]
  
  #### lambda
  for(iStrata in 1:nStrata){
    ## hazard
    if(nStrata==1){ # select only the time,lambda corresponding to the events and not censored observations
      timeStrata <- lambda0$time[lambda0$time %in% Ecpp[[1]]$Utime1]
      lambda0Strata <- lambda0$hazard[lambda0$time %in% Ecpp[[1]]$Utime1]
    }else{ # same within the strata
      index.strata <- which(lambda0$strata == object.levelStrata[iStrata])
      index.keep <- index.strata[lambda0$time[index.strata] %in% Ecpp[[iStrata]]$Utime1]
      
      timeStrata <- lambda0$time[index.keep]
      lambda0Strata <- lambda0$hazard[index.keep]
    }

    etime1.min[iStrata] <- timeStrata[1]
      
    ## tauHazard
    if(is.null(tauHazard)){
      tauHazard_strata <- object.time_strata[[iStrata]][object.status_strata[[iStrata]] == 1]
    }else if(is.list(tauHazard)){
      tauHazard_strata <- tauHazard[[nStrata]]
    }else{
      tauHazard_strata <- tauHazard
    }
    
    ## E
    nUtime1_strata <- length(Ecpp[[iStrata]]$Utime1)
    if(p>0){
      Etempo <- Ecpp[[iStrata]]$E[-NROW(Ecpp[[iStrata]]$E),,drop = FALSE]
    }else{
      Etempo <- matrix(0, ncol = 1, nrow = nUtime1_strata-1)
    }
    
    ## IF
    if(any(new.status_strata[[iStrata]]>0)){
      IFlambda_res <- IFlambda0_cpp(tau = tauHazard_strata,
                                    IFbeta = IFbeta,
                                    newT = new.time, neweXb = new.eXb, newStatus = new.status, newIndexJump = new.indexJump[[iStrata]], newStrata = as.numeric(new.strata),
                                    S01 = Ecpp[[iStrata]]$S0,
                                    E1 = Etempo,
                                    time1 = timeStrata, lastTime1 = Ecpp[[iStrata]]$Utime1[nUtime1_strata], # here lastTime1 will not correspond to timeStrata[length(timeStrata)] when there are censored observations
                                    lambda0 = lambda0Strata,
                                    p = p, strata = iStrata,
                                    exact = (store.iid!="approx"), minimalExport = (store.iid=="minimal")
                                    )      
    }else{
      if(length(tauHazard_strata)==0){tauHazard_strata <- max(object.time_strata[[iStrata]])}
      IFlambda_res <- list(hazard = matrix(0, ncol = length(tauHazard_strata), nrow = NROW(IFbeta)),
                           cumhazard = matrix(0, ncol = length(tauHazard_strata), nrow = NROW(IFbeta))
      )
      if(length(tauHazard_strata)==0){tauHazard_strata <- NA}
    }

      # output 
      ls.Utime1 <- c(ls.Utime1, list(tauHazard_strata))
      if(store.iid=="minimal"){
          calcIFhazard$delta_iS0 <- c(calcIFhazard$delta_iS0, list(IFlambda_res$delta_iS0))
          calcIFhazard$Elambda0 <- c(calcIFhazard$Elambda0, list(IFlambda_res$Elambda0))
          calcIFhazard$cumElambda0 <- c(calcIFhazard$cumElambda0, list(IFlambda_res$cumElambda0))
          calcIFhazard$lambda0_iS0 <- c(calcIFhazard$lambda0_iS0, list(IFlambda_res$lambda0_iS0))
          calcIFhazard$cumLambda0_iS0 <- c(calcIFhazard$cumLambda0_iS0, list(IFlambda_res$cumLambda0_iS0))
          calcIFhazard$time1 <- c(calcIFhazard$time1, list(timeStrata)) # event time by strata
      }else{
          if(keep.times){
              colnames(IFlambda_res$hazard) <- tauHazard_strata
              colnames(IFlambda_res$cumhazard) <- tauHazard_strata
          }
          IFhazard <- c(IFhazard, list(IFlambda_res$hazard))
          IFcumhazard <- c(IFcumhazard, list(IFlambda_res$cumhazard))
      }
  }

    #### export
    return(list(IFbeta = IFbeta,  
                IFhazard = IFhazard,
                IFcumhazard = IFcumhazard,
                calcIFhazard = calcIFhazard,
                time = ls.Utime1,  # time at which the IF is assessed
                etime1.min = etime1.min,
                etime.max = lambda0$lastEventTime,
                indexObs = new.order,
                store.iid = store.iid
                ))
}
