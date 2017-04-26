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
#' @param exact Logical. If \code{TRUE} then the exact influence function is computed.
#' Otherwise the influence function is estimated using the counting processes instead of the estimated martingales.
#' 
#' @details If there is no event in a strata, the influence function for the baseline hazard is set to 0.
#'
#' @examples
#' library(survival)
#' library(data.table)
#' set.seed(10)
#' d <- sampleData(100, outcome = "survival")[,.(eventtime,event,X1,X6)]
#' setkey(d, eventtime)
#' 
#' m.cox <- coxph(Surv(eventtime, event) ~ X1+X6, data = d, y = TRUE, x = TRUE)
#' system.time(IC.cox <- iidCox(m.cox))
#' system.time(IC.cox_approx <- iidCox(m.cox, exact = FALSE))
#'
#' 
#' IC.cox <- iidCox(m.cox, tauHazard = sort(unique(c(7,d$eventtime))))
#'  
#' 

#' @rdname iid
#' @export
iidCox <- function(object, newdata = NULL, tauHazard = NULL, 
                   keep.times = TRUE, exact = TRUE){
  
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
  
  #### Compute quantities of interest ####
  p <- NCOL(object.LPdata)
  
    ## baseline hazard
    lambda0 <- predictCox(object,
                          type = "hazard",
                          centered = FALSE,
                          keep.strata = TRUE)
  
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
  ICbeta <- NULL
  ICcumhazard <- NULL
  IChazard <- NULL
  ls.Utime1 <- NULL
  
  #### beta
  for(iStrata in 1:nStrata){
    
    new.indexJump_strata <- new.indexJump[[iStrata]][new.index_strata[[iStrata]][new.order_strata[[iStrata]]]]
    
      ## IF
      if(p>0){
          if(exact){
              ICbeta_tempo <- ICbeta_cpp(newT = new.time_strata[[iStrata]],
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
              ICbeta_tempo <- ICbetaApprox_cpp(newX = new.LPdata_strata[[iStrata]],
                                               newStatus = new.status_strata[[iStrata]],
                                               newIndexJump = new.indexJump_strata,  
                                               E1 = Ecpp[[iStrata]]$E,
                                               iInfo = iInfo,
                                               p = p)
          }
      }else{
          ICbeta_tempo <- matrix(NA, ncol = 1, nrow = length(new.index_strata[[iStrata]]))
      }
    
    ## output
    ICbeta <- rbind(ICbeta, ICbeta_tempo)
    
  }
  
  ## set original order
  ICbeta <- ICbeta[order(new.order),,drop=FALSE]
  
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
        
          IClambda_res <- IClambda0_cpp(tau = tauHazard_strata,
                                        ICbeta = ICbeta,
                                        newT = new.time, neweXb = new.eXb, newStatus = new.status, newIndexJump = new.indexJump[[iStrata]], newStrata = as.numeric(new.strata),
                                        S01 = Ecpp[[iStrata]]$S0,
                                        E1 = Etempo,
                                        time1 = timeStrata, lastTime1 = Ecpp[[iStrata]]$Utime1[nUtime1_strata], # here lastTime1 will not correspond to timeStrata[length(timeStrata)] when there are censored observations
                                        lambda0 = lambda0Strata,
                                        p = p, strata = iStrata, exact = exact)
      
      }else{
          if(length(tauHazard_strata)==0){tauHazard_strata <- max(object.time_strata[[iStrata]])}
          IClambda_res <- list(hazard = matrix(0, ncol = length(tauHazard_strata), nrow = NROW(ICbeta)),
                               cumhazard = matrix(0, ncol = length(tauHazard_strata), nrow = NROW(ICbeta))
                               )
          if(length(tauHazard_strata)==0){tauHazard_strata <- NA}
      }
    
        # output 
        ls.Utime1 <- c(ls.Utime1, list(tauHazard_strata))
        if(keep.times){
            colnames(IClambda_res$hazard) <- tauHazard_strata
            colnames(IClambda_res$cumhazard) <- tauHazard_strata
        }
        IChazard <- c(IChazard, list(IClambda_res$hazard))
        ICcumhazard <- c(ICcumhazard, list(IClambda_res$cumhazard))
    }
  
  #### export
  return(list(ICbeta = ICbeta,  
              IChazard = IChazard,
              ICcumhazard = ICcumhazard,
              time = ls.Utime1,
              indexObs = new.order
  ))
}
