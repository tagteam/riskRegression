#' @title Compute the average treatment effect using CSC
#'
#' Use the g-formula to estimate the average treatment effect
#' @param data dataset
#' @param treatVar name of the column containing the treatment variable
#' @param treatLevels the levels of treatment to be compared
#' @param ls.formula formula argument of CSC
#' @param times time argument of CSC
#' @param cause argument of CSC
#' @param return.model should the fitted CSC model be returned
#' @param addTime should the time be specified in the output
#' @param n.bootstrap the number of bootstrap replications used to compute the confidence interval
#' @param ncpus the number of CPU to use for the bootstrap
#' @param seed the seed to set at the beginning of the bootstrap
#' @param alpha the type 1 error
#' @param ... passed on
#' 
#' @details 
#' Require the snowfall package for the bootstrap.
#' 
#' @return A list with:
#' model: the CSC model
#' averageCIF: the predicted CIF averaged over patient for a given treatment modality (column) and time (row)
#' gFormula: A list (times) containing the value of the g formula comparing all pair of treatment modality, with the first one as a reference.
#' averageCIFCIsup: the lower bound of the CI for averageCIF
#' averageCIFCIinf: the upper bound of the CI for averageCIF
#' gFormulaCIsup: the lower bound of the CI for averageCIFCIsup
#' gFormulaCIinf: the upper bound of the CI for averageCIFCIinf
#' 
#' @examples 
#' #library(snowfall)
#' data <- SimCompRisk(1e3)
#' data$time <- round(data$time,1)
#' 
#' seqtimes <- sample(x = unique(sort(data$time)), size = 100) 
#' 
#' data$X1 <- factor(rbinom(1e3, prob = c(0.2,0.3,0.2) , size = 3), labels = paste0("T",0:3))
#' res.Gformula <- ateCSC(data = data, treatVar = "X1", treatLevels = NULL,
#'                          ls.formula = Hist(time,event)~ X1*X2, 
#'                          times = 1:10, cause = 1,
#'                          n.bootstrap = 10)
#'  paste0(res.Gformula$averageCIF[res.Gformula$averageCIF[,"time"] == 5,"T0"],
#'        " [",
#'         res.Gformula$averageCIFCIinf[res.Gformula$averageCIFCIinf[,"time"] == 5,"T0"],
#'         " ; ",
#'         res.Gformula$averageCIFCIsup[res.Gformula$averageCIFCIsup[,"time"] == 5,"T0"],
#'         "]")
#' @export
#' 
ateCSC <- function(data, treatVar, treatLevels = NULL,
                   ls.formula, times, cause, return.model = FALSE, addTime = TRUE,
                   n.bootstrap = 0, ncpus = 1, seed = NULL, alpha = 0.05, ...){
  
  #### Prepare
  if(treatVar %in% names(data) == FALSE){
    stop("Gformula: the column ",treatVar," has not been found in data \n",
         "names(data): ",paste(names(data), collapse = " "),"\n")
  }
  if(is.null(treatLevels)){
    treatLevels <- sort(unique(data[[treatVar]]))
  }
  test.factor <- is.factor(data[[treatVar]])
  n.treatLevels <- length(treatLevels)
  n.times <- length(times)
  n.obs <- nrow(data)
  
  #### calc G formula
  warperGformula <- function(data, treatVar, treatLevels, ls.formula, times, cause, return.model, ...){
    
    ## fit
    CSC.fit <- riskRegression::CSC(formula = ls.formula, data = data, cause = cause, ...)  
    
    ## prediction
    resPrediction <- matrix(NA, nrow = n.times, ncol = n.treatLevels, dimnames = list(NULL, treatLevels))
    for(iterLevel in  1:n.treatLevels){
      dataTempo <- data
      
      if(test.factor){
        dataTempo[[treatVar]] <- factor(treatLevels[iterLevel], levels = treatLevels)
      }else{
        dataTempo[[treatVar]] <- treatLevels[iterLevel]
      }
      
      resPrediction[,iterLevel]<- colMeans(predictEvent(CSC.fit, newdata = dataTempo, times = times, cause = cause))
    }
    
    ## postprocessing
    #Alldiff <- lapply(1:n.times, function(x){ dist(as.numeric(resPrediction[x,])) })
    ## alldiff <- simplify2array(lapply(1:n.times, function(x){ 
      ## res <- as.matrix(dist(as.numeric(resPrediction[x,])))
      ## dimnames(res)  <- list(treatLevels, treatLevels)
      ## return(res)
    ## }))
    
    ## output
    return(list(model = if(return.model){CSC.fit}else{NULL},
                averageCIF = resPrediction,
                gFormula = alldiff) 
    )
    
  }
  
  #### Punctual estimate
  res.export <- warperGformula(data, treatVar, treatLevels, ls.formula, times, cause, return.model, ...)
  
  #### Bootstrap
  if(n.bootstrap>0){
    if(!is.null(seed)){set.seed(seed)}
    
    if(ncpus > parallel::detectCores()){
        warning("Gformula: not enough available cores \n",
                "available: ",parallel::detectCores()," | requested: ",ncpus,"\n")
    }
    
    snowfall::sfInit(parallel = TRUE, cpus = ncpus)
    VarnamesToLoad <- c("n.obs", "data",  "treatVar", "treatLevels", "ls.formula", "times", "cause", "return.model", names(list(...)))
    snowfall::sfExport(list = as.list(VarnamesToLoad))
    LibraryToLoad <- c("riskRegression")
    sapply(LibraryToLoad, function(x){eval(call("sfLibrary", x, character.only = TRUE, verbose = FALSE))})
    
    res.boot <- snowfall::sfClusterApplyLB(1:n.bootstrap, function(iterN){
      dataBoot <- data[sample(1:n.obs, size = n.obs, replace = TRUE),]
      warperGformula(dataBoot, treatVar, treatLevels, ls.formula, times, cause, return.model = FALSE, ...)
    })
    snowfall::sfStop()
    
    ## confidence intervals
    res.export$averageCIFCIsup <- apply(simplify2array(lapply(res.boot, function(x){x[["averageCIF"]]})),
                                         1:2,
                                         quantile, prob = c(1-alpha/2)
    )
    res.export$averageCIFCIinf <- apply(simplify2array(lapply(res.boot, function(x){x[["averageCIF"]]})),
                                         1:2,
                                         quantile, prob = c(alpha/2)
    )
    res.export$gFormulaCIsup <- apply(simplify2array(lapply(res.boot, function(x){x[["gFormula"]]})),
                                       1:3,
                                       quantile, prob = c(1-alpha/2)
    )
    res.export$gFormulaCIinf <- apply(simplify2array(lapply(res.boot, function(x){x[["gFormula"]]})),
                                       1:3,
                                       quantile, prob = c(alpha/2)
    )
    
  }
  
  #### add time
  if(addTime){
    res.export$averageCIF <- cbind(time = times, res.export$averageCIF)
    res.export$averageCIFCIsup <- cbind(time = times, res.export$averageCIFCIsup)
    res.export$averageCIFCIinf <- cbind(time = times, res.export$averageCIFCIinf)
    
    dimnames(res.export$gFormula)[[3]] <- times
    dimnames(res.export$gFormulaCIsup)[[3]] <- times
    dimnames(res.export$gFormulaCIinf)[[3]] <- times
  }
  
  #### Export 
  return(res.export)
  
}
