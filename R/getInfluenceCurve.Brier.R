### getInfluenceCurve.Brier.R --- 
#----------------------------------------------------------------------
## Author: Thomas Alexander Gerds
## Created: Jun  4 2024 (11:46) 
## Version: 
## Last-Updated: Jun  4 2024 (11:46) 
##           By: Thomas Alexander Gerds
##     Update #: 1
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getInfluenceCurve.Brier <- function(t,
                                    time,
                                    IC0,
                                    residuals,
                                    IC.G,
                                    cens.model,
                                    nth.times=NULL,
                                    conservative,
                                    event){
  if (is.unsorted(time)){
    stop("Internal error. Time is not sorted in ascending order. ")
  }

  ##
  ## Compute influence function of Brier score estimator using weights of the reverse Cox model
  ## This function evaluates the part of influence function which is related to the IPCW weights
  ## The other part is IC0.
  ##
  ## \frac{1}{n}\sum_{i=1}^n
  ## m_{t,n}^{(1)}(X_i)
  ## [\frac{I_{T_i\leq t}\Delta_i}{G^2(T_i\vert Z_i)}IF_G(T_i,X_k; X_i)+\frac{I_{T_i>t}}{G^2(t|Z_i)}IF_G(t,X_k; X_i)]
  ## with
  ## IF_G(t,X_k; X_i)=-\exp(-\Lambda(t\vert Z_i))IF_{\Lambda}(t,X_k; X_i)
  ##
  ## IC_G(t,z;x_k) is an array with dimension (nlearn=N, gtimes, newdata)
  ## where gtimes = subject.times (Weights$IC$IC.subject) or times (Weights$IC$IC.times)
  ## and subject.times=Y[(((Y<=max(times))*status)==1)]
  ##
  ## don't square the weights because they will be multiplied with the
  ## residuals that are already weighted
  ##
  N <- length(residuals)
  if (conservative[[1]] || cens.model[[1]] == "none"){
    IC0
  }
  else if (cens.model[[1]]=="cox") {## Cox regression
    if (!IC.G$censoring.save.memory){
      ic.weights <- IC.G[[2]][[nth.times]]
      IF.Brier <- IC0+as.numeric(rowSumsCrossprod(as.matrix(residuals),ic.weights,0)) / N
    }
    else {
      wdata <- IC.G[[3]]
      fit <- IC.G[[2]]
      TiMinus <- IC.G[[4]]
      res1 <- (time <= t & event != 0)*residuals
      res2 <- (time > t)*residuals
      IF.Brier <- IC0 + (predictCoxWeights(fit,newdata = wdata,times = TiMinus,diag=TRUE,weights=res1,isBeforeTau = TRUE, tau=t) + predictCoxWeights(fit,newdata = wdata,times = t,diag=FALSE,weights=res2,isBeforeTau = FALSE, tau=t))/N
    }
    IF.Brier
  }else if (cens.model[[1]] == "KaplanMeier"){
    IC0 + getInfluenceFunctionBrierKMCensoringTerm(t,time,residuals,event)
  }
  else {
    stop("Non conservative options with cens.model not being a Cox model or KaplanMeier are not implemented. ")
  }
}


######################################################################
### getInfluenceCurve.Brier.R ends here
