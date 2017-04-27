### confBandCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 21 2017 (17:40) 
## Version: 
## last-updated: apr 27 2017 (14:54) 
##           By: Brice Ozenne
##     Update #: 21
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Compute quantiles of a gaussian process
#' @description Compute quantiles of a gaussian process
#' 
#' @param iid The iid decomposition of the estimator over time.
#' @param se The variance of the estimate over time.
#' @param times The time points covered by the confidence bands.
#' @param n.object the number of observations in dataset used to fit the model.
#' @param n.new the number of observations for which the quantiles should be computed. 
#' @param n.sim The number of simulations used to compute the quantiles.
#' @param conf.level Level of confidence.
#' 
confBandCox <- function(iid, se, times,
                        n.object, n.new,
                        n.sim = 500, conf.level = 0.95){  
    # NOTE
    # iid must be (n.new,n.times,n.object)
    #  se must be (n.new,n.times)

    new.quantile <- quantileProcess_cpp(nObject = n.object,
                                        nNew = n.new,
                                        nSim = n.sim,
                                        iid = aperm(iid, c(2,3,1)),
                                        se = t(se),
                                        confLevel = conf.level)
    
    return(new.quantile)
}


#----------------------------------------------------------------------
### confBandCox.R ends here
