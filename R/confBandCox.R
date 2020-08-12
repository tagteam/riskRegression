### confBandCox.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: apr 21 2017 (17:40) 
## Version: 
## last-updated: aug 11 2020 (12:51) 
##           By: Brice Ozenne
##     Update #: 37
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
#' @param n.sim The number of simulations used to compute the quantiles.
#' @param conf.level Level of confidence.
#' 
confBandCox <- function(iid, se, n.sim, conf.level){    
    # NOTE
    # iid must be (n.object,n.times,n.new)
    #  se must be (n.new,n.times)
    dimTempo <- dim(iid)
    new.quantile <- quantileProcess_cpp(nObject = dimTempo[1],
                                        nNew = dimTempo[3],
                                        nSim = n.sim,
                                        iid = aperm(iid, c(2,1,3)),
                                        se = t(se),
                                        confLevel = conf.level)
    
    return(new.quantile)
}


#----------------------------------------------------------------------
### confBandCox.R ends here
