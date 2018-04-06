### discreteRoot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2017 (13:39) 
## Version: 
## Last-Updated: apr  5 2018 (18:08) 
##           By: Brice Ozenne
##     Update #: 170
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * discreteRoot - Documentation
#' @title Dichotomic search for monotone function
#' @description Find the root of a monotone function on a discrete grid of value using dichotomic search
#' @name discreteRoot
#'
#' @param fn [function] objective function to minimize in absolute value.
#' @param grid [vector] possible minimizers.
#' @param increasing [logical] is the function fn increasing?
#' @param check [logical] should the program check that fn takes a different sign for the first vs. the last value of the grid?
#' @param tol [numeric] the absolute convergence tolerance.
#'
#' @examples
#'
#' ### find the position of a value in a vector
#' f <- function(x){abs(vec[x]-1)}
#' discreteRoot(function(x){x},grid = seq(-20,10,1))
#' 

## * discreteRoot
#' @rdname dicreteRoot
#' @export
discreteRoot <- function(fn, grid, increasing = TRUE, check = TRUE,
                         tol = .Machine$double.eps ^ 0.5) {

    n.grid <- length(grid)    
    value.grid <- rep(NA, n.grid)    
    iter <- 1
    ncv <- TRUE
    iSet <- 1:n.grid
    factor <- c(-1,1)[increasing+1]
    
### ** Check
    if(check){
        value.grid[1] <- fn(grid[1])
        value.grid[n.grid] <- fn(grid[n.grid])
        if(sign(value.grid[1])==value.grid[n.grid]){
            list(par = NA,
                 value = NA,
                 counts = 2,
                 cv = 1,
                 message = "Cannot find a solution because the function does not change sign \n")
        }
    }

    
### ** Expore the grid using dichotomic search
    while(iter <= n.grid && ncv && length(iSet)>0){
        iMiddle <- ceiling(length(iSet)/2)
        iIndexInSet <- iSet[iMiddle]
        if(check==FALSE || iIndexInSet %in% c(1,n.grid) == FALSE){
            ## if the current index we are looking at has not already been computed,
            ## then evaluate the objective function.
            ## this is only the case when check is TRUE and we look at the borders
            value.grid[iIndexInSet] <- fn(grid[iIndexInSet])
        }
        if(is.na(value.grid[iIndexInSet])){
            ## handle NA value by just removing the observation from the set of possibilities
            iSet <- setdiff(iSet,iMiddle)
            iter <- iter + 1
        }else if(factor*value.grid[iIndexInSet] > tol){
            ## look in subgrid corresponding to the lowest values (left part)
            iSet <- iSet[setdiff(1:iMiddle,iMiddle)]
            iter <- iter + 1
        }else if(factor*value.grid[iIndexInSet] < -tol){
            ## look in subgrid corresponding to the largest values (right part)
            iN.set <- length(iSet)
            iSet <- iSet[setdiff(iMiddle:iN.set,iMiddle)]
            iter <- iter + 1
        }else{
            ## convergence
            ncv <- FALSE
            solution <- grid[iIndexInSet]
            value <- value.grid[iIndexInSet]
        }
                
    }
    
### ** If did not find a value whose image matched tol, give the closest solution
    if(ncv){
        iIndexInSet <- which.min(abs(value.grid))

       ncv <- FALSE
       solution <- grid[iIndexInSet]
       value <- value.grid[iIndexInSet]
    }

    return(list(par = solution,
                value = value,
                ## grid = setNames(value.grid,grid),
                counts = iter,
                cv = ncv,
                message = NULL))
}

## * boot2pvalue - Documentation
#' @title Compute the p.value from the distribution under H1
#' @description Compute the p.value associated with the estimated statistic
#' using a bootstrap sample of its distribution under H1.
#' 
#' @param x [numeric vector] a vector of bootstrap estimates of the statistic.
#' @param estimate [numeric] the estimated statistic.
#' @param FUN.ci [function] the function used to compute the confidence interval.
#' Must take \code{x}, \code{alternative}, \code{conf.level} and \code{sign.estimate} as arguments
#' and only return the relevant limit (either upper or lower) of the confidence interval.
#' @param alternative [character] a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param tol [numeric] the absolute convergence tolerance.
#' @details
#' For test statistic close to 0, this function returns 1. \cr \cr
#' 
#' For positive test statistic, this function search the quantile alpha such that:
#'\itemize{
#' \item \code{quantile(x, probs = alpha)=0} when the argument alternative is set to \code{"greater"}.
#' \item \code{quantile(x, probs = 0.5*alpha)=0} when the argument alternative is set to \code{"two.sided"}.
#' }
#' If the argument alternative is set to \code{"less"}, it returns 1. \cr \cr
#' 
#' For negative test statistic, this function search the quantile alpha such that:
#' \itemize{
#' \item \code{quantile(x, probs = 1-alpha=0} when the argument alternative is set to \code{"less"}.
#' \item \code{quantile(x, probs = 1-0.5*alpha=0} when the argument alternative is set to \code{"two.sided"}.
#' }
#' If the argument alternative is set to \code{"greater"}, it returns 1.
#' 
#' @examples 
#' set.seed(10)
#' 
#' #### no effect ####
#' x <- rnorm(1e3) 
#' boot2pvalue(x, estimate = mean(x), alternative = "two.sided")
#' ## expected value of 1
#' boot2pvalue(x, estimate = mean(x), alternative = "greater")
#' ## expected value of 0.5
#' boot2pvalue(x, estimate = mean(x), alternative = "less")
#' ## expected value of 0.5
#' 
#' #### positive effect ####
#' x <- rnorm(1e3, mean = 1) 
#' boot2pvalue(x, estimate = 1, alternative = "two.sided")
#' ## expected value of 0.32 = 2*pnorm(q = 0, mean = -1) = 2*mean(x<=0)
#' boot2pvalue(x, estimate = 1, alternative = "greater")  
#' ## expected value of 0.16 = pnorm(q = 0, mean = 1) = mean(x<=0)
#' boot2pvalue(x, estimate = 1, alternative = "less")
#' ## expected value of 0.84 = 1-pnorm(q = 0, mean = 1) = mean(x>=0)
#'
#' #### negative effect ####
#' x <- rnorm(1e3, mean = -1) 
#' boot2pvalue(x, estimate = -1, alternative = "two.sided") 
#' ## expected value of 0.32 = 2*(1-pnorm(q = 0, mean = -1)) = 2*mean(x>=0)
#' boot2pvalue(x, estimate = -1, alternative = "greater")
#' ## expected value of 0.84 = pnorm(q = 0, mean = -1) = mean(x<=0)
#' boot2pvalue(x, estimate = -1, alternative = "less") # pnorm(q = 0, mean = -1)
#' ## expected value of 0.16 = 1-pnorm(q = 0, mean = -1) = mean(x>=0)

## * boot2pvalue
#' @rdname boot2pvalue
#' @export
boot2pvalue <- function(x, estimate = NULL, alternative = "two.sided",
                        FUN.ci = quantileCI,
                        tol = .Machine$double.eps ^ 0.5){ 
  
    x <- na.omit(x)
    n.x <- length(x)
    if(is.null(estimate)){
        estimate <- mean(x)
    }        
    sign.estimate <- estimate>=0
    if(sign(mean(x))!=sign(estimate)){
        warning("the estimate and the average bootstrap estimate do not have same sign \n")
    }

    if(abs(estimate) < tol){ ## too small test statistic
        p.value <- 0
    }else if(n.x < 10){ ## too few bootstrap samples
        p.value <- as.numeric(NA)
    }else if(all(x>0)){ ## clear p.value
        p.value <- switch(alternative,
                          "two.sided" = 0,
                          "less" = 1,
                          "greater" = 0)
    } else if(all(x<0)){ ## clear p.value 
        p.value <- switch(alternative,
                          "two.sided" = 0,
                          "less" = 0,
                          "greater" = 1)
    }else if(all(x==0)){ ## cannot compute quantiles 
        p.value <- switch(alternative,
                          "two.sided" = 0,
                          "less" = 0,
                          "greater" = 0)
    }else{ ## need search to obtain p.value
        ## when the p.value=1-coverage increases, does the quantile increases?
        monotone <- switch(alternative,
                           "two.sided" = sign.estimate,
                           "less" = FALSE,
                           "greater" = TRUE)
        ## grid of confidence level
        grid <- seq(0,by=1/n.x,length.out=n.x)
        ## search for critical confidence level
        
        resSearch <- discreteRoot(fn = function(p.value){
            CI <- do.call(FUN.ci, list(x,
                                       p.value = p.value,
                                       alternative = alternative,
                                       sign.estimate = sign.estimate))
            return(CI[1])
        },
        grid = grid,
        increasing = monotone,
        check = FALSE)

        if(is.na(resSearch$value) || length(resSearch$value)==0 || abs(resSearch$value)>10/n.x){
            warning("incorrect convergence of the algorithm finding the critical quantile \n",
                    "p-value may not be reliable \n")
        }
        p.value <- resSearch$par
    }
  
  return(p.value)
}

## * quantileCI
quantileCI <- function(x, alternative, p.value, sign.estimate, ...){
    probs <- switch(alternative,
                    "two.sided" = c(p.value/2,1-p.value/2)[2-sign.estimate], ## if positive p.value/2 otherwise 1-p.value/2
                    "less" = 1-p.value,
                    "greater" = p.value)

    return(quantile(x, probs = probs)[1])
}


##----------------------------------------------------------------------
### discreteRoot.R ends here
