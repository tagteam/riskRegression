### discreteRoot.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: nov 22 2017 (13:39) 
## Version: 
## Last-Updated: Jan 29 2019 (11:13) 
##           By: Thomas Alexander Gerds
##     Update #: 242
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
            return(list(par = NA,
                        value = NA,
                        counts = 0,
                        cv = 1,
                        message = "Cannot find a solution because the function does not change sign \n"))
        }
        if(increasing[[1]] && value.grid[[1]] > value.grid[[n.grid]]){
            return(list(par = NA,
                        value = NA,
                        counts = 0,
                        cv = 1,
                        message = "Cannot find a solution - argument \'increasing\' does not match the variations of the functions \n"))
        }
        if(!increasing[[1]] && value.grid[[1]] < value.grid[[n.grid]]){
            return(list(par = NA,
                        value = NA,
                        counts = 0,
                        cv = 1,
                        message = "Cannot find a solution - argument \'increasing\' does not match the variations of the functions \n"))
        }
    }

    
### ** Expore the grid using dichotomic search
    while(iter[[1]] <= n.grid[[1]] && ncv[[1]]==TRUE && length(iSet)>0){
        iMiddle <- ceiling(length(iSet)/2)
        iIndexInSet <- iSet[iMiddle]
        if(check[[1]]==FALSE || iIndexInSet %in% c(1,n.grid) == FALSE){
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
#' @param null [numeric] value of the statistic under the null hypothesis.
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
#' boot2pvalue(x, null = 0, estimate = mean(x), alternative = "two.sided")
#' ## expected value of 1
#' boot2pvalue(x, null = 0, estimate = mean(x), alternative = "greater")
#' ## expected value of 0.5
#' boot2pvalue(x, null = 0, estimate = mean(x), alternative = "less")
#' ## expected value of 0.5
#' 
#' #### positive effect ####
#' x <- rnorm(1e3, mean = 1) 
#' boot2pvalue(x, null = 0, estimate = 1, alternative = "two.sided")
#' ## expected value of 0.32 = 2*pnorm(q = 0, mean = -1) = 2*mean(x<=0)
#' boot2pvalue(x, null = 0, estimate = 1, alternative = "greater")  
#' ## expected value of 0.16 = pnorm(q = 0, mean = 1) = mean(x<=0)
#' boot2pvalue(x, null = 0, estimate = 1, alternative = "less")
#' ## expected value of 0.84 = 1-pnorm(q = 0, mean = 1) = mean(x>=0)
#'
#' #### negative effect ####
#' x <- rnorm(1e3, mean = -1) 
#' boot2pvalue(x, null = 0, estimate = -1, alternative = "two.sided") 
#' ## expected value of 0.32 = 2*(1-pnorm(q = 0, mean = -1)) = 2*mean(x>=0)
#' boot2pvalue(x, null = 0, estimate = -1, alternative = "greater")
#' ## expected value of 0.84 = pnorm(q = 0, mean = -1) = mean(x<=0)
#' boot2pvalue(x, null = 0, estimate = -1, alternative = "less") # pnorm(q = 0, mean = -1)
#' ## expected value of 0.16 = 1-pnorm(q = 0, mean = -1) = mean(x>=0)

## * boot2pvalue
#' @rdname boot2pvalue
#' @export
boot2pvalue <- function(x, null, estimate = NULL, alternative = "two.sided",
                        FUN.ci = quantileCI,
                        tol = .Machine$double.eps ^ 0.5){ 
  
    x.boot <- na.omit(x)
    n.boot <- length(x.boot)
    statistic.boot <- mean(x.boot) - null
    if(is.null(estimate)){
        statistic <- statistic.boot
    }else{
        statistic <- estimate - null
        if(sign(statistic.boot)!=sign(statistic)){
            warning("the estimate and the average bootstrap estimate do not have same sign \n")
        }
    }
    sign.statistic <- statistic>=0

    if(abs(statistic) < tol){ ## too small test statistic
        p.value <- 1
    }else if(n.boot < 10){ ## too few bootstrap samples
        p.value <- as.numeric(NA)
    }else if(all(x.boot>null)){ ## clear p.value
        p.value <- switch(alternative,
                          "two.sided" = 0,
                          "less" = 1,
                          "greater" = 0)
    } else if(all(x.boot<null)){ ## clear p.value 
        p.value <- switch(alternative,
                          "two.sided" = 0,
                          "less" = 0,
                          "greater" = 1)
    }else{ ## need search to obtain p.value
        ## when the p.value=1-coverage increases, does the quantile increases?
        increasing <- switch(alternative,
                             "two.sided" = sign.statistic,
                             "less" = FALSE,
                             "greater" = TRUE)
        ## grid of confidence level
        grid <- seq(0,by=1/n.boot,length.out=n.boot)
        
        ## search for critical confidence level
        resSearch <- discreteRoot(fn = function(p.value){
            CI <- FUN.ci(x = x.boot,
                         p.value = p.value,
                         alternative = alternative,
                         sign.estimate = sign.statistic)
            return(CI[1]-null)
        },
        grid = grid,
        increasing = increasing,
        check = FALSE)

        ## check change sign
        sign.before <- sign(FUN.ci(x = x.boot,
                                   p.value = max(0,resSearch$par-1/n.boot),
                                   alternative = alternative,
                                   sign.estimate = sign.statistic)-null)

        sign.after <- sign(FUN.ci(x = x.boot,
                                  p.value = min(1,resSearch$par+1/n.boot),
                                  alternative = alternative,
                                  sign.estimate = sign.statistic)-null)

        ##
        if(is.na(resSearch$value[[1]]) || length(resSearch$value)==0 || resSearch$par[[1]]<0 || resSearch$par[[1]]>1 || sign.before[[1]]==sign.after[[1]]){
            warning("incorrect convergence of the algorithm finding the critical quantile \n",
                    "p-value may not be reliable \n")

        }
        p.value <- resSearch$par
    }

    if(p.value %in% c(0,1)){
        message("Estimated p-value of ",p.value," - consider increasing the number of boostrap samples \n")
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
