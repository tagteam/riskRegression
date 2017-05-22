### influence2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 19 2017 (16:01) 
## Version: 
## last-updated: maj 19 2017 (16:59) 
##           By: Brice Ozenne
##     Update #: 6
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Influence test [Experimental!!]
#' @description Compare two estimates using their influence function
#' 
#' @param predict1 a list containing an estimate and its influence function 
#' @param predict2 same as predict1 but for another model
#' @param type the type of estimate
#' 
#' @examples
#' \dontrun{
#' warperNoStrata <- function(n, ...){
#'   newdata <- data.frame(X1=0,X6=0)
#'   time <- 1
#'   
#' m <- lvm()
#' regression(m) <- y ~ X6 + X1
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
#' distribution(m,~y) <- a <- coxWeibull.lvm(scale=0.3)
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#' setkey(d, eventtime)
#'   
#' m.cox <- coxph(Surv(eventtime, status) ~ X1+X6, 
#'                data = d, y = TRUE, x = TRUE)
#' survNoStrata <- predictCox(m.cox, type = "survival",
#'                            newdata = newdata, times = time, iid = TRUE, se = TRUE)
#'   
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1)+X6, 
#'                      data = d, y = TRUE, x = TRUE)
#' survStrata <- predictCox(mStrata.cox, type = "survival", 
#'                          newdata = newdata, times = time, iid = TRUE, se = TRUE)
#'   
#' res <- influenceTest(survStrata, survNoStrata, type = "survival")
#' return(res)
#' }
#' 
#' warperStrata <- function(n, ...){
#'   newdata <- data.frame(X1=0,X6=0)
#' time <- 1
#'   
#' m <- lvm()
#' regression(m) <- y ~ X6
#' regression(m) <- s ~ exp(-2*X1)
#' distribution(m,~X1) <- binomial.lvm()
#' distribution(m,~cens) <- coxWeibull.lvm(scale=1)
#' distribution(m,~y) <- a <- coxWeibull.lvm(scale=0.3,shape=~s)
#' eventTime(m) <- eventtime ~ min(y=1,cens=0)
#' d <- as.data.table(sim(m,n))
#'   setkey(d, eventtime)
#'   
#' m.cox <- coxph(Surv(eventtime, status) ~ X1+X6, data = d, y = TRUE, x = TRUE)
#' survNoStrata <- predictCox(m.cox, type = "survival", 
#'                            newdata = newdata, times = time, iid = TRUE, se = TRUE)
#'   
#' mStrata.cox <- coxph(Surv(eventtime, status) ~ strata(X1)+X6, data = d, y = TRUE, x = TRUE)
#' survStrata <- predictCox(mStrata.cox, type = "survival", 
#'                          newdata = newdata, times = time, iid = TRUE, se = TRUE)
#' res <- influenceTest(survStrata, survNoStrata, type = "survival")
#' return(res)
#' }
#' 
#' n <- 500
#' resNoStrata <- pbsapply(1:200, warperNoStrata, n = n)
#' 
#' sd(unlist(resNoStrata["delta",]))
#' quantile(unlist(resNoStrata["se.delta",]))
#' 
#' mean(unlist(resNoStrata["p.value",])<0.05)
#' 
#' resStrata <- pbsapply(1:100, warperStrata, n = n)
#' 
#' sd(unlist(resStrata["delta",]))
#' quantile(unlist(resStrata["se.delta",]))
#' 
#' mean(unlist(resStrata["p.value",])<0.05)
#' }
#' @export

influenceTest <- function(predict1,predict2, type){
  
  #### estimator
  if(type %in% names(predict1) == FALSE){
    stop("type ",type," is not in predict1 \n")
  }else{
    estimator1 <- predict1[[type]]
  }
  if(type %in% names(predict2) == FALSE){
    stop("type ",type," is not in predict2 \n")
  }else{
    estimator2 <- predict2[[type]]
  }
  
  #### influence function
  type.iid <- paste0(type,".iid")
  if(type.iid %in% names(predict1) == FALSE){
    stop("No iid decomposition for type ",type," in predict1 \n")
  }else{
    estimator1.iid <- predict1[[type.iid]]
  }
  if(type.iid %in% names(predict2) == FALSE){
    stop("No iid decomposition for type ",type," in predict2 \n")
  }else{
    estimator2.iid <- predict2[[type.iid]]
  }
  
  #### check sample size
  n <- length(estimator1.iid[1,1,])
  if(n!=length(estimator2.iid[1,1,])){
    stop("cannot compare estimates on different sample size \n")
  }
  
  #### test
  delta <- estimator1-estimator2
  se.delta <- sqrt(sum((estimator1.iid-estimator2.iid)^2))
  out <- data.frame(delta = delta, se.delta = se.delta, t.delta = delta/se.delta)
  out$p.value <- 2*(1-pnorm(abs(out$t.delta)))  
  
  return(out)
}

#----------------------------------------------------------------------
### influence2.R ends here
