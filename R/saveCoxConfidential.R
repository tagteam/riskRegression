#' library(survival)
#' library(lava)
#' set.seed(18)
#' trainSurv <- sampleData(300,outcome="survival")
#' testSurv <- sampleData(40,outcome="survival")
#' fit = coxph(Surv(time,event)~X1+X2+X3+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
#' u=saveCoxConfidential(fit,times=3)
#' #\dontrun{
#' # write object as plain text file
#' #sink("~/tmp/u.R")
#' #cat("U <- ")
#' #dput(u)
#' #sink(NULL)
#' ## reload object
#' #source("~/tmp/u.R")
#' #class(U) <- "CoxConfidential"
#' #}
#' predictRisk(U,newdata=testsurv)
#' cox1 = coxph(Surv(time,event)~strata(X1)+X2+X3+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
#' z<-saveCoxConfidential(cox1,c(2,5))
#' 
#' dput(z) ## get output to copy object
#' 
#' all.equal(predictRisk(z,newdata=testSurv),
#'           predictRisk(cox1,newdata=testSurv,times=c(2,5)))
#' 
#' cox2 = coxph(Surv(time,event)~X1+X2+X7+X9,data=trainSurv, y=TRUE, x = TRUE)
#' z<-saveCoxConfidential(cox2,c(2,5))
#' 
#' predictRisk(z,testSurv)-predictRisk(z,testSurv,c(2,5))
#'@export
saveCoxConfidential <- function(object, times){
  if (class(object) != "coxph"){
    stop("Object has to be a coxph object. ")
  }
  if (missing(times)) {
    stop("You have to specify the times. ")
  }
  times <- sort(times)
  fun.aggregate <- function(x, y){
    x.max <- max(x)
    x.min <- min(x)
    res <- rep(NA, length(times))
    ## handle times outside of the observed times for the specific strata
    res[times > x.min & times < x.max] <- y[sindex(x, times)][times > x.min & times < x.max] 
    res
  } 
  base.haz <- data.table::as.data.table(survival::basehaz(object))
  if (is.null(object$strata)){
    base.haz$strata <- rep(1,object$n)
    n.strata <- 1
  }
  else {
    n.strata <- length(levels(object$strata))
  }
  base.haz <- base.haz[,.(hazard0 = fun.aggregate(time,hazard)), by = strata]
  base.haz$times <- rep(times,n.strata)
  out <- list()
  out$Lambda0 <- as.data.frame(base.haz)
  out$times <- times
  out$coefficients <- object$coefficients
  out$formula <- object$formula
  out$means <- object$means
  class(out) <- "CoxConfidential"
  out
}

## saveTextCoxConfidential <- function(object){
  ## dput(object)
## }

## * predictRisk.CoxConfidential
##' @export
#' @rdname predictRisk
#' @method predictRisk CoxConfidential
predictRisk.CoxConfidential <- function(object,
                                        newdata,
                                        ...){
  if (missing(newdata)){
    stop("Argument newdata is missing. ")
  }
  formula <- object$formula
  formula[[2]] <- NULL
  frame <- Publish:::specialFrame(formula,
                                  data=newdata,
                                  specials="strata")
  design <- frame$design
  fit.times <- object$times
  n.times <- length(fit.times)
  n.newdata <- nrow(newdata)
  lp <- t(t(design)-object$means) %*% object$coefficients
  exp.lp <- exp(lp)
  if (!is.null(frame$strata[[1]])){
    dat.frame <- data.frame(eLP = rep(exp.lp, n.times), 
                            strata=rep(frame$strata[[1]], n.times), 
                            times=rep(fit.times,each=n.newdata), 
                            observation = rep(1:n.newdata, n.times))
  } 
  else {
    dat.frame <- data.frame(eLP = rep(exp.lp, n.times), 
                            times=rep(fit.times,each=n.newdata), 
                            observation = rep(1:n.newdata, n.times))  
  }
  Lambda0 <- object$Lambda0
  dat.frame2 <- data.frame(Lambda0)
  merged <- merge(dat.frame,dat.frame2)
  setorder(merged, times, observation)
  merged$risk <- 1-exp(-merged$eLP*merged$hazard0)
  matrix(merged$risk, ncol=n.times,nrow=n.newdata)
}
