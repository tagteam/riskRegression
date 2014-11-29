#' Formula interface for Fine-Gray regression competing risk models.
#' 
#' Formula interface for the function \code{crr} from the \code{cmprsk}
#' package.
#' 
#' The function \code{crr} allows to multiply some covariates by time before
#' they enter the linear predictor. This can be achieved with the formula
#' interface, however, the code becomes a little cumbersome. See the examples.
#'
#' @title Formula wrapper for crr from cmprsk 
#' @param formula A formula whose left hand side is a \code{Hist} object -- see
#' \code{\link{Hist}}.  The right hand side specifies (a linear combination of)
#' the covariates. See examples below.
#' @param data A data.frame in which all the variables of \code{formula} can be
#' interpreted.
#' @param cause The failure type of interest. Defaults to \code{1}.
#' @param y logical value: if \code{TRUE}, the response vector is returned in component \code{response}.
#' @param \dots ...
#' @return See \code{crr}.
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{riskRegression}}
#' @references Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute risk
#' regression for competing risks: interpretation, link functions and
#' prediction Research report 11/7. Department of Biostatistics, University of
#' Copenhagen
#' @keywords survival
##' @examples
##' 
##' library(riskRegression)
##' library(prodlim)
##' library(pec)
##' library(cmprsk)
##' library(lava)
##' d <- SimCompRisk(100)
##' f1 <- FGR(Hist(time,cause)~X1+X2,data=d)
##' print(f1)
##' 
##' ## crr allows that some covariates are multiplied by
##' ## a function of time (see argument tf of crr)
##' ## by FGR uses the identity matrix
##' f2 <- FGR(Hist(time,cause)~cov2(X1)+X2,data=d)
##' print(f2)
##' 
##' ## same thing, but more explicit:
##' f3 <- FGR(Hist(time,cause)~cov2(X1)+cov1(X2),data=d)
##' print(f3)
##' 
##' ## both variables can enter cov2:
##' f4 <- FGR(Hist(time,cause)~cov2(X1)+cov2(X2),data=d)
##' print(f4)
##' 
##' ## change the function of time
##' qFun <- function(x){x^2}
##' noFun <- function(x){x}
##' sqFun <- function(x){x^0.5}
##' 
##' ## multiply X1 by time^2 and X2 by time:
##' f5 <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2),data=d)
##' print(f5)
##' print(f5$crrFit)
##' ## same results as crr
##' with(d,crr(ftime=time,
##'            fstatus=cause,
##'            cov2=d[,c("X1","X2")],
##'            tf=function(time){cbind(qFun(time),time)}))
##' 
##' ## still same result, but more explicit
##' f5a <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=noFun),data=d)
##' f5a$crrFit
##' 
##' ## multiply X1 by time^2 and X2 by sqrt(time)
##' f5b <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=sqFun),data=d,cause=1)
##' 
##' ## additional arguments for crr
##' f6<- FGR(Hist(time,cause)~X1+X2,data=d, cause=1,gtol=1e-5)
#' @export
FGR <- function(formula,data,cause=1,y=TRUE,...){
    # {{{ read the data and the design
    cov1 <- function(x)x
    cov2 <- function(x)x
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=c("cov1","cov2"),
                                       specialsDesign=TRUE,
                                       stripSpecialNames=TRUE)
    EHF$cov1 <- cbind(EHF$cov1,EHF$design)
    EHF$design <- NULL
    theData <- data.frame(cbind(unclass(EHF$event.history),do.call("cbind",EHF[-1])))
    ## theData <- data.frame(EHF)
    # }}}
    # {{{ response
    response <- EHF$event.history
    cens.code <- attr(response,"cens.code")
    Y <- as.vector(response[,"time"])
    status <- as.vector(response[,"status"])
    if (match("event",colnames(response),nomatch=0)==0){
        event <- status
    }
    else{
        event <- as.numeric(prodlim::getEvent(response))
    }
    # }}}
    # {{{ cause of interest
    states <- prodlim::getStates(response)
    if (missing(cause)){
        cause <- 1
        message("Argument cause missing. Analyse cause: ",states[1])
    }
    else{
        if ((foundCause <- match(as.character(cause),states,nomatch=0))==0)
            stop(paste("Requested cause: ",cause," Available causes: ", states))
        else
            cause <- foundCause
    }  
    # }}}
    # {{{ covariate design matrices
    cov1 <- EHF$cov1
    cov2 <- EHF$cov2
    # }}}
    # {{{ call crr
    args <- list(ftime=Y,
                 fstatus=event,
                 cov1=cov1,
                 cov2=cov2,
                 failcode=cause,
                 cencode=length(states)+1,...)
    if (!is.null(cov2) && NCOL(cov2)>0){
        ## it is a bit silly to first remove special and then to add it again,
        ## but it seems to work.
        specialArguments <- sapply(paste("cov2(",colnames(cov2),")",sep=""),function(s)
                                   prodlim::parseSpecialNames(s,"cov2",list(cov2="tf")))
        this <- sapply(specialArguments,is.null)
        if (all(this)){
            args$tf <- bquote(function(x){matrix(rep(x,.(NCOL(cov2))),ncol=.(NCOL(cov2)),byrow=FALSE)})
        }
        else{
            tf.temp <- lapply(specialArguments,function(a)a$tf)
            if (any(this)){
                tf.temp[this] <- "I" 
            }
            names(tf.temp) <- NULL
            args$tf <- bquote(function(x){
                do.call("cbind",lapply(.(tf.temp),function(f){
                    do.call(f,list(x))
                }))
            })
        }
    }else{
        args$tf <- function(x){args$tf <- bquote(function(x){matrix(rep(x,.(NCOL(cov2))),
                                                                    ncol=.(NCOL(cov2)),
                                                                    byrow=FALSE)})}
    }
    args <- args[!sapply(args,is.null)]
    fit <- do.call(cmprsk::crr,args)
    fit$call <- NULL
    out <- list(crrFit=fit,
                cause=cause,
                "na.action"=attr(EHF,"na.action"))
    if (y==TRUE)
        out <- c(out,list(response=response))
    # }}}
    # {{{ clean up
    out$call <- match.call()
    if (is.null(out$call$cause)) out$call$cause <- cause
    out$terms <- terms(formula)
    class(out) <- "FGR"
    # }}}
    out
}
