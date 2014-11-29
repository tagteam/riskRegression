#' Formula interface for Fine-Gray regression competing risk models.
#' 
#' Formula interface for the function \code{crr} from the \code{cmprsk}
#' package.
#' 
#' The function \code{crr} allows to multiply some covariates by time before
#' they enter the linear predictor. This can be achieved with the formula
#' interface, however, the code becomes a little cumbersome. See the examples.
#' 
#' @param formula A formula whose left hand side is a \code{Hist} object -- see
#' \code{\link{Hist}}.  The right hand side specifies (a linear combination of)
#' the covariates. See examples below.
#' @param data A data.frame in which all the variables of \code{formula} can be
#' interpreted.
#' @param cause The failure type of interest. Defaults to \code{1}.
#' @param \dots ...
#' @return See \code{crr}.
#' @author Thomas Alexander Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{riskRegression}}
#' @references Gerds, TA and Scheike, T and Andersen, PK (2011) Absolute risk
#' regression for competing risks: interpretation, link functions and
#' prediction Research report 11/7. Department of Biostatistics, University of
#' Copenhagen
#' @keywords survival
#' @examples
#' 
#' library(riskRegression)
#' library(prodlim)
#' library(pec)
#' library(cmprsk)
#' d <- SimCompRisk(100)
#' f1 <- FGR(Hist(time,cause)~X1+X2,data=d)
#' print(f1)
#' 
#' ## crr allows that some covariates are multiplied by
#' ## a function of time (see argument tf of crr)
#' ## by FGR uses the identity matrix
#' f2 <- FGR(Hist(time,cause)~cov2(X1)+X2,data=d)
#' print(f2)
#' 
#' ## same thing, but more explicit:
#' f3 <- FGR(Hist(time,cause)~cov2(X1)+cov1(X2),data=d)
#' print(f3)
#' 
#' ## both variables can enter cov2:
#' f4 <- FGR(Hist(time,cause)~cov2(X1)+cov2(X2),data=d)
#' print(f4)
#' 
#' ## change the function of time
#' qFun <- function(x){x^2}
#' noFun <- function(x){x}
#' sqFun <- function(x){x^0.5}
#' 
#' ## multiply X1 by time^2 and X2 by time:
#' f5 <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2),data=d)
#' 
#' ## the same but more explicit
#' f5a <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=noFun),data=d)
#' 
#' ## multiply X1 by time^2 and X2 by sqrt(time)
#' f5b <- FGR(Hist(time,cause)~cov2(X1,tf=qFun)+cov2(X2,tf=sqFun),data=d,cause=1)
#' 
#'
#' @export
FGR <- function(formula,data,cause=1,...){
  # {{{ check if formula has the form Hist(time,event)~X1+X2+...

  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")
      ||
      (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
    stop("Invalid specification of formula. Perhaps forgotten right hand side?\nNote that any subsetting, ie data$var or data[,\"var\"], is invalid for this function.")
  }
  else
    if (!(formula.names[2] %in% c("Hist"))) stop("formula is NOT a proper event history formula,\nwhich must have a `Hist' object as response.")
  
  # }}}
  # {{{ read the data and the design
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  if (match("subset",names(call),nomatch=FALSE))
      stop("Subsetting of data is not possible.")
  m <- m[match(c("","formula","data","subset","na.action"),names(m),nomatch = 0)]
  m[[1]]  <-  as.name("model.frame")
  if (missing(data)) stop("Argument 'data' is missing")
  formList <- readFormula(formula,
                          specials=c("cov2","cov1"),
                          specialArgumentNames=list("cov2"="tf"),
                          unspecified="cov1")
  m$formula <- formList$allVars
  theData <- eval(m, parent.frame()) 
  if ((nMiss <- (NROW(data)-NROW(theData)))>0)
      warning(nMiss," lines have been removed from data due to missing values")
  if (NROW(theData) == 0) stop("No (non-missing) observations")
  # }}}
  # {{{ response
  response <- model.response(model.frame(formula=formList$Response,data=theData))
  cens.code <- attr(response,"cens.code")
  Y <- as.vector(response[,"time"])
  time  <- numeric(length(Y))
  status <- as.vector(response[,"status"])
  if (match("event",colnames(response),nomatch=0)==0){
      event <- status
  }
  else{
      event <- as.numeric(getEvent(response))
  }
  # }}}
  # {{{ cause of interest
  states <- getStates(response)
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
  cov1 <- modelMatrix(formula=formList$cov1$formula,
                      data=theData,
                      intercept=NULL)
  cov2 <- modelMatrix(formula=formList$cov2$formula,
                      data=theData,
                      intercept=NULL)
  if (!is.null(cov1))
      class(cov1) <- "matrix"
  if (!is.null(cov2))
      class(cov2) <- "matrix"
  # }}}
  # {{{ call crr
  args <- list(ftime=Y,
               fstatus=event,
               cov1=cov1,
               cov2=cov2,
               failcode=cause,
               cencode=length(states)+1,...)
  if (NCOL(cov2)>0){
      this <- sapply(formList$cov2$specialArguments,is.null)
      ## if (all(this <- sapply(formList$cov2$specialArguments,is.null))){
      ## args$tf <- function(x){matrix(x,ncol=NCOL(cov2),byrow=FALSE)}
      ## args$tf <- function(x){matrix(x,ncol=NCOL(cov2),byrow=FALSE)}
      ## }
      ## else{
      tf.temp <- lapply(formList$cov2$specialArguments,function(a)a$tf)
      if (any(this)){
          id <- function(x)x
          tf.temp[this] <- "id" 
      }
      args$tf <- function(x){
          do.call("cbind",lapply(tf.temp,function(f){
              do.call(f,list(x))
          }))
      }
  }
  else{
      args$tf <- function(x){matrix(x,ncol=NCOL(cov2),byrow=FALSE)}
  }
  args <- args[!sapply(args,is.null)]
  fit <- do.call("crr",args)
  fit$call <- NULL
  out <- list(crrFit=fit,response=response,cause=cause)
  # }}}
  # {{{ clean up
  out$call <- match.call()
  if (is.null(out$call$cause)) out$call$cause <- cause
  out$call$formula <- eval(out$call$formula)
  class(out) <- "FGR"
  # }}}
  out
}
