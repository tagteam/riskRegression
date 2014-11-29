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
