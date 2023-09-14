##' Parse hyperparameters for data splitting algorithm
##'
##' @title Input for data splitting algorithms
##' @param split.method A character string specifying the algorithm for data splitting:
##' \itemize{
##' \item{"loob"} leave one out bootstrap
##' \item{"bootcv"} bootstrap cross validation
##' \item{"cv5"} 5-fold cross validation
##' \item{"loocv"} leave one out cross validation aka N-1 fold cross validation
##' \item{"632plus"} Efron's .632+ bootstrap
##' }
##' @param B Number of repetitions of bootstrap or k-fold cross-validation
##' @param N Sample size
##' @param M Subsample size. Default is N (no subsampling).
##' @param seed Integer passed to set.seed. If not given or NA no seed is set.
##' @return A list with the following elements:
##' \itemize{
##' \item{split.methodName}: the print name of the algorithm
##' \item{split.method}: the internal name of the algorithm
##' \item{index}: the index for data splitting. For bootstrap splitting this
##' is a matrix with B columns and M rows identifying the in-bag subjects. For k-fold
##' cross-validation this is a matrix with B columns identifying the membership to the k groups.
##' \item{k}: the k of k-fold cross-validation
##' \item{N}: the sample size
##' \item{M}: the subsample size
##' }
##' @seealso Score
##' @examples
##' # 3-fold crossvalidation
##' getSplitMethod("cv3",B=4,N=37)
##'
##' # bootstrap with replacement
##' getSplitMethod("loob",B=4,N=37)
##'
##' # bootstrap without replacement
##' getSplitMethod("loob",B=4,N=37,M=20)
##'
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
getSplitMethod <- function(split.method,B,N,M,seed){
  if (!missing(seed) && !is.null(seed) && !is.na(seed[[1]])) set.seed(seed)
  if (missing(split.method)) split.method <- ""
  split.methodName <- NULL
  split.method <- tolower(split.method)
  k <- as.numeric(substring(grep("^cv[0-9]+$",split.method,value=TRUE,ignore.case=TRUE),3))
  if (length(k)==0) k <- NULL
  ## none
  if (split.method %in% c("","noplan","none","no data",NA,FALSE,0L)) {
    B <- 0
    split.method <- "noplan"
    split.methodName <- "no data splitting"
  }
  ## classical cross-validation
  if (!is.null(k)){ ## classical cross-validation
    split.method <- "crossval"
    split.methodName <- paste(k,"fold cross-validation",sep="-")
    if (missing(B)) B <- 1 # repeat k-fold CrossVal one time
  }
  else{
    if (length(grep("loocv",split.method,ignore.case=TRUE))>0){
      split.method <- "loocv"
      split.methodName <- "LeaveOneOutCV"
      k <- N-1
      B <- 1
    }
  }
  ## resample or subsample bootstrap
  if(length(grep("^boot",split.method,value=FALSE,ignore.case=TRUE))>0){
    split.method <- "BootCv"
    split.methodName <- "BootCv"
    if (missing(B)) B <- 100
  }
  if (length(grep("632",split.method,value=FALSE,ignore.case=TRUE))>0){
    if (length(grep("plus|\\+",split.method,value=FALSE,ignore.case=TRUE))>0){
      split.method <- "Boot632plus"
      split.methodName <- ".632+"
      if (missing(B)) B <- 100
    }
    else{
      split.method <- "Boot632"
      split.methodName <- ".632"
      if (missing(B)) B <- 100
    }
  }
  ## default is leave one out bootstrap
  ## if ((length(grep("looboot|loob|leaveoneoutb",split.method,value=FALSE,ignore.case=TRUE))>0) ||
  if (!(split.method %in% c("noplan","crossval","loocv","BootCv","Boot632","Boot632plus"))){
    split.method <- "LeaveOneOutBoot"
    split.methodName <- "LeaveOneOutBoot"
    if (missing(B)) B <- 100
  }
  if (missing(M)) M <- N
  stopifnot(M[[1]]>0 && M[[1]]<=N[[1]])
  subsampling <- M!=N
  if (M<1) M <- round(M*N)
  effective.seeds <- sample(1:2^15, B)
  if (split.method=="noplan"){
    ResampleIndex <- NULL
  }
  else {
    ResampleIndex <- function(b, getSeed=FALSE){
      if (b > B){
        stop("Cannot get index for b > B.")
      }
      res<-switch(split.method,
                  "loocv"={
                    matrix(1:N,ncol=1)
                  },
                  "crossval"={
                    set.seed(effective.seeds[b])
                    sample(rep(1:k,length.out=N))
                  },
                  { ## default is bootstrap
                    ## split.method <- "LeaveOneOutBoot"
                    ## split.methodName <- "LeaveOneOutBoot"
                    set.seed(effective.seeds[b])
                    sort(sample(1:N,size=M,replace=!subsampling))
                  })
      if (getSeed){
        list(indeces = res, seed = effective.seeds[b])
      }
      else {
        res
      }
    }
  }
  if (missing(B)) {
    B <- switch(split.method,"loocv"={1},"noplan"={0},{100})
  }
  else{
    stopifnot(B[[1]]<0 || B[[1]]==round(B[[1]]))
  }
  out <- list(name=split.methodName,
              internal.name=split.method,
              index=ResampleIndex,
              k=k,
              B=B,
              M=M,
              N=N)
  class(out) <- "split.method"
  out
}

