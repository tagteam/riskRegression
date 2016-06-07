getSplitMethod <- function(splitMethod,B,N,M){
  if (missing(splitMethod))
    splitMethod <- ""
  splitMethodName <- NULL
  k <- as.numeric(substring(grep("^cv[0-9]+$",splitMethod,value=TRUE,ignore.case=TRUE),3))
  if (length(k)==0) k <- NULL
  if (!is.null(k)){ ## classical cross-validation
    splitMethod <- "crossval"
    splitMethodName <- paste(k,"fold cross-validation",sep="-")
  }
  else{
    if (length(grep("loocv",splitMethod,ignore.case=TRUE))>0){
      splitMethodName <- "LeaveOneOutCV"
      k <- N-1
      B <- 1
    }
    else{
      ## some form of bootstrap
      match.BootCv <- length(grep("boot|outofbag",splitMethod,value=FALSE,ignore.case=TRUE))>0
      if (match.BootCv==FALSE){
        splitMethod <- "noPlan"
        splitMethodName <- "no data splitting"
      }
      else{
        match.632 <- length(grep("632",splitMethod,value=FALSE,ignore.case=TRUE))>0
        match.plus <- length(grep("plus|\\+",splitMethod,value=FALSE,ignore.case=TRUE))>0
        if (match.632==TRUE){
          if (match.plus==TRUE){
            splitMethod <- "Boot632plus"
            splitMethodName <- ".632+"
          }
          else{
            splitMethod <- "Boot632"
            splitMethodName <- ".632"
          }
        }
        else{
          splitMethod <- "BootCv"
          splitMethodName <- "BootCv"}
      }
    }
  }
  if (missing(M)) M <- N
  stopifnot(M>0 && M<=N) 
  subsampling <- M!=N
  ##   if (!subsampling && resampleTraining)
  ##     stop("Resampling the training data is only available for subsampling")
  if (splitMethod %in% c("","noPlan","none")) {
    B <- 0
    ##     resampleTraining <- FALSE
  }
  else{
    if (missing(B)){
      if (length(k)>0) B <- 1 # repeat k-fold CrossVal ones
      else B <- 100 # 
    }
    else if (B==0) stop("No. of resamples must be a positive integer.")
  }
  if (length(k)>0){
    if (splitMethod=="loocv")
      ResampleIndex <- data.frame(id=1:N)
    else
      ResampleIndex <- do.call("cbind",lapply(1:B,function(b){sample(rep(1:k,length.out=N))}))
  }
  else{
    if (splitMethod %in% c("Boot632plus","BootCv","Boot632")){
      ResampleIndex <- do.call("cbind",lapply(1:B,function(b){
        sort(sample(1:N,size=M,replace=!subsampling))
      }))
      colnames(ResampleIndex) <- paste("Train",1:B,sep=".")
    }
    else{
      ResampleIndex <- NULL
    }
  }
  ##   if (is.logical(resampleTraining)){
  ##     if (resampleTraining==TRUE)
  ##       resampleTrainingSize <- N
  ##   }
  ##   else{
  ##     stopifnot(resampleTraining>0 &&resampleTraining==round(resampleTraining))
  ##     resampleTrainingSize <- resampleTraining
  ##     resampleTraining <- TRUE
  ##   }
  ##   if (resampleTraining==TRUE){
  ##     if (subsampling==TRUE && resampleTrainingSize<=M)
  ##       stop("Size for resampling the training indices should exceed ",M)
  ##     ##     if (subsampling==FALSE)
  ##     ##       stop("Resampling the training indices is only allowed for subsampling")
  ##     ResampleIndex <- apply(ResampleIndex,2,function(x){
  ##       sort(c(x,sample(x,replace=TRUE,size=resampleTrainingSize-M)))
  ##     })
  ##   }
  out <- list(name=splitMethodName,
              internal.name=splitMethod,
              index=ResampleIndex,
              k=k,
              B=B,
              M=M,
              N=N)
  ##               resampleTraining=resampleTraining)
  class(out) <- "splitMethod"
  out
}

