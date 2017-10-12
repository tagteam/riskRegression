getSplitMethod <- function(split.method,B,N,M){
    if (missing(split.method))
        split.method <- ""
    split.methodName <- NULL
    k <- as.numeric(substring(grep("^cv[0-9]+$",split.method,value=TRUE,ignore.case=TRUE),3))
    if (length(k)==0) k <- NULL
    if (!is.null(k)){ ## classical cross-validation
        split.method <- "crossval"
        split.methodName <- paste(k,"fold cross-validation",sep="-")
    }
    else{
        if (length(grep("loocv",split.method,ignore.case=TRUE))>0){
            split.methodName <- "LeaveOneOutCV"
            k <- N-1
            B <- 1
        }
        else{
            ## some form of bootstrap
            match.BootCv <- length(grep("boot|outofbag",split.method,value=FALSE,ignore.case=TRUE))>0
            if (match.BootCv==FALSE){
                split.method <- "noPlan"
                split.methodName <- "no data splitting"
            }
            else{
                match.632 <- length(grep("632",split.method,value=FALSE,ignore.case=TRUE))>0
                match.plus <- length(grep("plus|\\+",split.method,value=FALSE,ignore.case=TRUE))>0
                match.bootloo <- length(grep("looboot",split.method,value=FALSE,ignore.case=TRUE))>0
                if (match.632==TRUE){
                    if (match.plus==TRUE){
                        split.method <- "Boot632plus"
                        split.methodName <- ".632+"
                    }
                    else{
                        split.method <- "Boot632"
                        split.methodName <- ".632"
                    }
                }
                else{
                    if (match.bootloo==TRUE){  
                        split.method <- "LeaveOneOutBoot"
                        split.methodName <- "LeaveOneOutBoot"
                    }
                    else {
                        split.method <- "BootCv"
                        split.methodName <- "BootCv"
                    }
                }
            }
        }
        if (missing(M)) M <- N
        stopifnot(M>0 && M<=N) 
        subsampling <- M!=N
        ##   if (!subsampling && resampleTraining)
        ##     stop("Resampling the training data is only available for subsampling")
        if (split.method %in% c("","noPlan","none")) {
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
            if (split.method=="loocv")
                ResampleIndex <- data.frame(id=1:N)
            else
                ResampleIndex <- do.call("cbind",lapply(1:B,function(b){sample(rep(1:k,length.out=N))}))
        }
        else{
            if (split.method %in% c("Boot632plus","BootCv","Boot632","LeaveOneOutBoot")){
                ResampleIndex <- do.call("cbind",lapply(1:B,function(b){
                    sort(sample(1:N,size=M,replace=!subsampling))
                }))
                colnames(ResampleIndex) <- paste("Train",1:B,sep=".")
            }
            else{
                ResampleIndex <- NULL
            }
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
    out <- list(name=split.methodName,
                internal.name=split.method,
                index=ResampleIndex,
                k=k,
                B=B,
                M=M,
                N=N)
    ##               resampleTraining=resampleTraining)
    class(out) <- "split.method"
    out
}

