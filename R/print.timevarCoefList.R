print.timevarCoefList <- function(x,digits=3,eps=10^-4,...){
  times <- names(x)
  lapply(1:length(x),function(i){
    cat("\ntime:",signif(as.numeric(times[i]),digits),"\n")
    xx <- data.frame(x[[i]])
    uu <- do.call("cbind",lapply(names(xx),function(a){
      switch(a,"Factor"={xx[,a]},
             "Pvalue"={
               vv <- format.pval(as.numeric(unlist(xx[,a,drop=FALSE])),digits=digits,eps=eps)
               vv[vv=="NA"] <- "--"
               vv
             },
             {
               vv <- as.character(signif(as.numeric(unlist(xx[,a,drop=FALSE])),digits=digits))
               vv[is.na(vv)] <- "--"
               vv
             })
    }))
    colnames(uu) <- colnames(x[[1]])
    print(uu,quote=FALSE)
  })
}
