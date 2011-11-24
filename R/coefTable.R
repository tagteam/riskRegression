coefTable <- function(formula,
                      coef,
                      se,
                      levels,
                      ref,
                      trans="exp",
                      test="wald"){
  vars <- all.vars(formula)
  if (length(vars)==0){
    NULL
  }
  else{
    if (test=="wald"){
      testStatistic <- coef/se
      pvalue <- (1 - pnorm(abs(testStatistic))) * 2
      ## pvalue <- format.pval(pvalue,digits=digits,eps=eps)
      names(pvalue) <- names(pvalue)
      pvalue[se==0] <- NA
      if (any(se==0))
        warning("Some standard errors are zero. It seems that the model did not converge")
    }
    else{
      stop(paste("test ",test," not availabe."))
    }
    coefMat <- do.call("rbind",lapply(vars,function(v){
      covname <- strsplit(v,":")[[1]][[1]]
      if (is.null(levels[[covname]])){
        out <- data.frame(v)
            out <- cbind(v,coef[v],exp(coef[v]),se[v],testStatistic[v],pvalue[v])
      }
      else{
        rlev <- ref[[covname]]
        out <- do.call("rbind",lapply(levels[[covname]],function(l){
          V <- paste(covname,l,sep=":")
          if (match(V,paste(covname,rlev,sep=":"),nomatch=FALSE)){
            dV <- data.frame(paste(covname,rlev,sep=":"))
            lout <- cbind(dV,as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA),as.numeric(NA))
          }
          else{
            dV <- data.frame(V)
            lout <- cbind(dV,coef[V],exp(coef[V]),se[V],testStatistic[V],pvalue[V])
          }
          colnames(lout) <- c("Factor","Coef","exp(Coef)","StandardError","testStatistic","Pvalue")
          lout
        }))
      }
      colnames(out) <- c("Factor","Coef","exp(Coef)","StandardError","testStatistic","Pvalue")
      out
    }))
    ## colnames(coefMat) <- c("Factor","Coef","exp(Coef)","StandardError","testStatistic","Pvalue")
    ## rownames(coefMat) <- rep("",NROW(coefMat))
    coefMat
  }
}
