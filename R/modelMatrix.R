modelMatrix <- function(formula,
                        data,
                        factorLevels=NULL,
                        intercept){
  if (is.null(formula)){
    if (!missing(intercept))
      cbind(Intercept=rep(intercept,NROW(data)))
    else
      NULL
  }
  else{
    m <- model.frame(formula,data)
    varnames <- names(m)
    res <- lapply(varnames,function(v){
      ## if (match(v,names(data),nomatch=0)==0){
      ## stop("modelMatrix: Variable ",v," not found in data",call.=FALSE)
      ## }
      ## if (match(v,names(factorLevels),nomatch=0)==0){
      ## stop("in modelMatrix: Variable ",v," not found in factorLevels",call.=FALSE)
      ## }
      if (match(v,names(factorLevels),nomatch=0)==0){
        ff <- data[,v,drop=TRUE]
        if (is.factor(ff)){
          fflevels <- levels(ff)
          x <- contrast(ff,name=v,levels=fflevels)
        }
        else{
          ## DOES NOT MAKE SENSE
          ## BECAUSE WE DO NOT KNOW THE REF LEVEL
          ## OF ORIGINAL DATA WHEN LOOKING AT
          ## NEW DATA (which does not have to include all levels)
          ## if (length(unique(ff))<3){
          ## fflevels <- unique(ff)
          ## x <- contrast.factor(ff,name=v,levels=fflevels)
          ## } else{
          x <- contrast(ff,name=v)
          ## }
        }
      }
      else{
        x <- contrast(data[,v,drop=TRUE],name=v,levels=factorLevels[[v]])
      }
    })
    names(res) <- varnames
    out <- do.call("cbind",res)
    if (!missing(intercept))
      if (NROW(out)==0)
        out <- cbind("Intercept"=rep(intercept,NROW(data)))
      else
        out <- cbind("Intercept"=rep(intercept,NROW(out)),out)
    attr(out,"variableNames") <- varnames
    factorLevels <- lapply(res,function(r){attr(r,"levels")})
    factorLevels <- factorLevels[!sapply(factorLevels,is.null)]
    attr(out,"factorLevels") <- factorLevels
    refLevels <- lapply(res,function(r){attr(r,"ref")})
    refLevels <- refLevels[!sapply(refLevels,is.null)]
    attr(out,"refLevels") <- refLevels
    class(out) <- "modelMatrix"
    if (NROW(out)==0 && length(all.vars(formula))>0)
      stop("No data in model matrix.")
    out
  }
}

                
