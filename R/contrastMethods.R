contrast <- function(x,...){
  ## this leads to problems with small (1-2 row) newdata
  ##   if (length(unique(x))<=2){
  ##     contrast.character(as.character(x),...)
  ##   }
  ##   else{
  UseMethod("contrast",x)
  ##   }
}
contrast.default <- function(x,...){
  stop("No contrast method for class",class(x))
}
contrast.matrix <- function(x,...){
  stop("No contrast method for a matrix")
}
contrast.numeric <- function(x,name,...){
  out <- matrix(x,ncol=1)
  colnames(out) <- name
  out
}
contrast.integer <- function(x,name,...){
  contrast.numeric(x=x,name=name,...)
}
contrast.factor <- function(x,name,levels,sep=":",...){
  if (missing(levels)|is.null(levels)){
    levels <- levels(x)
  }
  else{
    if (any(wo <- match(levels(x),levels,nomatch=0)==0))
      stop("Level(s): '",paste(levels(x)[wo],collapse=","),"' unknown.")
  }
  ref <- levels[1]
  ## fix for factors with one level
  ## used e.g. in predicting based
  ## expand.grid values
  if (length(levels)==1) {
    out <- matrix(as.numeric(x==levels),ncol=1)    
  }
  else{
    if (length(levels)==2)
      out <- matrix(as.numeric(x==levels[-1]),ncol=1)
    else{
      out <- do.call("cbind",lapply(levels[-1],function(l){
        as.numeric(x==l)
      }))
    }
  }
  if (length(levels)==1){
    colnames(out) <- paste(name,levels,sep=sep)
  }else{
    colnames(out) <- sapply(levels[-1],function(l){
      paste(name,l,sep=sep)
    })}
  attr(out,"levels") <- levels
  attr(out,"ref") <- ref
  out
}
contrast.character <- function(x,name,...){
  contrast.factor(factor(x),name=name,...)
}
contrast.logical <- function(x,name,...){
  contrast.factor(factor(x),name=name,...)
}

