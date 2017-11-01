### getComparisons.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Jan  3 2016 (13:30) 
## Version: 
## last-updated: Oct 21 2017 (10:09) 
##           By: Thomas Alexander Gerds
##     Update #: 42
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
getComparisons <- function(dt,NF,N,alpha,dolist=NF:1,se.fit){
    ## IMPORTANT: this function assumes that the
    ##            data are ordered according to model,times
    x=model=IF=NULL
    if (length(dolist)>0){
        Qnorm <- qnorm(1 - alpha/2)
        data.table::rbindlist(lapply(dolist,function(g){
            theta <- dt[,list(x=x[1]),by=model]
            delta <- theta[model%in%g[-1]][["x"]]-theta[model==g[1]][["x"]]
            se.delta <- dt[model%in%g[-1],list(se=sd(dt[model==g[1]][["IF"]]-IF)/sqrt(N)),by=model][["se"]]
            p <-2*pnorm(abs(delta)/se.delta,lower.tail=FALSE)
            if (se.fit==1L){
                lower <- delta - Qnorm * se.delta
                upper <- delta + Qnorm * se.delta
                data.table(model=theta[model%in%g[-1]][["model"]],
                           reference=g[1],
                           delta=delta,
                           se=se.delta,
                           lower=lower,
                           upper=upper,
                           p=p)
            }else{ ## only multisplit test
                data.table(model=theta[model%in%g[-1]][["model"]],
                           reference=g[1],
                           delta=delta,
                           p=p)
            }
        }))
    } else {
        NULL
    }
}
#----------------------------------------------------------------------
### getComparisons.R ends here
