### boxplot.Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Aug 15 2016 (09:45) 
## Version: 
## last-updated: Jan 27 2020 (12:26) 
##           By: Thomas Alexander Gerds
##     Update #: 141
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:
#' Boxplot risk quantiles
#' 
#' Retrospective boxplots of risk quantiles conditional on outcome
#' @param x Score object obtained by calling function \code{Score}.
#' @param model Choice of risk prediction model
#' @param reference Choice of reference risk prediction model for
#'     calculation of risk differences.
#' @param type Either \code{"risk"} for predicted risks or
#'     \code{"diff"} for differences between predicted risks.
#' @param timepoint time point specifying the prediction horizon
#' @param lwd line width
#' @param xlim x-axis limits
#' @param xlab x-axis label
#' @param main title of plot
#' @param outcome.label Title label for column which shows the outcome
#'     status
#' @param outcome.label.offset Vertical offset for outcome.label 
#' @param event.labels Labels for the different events (causes).
#' @param refline Logical, for \code{type="diff"} only. If \code{TRUE}
#'     draw a red vertical line at \code{0}.
#' @param overall Logical. Tag to be documented.
#' @param add Logical. Tag to be documented.
#' @param ... not used
##' @examples
##' # binary outcome
##' library(data.table)
##' library(prodlim)
##' db=sampleData(40,outcome="binary")
##' fitconv=glm(Y~X3+X5,data=db,family=binomial)
##' fitnew=glm(Y~X1+X3+X5+X6+X7,data=db,family=binomial)
##' x=Score(list(new=fitnew,conv=fitconv),
##'         formula=Y~1,contrasts=list(c(2,1)),
##'                data=db,summary="riskQuantile",null.model=FALSE)
##' boxplot(x)
##'
##' # survival outcome
##' library(survival)
##' ds=sampleData(40,outcome="survival")
##' fitconv=coxph(Surv(time,event)~X6,data=ds,x=TRUE,y=TRUE)
##' fitnew=coxph(Surv(time,event)~X6+X9,data=ds,x=TRUE,y=TRUE)
##' \dontrun{ 
##' scoreobj=Score(list("conventional model"=fitconv,"new model"=fitnew),
##'                 formula=Hist(time,event)~1, data=ds,
##'                 summary="riskQuantile",metrics=NULL, plots=NULL,
##'                 c(0,0.25,0.5,0.75,1),
##'                 times=5,null.model=FALSE)
##' boxplot(scoreobj)
##' 
##' scoreobj1=Score(list("conventional model"=fitconv,"new model"=fitnew),
##'                 formula=Hist(time,event)~1, data=ds,
##'                 summary="riskQuantile",metrics=NULL, plots=NULL,
##'                 times=5,null.model=FALSE,compare=list(c(2,1)))
##' boxplot(scoreobj1)
##' }
##' 
##' # competing risks outcome
##' library(survival)
##' data(Melanoma, package = "riskRegression")
##' fitconv = CSC(Hist(time,status)~invasion+age+sex,data=Melanoma)
##' fitnew = CSC(Hist(time,status)~invasion+age+sex+logthick,data=Melanoma)
##' scoreobj=Score(list("Conventional model"=fitconv,"New model"=fitnew),
##'                formula=Hist(time,status)~1,
##'                data=Melanoma,metrics=NULL,summary="riskQuantile",times=5*365.25,null.model=FALSE)
##' boxplot(scoreobj)
##'
##' # more than 2 competing risks
##' m=lava::lvm(~X1+X2+X3)
##' lava::distribution(m, "eventtime1") <- lava::coxWeibull.lvm(scale = 1/100)
##' lava::distribution(m, "eventtime2") <- lava::coxWeibull.lvm(scale = 1/100)
##' lava::distribution(m, "eventtime3") <- lava::coxWeibull.lvm(scale = 1/100)
##' lava::distribution(m, "censtime") <- lava::coxWeibull.lvm(scale = 1/100)
##' lava::regression(m,eventtime2~X3)=1.3
##' m <- lava::eventTime(m,
##' time ~ min(eventtime1 = 1, eventtime2 = 2, eventtime3 = 3, censtime = 0), "event")
##' set.seed(101)
##' dcr=as.data.table(lava::sim(m,101))
##' fitOld = CSC(Hist(time,event)~X1+X2,data=dcr)
##' fitNew = CSC(Hist(time,event)~X1+X2+X3,data=dcr)
##' scoreobj=Score(list("Conventional model"=fitOld,"New model"=fitNew),
##'                formula=Hist(time,event)~1,
##'                data=dcr,summary="riskQuantile",times=5,null.model=FALSE)
##' boxplot(scoreobj)
##' 
##' 
#' @export
boxplot.Score <- function(x,
                          model,
                          reference,
                          type="risk",
                          timepoint,
                          overall=1L,
                          lwd=3,
                          xlim,
                          xlab="",
                          main,
                          outcome.label,
                          outcome.label.offset=0,
                          event.labels,
                          refline=(type!="risk"),
                          add=FALSE,
                          ...){
    if (is.null(x$riskQuantile))
        stop("Object has no information for reclassification boxplots.\nYou should call the function \"riskRegression::Score\" with summary=\"riskQuantile\".")
    times=cause=models=NULL
    fitted <- x$models
    models <- names(x$models)
    if (missing(type)) { type = "risk"}
    ## if (length(models)==1) 
    ## type="risk"
    ## else
    ## type=ifelse(NROW(x$riskQuantile$contrasts)>0,"diff","risk")
    ## }

    if (type=="diff"){
        warning("Boxplots of differences of predicted risk are non-proper\nin the sense that a miscalibrated model can apparently look good.")
        pframe <- x$riskQuantile$contrasts[model!="Null model"]
    } else{
        pframe <- x$riskQuantile$score[model!="Null model"]
    }
    if (x$response.type!='binary'){
        if (missing(timepoint))
            timepoint <- max(pframe[["times"]])
        else ## can only do one timepoint
            timepoint <- timepoint[[1]]
        pframe <- pframe[times==timepoint]
    }
    if(missing(model)) mod <- pframe[,model[1]] else mod <- model
    if (!(mod %in% unique(pframe[["model"]])))
        stop(paste("Requested model ",mod, "not fitted in object."))
    if (type=="diff"){
        if(missing(reference)) ref <- pframe[,reference[1]] else ref <- reference
        pframe <- pframe[model==mod & reference==ref]
    }else{
        pframe <- pframe[model==mod ]
    }
    qq.pos <- grep("^Q",colnames(pframe))
    causes <- pframe[,unique(cause)]
    if (overall==FALSE) causes <- causes[causes!="overall"]
    if (missing(xlim))
        if (type=="risk"){xlim=c(0,1) 
        } else {
            max <- ceiling(max(abs(100*(pframe[,qq.pos,with=FALSE]))))/100
            xlim=c(-max,max)
        }
    if (missing(main))
        if (type=="risk") main=mod else main=switch(x$response.type,
                                                    "competing.risks"={paste("Difference in predicted risks\nof an event of type",x$cause,"\nbefore time",timepoint)},
                                                    "survival"={paste("Difference in predicted risks\nof an event before time",timepoint)},
                                                    "binary"={"Difference in predicted risks"})
    if (missing(outcome.label)) outcome.label <- switch(x$response.type,
                                                        "competing.risks"={paste("Event status\nat time",timepoint)},
                                                        "survival"={paste("Event status\nat time",timepoint)},
                                                        "binary"={"Event status"})
    plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = c(0,length(causes)),xlim = xlim,ylab="",xlab=xlab)
    axis1.DefaultArgs <- list(side=1,
                              las=1,
                              at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),
                              paste(100*seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),"%"))
    abline.DefaultArgs <- list(v=0,col=2,lwd=2)
    bxp.DefaultArgs <- list(border="black",add=TRUE,horizontal=TRUE,axes=FALSE)
    control <- prodlim::SmartControl(call= list(...),
                                     keys=c("plot","bxp","axis1","abline"),
                                     ignore=NULL,
                                     ignore.case=TRUE,
                                     defaults=list("plot"=plot.DefaultArgs,
                                                   "bxp"=bxp.DefaultArgs,
                                                   "abline"=abline.DefaultArgs,
                                                   "axis1"=axis1.DefaultArgs),
                                     forced=list("plot"=list(axes=FALSE),
                                                 "axis1"=list(side=1)),
                                     verbose=TRUE)
    
    if (missing(xlab))
        if (type=="risk") xlab="Predicted risk" else xlab=""
    if (add==0L) do.call("plot",control$plot)
    if (missing(event.labels)) event.labels <- causes
    else if (length(event.labels)!=length(causes))
        stop(paste0("Argument event.labels has wrong length: ",length(event.labels),"\nShould be a character vector of length: ",length(causes),"\none for each cause: ",paste(causes,collapse=",")))
    ypos <- c((1:(length(causes)))-0.5,length(causes))
    text(x=xlim[1],
         y=ypos + c(rep(0,length(event.labels)),outcome.label.offset),
         labels=c(event.labels,outcome.label),
         pos=2,
         xpd=NA)
    if (type=="diff"){
        mtext(paste(ref,"higher risk"),side=1,adj=0,line=par()$mgp[1])
        mtext(paste(mod,"higher risk"),side=1,adj=1,line=par()$mgp[1])
    }
    for (i in 1:length(causes)){
        cc <- causes[[i]]
        args.i <- control$bxp
        args.i$z <- list(stats=t(pframe[cause==cc,qq.pos,with=FALSE,drop=FALSE]),n=10)
        args.i$at=ypos[i]
        do.call("bxp",args.i)
        ## bxp(,add=TRUE,horizontal=TRUE,at=ypos[i],axes=FALSE)
    }
    if (refline==TRUE)
        do.call("abline",control$abline)
    do.call("axis",control$axis1)
    invisible(x)
}



#----------------------------------------------------------------------
### boxplot.Score.R ends here
