
### boxplot.Score.R --- 
#----------------------------------------------------------------------
## author: Thomas Alexander Gerds
## created: Aug 15 2016 (09:45) 
## Version: 
## last-updated: Dec 12 2016 (07:03) 
##           By: Thomas Alexander Gerds
##     Update #: 50
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
#' @param outcomeLabel Title label for column which shows the outcome status 
#' @param refline Logical, for \code{type="diff"} only. If \code{TRUE} draw a red vertical line at \code{0}.
#' @param ... not used
##' @examples
##' # binary outcome
##' db=sampleData(100,outcome="binary")
##' fitconv=glm(Y~X3+X5,data=db,family=binomial)
##' fitnew=glm(Y~X1+X3+X5+X6+X7,data=db,family=binomial)
##' scoreobj=Score(list(new=fitnew,conv=fitconv),
##'         formula=Y~1,contrasts=list(c(2,1)),
##'                data=db,summary="riskQuantile",nullModel=FALSE)
##' boxplot(scoreobj)
##'
##' # survival outcome
##' library(survival)
##' ds=sampleData(100,outcome="survival")
##' fitconv=coxph(Surv(time,event)~X6,data=ds,x=TRUE,y=TRUE)
##' fitnew=coxph(Surv(time,event)~X6+X9,data=ds,x=TRUE,y=TRUE)
##' scoreobj=Score(list("conventional model"=fitconv,"new model"=fitnew),
##'                 formula=Hist(time,event)~1, data=ds,
##'                 summary="riskQuantile",metrics=NULL, plots=NULL,
##'                 c(0,0.25,0.5,0.75,1),
##'                 times=5,nullModel=FALSE)
##' boxplot(scoreobj)
##'
##' scoreobj1=Score(list("conventional model"=fitconv,"new model"=fitnew),
##'                 formula=Hist(time,event)~1, data=ds,
##'                 summary="riskQuantile",metrics=NULL, plots=NULL,
##'                 times=5,nullModel=FALSE,compare=list(c(2,1)))
##' boxplot(scoreobj1)
##'
##' # competing risks outcome
##' data(Melanoma)
##' fitconv = CSC(Hist(time,status)~invasion+age+sex,data=Melanoma)
##' fitnew = CSC(Hist(time,status)~invasion+age+sex+logthick,data=Melanoma)
##' scoreobj=Score(list("Conventional model"=fitconv,"New model"=fitnew),
##'                formula=Hist(time,status)~1,
##'                data=Melanoma,summary="riskQuantile",times=5*365.25,nullModel=FALSE)
##' boxplot(scoreobj)
##' 
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
##'                data=dcr,summary="riskQuantile",times=5,nullModel=FALSE)
##' boxplot(scoreobj)
##' 
##' 
#' @export
boxplot.Score <- function(x,model,reference,type,timepoint,lwd=3,xlim,xlab,main,outcomeLabel,eventLabels,refline=TRUE,...){
    times=cause=models=NULL
    fitted <- x$models
    models <- names(x$models)
    if (missing(type)) {
        if (length(models)==1) 
            type="risk"
        else
            type=ifelse(NROW(x$riskQuantile$contrasts)>0,"diff","risk")
    }
    if (type=="diff"){
        pframe <- x$riskQuantile$contrasts
    } else{
        pframe <- x$riskQuantile$score
    }
    if (x$responseType!='binary'){
        if (missing(timepoint))
            timepoint <- max(pframe[["times"]])
        else ## can only do one timepoint
            timepoint <- timepoint[[1]]
        pframe <- pframe[times==timepoint]
    }
    if(missing(model)) mod <- pframe[,model[1]] else mod <- model
    if (type=="diff"){
        if(missing(reference)) ref <- pframe[,reference[1]] else ref <- reference
        pframe <- pframe[model==mod & reference==ref]
    }else{
        pframe <- pframe[model==mod ]
    }
    qq.pos <- grep("^Q",colnames(pframe))
    if (missing(xlim))
        if (type=="risk"){xlim=c(0,100) 
        } else {
            max <- ceiling(max(abs(100*(pframe[,qq.pos,with=FALSE]))))
            xlim=c(-max,max)
        }
    if (missing(main))
        if (type=="risk") main=mod else main=switch(x$responseType,
                                                    "competing.risks"={paste("Difference in predicted risks\nof an event of type",x$cause,"\nbefore time",timepoint)},
                                                    "survival"={paste("Difference in predicted risks\nof an event before time",timepoint)},
                                                    "binary"={"Difference in predicted risks"})
    if (missing(outcomeLabel)) outcomeLabel <- switch(x$responseType,
                                                      "competing.risks"={paste("Event status\nat time",timepoint)},
                                                      "survival"={paste("Event status\nat time",timepoint)},
                                                      "binary"={"Event status"})
    if (missing(xlab))
        if (type=="risk") xlab="Predicted risk" else xlab=""
    if (x$responseType!="competing.risks"){
        plot(0,0,type="n",
             main=main,
             xlim = xlim,
             ylim = c(0,NROW(pframe)),
             axes=FALSE,
             xlab = xlab,
             ylab = "")
        if (refline==TRUE)
            abline(v=0,col=2,lwd=2)
        axis(1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),labels=paste(seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),"%"))
        causes <- pframe[,unique(cause)]
        if (missing(eventLabels)) eventLabels <- causes
        text(x=xlim[1],y=c(0.5,1.5,2.5,3),labels=c(eventLabels,outcomeLabel),pos=2,xpd=NA)
        if (type=="diff"){
            mtext(paste(ref,"higher risk"),side=1,adj=0,line=par()$mgp[1])
            mtext(paste(mod,"higher risk"),side=1,adj=1,line=par()$mgp[1])
        }
        bxp(list(stats=t(100*pframe[cause=="overall",qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=0.5,axes=FALSE)
        bxp(list(stats=t(100*pframe[cause=="event",qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=1.5,axes=FALSE)
        bxp(list(stats=t(100*pframe[cause=="event-free",qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=2.5,axes=FALSE)
    }else{
        plot(0,0,type="n",
             main=main,
             xlim = xlim,
             ylim = c(0,NROW(pframe)),
             axes=FALSE,
             xlab = xlab,
             ylab = "")
        if (refline==TRUE)
            abline(v=0,col=2,lwd=2)
        ## axis(1,at=seq(0,100,25),labels=paste(seq(0,100,25),"%"))
        axis(1,at=seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),labels=paste(seq(xlim[1],xlim[2],(xlim[2]-xlim[1])/4),"%"))
        causes <- pframe[,unique(cause)]
        if (missing(eventLabels)) eventLabels <- causes
        ypos <- c((1:(length(causes)))-0.5,length(causes))
        text(x=xlim[1],
             y=ypos,
             labels=c(eventLabels,outcomeLabel),
             pos=2,
             xpd=NA)
        if (type=="diff"){
            mtext(paste(ref,"higher risk"),side=1,adj=0,line=par()$mgp[1])
            mtext(paste(mod,"higher risk"),side=1,adj=1,line=par()$mgp[1])
        }
        for (i in 1:length(causes)){
            cc <- causes[[i]]
            bxp(list(stats=t(100*pframe[cause==cc,qq.pos,with=FALSE,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=ypos[i],axes=FALSE)
        }
    }
    invisible(x)
}

## plot.riskQuantile <- function(x,text.title="Outcome after 10 years",xlab="",text=rownames(x),...){
    ## plot(0,0,type="n",axes=FALSE,xlim=c(-10,10),ylim=c(1,5),xlab=xlab,ylab="")
    ## axis(1,at=c(-10,-5,-2.5,0,2.5,5,10))
    ## bxp(list(stats=t(100*x[4,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=1.5,axes=FALSE)
    ## bxp(list(stats=t(100*x[3,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=2.5,axes=FALSE)
    ## bxp(list(stats=t(100*x[2,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=3.5,axes=FALSE)
    ## bxp(list(stats=t(100*x[1,,drop=FALSE]),n=10),add=TRUE,horizontal=TRUE,at=4.5,axes=FALSE)
    ## text(x=-10,y=c(1.5,2.5,3.5,4.5),labels=rev(text),xpd=NA)
    ## text(x=-10,y=5,labels=text.title,xpd=NA)
## }


#----------------------------------------------------------------------
### boxplot.Score.R ends here
