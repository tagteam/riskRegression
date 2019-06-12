##' Level plots for predicted risks
##'
##' Level plots for predicted risks
##' @title Level plots for risk prediction models 
##' @param object risk prediction model object
##' @param formula formula 
##' @param data data
##' @param times prediction horizon time point for survival models
##' @param cause cause of interest for competing risk models
##' @param ... passed to lattice::levelplot
##' @examples
##' 
##' # ---------- logistic regression --------------------
##' expit <- function(x){exp(x)/(1+exp(x))}
##' partyData <- function(N){
##'   Age <- runif(N,.5,15)
##'   Parasites <- rnorm(N,mean=3.5-0.03*Age)
##'   Fever <- factor(rbinom(N,1,expit(-3.5-.3*Age+.55*Parasites+0.15*Age*Parasites)))
##'   data.frame(Fever,Age,Parasites)
##' }
##' d <- partyData(100)
##' f <- glm(Fever~Age+Parasites,data=d,family="binomial")
##' riskLevelPlot(f,Fever~Age+Parasites,data=d)
##' library(party)
##' ct <- ctree(Fever~Age+Parasites,data=d)
##' riskLevelPlot(ct,Fever~Age+Parasites,d)
##' library(randomForest)
##' rf <- randomForest(Fever~Age+Parasites,data=d)
##' riskLevelPlot(f,Fever~Age+Parasites,d)
##' riskLevelPlot(rf,Fever~Age+Parasites,d)
##'
##' # ---------- survival analysis --------------------
##'
##' # --simulate an artificial data frame
##' # with survival response and three predictors
##' 
##' library(survival)
##' library(prodlim)
##' survData <- function( N ){
##'   sdata <- SimSurv( N )
##'   sdata$X3 <- rnorm(N)
##'   sdata
##' }
##' set.seed(140515)
##' sdat <- survData(100)
##' sdat <- sampleData(100,outcome="survival")
##' # -- fit a Cox regression model 
##' survForm = Surv(time,event) ~ X7 + X8
##' cox <- coxph(survForm, data = sdat,x=TRUE)
##'
##' # --choose a time horizon for the predictions and plot the risks
##' timeHorizon <- floor(median(sdat$time))
##' riskLevelPlot(cox, survForm, data = sdat, horizon = timeHorizon)
##'
##' survForm2 = Surv(time,status) ~ strata(X1) + X2 + X3
##' cox2 <- coxph(survForm2, data = sdat)
##' riskLevelPlot(cox2, survForm2, data = sdat, horizon = timeHorizon)
##'
##' # ---------- competing risks --------------------
##'
##' # -- simulate an artificial data frame
##' # with competing cause response and three predictors
##' library(cmprsk)
##' library(riskRegression)
##' compRiskData <- function( N ){
##'   crdata <- SimCompRisk( N )
##'   crdata$X3 <- rnorm(N)
##'   crdata
##' }
##' set.seed(140515)
##' crdat <- compRiskData(1000)
##'
##' # -- fit a cause-specific Cox regression model
##' crForm <- Hist(time,event)~X2+X3
##' csCox  <- CSC(crForm, data=crdat)
##'
##' # -- choose a time horizon and plot the risk for a given cause
##' timeHorizon <- floor(median(crdat$time))
##' riskLevelPlot(csCox, crForm, data = crdat, horizon = timeHorizon, cause = 1)
##' 
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
riskLevelPlot <- function(object,
                          formula,
                          data = parent.frame(),
                          horizon=NULL,
                          cause=1,
                          ## colorkey=TRUE, Add to option? default now =TRUE
                          ...){
                                        # {{{ Call
    require(lattice)
    call <- match.call()
    formula.names <- try(all.names(formula),silent=TRUE)
                                        # }}}
                                        # {{{ Data 
    if(is.null(data)){
        try(data <- object$data)
        if(is.null(data)){
            stop("The data object is missing")
        }
    }
                                        # }}}
                                        # {{{ Response type 
    form <- formula
    if (formula.names[2] == as.name("Surv")){ # Use Hist() if censored data
        form[[2]][[1]] <- as.name("Hist")
    }
    m <- model.frame(form,data,na.action=na.fail)
    
    response <- stats::model.response(m)
                                        # }}}
                                        # {{{ Risk factors and grid sides
    formTerms <- terms(formula)
    riskfactors <- attr(formTerms,"term.labels")
    ## print(riskfactors)
    if (length(riskfactors) != 2)
        stop("Can currently only handle exactly 2 risk factors") 
    gridSides <- lapply(riskfactors,function(x){
        xVal <- data[[x]]
        sideValues <- getSeqForGrid(xVal,gridSize=25)
        ## if (is.factor(xVal) || unique(xVal) < 3) stop("Currently both risk factors must be continuous")
    })
                                        # }}}
                                        # {{{ Create newdata grid
    ## browser()
    newData <- expand.grid(gridSides)  
    names(newData) <- riskfactors
                                        # }}}
                                        # {{{ Predict grid probabilities z
    pr.args <- list(object=object,newdata=newData)
    ## if (horizon) times=horizon,cause=cause)
    z <- do.call(predictRisk,pr.args)
                                        # }}}  
                                        # ----------------------------copy/pasted----------------------------
                                        # {{{ Find the limits of z
    zLim <- NA # added
    if(is.na(zLim)){
        if(min(z) >= 0 && max(z) <= 1){
            zLim <- c(0,1)
        } else {
            zLim <- range(z)
        }
    }
    if(min(z) < min(zLim) || max(z) > max(zLim)){
        warning("The range of z (min=", min(z), ", max=", max(z), ") exceeds the range of the zLim argument provided (min=", min(zLim), ", max=", max(zLim), ")")
    }
    nColourValues <- 200 #added
    atVector <- seq(from=min(zLim),to=max(zLim),length.out=nColourValues+2)[-c(1,nColourValues+2)]
                                        # }}}
                                        # {{{ truncate extreme values to avoid missing colors
    if (any(z>0.995)){
        ## warning("predicted values above 0.995 are set to 0.995")
        z <- ifelse(z>0.995, 0.995, z)
    }
    if (any(z<0.005)){
        ## warning("predicted values below 0.005 are set to 0.005")
        z <- ifelse(z<0.005, 0.005, z)
    }
                                        # }}}
                                        # {{{ Include predicted values in newData
    newData$z <- z
                                        # }}}
                                        # {{{ Define default colors: rainbow color scheme
    colorVector <- rev(rainbow(nColourValues+1,end=0.6))
    colorkey <- TRUE
    colorkeyList <- colorkey
    if(colorkey == TRUE){
        if(all(zLim == c(0,1))){
            colorkeyList <- list(at = seq(from=min(zLim),to = max(zLim),length.out = nColourValues+2),
                                 col = colorVector,
                                 labels = list(at=(0:10)/10,labels = paste((0:10)*10,'%',sep='')))
        }
        else {
            colorkeyList <- list(at = seq(from=min(zLim),
                                          to=max(zLim),
                                          length.out=nColourValues+2),
                                 col = colorVector,
                                 labels = list(at=seq(from=min(zLim),
                                                      to=max(zLim),
                                                      length.out=10)))
        }
    }
                                        # }}}
                                        # --------------------------------------------------------------------
                                        # {{{ Heat map 
    plotObs <- TRUE
    levelplotForm <- as.formula(paste("z~",paste(riskfactors,collapse="+")))
    trellisObject <- lattice::levelplot(levelplotForm,
                                        newData,
                                        colorkey = colorkeyList,
                                        col.regions = colorVector,
                                        at=atVector,
                                        ...)
                                        # }}}
    return(trellisObject)
}

# {{{ getSeqForGrid
getSeqForGrid <- function(x,gridSize=25){
    if(is.factor(x)){  # Remember jitter (=noise) later
        xlev <- levels(x)
        sideLength <- length(xlev) 
        out <- factor(xlev)
    }
    else if(is.logical(x)) { # Remember jitter (=noise) later
        out <- unique(x)    
        sideLength <- length(out)    
    }
    else {
        sideLength <- gridSize
        x.range <- range(x,finite=TRUE)
        step <- (max(x.range)-min(x.range))/sideLength
        out <- seq(min(x.range),max(x.range),step)
    }
    return(out)
}
# }}}


# {{{ Not included yet-- need to be? : panelFunction

panelFunction <- function(...){
    ## browser()
    panel.levelplot(...)
    currentRow <- current.row()
    currentColumn <- current.column()
    if(list(...)$passMeToPanelFunction$plotObservations == TRUE){
        ## browser()
        data <- list(...)$passMeToPanelFunction$data
        ## innerX <- list(...)$passMeToPanelFunction$innerX
        ## innerY <- list(...)$passMeToPanelFunction$innerY
        ## outerX <- list(...)$passMeToPanelFunction$outerX
        ## outerY <- list(...)$passMeToPanelFunction$outerY
        ## pchObs <- list(...)$passMeToPanelFunction$pchObs
        ## pchCol <- list(...)$passMeToPanelFunction$pchCol
        ## xLabels <- list(...)$passMeToPanelFunction$xLabels
        ## yLabels <- list(...)$passMeToPanelFunction$yLabels
        ## jitterAmounts <- list(...)$passMeToPanelFunction$jitterAmounts
        ## jitterHow <- list(...)$passMeToPanelFunction$jitterHow
        groups <- list(...)$groups
        subscripts <- list(...)$subscripts

    #
    ## if(is.na(outerX) && is.na(outerY)){
      ## whichValues <- 1:nrow(data)
    ## } else if(!is.na(outerX) && is.na(outerY)){
      ## whichValues <- which(data[[outerX]] == sapply(strsplit(as.character(groups[subscripts[1]]),' = '),function(x){x[2]}))
    ## } else if(is.na(outerX) && !is.na(outerY)){
      ## whichValues <- which(data[[outerY]] == sapply(strsplit(as.character(groups[subscripts[1]]),' = '),function(x){x[2]}))
    ## } else {
      ## whichValues <- which(as.character(data[[outerX]]) == sort(levels(data[[outerX]]))[currentColumn] &
                           ## as.character(data[[outerY]]) == sort(levels(data[[outerY]]))[currentRow])
    ## }
    ## if(length(pchObs) > 1){
      ## pchObs <- pchObs[whichValues]
    ## } else if (length(pchObs) == 0){
      ## pchObs <- 1
    ## }
#
#    
    ## if(jitterHow == ''){
      ## lpoints(data[[innerX]][whichValues],data[[innerY]][whichValues],pch=pchObs,col=pchCol)
    ## } else if(jitterHow == 'x') {
      ## lpoints(jitter(match(as.character(data[[innerX]][whichValues]),xLabels),amount = jitterAmounts[1]),
              ## data[[innerY]][whichValues],pch=pchObs,col=pchCol)
    ## } else if(jitterHow == 'y') {
      ## lpoints(data[[innerX]][whichValues],
              ## jitter(match(as.character(data[[innerY]][whichValues]),yLabels),amount = jitterAmounts[2]),pch=pchObs,col=pchCol)
    ## } else {
      ## lpoints(jitter(match(as.character(data[[innerX]][whichValues]),xLabels),amount = jitterAmounts[1]),
              ## jitter(match(as.character(data[[innerY]][whichValues]),yLabels),amount = jitterAmounts[2]),
              ## pch=pchObs,col=pchCol)
    ## }
  }
  
}
# }}}
