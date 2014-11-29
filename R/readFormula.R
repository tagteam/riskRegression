readFormula <- function(formula,
                        specials,
                        specialArgumentNames,
                        alias,
                        unspecified="unSpec"){
    # {{{ check formula and convert to character
    ff <- as.character(formula)
    # }}}
    # {{{ extract the response = lhs
    lhs <- ff[[2]]
    browser()
    if (length(grep(":|\\*",lhs))>0)
        stop("Interaction terms are not allowed")
    rForm <- formula(paste(lhs,"~1",sep=""))
    formList <- list(Response=rForm)
    # }}}
    # {{{ split the rhs into terms
    rhs <- ff[[3]]
    terms <- strsplit(rhs,"\\+|\\-")[[1]]
    terms <- sapply(terms,function(tt){ ## remove whitespace
        gsub(" ","",tt)
    })
    # }}}
    # {{{ extract the intercept
    numterms <- grep("^[0-9]+$",terms,value=TRUE)
    switch(as.character(length(numterms)),
           "0"={
               intercept <- 1
           },
           "1"={
               sign <- ifelse(length(grep(paste("\\-[ ]*",numterms,sep=""),rhs))>0,-1,1)
               intercept <- sign*as.numeric(numterms)
           },
           "2"={
               stop("Error in readFormula:\nMultiple intercept terms in formula: ",rhs,call.=FALSE)
           })
    formList <- c(formList,list(Intercept=intercept))
    # }}}
    # {{{ find and check special terms
    names(specials) <- specials
    listTerms <- strsplit(terms,"[()]")
    ##   listTerms <- listTerms[-grep("^[0-9]+$",terms,value=FALSE)]
    isSpecial <- sapply(listTerms,length)
    specialNames <- sapply(listTerms[isSpecial==2],function(x){x[[1]]})
    # {{{ resolve alias names
    if (!missing(alias)){
        isAlias <- match(specialNames,names(alias),nomatch=0)
        if (any(isAlias>0)){
            specialNames[isAlias!=0] <- alias[isAlias]
        }}
    # }}}
    if (any(isSpecial>2))
        stop("Mispecified formula: constructions like 'const(factor(sex))' are not supported")
    if (any(notFound <- match(specialNames,specials,nomatch=0)==0))
        stop("Mispecified formula: special(s) '",
             paste(specialNames[notFound],collapse=" and "),
             "' ",
             ifelse(length(notFound)>1,"are","is"),
             " not supported")
    # }}}
    # {{{ find special variable names and extra arguments
    listTermsWithArguments <- unlist(lapply(listTerms[isSpecial==2],function(x)strsplit(x[[2]],",")),recursive=FALSE)
    specialVarnames <- sapply(listTermsWithArguments,function(x){x[[1]]})
    if (length(who <- grep("=",specialVarnames,value=TRUE))>0)
        stop("Problematic variable name '",who,"'")
    specialArguments <- lapply(listTermsWithArguments,function(x){
        if (length(x)==1) out <- NULL else out <- x[2:length(x)]
    })
    names(specialArguments) <- specialVarnames
    if (missing(specialArgumentNames)) {
        specialArgumentNames <- lapply(specials,function(x)NULL)
    }
    else {
        if (!(all(match(names(specialArgumentNames),specials,nomatch=0))))
            stop("Mispecified argument specialArgumentNames")
    }
    if (length(specialArguments)>0){
        specialArgumentList <- lapply(1:length(specialArguments),function(i){
            args <- specialArguments[[i]]
            if (!is.null(args)){
                fullvalue <- strsplit(args,"=")
                fullvalue <- lapply(fullvalue,function(x){ ## remove whitespace
                    gsub(" ","",x)
                })
                givennames <- sapply(fullvalue,function(x){
                    if (length(x)==1)
                        ""
                    else
                        x[[1]]
                })
                values <- lapply(fullvalue,function(x){
                    if (length(x)==1)
                        x[[1]]
                    else
                        x[[2]]
                })
                specName <- specialNames[[i]]
                if (is.null(specialArgumentNames[[specName]]))
                    wantednames <- paste("Arg",1:length(args),sep=".")
                else{
                    wantednames <- specialArgumentNames[[specName]]
                    if(length(wantednames)<length(args)) stop("Too many arguments for special function ",specName)
                }
                realnames <- givennames[givennames!=""]
                thismatch <- match(realnames,wantednames,nomatch=0)
                if (length(realnames)>0)
                    if (!all(thismatch))
                        stop("Argument(s) ",realnames," is not an argument of  ",specName, ". Allowed are\n\n",paste(wantednames,collapse=", "))
                names(values) <- givennames
                nadd <- length(wantednames)-length(values)
                if (nadd>0){
                    values <- c(values,rep(NA,nadd))
                }
                thatmatch <- match(wantednames,names(values),nomatch=0)
                names(values)[names(values)==""] <- wantednames[thatmatch==0]
                values
            }
            else NULL
        })
        names(specialArgumentList) <- names(specialArguments)
    }
    # }}}
    # {{{ make new formulas for all specials
    specialList <- lapply(specials,function(spec){
        found <- specialNames==spec
        if (any(found)){
            sform <- formula(paste("~",paste(specialVarnames[found],collapse="+")))
            if (length(specialArguments[found])>0){
                list(formula=sform,specialArguments=specialArgumentList[found])
            }
            else{
                list(formula=sform)
            }
        }
        else NULL
    })
    formList <- c(formList,specialList)

  # }}}
  # {{{ collect the unspecial variables
  rhsVars <- specialVarnames
  if (any(isSpecial==1)){
    unSpecVars <- sapply(listTerms[isSpecial==1],function(x){x[[1]]})
    rhsVars <- c(rhsVars,unSpecVars)
    if (length(formList[[unspecified]])>0){
      oldform <- formList[[unspecified]]$formula
      formList[[unspecified]]$formula <- update.formula(oldform,paste("~ . + ",paste(unSpecVars,collapse="+")))
    }
    else{
      unSpecForm <- formula(paste("~",paste(unSpecVars,collapse="+")))
      unSpecEl <- list(list(formula=unSpecForm,specialArguments=NULL))
      names(unSpecEl) <- unspecified
      formList <- c(formList[-match(unspecified,names(formList),nomatch=length(formList)+1)],unSpecEl)
    }
  }
  # }}}  
  # {{{ collect all variables
  allForm <- formula(paste("~",paste(c(all.vars(rForm),rhsVars),collapse="+")))
  formList <- c(formList,list(allVars=allForm))
  # }}}
  formList
}
