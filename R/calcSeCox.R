### calcSeCox.R ---
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: maj 27 2017 (11:46)
## Version:
## last-updated: okt 22 2021 (18:26)
##           By: Brice Ozenne
##     Update #: 860
#----------------------------------------------------------------------
##
### Commentary:
##
### Change Log:
#----------------------------------------------------------------------
##
### Code:


## * calcSeCox (documentation)
#' @title Computation of standard errors for predictions
#' @description Compute the standard error associated to the predictions from Cox regression model
#' using a first order von Mises expansion of the functional (cumulative hazard or survival).
#' @name calcSeCox
#'
#' @param object The fitted Cox regression model object either
#'     obtained with \code{coxph} (survival package) or \code{cph}
#'     (rms package).
#' @param times Vector of times at which to return the estimated
#'      hazard/survival.
#' @param nTimes the length of the argument \code{times}.
#' @param type One or several strings that match (either in lower or upper case or mixtures) one
#' or several of the strings \code{"hazard"},\code{"cumhazard"}, \code{"survival"}.
#' @param diag [logical] when \code{FALSE} the hazard/cumlative hazard/survival for all observations at all times is computed,
#' otherwise it is only computed for the i-th observation at the i-th time.
#' @param Lambda0 the baseline hazard estimate returned by \code{BaseHazStrata_cpp}.
#' @param object.n the number of observations in the dataset used to estimate the object.
#' @param object.time the time to event of the observations used to estimate the object.
#' @param object.eXb the exponential of the linear predictor relative to the observations used to estimate the object.
#' @param object.strata the strata index of the observations used to estimate the object.
#' @param nStrata the number of strata.
#' @param new.n the number of observations for which the prediction was performed.
#' @param new.eXb the linear predictor evaluated for the new observations.
#' @param new.LPdata the variables involved in the linear predictor for the new observations.
#' @param new.strata the strata indicator for the new observations.
#' @param new.survival the survival evaluated for the new observations.
#' @param nVar.lp the number of variables that form the linear predictor.
#' @param export can be "iid" to return the value of the influence function for each observation.
#'                      "se" to return the standard error for a given timepoint.
#'                      "average.iid" to return the value of the average influence function over the observations for which the prediction was performed.
#'
#' @param store.iid Implementation used to estimate the influence function and the standard error.
#' Can be \code{"full"} or \code{"minimal"}. See the details section.
#'
#' @details \code{store.iid="full"} compute the influence function for each observation at each time in the argument \code{times}
#' before computing the standard error / influence functions.
#' \code{store.iid="minimal"} recompute for each subject specific prediction the influence function for the baseline hazard.
#' This avoid to store all the influence functions but may lead to repeated evaluation of the influence function.
#' This solution is therefore more efficient in memory usage but may not be in terms of computation time.
#'
# #' @inheritParams  predict.CauseSpecificCox
#'
#' @author Brice Ozenne broz@@sund.ku.dk, Thomas A. Gerds tag@@biostat.ku.dk
#'
#' @return A list optionally containing the standard error for the survival, cumulative hazard and hazard.
#'

## * calcSeCox (code)
#' @rdname calcSeCox
calcSeCox <- function(object, times, nTimes, type, diag,
                      Lambda0, object.n, object.time, object.eXb, object.strata, nStrata,
                      new.n, new.eXb, new.LPdata, new.strata, new.survival,
                      nVar.lp, export, store.iid) {
  ## ** Computation of the influence function
  if (is.iidCox(object)) {
    store.iid <- object$iid$store.iid
    iid.object <- selectJump(object$iid, times = times, type = type)
  } else {
    iid.object <- iidCox(object, tau.hazard = times, store.iid = store.iid, return.object = FALSE)
  }

  ## ** Prepare arguments
  if (diag) {
    nTimes <- 1
  }
  new.strata <- as.numeric(new.strata)

  Lambda0$strata <- as.numeric(Lambda0$strata)
  if ("hazard" %in% type) {
    Lambda0$hazard <- lapply(1:nStrata, function(s) {
      Lambda0$hazard[Lambda0$strata == s][Lambda0$oorder.times]
    })
  }
  if ("cumhazard" %in% type || "survival" %in% type) {
    Lambda0$cumhazard <- lapply(1:nStrata, function(s) {
      Lambda0$cumhazard[Lambda0$strata == s][Lambda0$oorder.times]
    })
  }

  if (is.null(attr(export, "factor"))) {
    rm.list <- TRUE
    factor <- list(matrix(1, nrow = new.n, ncol = nTimes))
  } else {
    rm.list <- FALSE
    factor <- attr(export, "factor")
  }
  out <- list()
  if ("se" %in% export) {
    if ("cumhazard" %in% type) {
      out$cumhazard.se <- matrix(NA, nrow = new.n, ncol = nTimes)
    }
    if ("survival" %in% type) {
      out$survival.se <- matrix(NA, nrow = new.n, ncol = nTimes)
    }
  }
  if ("iid" %in% export) {
    if ("hazard" %in% type) {
      out$hazard.iid <- array(NA, dim = c(object.n, nTimes, new.n))
    }
    if ("cumhazard" %in% type) {
      out$cumhazard.iid <- array(NA, dim = c(object.n, nTimes, new.n))
    }
    if ("survival" %in% type) {
      out$survival.iid <- array(NA, dim = c(object.n, nTimes, new.n))
    }
  }
  if ("average.iid" %in% export) {
    if ("cumhazard" %in% type) {
      out$cumhazard.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)
    }
    if ("survival" %in% type) {
      out$survival.average.iid <- matrix(0, nrow = object.n, ncol = nTimes)
    }
  }

  ## ** Linear predictor
  if ("lp" %in% type || (store.iid[[1]] == "full" && nVar.lp > 0 && ("iid" %in% export || "se" %in% export))) {
    X_IFbeta_mat <- tcrossprod(iid.object$IFbeta, new.LPdata)
  }

  if ("lp" %in% type) {
    if ("iid" %in% export) {
      if (nVar.lp > 0) {
        out$lp.iid <- X_IFbeta_mat
      } else {
        out$lp.iid <- matrix(0, nrow = object.n, ncol = new.n)
      }
    }
    if ("se" %in% export) {
      if (nVar.lp > 0) {
        out$lp.se <- cbind(colSums(X_IFbeta_mat^2))
      } else {
        out$lp.se <- matrix(0, nrow = length(new.n), ncol = 1)
      }
    }
    if ("average.iid" %in% export) {
      if (nVar.lp > 0) {
        out$lp.average.iid <- rowMeans(X_IFbeta_mat)
      } else {
        out$lp.average.iid <- rep(0, times = length(object.n))
      }
    }

    if (length(type) == 1) { ## do not look at anything else like survival or hazard
      return(out)
    }
  }

  ## ** hazard / cumulative hazard / survival

  if (store.iid[[1]] == "minimal") {
    ## *** method 1: minimal storage of the influence function
    resCpp <- calcSeMinimalCox_cpp(
      seqTau = times,
      newSurvival = if ("survival" %in% type) {
        new.survival
      } else {
        new.survival <- matrix(NA)
      },
      hazard0 = if ("hazard" %in% type) {
        Lambda0$hazard
      } else {
        list(NA)
      },
      cumhazard0 = if ("cumhazard" %in% type || "survival" %in% type) {
        Lambda0$cumhazard
      } else {
        list(NA)
      },
      newX = new.LPdata,
      neweXb = new.eXb,
      IFbeta = iid.object$IFbeta,
      Ehazard0 = iid.object$calcIFhazard$Elambda0,
      cumEhazard0 = iid.object$calcIFhazard$cumElambda0,
      hazard_iS0 = iid.object$calcIFhazard$lambda0_iS0,
      cumhazard_iS0 = iid.object$calcIFhazard$cumLambda0_iS0,
      delta_iS0 = iid.object$calcIFhazard$delta_iS0,
      sample_eXb = iid.object$calcIFhazard$eXb,
      sample_time = iid.object$obstime,
      indexJumpSample_time = lapply(
        iid.object$calcIFhazard$time1,
        function(iTime) {
          prodlim::sindex(jump.times = iTime, eval.times = iid.object$obstime) - 1
        }
      ),
      jump_time = iid.object$calcIFhazard$time1,
      indexJumpTau = lapply(
        iid.object$calcIFhazard$time1,
        function(iTime) {
          prodlim::sindex(jump.times = iTime, eval.times = times) - 1
        }
      ),
      lastSampleTime = iid.object$etime.max,
      newdata_index = lapply(
        1:nStrata,
        function(iS) {
          which(new.strata == iS) - 1
        }
      ),
      factor = factor,
      nTau = nTimes, nNewObs = new.n, nSample = object.n, nStrata = nStrata, p = nVar.lp,
      diag = diag, exportSE = "se" %in% export, exportIF = "iid" %in% export, exportIFmean = "average.iid" %in% export,
      exportHazard = "hazard" %in% type, exportCumhazard = "cumhazard" %in% type, exportSurvival = "survival" %in% type,
      debug = 0
    )

    if ("iid" %in% export) {
      if ("hazard" %in% type) {
        out$hazard.iid <- aperm(resCpp$IF_hazard, perm = c(1, 3, 2))
      }
      if ("cumhazard" %in% type) {
        out$cumhazard.iid <- aperm(resCpp$IF_cumhazard, perm = c(1, 3, 2))
      }
      if ("survival" %in% type) {
        out$survival.iid <- aperm(resCpp$IF_survival, perm = c(1, 3, 2))
      }
    }
    if ("se" %in% export) {
      if ("cumhazard" %in% type) {
        out$cumhazard.se <- resCpp$SE_cumhazard
      }
      if ("survival" %in% type) {
        out$survival.se <- resCpp$SE_survival
      }
    }
    if ("average.iid" %in% export) { ## average over strata
      if (rm.list) {
        if ("hazard" %in% type) {
          out$hazard.average.iid <- matrix(resCpp$IFmean_hazard[[1]], nrow = object.n, ncol = nTimes)
        }
        if ("cumhazard" %in% type) {
          out$cumhazard.average.iid <- matrix(resCpp$IFmean_cumhazard[[1]], nrow = object.n, ncol = nTimes)
        }
        if ("survival" %in% type) {
          out$survival.average.iid <- matrix(resCpp$IFmean_survival[[1]], nrow = object.n, ncol = nTimes)
        }
      } else {
        if ("hazard" %in% type) {
          out$hazard.average.iid <- lapply(resCpp$IFmean_hazard, function(iVec) {
            matrix(iVec, nrow = object.n, ncol = nTimes)
          })
        }
        if ("cumhazard" %in% type) {
          out$cumhazard.average.iid <- lapply(resCpp$IFmean_cumhazard, function(iVec) {
            matrix(iVec, nrow = object.n, ncol = nTimes)
          })
        }
        if ("survival" %in% type) {
          out$survival.average.iid <- lapply(resCpp$IFmean_survival, function(iVec) {
            matrix(iVec, nrow = object.n, ncol = nTimes)
          })
        }
      }
    }

    # }}}
  } else if ("iid" %in% export || "se" %in% export) {
    ## *** method 2: using the influence function of the baseline hazard/baseline cumulative hazard
    if (diag || (nVar.lp == 0)) {
      for (iStrata in 1:nStrata) { ## iStrata <- 1
        indexStrata <- which(new.strata == iStrata)
        if (length(indexStrata) == 0) {
          next
        }
        iPrevalence <- length(indexStrata) / new.n

        ## compute iid
        if ("hazard" %in% type) {
          if (diag) {
            if (nVar.lp == 0) {
              iIFhazard <- iid.object$IFhazard[[iStrata]][, indexStrata, drop = FALSE]
            } else {
              iIFhazard <- rowMultiply_cpp(
                iid.object$IFhazard[[iStrata]][, indexStrata, drop = FALSE] + rowMultiply_cpp(X_IFbeta_mat[, indexStrata, drop = FALSE],
                  scale = Lambda0$hazard[[iStrata]][indexStrata]
                ),
                scale = new.eXb[indexStrata]
              )
            }
          } else {
            ## nVar.lp==0
            iIFhazard <- iid.object$IFhazard[[iStrata]]
            ## tiIFhazard <- t(iIFhazard)
          }
        }
        if ("cumhazard" %in% type || "survival" %in% type) {
          if (diag) {
            if (nVar.lp == 0) {
              iIFcumhazard <- iid.object$IFcumhazard[[iStrata]][, indexStrata, drop = FALSE]
            } else {
              iIFcumhazard <- rowMultiply_cpp(iid.object$IFcumhazard[[iStrata]][, indexStrata, drop = FALSE] + rowMultiply_cpp(X_IFbeta_mat[, indexStrata, drop = FALSE], scale = Lambda0$cumhazard[[iStrata]][indexStrata]),
                scale = new.eXb[indexStrata]
              )
            }

            if ("survival" %in% type && ("iid" %in% export || "average.iid" %in% export)) {
              iIFsurvival <- rowMultiply_cpp(-iIFcumhazard, scale = new.survival[indexStrata, ])
            }
          } else {
            ## nVar.lp == 0
            iIFcumhazard <- iid.object$IFcumhazard[[iStrata]]
            ## tiIFcumhazard <- t(iIFcumhazard)

            if ("survival" %in% type && ("iid" %in% export || "average.iid" %in% export)) {
              iIFsurvival <- rowMultiply_cpp(-iIFcumhazard, scale = new.survival[indexStrata[1], ]) ## everybody has the same survival
              ## tiIFsurvival <- t(iIFsurvival)
            }
          }
        }

        ## export
        if (diag) {
          if ("iid" %in% export) {
            if ("hazard" %in% type) {
              out$hazard.iid[, 1, indexStrata] <- iIFhazard
            }
            if ("cumhazard" %in% type) {
              out$cumhazard.iid[, 1, indexStrata] <- iIFcumhazard
            }
            if ("survival" %in% type) {
              out$survival.iid[, 1, indexStrata] <- iIFsurvival
            }
          }
          if ("se" %in% export) {
            iSEcumhazard <- sqrt(colSums(iIFcumhazard^2))
            if ("cumhazard" %in% type) {
              out$cumhazard.se[indexStrata, 1] <- iSEcumhazard
            }
            if ("survival" %in% type) {
              out$survival.se[indexStrata, 1] <- iSEcumhazard * new.survival[indexStrata, 1]
            }
          }
          if ("average.iid" %in% export) {
            if ("hazard" %in% type) {
              out$hazard.average.iid[, 1] <- out$hazard.average.iid[, 1] + rowSums(iIFhazard) / new.n
            }
            if ("cumhazard" %in% type) {
              out$cumhazard.average.iid[, 1] <- out$cumhazard.average.iid[, 1] + rowSums(iIFcumhazard) / new.n
            }
            if ("survival" %in% type) {
              out$survival.average.iid[, 1] <- out$survival.average.iid[, 1] + rowSums(iIFsurvival) / new.n
            }
          }
        } else { ## nVar.lp==0

          if ("se" %in% export) {
            iSEcumhazard <- sqrt(colSums(iIFcumhazard^2))
            if ("survival" %in% type) {
              iSEsurvival <- iSEcumhazard * new.survival[indexStrata[1], ] ## everybody has the same survival
            }
          }

          for (iObs in indexStrata) { ## iObs <- 1
            if ("iid" %in% export) {
              if ("hazard" %in% type) {
                out$hazard.iid[, , iObs] <- iIFhazard
              }
              if ("cumhazard" %in% type) {
                out$cumhazard.iid[, , iObs] <- iIFcumhazard
              }
              if ("survival" %in% type) {
                out$survival.iid[, , iObs] <- iIFsurvival
              }
            }
            if ("se" %in% export) {
              if ("cumhazard" %in% type) {
                out$cumhazard.se[iObs, ] <- iSEcumhazard
              }
              if ("survival" %in% type) {
                out$survival.se[iObs, ] <- iSEsurvival
              }
            }
          }
          if ("average.iid" %in% export) {
            if ("hazard" %in% type) {
              out$hazard.average.iid <- out$hazard.average.iid + iIFhazard * iPrevalence
            }
            if ("cumhazard" %in% type) {
              out$cumhazard.average.iid <- out$cumhazard.average.iid + iIFcumhazard * iPrevalence
            }
            if ("survival" %in% type) {
              out$survival.average.iid <- out$survival.average.iid + iIFsurvival * iPrevalence
            }
          }
        }
      }
    } else { ## nVar.lp > 0
      for (iObs in 1:new.n) { ## iObs <- 1
        iObs.strata <- new.strata[iObs]

        ## compute iid
        if ("hazard" %in% type) {
          iIFhazard <- (new.eXb[iObs] * (iid.object$IFhazard[[iObs.strata]] + crossprod(t(X_IFbeta_mat[, iObs, drop = FALSE]), Lambda0$hazard[[iObs.strata]])))
        }
        if ("cumhazard" %in% type || "survival" %in% type) {
          iIFcumhazard <- new.eXb[iObs] * (iid.object$IFcumhazard[[iObs.strata]] + crossprod(t(X_IFbeta_mat[, iObs, drop = FALSE]), Lambda0$cumhazard[[iObs.strata]]))
        }
        if ("survival" %in% type && ("iid" %in% export || "average.iid" %in% export)) {
          iIFsurvival <- rowMultiply_cpp(-iIFcumhazard, scale = new.survival[iObs, ])
        }

        ## export
        if ("iid" %in% export) {
          if ("hazard" %in% type) {
            out$hazard.iid[, , iObs] <- iIFhazard
          }
          if ("cumhazard" %in% type) {
            out$cumhazard.iid[, , iObs] <- iIFcumhazard
          }
          if ("survival" %in% type) {
            out$survival.iid[, , iObs] <- iIFsurvival
          }
        }
        if ("se" %in% export) {
          iSEcumhazard <- sqrt(colSums(iIFcumhazard^2))
          if ("cumhazard" %in% type) {
            out$cumhazard.se[iObs, ] <- iSEcumhazard
          }
          if ("survival" %in% type) {
            out$survival.se[iObs, ] <- iSEcumhazard * new.survival[iObs, , drop = FALSE]
          }
        }
        if ("average.iid" %in% export) { ## average over observations
          if ("hazard" %in% type) {
            out$hazard.average.iid <- out$hazard.average.iid + iIFhazard / new.n
          }
          if ("cumhazard" %in% type) {
            out$cumhazard.average.iid <- out$cumhazard.average.iid + iIFcumhazard / new.n
          }
          if ("survival" %in% type) {
            out$survival.average.iid <- out$survival.average.iid + iIFsurvival / new.n
          }
        }
      }
    }
  } else if ("average.iid" %in% export) { ## fast average over observations
    ## *** method 3: computation of the average influence function

    ## prepare strata
    new.Ustrata <- sort(unique(new.strata))
    new.nStrata <- length(new.Ustrata)

    new.indexStrata <- lapply(new.Ustrata, function(iStrata) {
      which(new.strata == iStrata) - 1
    })
    new.prevStrata <- sapply(new.indexStrata, length) / new.n

    ## normalize arguments for C++
    attr(new.LPdata, "levels") <- NULL
    if (is.null(new.survival)) {
      new.survival <- matrix()
    }

    ## C++
    if ("hazard" %in% type) {
      outRcpp.hazard <- calcAIFsurv_cpp(
        ls_IFcumhazard = iid.object$IFhazard[new.Ustrata],
        IFbeta = iid.object$IFbeta,
        cumhazard0 = Lambda0$hazard[new.Ustrata],
        survival = matrix(0),
        eXb = new.eXb,
        X = new.LPdata,
        prevStrata = new.prevStrata,
        ls_indexStrata = new.indexStrata,
        factor = factor,
        nTimes = nTimes,
        nObs = object.n,
        nStrata = new.nStrata,
        nVar = nVar.lp,
        diag = diag,
        exportCumHazard = TRUE,
        exportSurvival = FALSE
      )
    }
    if (("cumhazard" %in% type) || ("survival" %in% type)) {
      outRcpp.cumhazard <- calcAIFsurv_cpp(
        ls_IFcumhazard = iid.object$IFcumhazard[new.Ustrata],
        IFbeta = iid.object$IFbeta,
        cumhazard0 = Lambda0$cumhazard[new.Ustrata],
        survival = new.survival,
        eXb = new.eXb,
        X = new.LPdata,
        prevStrata = new.prevStrata,
        ls_indexStrata = new.indexStrata,
        factor = factor,
        nTimes = nTimes,
        nObs = object.n,
        nStrata = new.nStrata,
        nVar = nVar.lp,
        diag = diag,
        exportCumHazard = ("cumhazard" %in% type),
        exportSurvival = ("survival" %in% type)
      )
    }

    ## reshape
    if ("hazard" %in% type) {
      if (rm.list) {
        out$hazard.average.iid <- matrix(outRcpp.hazard[[1]][[1]], nrow = object.n, ncol = nTimes)
      } else {
        out$hazard.average.iid <- lapply(outRcpp.hazard[[1]], function(iMat) {
          matrix(iMat, nrow = object.n, ncol = nTimes)
        })
      }
    }
    if ("cumhazard" %in% type) {
      if (rm.list) {
        out$cumhazard.average.iid <- matrix(outRcpp.cumhazard[[1]][[1]], nrow = object.n, ncol = nTimes)
      } else {
        out$cumhazard.average.iid <- lapply(outRcpp.cumhazard[[1]], function(iMat) {
          matrix(iMat, nrow = object.n, ncol = nTimes)
        })
      }
    }
    if ("survival" %in% type) {
      if (rm.list) {
        out$survival.average.iid <- matrix(outRcpp.cumhazard[[2]][[1]], nrow = object.n, ncol = nTimes)
      } else {
        out$survival.average.iid <- lapply(outRcpp.cumhazard[[2]], function(iMat) {
          matrix(iMat, nrow = object.n, ncol = nTimes)
        })
      }
    }
  }

  ## export
  return(out)
}

## * selectJump
#' @title Evaluate the influence function at selected times
#'
#' @description Evaluate the influence function at selected times
#' @param IF influence function returned by iidCox
#' @param times the times at which the influence function should be assessed
#' @param type can be \code{"hazard"} or/and \code{"cumhazard"}.
#'
#' @author Brice Ozenne broz@@sund.ku.dk
#'
#' @return An object with the same dimensions as IF
#'
selectJump <- function(IF, times, type) {
  if (any(times < 0)) {
    warning("selectJump may not handle correctly negative times")
  }

  nStrata <- length(IF$time)
  for (iStrata in 1:nStrata) {
    if (IF$store.iid == "minimal") {
      isJump <- times %in% IF$time[[iStrata]]
      indexJump <- prodlim::sindex(jump.times = c(0, IF$time[[iStrata]]), eval.times = times)
      if (NROW(IF$calcIFhazard$Elambda0[[iStrata]]) > 0) {
        IF$calcIFhazard$Elambda0[[iStrata]] <- rowMultiply_cpp(cbind(0, IF$calcIFhazard$Elambda0[[iStrata]])[, indexJump, drop = FALSE], scale = isJump)
      } else {
        IF$calcIFhazard$Elambda0[[iStrata]] <- matrix(NA, nrow = 0, ncol = length(isJump))
      }
      if (NROW(IF$calcIFhazard$cumElambda0[[iStrata]]) > 0) {
        IF$calcIFhazard$cumElambda0[[iStrata]] <- cbind(0, IF$calcIFhazard$cumElambda0[[iStrata]])[, indexJump, drop = FALSE]
      } else {
        IF$calcIFhazard$cumElambda0[[iStrata]] <- matrix(NA, nrow = 0, ncol = length(isJump))
      }
      IF$calcIFhazard$lambda0_iS0[[iStrata]] <- IF$calcIFhazard$lambda0_iS0[[iStrata]] * (IF$calcIFhazard$time1[[iStrata]] <= max(times))
      IF$calcIFhazard$cumLambda0_iS0[[iStrata]] <- IF$calcIFhazard$cumLambda0_iS0[[iStrata]] * (IF$calcIFhazard$time1[[iStrata]] <= max(times))
    } else {
      if ("hazard" %in% type) {
        match.times <- match(times, table = IF$time[[iStrata]])
        match.times[is.na(match.times)] <- 0
        if (any(times > IF$etime.max[[iStrata]])) {
          match.times[times > IF$etime.max[[iStrata]]] <- NA
        }
        IF$IFhazard[[iStrata]] <- subsetIndex(IF$IFhazard[[iStrata]], index = match.times, default = 0, col = TRUE)
      }
      if ("cumhazard" %in% type || "survival" %in% type) {
        indexJump <- prodlim::sindex(jump.times = IF$time[[iStrata]], eval.times = times)
        if (any(times > IF$etime.max[[iStrata]])) {
          indexJump[times > IF$etime.max[[iStrata]]] <- NA
        }
        IF$IFcumhazard[[iStrata]] <- subsetIndex(IF$IFcumhazard[[iStrata]], index = indexJump, default = 0, col = TRUE)
      }
    }
    IF$time[[iStrata]] <- times
  }

  return(IF)
}

#----------------------------------------------------------------------
### calcSeCox.R ends here
