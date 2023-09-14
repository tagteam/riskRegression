getNullModel <- function(formula, data, response.type) {
  nullform <- stats::reformulate("1", response = formula[[2]])
  nullfit <- switch(response.type,
    "binary" = {
      list("Null model" = stats::glm(nullform, data, family = "binomial"))
    },
    "continuous" = {
      stop("FIXME: Need to adapt a function for the empirical distribution function or glm.")
    },
    "survival" = {
      list("Null model" = prodlim::prodlim(nullform, data))
    },
    "competing.risks" = {
      list("Null model" = prodlim::prodlim(nullform, data))
    }
  )
  nullfit[[1]]$call$data <- NULL
  nullfit[[1]]$formula <- NULL
  nullfit[[1]]$call$formula <- nullform
  nullfit
}
