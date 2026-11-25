library(survival)
library(prodlim)

test_that("default IBS preserves Dint endpoint and missing-value behavior", {
    integrate <- riskRegression:::integratedBrier
    expect_equal(integrate(c(1, 3), c(1/4, 5/16)), c(0, 1/6))
    expect_equal(integrate(c(0, 1, 3), c(1/16, 1/4, 5/16)),
                 c(0, 1/16, 3/16))
    expect_equal(integrate(c(1, 3, 4), c(1/4, NA, 1/2)), c(0, 1/6, NA))
    expect_equal(integrate(c(0, 1), c(NA, 1/4)), c(0, NA))
    expect_equal(integrate(0, NA_real_), 0)
    expect_equal(integrate(2, 1/4), 0)
    expect_equal(integrate(2, 1/4, "trapezoidal"), 1/8)
    expect_equal(integrate(c(1, 3, 4), c(1/4, NA, 1/2), "trapezoidal"),
                 c(1/8, NA, NA))
})

test_that("Score defaults to the original left-endpoint integration", {
    d <- data.frame(time = c(1, 2, 4, 5), event = c(1, 1, 1, 0))
    risks <- matrix(rep(c(0.25, 0.5, 0.75), each = nrow(d)), nrow = nrow(d))
    args <- list(object = list(model = risks), formula = Surv(time, event) ~ 1,
                 data = d, times = c(0, 1, 3), metrics = "brier", summary = "ibs",
                 null.model = FALSE, se.fit = FALSE, verbose = FALSE)
    default <- do.call(Score, args)$Brier$score
    left <- do.call(Score, c(args, list(ibs.method = "left")))$Brier$score
    trap <- do.call(Score, c(args, list(ibs.method = "trapezoidal")))$Brier$score
    expect_equal(default$IBS, c(0, 1/16, 3/16))
    expect_equal(default$IBS, left$IBS)
    expect_equal(trap$IBS, c(0, 5/32, 23/96))
    expect_equal(default$Brier, trap$Brier)
    expect_error(do.call(Score, c(args, list(ibs.method = "unknown"))), "arg.*should be")
})

test_that("sparse-grid IBS warnings respect method, horizon, and verbosity", {
    d <- data.frame(time = c(1, 2, 4, 5), event = c(1, 1, 1, 0))
    score <- function(times, data = d, ...) {
        Score(list(model = matrix(0.2, nrow(data), length(times))),
              formula = Surv(time, event) ~ 1, data = data, times = times,
              metrics = "brier", null.model = FALSE, se.fit = FALSE, ...)
    }
    expect_warning(score(c(1, 3), summary = "ibs"), 'ibs.method = "trapezoidal"')
    expect_warning(score(3, summary = "ibs"), "Dint")
    expect_warning(score(c(1, 3), summary = "crps"), "Dint")
    for (quiet in list(FALSE, 0, -1)) {
        expect_warning(score(c(1, 3), summary = "ibs", verbose = quiet), NA)
    }
    expect_warning(score(c(1, 3), summary = "ibs", ibs.method = "trapezoidal"), NA)
    expect_warning(score(c(1, 3)), NA)
    # Later follow-up times are irrelevant to integration up to time 3.
    expect_warning(score(c(1, 2, 3), summary = "ibs"), NA)
    expect_warning(score(0.5, summary = "ibs"), NA)
    # Omitting a censoring time is also a potential missed Brier score change.
    censored <- transform(d, event = c(1, 0, 1, 0))
    expect_warning(score(c(1, 3), data = censored, summary = "ibs"), "censoring times")
    expect_warning(score(c(1, 2, 3), data = censored, summary = "ibs"), NA)

    fit <- prodlim::prodlim(Hist(time, event) ~ 1, data = d)
    expect_warning(Score(list(KM = fit), formula = Surv(time, event) ~ 1,
                         data = d, use.event.times = TRUE, summary = "ibs",
                         metrics = "brier", null.model = FALSE, se.fit = FALSE), NA)
})

test_that("IBS options also apply to competing risks and leave binary scores alone", {
    d <- data.frame(time = c(1, 2, 4, 5), event = c(1, 2, 1, 0))
    args <- list(object = list(model = matrix(0.2, 4, 2)),
                 formula = Hist(time, event) ~ 1, data = d, cause = 1,
                 times = c(1, 3), metrics = "brier", summary = "ibs",
                 null.model = FALSE, se.fit = FALSE)
    expect_warning(do.call(Score, args), "Dint")
    for (method in c("left", "trapezoidal")) {
        scores <- do.call(Score, c(args, list(ibs.method = method, verbose = FALSE)))$Brier$score
        expected <- if (method == "left") c(0, 2 * scores$Brier[1] / 3) else
            c(scores$Brier[1] / 2, (1.5 * scores$Brier[1] + scores$Brier[2]) / 3)
        expect_equal(scores$IBS, expected)
        binary <- data.frame(y = c(0, 1, 0, 1))
        expect_warning(out <- Score(list(model = rep(0.2, 4)), formula = y ~ 1,
                                    data = binary, metrics = "brier", summary = "ibs",
                                    ibs.method = method, null.model = FALSE, se.fit = FALSE), NA)
        expect_false("IBS" %in% names(out$Brier$score))
    }
})

test_that("resampling propagates the IBS method to both integration paths", {
    set.seed(471)
    d <- sampleData(60, outcome = "survival")
    fit <- survival::coxph(Surv(time, event) ~ X1 + X2, data = d, x = TRUE)
    settings <- list(list(split.method = "cv3", B = 1),
                     list(split.method = "cv3", B = 2),
                     list(split.method = "loob", B = 20),
                     list(split.method = "bootcv", B = 3))
    for (setting in settings) {
        for (method in c("left", "trapezoidal")) {
            out <- do.call(Score, c(list(object = list(Cox = fit),
                                        formula = Surv(time, event) ~ 1, data = d,
                                        times = c(1, 3), metrics = c("auc", "brier"), summary = "ibs",
                                        ibs.method = method, null.model = FALSE, se.fit = FALSE,
                                        verbose = FALSE, progress.bar = NULL, seed = 29, keep = "cv"),
                                   setting))
            scores <- if (setting$split.method == "bootcv") out$Brier$cv.score else out$Brier$score
            first <- scores[scores$times == 1, ]
            last <- scores[scores$times == 3, ]
            if (method == "left") {
                expect_equal(first$IBS, rep(0, nrow(first)))
                expect_equal(last$IBS, 2 * first$Brier / 3)
            } else {
                expect_equal(first$IBS, first$Brier / 2)
                expect_equal(last$IBS, (1.5 * first$Brier + last$Brier) / 3)
            }
        }
    }
})
