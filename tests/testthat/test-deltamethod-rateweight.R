## test-deltamethod-rateweight.R
##
## Guard test: rateWeight specs must force finite-difference Jacobian
## in deltaMethodUncertainty (analytical path is not yet corrected).

context("deltaMethodUncertainty — rateWeight forces FD")

make_fd_only_estimator <- function() {
  function(theta, perturbations = NULL,
           useChangeContributions = FALSE, mode = "outcome") {
    if (identical(mode, "jacobian"))
      stop("analytical jacobian path should not be called for rateWeight")

    eval_at <- function(th) {
      list(spec_rw = data.frame(outcome = th[["t1"]] + 2 * th[["t2"]]))
    }

    hat <- eval_at(theta)
    if (is.null(perturbations)) return(hat)
    perts_out <- lapply(perturbations, eval_at)
    list(hat = hat, perturbations = perts_out)
  }
}

test_that("rateWeight spec uses FD and skips analytical mode", {
  est <- make_fd_only_estimator()

  thetaHat <- c(t1 = 1, t2 = 2)
  covTheta <- diag(c(0.5, 0.25))

  specs <- list(
    spec_rw = list(
      outcomeName = "outcome",
      accumulated = FALSE,
      rateWeight = TRUE
    )
  )

  ## The FD fallback is a correct-but-slower path, not a defect: it must be
  ## silent by default and only reported under verbose. See
  ## docs/design/postestimation/postestimate_api_redesign.md (Sec. 1.3).
  expect_no_warning(
    res <- deltaMethodUncertainty(
      wide = NULL,
      estimator = est,
      ssc_sum = NULL,
      thetaHat = thetaHat,
      covTheta = covTheta,
      specs = specs,
      type = "changeProb",
      fullMode = FALSE
    )
  )

  ## ...but the diagnostic is still available when asked for.
  expect_message(
    deltaMethodUncertainty(
      wide = NULL,
      estimator = est,
      ssc_sum = NULL,
      thetaHat = thetaHat,
      covTheta = covTheta,
      specs = specs,
      type = "changeProb",
      fullMode = FALSE,
      verbose = 1
    ),
    regexp = "rateWeight detected: forcing finite-difference Jacobian"
  )

  ## Linear fixture: true gradient is [1, 2].
  expect_equal(as.numeric(res$spec_rw$J_cond[1, ]), c(1, 2), tolerance = 1e-6)
})
