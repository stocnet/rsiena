# Accumulated + deltaFull: parity against the bootstrap reference.
#
# The bucketed REINFORCE path for accumulated specs (see
# dev-notes/accumulated-deltafull-patch-plan.md) is exercised structurally by
# test-deltamethod-buckets.R on synthetic fixtures.  Those confirm the plumbing
# but cannot confirm the numbers: their scores are uncorrelated rnorm draws, so
# the path-distribution term is noise by construction.
#
# This file is the numerical check, on a real (small) fitted model:
#   uncertaintyMode = "deltaFull"  vs  uncertaintyMode = "bootstrap"
# with bootstrap as the reference, since it re-simulates chains per draw and so
# captures the path-distribution channel by construction.
#
# Cost: ~30 s (bootstrap dominates; deltaFull is ~2 s).  Slow-gated.

test_that("accumulated deltaFull SE tracks the bootstrap reference", {
  skip_slow()
  skip_on_cran()

  ans        <- load_fixture("ans")
  mydata     <- load_fixture("mydata")
  mymodel    <- load_fixture("mymodel")
  mycontrols <- load_fixture("mycontrols")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel) ||
          is.null(mycontrols), "base fixtures unavailable")

  common <- list(
    object = ans, data = mydata,
    effectName1 = "transTrip", diff1 = 1,
    effects = mymodel, algorithm = mycontrols,
    dynamic = TRUE, n3 = 20,
    type = "tieProb", accumulated = TRUE, level = "period",
    uncertainty = TRUE, verbose = FALSE
  )

  ## deltaFull must NOT warn: a warning here means the bucketed REINFORCE path
  ## was unavailable and we silently got the conditional SE instead.
  set.seed(11L)
  expect_no_warning(
    del <- do.call(marginalEffects,
                   c(common, list(uncertaintyMode = "deltaFull")))
  )

  set.seed(11L)
  boot <- do.call(marginalEffects,
                  c(common, list(uncertaintyMode = "bootstrap", nsim = 40)))

  se_cond <- del$delta_se
  se_full <- del$delta_full_se
  se_boot <- boot$SE

  expect_true(all(is.finite(se_cond)), info = "conditional SE must be finite")
  expect_true(all(is.finite(se_full)), info = "deltaFull SE must be finite")
  expect_true(all(is.finite(se_boot)), info = "bootstrap SE must be finite")
  expect_equal(length(se_full), length(se_boot))

  ## Emit the comparison.  testthat's reporters swallow message() during a
  ## suite run, so this is only visible when the test is executed directly;
  ## the measured values are recorded in
  ## dev-notes/accumulated-deltafull-patch-plan.md Sec. 6.
  message(sprintf(
    "\n  accumulated deltaFull parity (n3=20, nsim=40):\n%s",
    paste(sprintf(
      "    period %s: bootstrap=%.5f  deltaFull=%.5f (%+.0f%%)  conditional=%.5f (%+.0f%%)",
      seq_along(se_boot), se_boot,
      se_full, 100 * (se_full / se_boot - 1),
      se_cond, 100 * (se_cond / se_boot - 1)), collapse = "\n")))

  ## (1) Order-of-magnitude agreement with the reference.  Loose on purpose:
  ##     bootstrap at nsim=40 carries roughly 10% Monte-Carlo error on the SE
  ##     itself, so a tight bound would flake rather than inform.
  expect_true(all(se_full / se_boot > 0.4 & se_full / se_boot < 2.5),
    info = "deltaFull SE must be within a factor ~2 of the bootstrap SE")

  ## (2) The substantive claim: adding the path-distribution channel moves the
  ##     SE toward the bootstrap reference.  Compared on mean absolute relative
  ##     error across periods rather than per period, so a single noisy cell
  ##     cannot flip the result.
  err_full <- mean(abs(se_full / se_boot - 1))
  err_cond <- mean(abs(se_cond / se_boot - 1))
  expect_lt(err_full, err_cond)
})
