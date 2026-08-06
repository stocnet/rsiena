## test-bootstrap-draws.R
##
## Regression tests for the parametric-bootstrap draws in drawSimBatch.
##
## The bug these pin: theta used to be drawn INSIDE the draw loop, one
## mvrnorm call per draw.  Every dynamic estimator reaches siena07, which
## does set.seed(algorithm$randomSeed) (and, when no seed is set, removes
## .Random.seed and reseeds from entropy).  Either way the global RNG stream
## is reset underneath the loop, so each theta after the first was a
## deterministic function of its predecessor rather than an independent
## draw.  The bootstrap SE then came out biased LOW by an amount that GREW
## with nsim -- a variance estimate that decays instead of converging.
##
## No simulation is needed to test this: an estimator that resets the seed
## the way siena07 does reproduces the pathology exactly.

context("bootstrap draws — independence under a seed-resetting estimator")

## An estimator that (a) records the theta it was handed and (b) resets the
## global RNG, exactly as a dynamic RSiena estimator does via siena07.
make_recording_estimator <- function(store, reset_seed = TRUE) {
  est <- function(theta, ...) {
    store$thetas[[length(store$thetas) + 1L]] <- theta
    if (reset_seed) set.seed(11)
    # consume a theta-dependent number of variates, as a simulation would
    invisible(runif(5L + as.integer(abs(theta[[1L]]) * 37) %% 23L))
    list(eff = data.frame(period = 1L, value = unname(theta[[1L]])))
  }
  attr(est, "eff_names") <- "eff"
  est
}

test_that("theta draws stay iid when the estimator resets the RNG seed", {
  thetaHat <- c(a = 0, b = 0)
  covTheta <- diag(c(1, 4))

  store <- new.env(); store$thetas <- list()
  est   <- make_recording_estimator(store, reset_seed = TRUE)

  set.seed(1)
  invisible(drawSimBatch(estimator = est, thetaHat = thetaHat,
                         covTheta = covTheta, nsim = 400L,
                         batchSize = 100L, verbose = FALSE))

  th <- do.call(rbind, store$thetas)
  expect_equal(nrow(th), 400L)

  # All draws distinct: the degenerate sequence collapsed onto a short cycle.
  expect_equal(nrow(unique(th)), 400L)

  # Spread must recover the specified covariance, not a fraction of it.
  # With 400 draws the sd of an sd is ~1/sqrt(2*399) ~ 3.5%, so 15% is loose
  # enough to never flake and tight enough to catch the old ~55% deficit.
  expect_equal(unname(sd(th[, 1])), 1, tolerance = 0.15)
  expect_equal(unname(sd(th[, 2])), 2, tolerance = 0.15)
  expect_equal(unname(mean(th[, 1])), 0, tolerance = 0.20)
})

test_that("draw spread does not decay as nsim grows", {
  thetaHat <- c(a = 0)
  covTheta <- matrix(1, 1, 1, dimnames = list("a", "a"))

  sds <- vapply(c(50L, 100L, 200L, 400L), function(n) {
    store <- new.env(); store$thetas <- list()
    est   <- make_recording_estimator(store, reset_seed = TRUE)
    set.seed(7)
    invisible(drawSimBatch(estimator = est, thetaHat = thetaHat,
                           covTheta = covTheta, nsim = n,
                           batchSize = 100L, verbose = FALSE))
    sd(do.call(rbind, store$thetas)[, 1])
  }, numeric(1L))

  # The old code produced a monotone decay (0.48, 0.45, 0.43, 0.42 against a
  # true sd of 1).  Every value must now sit near 1 with no downward trend.
  expect_true(all(abs(sds - 1) < 0.25),
              info = paste("sds:", paste(round(sds, 4), collapse = ", ")))
  expect_gt(sds[[4L]], sds[[1L]] * 0.7)
})

test_that("draws are reproducible from the caller's seed and ignore batchSize", {
  thetaHat <- c(a = 0, b = 0)
  covTheta <- diag(c(1, 4))

  run <- function(bs) {
    store <- new.env(); store$thetas <- list()
    est   <- make_recording_estimator(store, reset_seed = TRUE)
    set.seed(123)
    invisible(drawSimBatch(estimator = est, thetaHat = thetaHat,
                           covTheta = covTheta, nsim = 200L,
                           batchSize = bs, verbose = FALSE))
    do.call(rbind, store$thetas)
  }

  # Same seed, same draws; and batching must not change them.
  expect_equal(run(100L), run(100L))
  expect_equal(run(100L), run(50L))
})

test_that("a stale or foreign batch file errors instead of corrupting draws", {
  thetaHat <- c(a = 0)
  covTheta <- matrix(1, 1, 1, dimnames = list("a", "a"))

  bdir <- file.path(tempdir(), "stalebatch")
  dir.create(bdir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(bdir, recursive = TRUE), add = TRUE)

  # A leftover batch from a run with a different nsim: right name, wrong size.
  saveRDS(list(vector("list", 7L)),
          file.path(bdir, "simBatchB_b001.rds"))

  store <- new.env(); store$thetas <- list()
  est   <- make_recording_estimator(store, reset_seed = TRUE)

  set.seed(3)
  expect_error(
    drawSimBatch(estimator = est, thetaHat = thetaHat, covTheta = covTheta,
                 nsim = 200L, batchSize = 100L, batchDir = bdir,
                 verbose = FALSE),
    "different run|Stale or foreign"
  )
})

test_that("nsim = 1 still yields a well-formed single draw", {
  thetaHat <- c(a = 1, b = 2)
  covTheta <- diag(c(1, 1))

  store <- new.env(); store$thetas <- list()
  est   <- make_recording_estimator(store, reset_seed = TRUE)

  set.seed(5)
  out <- drawSimBatch(estimator = est, thetaHat = thetaHat,
                      covTheta = covTheta, nsim = 1L,
                      batchSize = 100L, verbose = FALSE)

  expect_equal(length(store$thetas), 1L)
  expect_equal(names(store$thetas[[1L]]), c("a", "b"))
  expect_equal(nrow(out$eff), 1L)
})

context("postestimation seed resolution")

test_that("control seed wins, then algorithm seed, then nothing", {
  alg <- list(randomSeed = 11)

  # 1. explicit control seed overrides the algorithm object
  r <- .postestSeeds(alg, 77)
  expect_equal(r$simSeed, 77L)
  expect_equal(r$algorithm$randomSeed, 77L)
  # ... on a COPY: the caller's object is untouched
  expect_equal(alg$randomSeed, 11)

  # 2. no control seed -> the algorithm's own seed
  r <- .postestSeeds(alg, NULL)
  expect_equal(r$simSeed, 11L)

  # 3. neither -> unreproducible, as in siena07
  r <- .postestSeeds(list(randomSeed = NULL), NULL)
  expect_null(r$simSeed)
  expect_null(r$drawSeed)
  expect_null(r$chainSeed)

  # 4. no algorithm at all (static paths) still gets a draw seed
  r <- .postestSeeds(NULL, 5)
  expect_equal(r$simSeed, 5L)
  expect_false(is.null(r$drawSeed))
})

test_that("derived seeds are drawn, not offset, and never collide", {
  r <- .postestSeeds(list(randomSeed = 42), NULL)
  expect_equal(r$simSeed, 42L)

  # NOT simSeed + 1: arithmetic offsets made the draw seed equal chain batch
  # 2's seed (both simSeed + 1), coupling the theta draws to those chains.
  expect_false(identical(r$drawSeed, r$simSeed + 1L))
  expect_false(identical(r$chainSeed, r$simSeed + 1L))
  expect_false(identical(r$drawSeed, r$chainSeed))
  expect_false(identical(r$drawSeed, r$simSeed))

  # Per-batch chain seeds are likewise drawn, and distinct from each other
  # and from the draw seed.
  bs <- .withSeed(r$chainSeed, sample.int(.Machine$integer.max, 8L))
  expect_equal(length(unique(bs)), 8L)
  expect_false(r$drawSeed %in% bs)
  expect_false(r$simSeed %in% bs)

  # Still deterministic: same master seed, same derived seeds.
  r2 <- .postestSeeds(list(randomSeed = 42), NULL)
  expect_equal(r$drawSeed,  r2$drawSeed)
  expect_equal(r$chainSeed, r2$chainSeed)
})

test_that(".withSeed restores the caller's stream", {
  set.seed(321)
  before <- get(".Random.seed", envir = globalenv())
  x <- .withSeed(999, sample.int(1000, 3))
  # base identical(): .Random.seed holds large integers and testthat's
  # numeric comparison overflows on them.
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))
  # and it is deterministic in its own right
  expect_equal(x, .withSeed(999, sample.int(1000, 3)))
})

test_that("a non-numeric algorithm seed warns rather than running silently", {
  expect_warning(r <- .postestSeeds(list(randomSeed = "1293"), NULL),
                 "not numeric")
  expect_null(r$simSeed)
})

test_that("drawSeed makes theta draws reproducible and restores the stream", {
  thetaHat <- c(a = 0, b = 0)
  covTheta <- diag(c(1, 4))

  grab <- function(seed) {
    store <- new.env(); store$thetas <- list()
    est   <- make_recording_estimator(store, reset_seed = TRUE)
    drawSimBatch(estimator = est, thetaHat = thetaHat, covTheta = covTheta,
                 nsim = 30L, batchSize = 100L, drawSeed = seed,
                 verbose = FALSE)
    do.call(rbind, store$thetas)
  }

  # Same drawSeed -> same draws, regardless of the ambient stream.
  set.seed(1);   a <- grab(2024L)
  set.seed(999); b <- grab(2024L)
  expect_equal(a, b)
  # Different drawSeed -> different draws.
  set.seed(1);   c <- grab(7L)
  expect_false(isTRUE(all.equal(a, c)))
})
