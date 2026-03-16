# Unit tests for mlogit_update_r (R wrapper), mlogit_update (Rcpp),
# resolvePerturbType, and perturbTypeToInt.

# ── mlogit_update_r: alter (one-alternative) mode ────────────────────────────

test_that("mlogit_update_r alter: identity when delta_u = 0", {
  p <- c(0.3, 0.5, 0.2)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, c(0, 0, 0), gid, "alter")
  expect_equal(res, p, tolerance = 1e-12)
})

test_that("mlogit_update_r alter: matches closed-form p*exp(d)/(1-p+p*exp(d))", {
  p <- c(0.25, 0.50, 0.10)
  d <- c(1.0, -0.5, 2.0)
  gid <- c(1L, 1L, 1L)
  expected <- p * exp(d) / (1 - p + p * exp(d))
  res <- mlogit_update_r(p, d, gid, "alter")
  expect_equal(res, expected, tolerance = 1e-12)
})

test_that("mlogit_update_r alter: result does NOT sum to 1 (per group)", {
  p <- c(0.3, 0.5, 0.2)
  d <- c(1, 0, 0)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, d, gid, "alter")
  # One-alternative update only changes the focal; doesn't renormalize
  expect_true(abs(sum(res) - 1) > 0.01)
})

test_that("mlogit_update_r alter: NA in delta_u propagates to NA", {
  p <- c(0.3, 0.5, 0.2)
  d <- c(1, NA, 0)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, d, gid, "alter")
  expect_true(is.na(res[2]))
  expect_false(is.na(res[1]))
  expect_false(is.na(res[3]))
})

test_that("mlogit_update_r: invalid perturbType errors", {
  p <- c(0.3, 0.5, 0.2)
  gid <- c(1L, 1L, 1L)
  expect_error(mlogit_update_r(p, c(0, 0, 0), gid, "invalid"))
})

# ── mlogit_update_r: ego (ego-wide) mode ─────────────────────────────────────

test_that("mlogit_update_r ego: identity when delta_u = 0", {
  p <- c(0.3, 0.5, 0.2)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, c(0, 0, 0), gid, "ego")
  expect_equal(res, p, tolerance = 1e-12)
})

test_that("mlogit_update_r ego: result sums to 1 within group", {
  p <- c(0.3, 0.5, 0.2)
  d <- c(1, 0, 0)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, d, gid, "ego")
  expect_equal(sum(res), 1.0, tolerance = 1e-12)
})

test_that("mlogit_update_r ego: matches p*exp(d)/sum(p*exp(d))", {
  p <- c(0.3, 0.5, 0.2)
  d <- c(1.5, -0.3, 0.7)
  gid <- c(1L, 1L, 1L)
  weighted <- p * exp(d)
  expected <- weighted / sum(weighted)
  res <- mlogit_update_r(p, d, gid, "ego")
  expect_equal(res, expected, tolerance = 1e-12)
})

test_that("mlogit_update_r ego: respects group boundaries", {
  p <- c(0.3, 0.5, 0.2, 0.6, 0.4)
  d <- c(1, 0, 0, 0.5, -0.5)
  gid <- c(1L, 1L, 1L, 2L, 2L)
  res <- mlogit_update_r(p, d, gid, "ego")

  # Group 1 sums to 1
  expect_equal(sum(res[1:3]), 1.0, tolerance = 1e-12)
  # Group 2 sums to 1
  expect_equal(sum(res[4:5]), 1.0, tolerance = 1e-12)

  # Verify group 2 independently
  w2 <- p[4:5] * exp(d[4:5])
  expect_equal(res[4:5], w2 / sum(w2), tolerance = 1e-12)
})

test_that("mlogit_update_r ego: NA propagates but doesn't break group sum", {
  p <- c(0.3, 0.5, 0.2)
  d <- c(1, NA, 0)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, d, gid, "ego")
  expect_true(is.na(res[2]))
  expect_false(is.na(res[1]))
  expect_false(is.na(res[3]))
})

test_that("mlogit_update_r ego: constant shift across group has no effect", {
  p <- c(0.3, 0.5, 0.2)
  d <- c(2, 2, 2)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, d, gid, "ego")
  # p * exp(c) / sum(p * exp(c)) = p * exp(c) / (exp(c) * sum(p)) = p
  expect_equal(res, p, tolerance = 1e-12)
})

test_that("mlogit_update_r: returns plain vector, not matrix", {
  p <- c(0.3, 0.5, 0.2)
  gid <- c(1L, 1L, 1L)
  res <- mlogit_update_r(p, c(1, 0, 0), gid, "alter")
  expect_null(dim(res))
  expect_true(is.numeric(res))
})

# ── resolvePerturbType ───────────────────────────────────────────────────────

test_that("resolvePerturbType: override takes precedence", {
  iTypes <- c(mynet_egoX_eval = "ego")
  expect_equal(resolvePerturbType("mynet_egoX_eval", iTypes, "alter"), "alter")
  expect_equal(resolvePerturbType("mynet_egoX_eval", iTypes, "ego"), "ego")
})

test_that("resolvePerturbType: NULL override uses auto-detection", {
  iTypes <- c(mynet_egoX_eval = "ego", mynet_density_eval = "dyadic",
              mynet_transTrip_eval = "")
  expect_equal(resolvePerturbType("mynet_egoX_eval", iTypes, NULL), "ego")
  expect_equal(resolvePerturbType("mynet_density_eval", iTypes, NULL), "alter")
  expect_equal(resolvePerturbType("mynet_transTrip_eval", iTypes, NULL), "alter")
})

test_that("resolvePerturbType: NULL iTypes defaults to alter", {
  expect_equal(resolvePerturbType("anything", NULL, NULL), "alter")
})

test_that("resolvePerturbType: unknown effect defaults to alter", {
  iTypes <- c(mynet_density_eval = "dyadic")
  expect_equal(resolvePerturbType("missing_effect", iTypes, NULL), "alter")
})

# ── perturbTypeToInt ─────────────────────────────────────────────────────────

test_that("perturbTypeToInt: correct mapping", {
  expect_equal(perturbTypeToInt("alter"), 0L)
  expect_equal(perturbTypeToInt("ego"), 1L)
})

test_that("perturbTypeToInt: invalid input errors", {
  expect_error(perturbTypeToInt("invalid"))
})

# ── computeMassContrasts ─────────────────────────────────────────────────────

test_that("computeMassContrasts: basic creation/dissolution sums", {
  # Two egos, one period, one group.
  # Ego 1: 2 creation rows, 1 dissolution row
  # Ego 2: 1 creation row, 2 dissolution rows
  fd      <- c(0.1, 0.2, -0.05,     0.3, -0.1, -0.15)
  density <- c(  1,   1,     -1,       1,   -1,    -1)
  ego     <- c(  1,   1,      1,       2,    2,     2)
  period  <- c(  1,   1,      1,       1,    1,     1)
  group   <- c(  1,   1,      1,       1,    1,     1)

  mc <- computeMassContrasts(fd, density, ego, period, group, type = "changeProb")
  # Ego 1: massCreation = 0.1+0.2 = 0.3, massDissolution = -0.05
  expect_equal(mc$massCreation[1], 0.3, tolerance = 1e-12)
  expect_equal(mc$massDissolution[1], -0.05, tolerance = 1e-12)
  # Same values broadcast to all ego-1 rows
  expect_equal(mc$massCreation[3], 0.3, tolerance = 1e-12)
  expect_equal(mc$massDissolution[3], -0.05, tolerance = 1e-12)
  # Ego 2: massCreation = 0.3, massDissolution = -0.1 + -0.15 = -0.25
  expect_equal(mc$massCreation[4], 0.3, tolerance = 1e-12)
  expect_equal(mc$massDissolution[4], -0.25, tolerance = 1e-12)
})

test_that("computeMassContrasts: tieProb mode flips dissolution", {
  # In tieProb mode, dissolution firstDiff = -(changeProbDiff)
  # so we need to flip sign for dissolution rows
  fd      <- c(0.1, 0.05)  # creation, dissolution (tieProb scale)
  density <- c(  1,   -1)
  ego     <- c(  1,    1)
  period  <- c(  1,    1)
  group   <- c(  1,    1)

  mc <- computeMassContrasts(fd, density, ego, period, group, type = "tieProb")
  # massCreation = 0.1 (no flip for creation)
  expect_equal(mc$massCreation[1], 0.1, tolerance = 1e-12)
  # massDissolution: changeProbDiff = -(0.05) = -0.05 (flip back from tieProb)
  expect_equal(mc$massDissolution[1], -0.05, tolerance = 1e-12)
})

test_that("computeMassContrasts: multiple periods kept separate", {
  fd      <- c(0.1, 0.2, 0.3, 0.4)
  density <- c(  1,   1,   1,   1)
  ego     <- c(  1,   1,   1,   1)
  period  <- c(  1,   1,   2,   2)
  group   <- c(  1,   1,   1,   1)

  mc <- computeMassContrasts(fd, density, ego, period, group, type = "changeProb")
  expect_equal(mc$massCreation[1], 0.3, tolerance = 1e-12)   # period 1
  expect_equal(mc$massCreation[3], 0.7, tolerance = 1e-12)   # period 2
  expect_equal(mc$massDissolution[1], 0, tolerance = 1e-12)   # no dissolution
})

test_that("computeMassContrasts: ego with only creation has massDissolution = 0", {
  fd      <- c(0.1, 0.2)
  density <- c(  1,   1)
  ego     <- c(  1,   1)
  period  <- c(  1,   1)
  group   <- c(  1,   1)

  mc <- computeMassContrasts(fd, density, ego, period, group, type = "changeProb")
  expect_equal(mc$massCreation[1], 0.3, tolerance = 1e-12)
  expect_equal(mc$massDissolution[1], 0, tolerance = 1e-12)
})
