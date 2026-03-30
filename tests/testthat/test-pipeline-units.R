# testthat::skip_on_cran()

# Unit tests for pipeline sub-functions in predict.R, sienaMargins.r, postestimate.R.
# These tests use synthetic data: no siena07() call required.
# Goal: catch errors early in the computational pipeline before they surface
# at the integration level (marginalEffects / predict.sienaFit).

# ── Synthetic helper ──────────────────────────────────────────────────────────

# Minimal static changeContributions struct (2 effects, 3 egos × 3 choices = 9 rows)
make_cc2 <- function(seed = 42) {
  set.seed(seed)
  # density in {1=add, -1=drop, 0=no-change}; one of each per ego
  density <- rep(c(1L, -1L, 0L), times = 3)
  recip   <- c(0L, 1L, 0L,  1L, 0L, 1L,  0L, 0L, 1L)
  mat <- matrix(c(density, recip), nrow = 9,
                dimnames = list(NULL, c("density", "recip")))
  list(
    contribMat = mat,
    group_id    = rep(1L:3L, each = 3L),
    group       = 1L,
    period      = rep(c(1L, 2L, 1L), each = 3L),
    ego         = rep(c(1L, 2L, 3L), each = 3L),
    choice      = density,
    effectNames = c("density", "recip")
  )
}

# Same structure but with 3 effects
make_cc3 <- function(seed = 42) {
  set.seed(seed)
  density   <- rep(c(1L, -1L, 0L), times = 3)
  recip     <- c(0L, 1L, 0L,  1L, 0L, 1L,  0L, 0L, 1L)
  transTrip <- c(2L, 1L, 0L,  1L, 2L, 0L,  0L, 1L, 3L)
  mat <- matrix(c(density, recip, transTrip), nrow = 9,
                dimnames = list(NULL, c("density", "recip", "transTrip")))
  list(
    contribMat = mat,
    group_id    = rep(1L:3L, each = 3L),
    group       = 1L,
    period      = rep(c(1L, 2L, 1L), each = 3L),
    ego         = rep(c(1L, 2L, 3L), each = 3L),
    choice      = density,
    effectNames = c("density", "recip", "transTrip")
  )
}

# ── calculateUtility ─────────────────────────────────────────────────────────

test_that("calculateUtility: matrix-vector product (numeric)", {
  contribMat   <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  theta <- c(0.5, -1)
  expect_equal(calculateUtility(contribMat, theta), as.numeric(contribMat %*% theta))
})

test_that("calculateUtility: single-column matrix", {
  contribMat   <- matrix(c(1, -1, 0), nrow = 3)
  theta <- c(-2)
  expect_equal(calculateUtility(contribMat, theta), c(-2, 2, 0))
})

test_that("calculateUtility: requires matrix input", {
  expect_error(calculateUtility(c(1, 2, 3), c(1)))
})

# ── calculateTieProb ───────────────────────────────────────────────────────────


test_that("calculateTieProb: density==-1 flips, density==0 → NA, density==1 unchanged", {
  prob        <- c(0.3, 0.7, 0.5)
  density_col <- c(1L, -1L, 0L)
  result <- calculateTieProb(prob, density_col)
  expect_equal(result[1], 0.3)           # density==1:  no flip
  expect_equal(result[2], 1 - 0.7)      # density==-1: flipped
  expect_true(is.na(result[3]))          # density==0:  NA
})

test_that("calculateTieProb: all density==1 → no change", {
  prob <- runif(5)
  result <- calculateTieProb(prob, rep(1L, 5))
  expect_equal(result, prob)
})

# ── calculateUtilityDiff ─────────────────────────────────────────────────────

test_that("calculateUtilityDiff: simple additive case", {
  effectNames <- c("density", "recip", "transTrip")
  theta       <- c(density = -2, recip = 1.5, transTrip = 0.5)
  densityValue <- c(1L, -1L, 0L)
  result <- calculateUtilityDiff(
    effectName = "recip", diff = 1,
    theta = theta, densityValue = densityValue,
    effectNames = effectNames
  )
  # util_diff = densityValue * diff * theta["recip"]
  expect_equal(result, densityValue * 1 * 1.5)
})

test_that("calculateUtilityDiff: interaction adds moderator term", {
  effectNames  <- c("density", "recip", "transTrip", "unspInt")
  theta        <- c(density = -2, recip = 1.0, transTrip = 0.5, unspInt = 0.3)
  densityValue <- c(1L, 1L, 1L)
  modContrib   <- c(0L, 2L, 1L)
  result <- calculateUtilityDiff(
    effectName = "recip", diff = 1,
    theta = theta, densityValue = densityValue,
    interaction = TRUE,
    intEffectNames = "unspInt", modEffectNames = "transTrip",
    modContribution = modContrib,
    effectNames = effectNames
  )
  # CS-space formula: d * diff * (theta[recip] + modContrib * theta[unspInt])
  expected <- densityValue * (1 * theta["recip"] + modContrib * theta["unspInt"])
  expect_equal(result, unname(expected))
})

# ── calculateFirstDiff ────────────────────────────────────────────────────────

test_that("calculateFirstDiff: returns named list, density==0 → NA", {
  n <- 9
  densityValue  <- rep(c(1L, -1L, 0L), 3)
  changeProb    <- rep(0.5, n)
  changeUtil    <- rep(0, n)
  effectContrib <- rep(1L, n)
  theta         <- c(density = -2, recip = 1.5)
  effectNames   <- c("density", "recip")

  result <- calculateFirstDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil, effectName = "recip",
    effectContribution = effectContrib, diff = 1,
    effectNames = effectNames, theta = theta, type = "changeProb",
    details = FALSE
  )
  expect_type(result, "list")
  expect_named(result, "firstDiff")
  expect_equal(length(result$firstDiff), n)
  expect_true(all(is.na(result$firstDiff[densityValue == 0L])))
  expect_true(all(!is.na(result$firstDiff[densityValue != 0L])))
})

test_that("calculateFirstDiff: positive density → positive diff, symmetry", {
  densityValue  <- c(1L, -1L, 0L)
  changeProb    <- rep(0.5, 3)
  changeUtil    <- rep(0, 3)
  effectContrib <- rep(1L, 3)
  theta         <- c(density = -2, recip = 1.5)

  result <- calculateFirstDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil, effectName = "recip",
    effectContribution = effectContrib, diff = 1,
    effectNames = c("density", "recip"), theta = theta,
    type = "changeProb", details = FALSE
  )
  fd <- result$firstDiff
  expect_gt(fd[1], 0)    # positive density: tie probability increases
  expect_lt(fd[2], 0)    # negative density: tie probability decreases
  expect_equal(fd[1], -fd[2], tolerance = 1e-10)  # symmetric around 0.5 base prob
})

test_that("calculateFirstDiff: details=TRUE returns data.frame with named columns", {
  densityValue  <- c(1L, -1L, 0L)
  changeProb    <- rep(0.5, 3)
  changeUtil    <- rep(0, 3)
  effectContrib <- rep(1L, 3)
  theta         <- c(density = -2, recip = 1.5)

  result <- calculateFirstDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil, effectName = "recip",
    effectContribution = effectContrib, diff = 1,
    effectNames = c("density", "recip"), theta = theta,
    type = "changeProb", details = TRUE
  )
  expect_true(is.data.frame(result))
  expect_true(all(c("firstDiff", "utilDiff", "newChangeProb", "oldChangeProb") %in% names(result)))
  expect_equal(nrow(result), 3)
})

test_that("calculateFirstDiff: tieProb type returns correct tieProb columns with details", {
  densityValue  <- c(1L, -1L, 1L)
  changeProb    <- rep(0.5, 3)
  changeUtil    <- rep(0, 3)
  effectContrib <- rep(1L, 3)
  theta         <- c(density = -2, recip = 1.5)
  tieProb_in    <- calculateTieProb(changeProb, densityValue)

  result <- calculateFirstDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil, effectName = "recip",
    effectContribution = effectContrib, diff = 1,
    effectNames = c("density", "recip"), theta = theta,
    type = "tieProb", tieProb = tieProb_in, details = TRUE
  )
  expect_true(all(c("newTieProb", "oldTieProb") %in% names(result)))
})

test_that("calculateFirstDiff: riskRatio path returns firstRiskRatio", {
  densityValue  <- c(1L, 1L, 1L)
  changeProb    <- rep(0.4, 3)
  changeUtil    <- rep(0.5, 3)
  effectContrib <- rep(1L, 3)
  theta         <- c(density = -2, recip = 1.0)

  result <- calculateFirstDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil, effectName = "recip",
    effectContribution = effectContrib, diff = 1,
    effectNames = c("density", "recip"), theta = theta,
    type = "changeProb", details = FALSE, mainEffect = "riskRatio"
  )
  expect_named(result, "firstRiskRatio")
  expect_true(all(result$firstRiskRatio > 0))  # ratios always positive
})

# ── calculateSecondDiff ───────────────────────────────────────────────────────

test_that("calculateSecondDiff: returns named list, density==0 → NA", {
  n <- 9
  densityValue  <- rep(c(1L, -1L, 0L), 3)
  changeProb    <- rep(0.5, n)
  changeUtil    <- rep(0, n)
  theta         <- c(density = -2, recip = 1.5, transTrip = 0.5)
  effectNames   <- c("density", "recip", "transTrip")

  result <- calculateSecondDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil,
    effectName1 = "recip",   effectContribution1 = rep(1L, n), diff1 = 1,
    effectName2 = "transTrip", effectContribution2 = rep(2L, n), diff2 = 1,
    effectNames = effectNames, theta = theta,
    type = "changeProb", details = FALSE
  )
  expect_type(result, "list")
  expect_named(result, "secondDiff")
  expect_equal(length(result$secondDiff), n)
  expect_true(all(is.na(result$secondDiff[densityValue == 0L])))
})

test_that("calculateSecondDiff: details=TRUE returns broad data.frame", {
  densityValue  <- c(1L, -1L, 1L)
  changeProb    <- rep(0.5, 3)
  changeUtil    <- rep(0, 3)
  theta         <- c(density = -2, recip = 1.5, transTrip = 0.5)
  effectNames   <- c("density", "recip", "transTrip")

  result <- calculateSecondDiff(
    densityValue = densityValue, changeProb = changeProb,
    changeUtil = changeUtil,
    effectName1 = "recip",     effectContribution1 = rep(1L, 3), diff1 = 1,
    effectName2 = "transTrip", effectContribution2 = rep(1L, 3), diff2 = 1,
    effectNames = effectNames, theta = theta,
    type = "changeProb", details = TRUE
  )
  expect_true(is.data.frame(result))
  expect_true(all(c("changeProb_base", "changeProb_main", "changeProb_mod",
                    "changeProb_both", "firstDiff1", "firstDiff2", "secondDiff") %in% names(result)))
})

# ── attachContributions ─────────────────────────────────────────────────────

test_that("attachContributions: density==-1 flips non-density columns only", {
  density <- c(1L, -1L, 0L)
  recip   <- c(2L, 3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density", "recip")))
  df <- data.frame(x = 1:3)

  result <- attachContributions(df, c("density", "recip"), contrib, flip = TRUE)

  expect_equal(result$density, c(1L, -1L, 0L))      # density keeps ±1 sign
  expect_equal(result$recip,   c(2L, -3L, 4L))      # density==-1 row → negated
})

test_that("attachContributions: flip=FALSE leaves all columns unchanged", {
  density <- c(1L, -1L, -1L)
  recip   <- c(2L, 3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density", "recip")))
  df <- data.frame(x = 1:3)

  result <- attachContributions(df, c("density", "recip"), contrib, flip = FALSE)

  expect_equal(result$recip, c(2L, 3L, 4L))        # no flipping
})

# ── groupColsList ───────────────────────────────────────────────────────────

test_that("groupColsList: static (no chain) returns correct names and values", {
  pb <- list(
    group = 1L, period = c(1L, 2L, 3L),
    ego = c(5L, 6L, 7L), choice = c(1L, 0L, 1L)
  )
  result <- groupColsList(pb)
  expect_named(result, c("group", "period", "ego", "choice"))
  expect_equal(result$ego, c(5L, 6L, 7L))
})

test_that("groupColsList: static with keep subsets correctly", {
  pb <- list(
    group = 1L, period = c(1L, 2L, 3L),
    ego = c(5L, 6L, 7L), choice = c(1L, 0L, 1L)
  )
  result <- groupColsList(pb, keep = c(TRUE, FALSE, TRUE))
  expect_equal(result$ego,    c(5L, 7L))
  expect_equal(result$period, c(1L, 3L))
})

test_that("groupColsList: dynamic (with chain) returns correct names", {
  pb <- list(
    chain = c(1L, 1L, 2L), group = c(1L, 1L, 1L),
    period = c(1L, 2L, 1L), ministep = c(10L, 20L, 30L),
    choice = c(0L, 1L, 1L)
  )
  result <- groupColsList(pb)
  expect_named(result, c("chain", "group", "period", "ministep", "choice"))
  expect_equal(result$ministep, c(10L, 20L, 30L))
})

# ── alignThetaNoRate ─────────────────────────────────────────────────────────

test_that("alignThetaNoRate: named theta — selects by name", {
  theta <- c(density = -2, recip = 1.5, transTrip = 0.5)
  result <- alignThetaNoRate(theta, c("recip", "transTrip"))
  expect_named(result, c("recip", "transTrip"))
  expect_equal(unname(result["recip"]), 1.5)
})

test_that("nameThetaFromEffects + alignThetaNoRate: cond=FALSE full-theta alignment", {
  # Simulate a cond=FALSE model where theta includes rate parameters.
  # nameThetaFromEffects warns and names unnamed theta; alignThetaNoRate then
  # selects the eval params.
  theta_raw <- c(-2.38, 5.1, 4.9, 3.06, -0.07)  # density, rate1, rate2, recip, unspInt
  effs <- data.frame(
    shortName = c("density", "rateX", "rateX", "recip", "unspInt"),
    type      = c("eval",    "rate",  "rate",  "eval",  "eval"),
    stringsAsFactors = FALSE
  )
  theta_named <- expect_warning(
    nameThetaFromEffects(theta_raw, effs),
    "theta has no names"
  )
  effectNames <- c("density_eval", "recip_eval", "unspInt_eval")
  result <- alignThetaNoRate(theta_named, effectNames)
  expect_named(result, effectNames)
  expect_equal(unname(result[["density_eval"]]), -2.38)
  expect_equal(unname(result[["recip_eval"]]),    3.06)
  expect_equal(unname(result[["unspInt_eval"]]), -0.07)
})

test_that("nameThetaFromEffects: returns named theta unchanged; warns for unnamed", {
  # Already-named theta is returned as-is.
  theta_named <- c(mynet_density_eval = -2, mynet_recip_eval = 1.5, mynet_transTrip_eval = 0.5)
  effs <- data.frame(
    shortName = c("density", "recip", "transTrip"),
    type      = c("eval",    "eval",  "eval"),
    stringsAsFactors = FALSE
  )
  result <- nameThetaFromEffects(theta_named, effs)
  expect_identical(result, theta_named)  # returned as-is

  # Unnamed theta triggers warning and gets names from getNamesFromEffects fallback.
  theta_unnamed <- c(-2, 1.5, 0.5)
  result2 <- expect_warning(
    nameThetaFromEffects(theta_unnamed, effs),
    "theta has no names"
  )
  expect_named(result2, c("density_eval", "recip_eval", "transTrip_eval"))
  expect_equal(unname(result2), c(-2, 1.5, 0.5))
})

# ── planBatch ───────────────────────────────────────────────────────────────

test_that("planBatch: result is always >= 1 and <= nsim", {
  dv <- array(0, dim = c(50, 50, 3))
  mock_data <- list(depvars = list(mynet = dv))
  result <- planBatch(mock_data, "mynet", nsim = 100, dynamic = FALSE)
  expect_gte(result, 1L)
  expect_lte(result, 100L)
})

test_that("planBatch: nsim=1 always returns 1", {
  dv <- array(0, dim = c(50, 50, 3))
  mock_data <- list(depvars = list(mynet = dv))
  expect_equal(planBatch(mock_data, "mynet", nsim = 1L), 1L)
})

test_that("planBatch: returns positive integer within nsim", {
  dv <- array(0, dim = c(10, 10, 3))
  mock_data <- list(depvars = list(net = dv))
  result <- planBatch(mock_data, "net", nsim = 50, dynamic = FALSE)
  expect_true(is.integer(result))
  expect_gte(result, 1L)
  expect_lte(result, 50L)
})

test_that("planBatch: memoryScale reduces batch size", {
  dv <- array(0, dim = c(50, 50, 3))
  mock_data <- list(depvars = list(mynet = dv))
  b1 <- planBatch(mock_data, "mynet", nsim = 200, dynamic = FALSE)
  b2 <- planBatch(mock_data, "mynet", nsim = 200, dynamic = FALSE, memoryScale = 4L)
  expect_lte(b2, b1)
})

# ── predictFirstDiff ─────────────────────────────────────────────────────────

test_that("predictFirstDiff: density==0 rows excluded; contrib cols present", {
  cc         <- make_cc2()
  theta_use  <- c(density = -2, recip = 1.5)
  n_valid    <- sum(cc$contribMat[, "density"] != 0L)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_equal(nrow(result), n_valid)
  expect_true("firstDiff" %in% names(result))
  expect_true("density"   %in% names(result))   # contrib columns attached
  expect_true("recip"     %in% names(result))
  expect_true(all(!is.na(result$firstDiff)))     # all remaining rows are non-NA
})

test_that("predictFirstDiff: details=TRUE attaches utilDiff, changeProb", {
  cc        <- make_cc2()
  theta_use <- c(density = -2, recip = 1.5)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = TRUE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true(all(c("changeUtil", "changeProb") %in% names(result)))
  expect_false("tieProb" %in% names(result))    # tieProb only for tieProb type
})

test_that("predictFirstDiff: tieProb type attaches tieProb col when details=TRUE", {
  cc        <- make_cc2()
  theta_use <- c(density = -2, recip = 1.5)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "tieProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = TRUE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true("tieProb" %in% names(result))
})

# data.table removed — output is always data.frame now
test_that("predictFirstDiff: output is always data.frame", {
  # skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  # library(data.table)
  cc        <- make_cc2()
  theta_use <- c(density = -2, recip = 1.5)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true(is.data.frame(result))
  expect_false("data.table" %in% class(result))
})

test_that("predictFirstDiff: base R output is plain data.frame", {
  # with_mocked_bindings no longer needed — data.table code paths removed
      cc        <- make_cc2()
      theta_use <- c(density = -2, recip = 1.5)

      result <- predictFirstDiff(
        changeContributions = cc, theta_use = theta_use, type = "changeProb",
        effectName = "recip", diff = 1, contrast = NULL,
        interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
        details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
      )
      expect_true(is.data.frame(result))
      expect_false("data.table" %in% class(result))
})

# ── predictSecondDiff ────────────────────────────────────────────────────────

test_that("predictSecondDiff: density==0 rows excluded; secondDiff col present", {
  cc         <- make_cc3()
  theta_use  <- c(density = -2, recip = 1.5, transTrip = 0.5)
  n_valid    <- sum(cc$contribMat[, "density"] != 0L)

  result <- predictSecondDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName1 = "recip",     diff1 = 1, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = "transTrip", diff2 = 1, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_equal(nrow(result), n_valid)
  expect_true("secondDiff" %in% names(result))
  expect_true(all(!is.na(result$secondDiff)))
})

test_that("predictSecondDiff: details=TRUE attaches utility and prob columns", {
  cc        <- make_cc3()
  theta_use <- c(density = -2, recip = 1.5, transTrip = 0.5)

  result <- predictSecondDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName1 = "recip",     diff1 = 1, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = "transTrip", diff2 = 1, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    details = TRUE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true(all(c("changeUtil", "changeProb") %in% names(result)))
})

# ── getGroupVars ─────────────────────────────────────────────────────────────

test_that("getGroupVars: composite condition appended to level vars", {
  expect_equal(getGroupVars("period", "transTrip_eval"),
               c("period", "transTrip_eval"))
})

test_that("getGroupVars: multiple composite conditions appended", {
  expect_equal(getGroupVars("ego", c("density_eval", "recip_creation")),
               c("period", "ego", "density_eval", "recip_creation"))
})

test_that("getGroupVars: NULL condition returns only level vars", {
  expect_equal(getGroupVars("period"),      "period")
  expect_equal(getGroupVars("none"),        character(0))
})

# ── makeGroupKey (removed — streaming path deleted) ──────────────────────────
#
# test_that("makeGroupKey: composite column name indexes correctly", { ... })
# test_that("makeGroupKey: empty group_vars returns __all__", { ... })

# ── agg: bare-name resolution ────────────────────────────────────────────────

test_that("agg: bare condition resolves to _eval variant column", {
  df <- data.frame(period = c(1L, 1L, 2L),
                   `transTrip_eval` = c(0, 1, 0),
                   changeProb = c(0.3, 0.7, 0.4),
                   check.names = FALSE)
  result <- with_mocked_bindings(
    agg("changeProb", df, level = "period", condition = "transTrip"),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
  expect_true("transTrip_eval" %in% names(result))
  expect_false("transTrip"     %in% names(result))
})

test_that("agg: composite condition passes through unchanged", {
  df <- data.frame(period = c(1L, 1L, 2L),
                   `transTrip_eval` = c(0, 1, 0),
                   changeProb = c(0.3, 0.7, 0.4),
                   check.names = FALSE)
  result <-
    agg("changeProb", df, level = "period", condition = "transTrip_eval")
  expect_true("transTrip_eval" %in% names(result))
  expect_equal(nrow(result), 3L)   # period + transTrip_eval + changeProb value
})

test_that("agg: creation/endow composite condition aggregates correctly", {
  df <- data.frame(period = c(1L, 1L, 2L, 2L),
                   `recip_creation` = c(0L, 1L, 0L, 1L),
                   changeProb = c(0.2, 0.8, 0.3, 0.7),
                   check.names = FALSE)
  result <- 
    agg("changeProb", df, level = "period", condition = "recip_creation")
  expect_true("recip_creation" %in% names(result))
  expect_equal(nrow(result), 4L)   # 2 periods × 2 condition values
})

test_that("agg: co-evolution — separate depvar data aggregates independently", {
  # Each depvar's prediction data has its own "density_eval" column;
  # agg() is called separately on each → results must differ.
  df_a <- data.frame(period = c(1L, 1L), `density_eval` = c(1L, -1L),
                     changeProb = c(0.8, 0.2), check.names = FALSE)
  df_b <- data.frame(period = c(1L, 1L), `density_eval` = c(1L, -1L),
                     changeProb = c(0.6, 0.4), check.names = FALSE)
  agg_a <-
    agg("changeProb", df_a, level = "period", condition = "density_eval")

  agg_b <-
    agg("changeProb", df_b, level = "period", condition = "density_eval")
  expect_false(identical(agg_a$changeProb, agg_b$changeProb))
})

# ── updateStream / finalizeStream (removed — streaming path deleted) ──────────
#
# test_that("updateStream + finalizeStream: composite condition column works end-to-end", { ... })
# test_that("updateStream + finalizeStream: creation/endow condition column works", { ... })

# ── mergeEstimates ────────────────────────────────────────────────────────────

test_that("mergeEstimates: composite condition merges point estimate with uncertainty", {
  df1 <- data.frame(period = 1L:2L, `transTrip_eval` = c(0, 1),
                    changeProb = c(0.3, 0.7), check.names = FALSE)
  df2 <- data.frame(period = 1L:2L, `transTrip_eval` = c(0, 1),
                    Mean = c(0.28, 0.72),     check.names = FALSE)
  result <- with_mocked_bindings(
    mergeEstimates(df1, df2, level = "period", condition = "transTrip_eval"),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
  expect_true("changeProb"     %in% names(result))
  expect_true("Mean"           %in% names(result))
  expect_true("transTrip_eval" %in% names(result))
  expect_equal(nrow(result), 2L)
})

# ── attachContributions: canonical output with flip ───────────────────────────

test_that("attachContributions: flip=TRUE produces canonical columns with sign flip", {
  density <- c(1L, -1L, 0L)
  recip   <- c(2L,  3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density_eval", "recip_eval")))
  result <- attachContributions(data.frame(x = 1:3),
                                 c("density_eval", "recip_eval"),
                                 contrib, flip = TRUE)
  # Canonical names: density_eval → density, recip_eval → recip
  expect_equal(result[["density"]], c(1L, -1L, 0L))  # density keeps ±1 sign
  expect_equal(result[["recip"]],   c(2L, -3L, 4L)) # density==-1 row negated
})

test_that("attachContributions: eval+creation collapse to single canonical column", {
  density  <- c(1L, -1L, 1L)
  recip_e  <- c(1L,  1L, 1L)
  recip_c  <- c(1L,  0L, 1L)   # 0 on deletion rows by construction
  contrib  <- matrix(c(density, recip_e, recip_c), nrow = 3,
                     dimnames = list(NULL,
                       c("density_eval", "recip_eval", "recip_creation")))
  result <- attachContributions(data.frame(x = 1:3),
                                 colnames(contrib), contrib, flip = TRUE)
  # Canonical: recip_eval + recip_creation → single "recip" column
  # Eval column used when hasEval is TRUE (eval-space change statistic)
  # density == -1 row: eval column flipped → -1
  expect_equal(result[["recip"]],   c(1L, -1L, 1L))
  expect_null(result[["recip_eval"]])       # separate columns gone
  expect_null(result[["recip_creation"]])   # separate columns gone
})

# ── creation/endow: utility calculation and theta alignment ──────────────────
#
# These tests guard three distinct failure modes:
#
#   (A) theta name-mismatch: if effectNames contain "_creation" / "_endow"
#       but theta is named without the type suffix, theta[effectNames] returns
#       NAs that silently corrupt the softmax — confirmed here so any future
#       change that reintroduces position-based indexing is caught.
#
#   (B) endow theta is a DISTINCT parameter from eval theta: the utility for a
#       dissolution choice is c_eval*theta_eval + c_endow*theta_endow, NOT
#       (c_eval + c_endow)*theta_eval.  When endow = −eval, the latter collapses
#       to zero regardless of theta_eval, masking the endowment effect entirely.
#
#   (C) predictProbability end-to-end with creation/endow composite effectNames:
#       probabilities must sum to 1 per ego, utility must match manual product,
#       and no NaN/NA should appear.

# Synthetic creation/endow contributions struct (2 egos × 3 choices = 6 rows)
# Effects: density (eval), recip (eval + creation), transTrip (eval + endow)
# Convention: depvarName_shortName_type  (mirrors getStaticChangeContributions returnWide)
#
# C++ convention: on dissolution rows (density=-1), BOTH eval and endow
# contributions = -calculateContribution (same sign, both negative for cc>0).
# So c_endow = c_eval on dissolution rows.
make_cc_cm <- function() {
  eff_names <- c("net_density_eval", "net_recip_eval", "net_recip_creation",
                 "net_transTrip_eval", "net_transTrip_endow")
  # choice: 1=add, -1=drop, 0=no-change  (one of each per ego, 2 egos)
  density_col    <- c( 1L, -1L,  0L,   1L, -1L,  0L)
  recip_eval     <- c( 0L, -1L,  0L,   1L,  0L,  0L)  # dissolution: -cc
  recip_creat    <- c( 0L,  0L,  0L,   1L,  0L,  0L)  # only new-tie rows
  transTrip_eval <- c( 2L, -1L,  0L,   1L, -3L,  0L)  # dissolution: -cc
  transTrip_endow <- c(0L, -1L,  0L,   0L, -3L,  0L)  # endow = eval on dissolution
  mat <- matrix(
    c(density_col, recip_eval, recip_creat, transTrip_eval, transTrip_endow),
    nrow = 6, dimnames = list(NULL, eff_names)
  )
  list(
    contribMat  = mat,
    effectNames = eff_names,
    group_id    = c(1L, 1L, 1L, 2L, 2L, 2L),
    group       = 1L,
    period      = c(1L, 1L, 1L, 1L, 1L, 1L),
    ego         = c(1L, 1L, 1L, 2L, 2L, 2L),
    choice      = density_col
  )
}

theta_cm <- function()
  c(net_density_eval = -2.0, net_recip_eval = 1.5, net_recip_creation = 0.8,
    net_transTrip_eval = 0.5, net_transTrip_endow = -0.3)

# (A) Mismatched theta names → NAs (regression guard: change must be explicit)
test_that("theta name mismatch on creation/endow names returns NAs — not silent wrong values", {
  cc          <- make_cc_cm()
  theta_wrong <- c(net_density = -2, net_recip = 1.5,
                   net_transTrip = 0.5)         # old-style: no _eval / _creation suffix
  theta_use   <- theta_wrong[cc$effectNames]
  expect_true(any(is.na(theta_use)),
    info = "mismatched names must produce NAs, not silently wrong values")
})

# (B) endow theta is separate: endow contribution does NOT cancel eval contribution
test_that("calculateUtility: endow theta is distinct from eval theta — utility is not zero", {
  # One row: choice == -1 (dissolution), endow = -eval by the invariant
  c_eval  <-  2.0
  c_endow <- -2.0
  theta_eval  <-  0.5
  theta_endow <- -0.3   # distinct — if accidentally shared, net would be 0

  eff <- c("net_density_eval", "net_tt_eval", "net_tt_endow")
  mat <- matrix(c(-1L, c_eval, c_endow), nrow = 1,
                dimnames = list(NULL, eff))
  theta <- c(net_density_eval = -2, net_tt_eval = theta_eval,
             net_tt_endow = theta_endow)

  util <- calculateUtility(mat, theta)

  # Correct: -2*(-1) + c_eval*theta_eval + c_endow*theta_endow = 2 + 1.0 - 0.6 = 2.4
  expected <- -2 * (-1L) + c_eval * theta_eval + c_endow * theta_endow
  expect_equal(util, expected, tolerance = 1e-14)

  # What you'd get if theta_endow were accidentally set to theta_eval:
  # c_eval*theta_eval + c_endow*theta_eval = (c_eval + c_endow)*theta_eval = 0 → net = 2
  wrong_util <- as.numeric(mat %*% c(-2, theta_eval, theta_eval))
  expect_false(isTRUE(all.equal(util, wrong_util)),
    info = "endow and eval thetas must produce different utilities")
})

# (C) predictProbability end-to-end with creation/endow composite effectNames
test_that("predictProbability: creation/endow effectNames — probabilities sum to 1, no NaN", {
  cc    <- make_cc_cm()
  theta <- theta_cm()

  result <- predictProbability(cc, theta)

  # Probabilities per ego (softmax group) must sum to 1
  expect_equal(sum(result$changeProb[cc$group_id == 1L]), 1, tolerance = 1e-10)
  expect_equal(sum(result$changeProb[cc$group_id == 2L]), 1, tolerance = 1e-10)

  # Utility must match raw matrix product (the canonical transform
  # with d× is equivalent to the raw c_raw %*% theta)
  expected_util <- as.numeric(cc$contribMat %*% theta[cc$effectNames])
  expect_equal(result$changeUtil, expected_util, tolerance = 1e-12)

  # No NaN or NA
  expect_false(any(is.nan(result$changeProb)))
  expect_false(any(is.na(result$changeProb)))
})

test_that("predictProbability: creation/endow — tieProb type produces NA for no-change rows", {
  cc    <- make_cc_cm()
  theta <- theta_cm()

  result <- predictProbability(cc, theta, type = "tieProb")

  # choice == 0L rows → tieProb must be NA
  no_change <- cc$contribMat[, "net_density_eval"] == 0L
  expect_true(all(is.na(result$tieProb[no_change])))
  # choice == -1L rows → tieProb = 1 - changeProb
  drop_rows <- cc$contribMat[, "net_density_eval"] == -1L
  expect_equal(result$tieProb[drop_rows],
               1 - result$changeProb[drop_rows], tolerance = 1e-14)
})

# ── egoNormalize / preAggEgo / detectEgoUnit ─────────────────────────────────

# Shared fixture: 3 egos with unequal alter counts per condition cell
make_ego_norm_df <- function() {
  data.frame(
    period  = rep(1L, 8),
    ego     = c(1L,1L,1L, 2L,2L,2L, 3L,3L),
    choice  = c(2L,3L,4L, 1L,3L,4L, 1L,2L),
    cond    = c("A","A","B", "A","B","B", "A","B"),
    tieProb = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8),
    stringsAsFactors = FALSE
  )
}

test_that("detectEgoUnit: static data returns group/period/ego", {
  df <- make_ego_norm_df()
  expect_equal(detectEgoUnit(df), c("period", "ego"))
})

test_that("detectEgoUnit: dynamic data includes chain/ministep", {
  df <- make_ego_norm_df()
  df$chain <- 1L
  df$ministep <- seq_len(nrow(df))
  expect_equal(detectEgoUnit(df),
               c("chain", "period", "ego", "ministep"))
})

test_that("agg with egoNormalize=TRUE gives ego-first means", {
  df <- make_ego_norm_df()
  # Ego 1 in A: mean(0.1, 0.2) = 0.15
  # Ego 2 in A: 0.4
  # Ego 3 in A: 0.7
  # Ego-first mean for A: mean(0.15, 0.4, 0.7) = 0.4167
  expected_A <- mean(c(mean(c(0.1, 0.2)), 0.4, 0.7))
  r <- agg("tieProb", df, level = "none", condition = "cond",
           egoNormalize = TRUE)
  actual_A <- r$tieProb[r$cond == "A"]
  expect_equal(actual_A, expected_A, tolerance = 1e-10)
})

test_that("agg egoNormalize=TRUE differs from flat mean when egos have unequal alters", {
  df <- make_ego_norm_df()
  flat    <- agg("tieProb", df, level = "none", condition = "cond",
    egoNormalize = FALSE)
  egofirst <- agg("tieProb", df, level = "none", condition = "cond",
                   egoNormalize = TRUE)
  # Cell A: flat = mean(0.1,0.2,0.4,0.7) = 0.35 vs ego-first = 0.4167
  expect_false(isTRUE(all.equal(flat$tieProb[flat$cond == "A"],
                                 egofirst$tieProb[egofirst$cond == "A"])))
})


test_that("agg egoNormalize works with complex sum_fun (summarizeValue)", {
  df <- make_ego_norm_df()
  r <- agg("tieProb", df, level = "none", condition = "cond",
           sum_fun = summarizeValue, egoNormalize = TRUE)
  expect_true("Mean" %in% names(r))
  expect_true("cases" %in% names(r))
})

test_that("agg egoNormalize: Rcpp grouped_agg_cpp path works correctly", {
  df <- make_ego_norm_df()
  r_result <- agg("tieProb", df, level = "none", condition = "cond",
                   egoNormalize = TRUE)
  expect_true(is.data.frame(r_result))
  expect_true("tieProb" %in% names(r_result))
  expect_true("cond"    %in% names(r_result))
  # Ego-normalize: each ego contributes equally, so within-ego means are
  # averaged across egos → result should differ from raw mean of all rows.
  raw <- agg("tieProb", df, level = "none", condition = "cond",
             egoNormalize = FALSE)
  # At least one condition value should differ
  r_sorted <- r_result[order(r_result$cond), ]
  raw_sorted <- raw[order(raw$cond), ]
  # egoNormalize with extra ego columns should produce different results than raw
  expect_false(identical(r_sorted$tieProb, raw_sorted$tieProb))
})

# ── Two user-specified interactions (unspInt numbering) ───────────────────────
# The model has two unspInt effects; getNamesFromEffects() should number them
# unspInt1 and unspInt2 so they are uniquely addressable.

test_that("resolveEffectName: canonical and alias mapping", {
  effectNames <- c("mynet_transTrip_eval", "mynet_cm_recip_eval")
  expect_equal(resolveEffectName("transTrip_eval", effectNames), "mynet_transTrip_eval")
  expect_equal(resolveEffectName("transTrip", effectNames), "mynet_transTrip_eval")
  expect_equal(resolveEffectName("mynet_cm_recip_eval", effectNames), "mynet_cm_recip_eval")

  effectNames <- c("mynet_transTrip1_eval")
  expect_error(resolveEffectName("transTrip_eval", effectNames),
               "not found")
  expect_equal(resolveEffectName("transTrip1_eval", effectNames), "mynet_transTrip1_eval")

  effectNames <- c("mynet_transTrip_eval")
  expect_equal(resolveEffectName("transTrip", effectNames), "mynet_transTrip_eval")

  effectNames <- c("mynet_transTrip1_eval", "mynet_transTrip2_eval")
  expect_error(resolveEffectName("transTrip_eval", effectNames),
               "not found|ambiguous")
})

# ── getCovCenteringMean ──────────────────────────────────────────────────────

test_that("getCovCenteringMean: returns mean for centered covar", {
  # Build a minimal effects-like data frame and data object
  eff <- data.frame(
    name = "mynet", shortName = c("density", "egoX"),
    interaction1 = c("", "mybeh"), interaction2 = c("", ""),
    type = c("eval", "eval"), include = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  # Simulate a centered covariate
  cov <- 1:10 - 5.5
  attr(cov, "centered") <- TRUE
  attr(cov, "mean") <- 5.5
  dat <- list(cCovars = list(mybeh = cov), vCovars = list(),
              dycCovars = list(), dyvCovars = list())

  resolved <- resolveEffectName("egoX", getEffectNamesNoRate(eff, "mynet"))
  expect_equal(getCovCenteringMean(resolved, eff, dat, "mynet"), 5.5)
})

test_that("getCovCenteringMean: returns 0 for structural effect", {
  eff <- data.frame(
    name = "mynet", shortName = c("density", "recip"),
    interaction1 = c("", ""), interaction2 = c("", ""),
    type = c("eval", "eval"), include = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  dat <- list(cCovars = list(), vCovars = list(),
              dycCovars = list(), dyvCovars = list())

  resolved <- resolveEffectName("recip", getEffectNamesNoRate(eff, "mynet"))
  expect_equal(getCovCenteringMean(resolved, eff, dat, "mynet"), 0)
})

test_that("getCovCenteringMean: returns 0 for uncentered covar", {
  eff <- data.frame(
    name = "mynet", shortName = c("density", "egoX"),
    interaction1 = c("", "mybeh"), interaction2 = c("", ""),
    type = c("eval", "eval"), include = c(TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  cov <- 1:10
  attr(cov, "centered") <- FALSE
  attr(cov, "mean") <- 5.5
  dat <- list(cCovars = list(mybeh = cov), vCovars = list(),
              dycCovars = list(), dyvCovars = list())

  resolved <- resolveEffectName("egoX", getEffectNamesNoRate(eff, "mynet"))
  expect_equal(getCovCenteringMean(resolved, eff, dat, "mynet"), 0)
})