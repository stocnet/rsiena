# testthat::skip_on_cran()

# Unit tests for pipeline sub-functions in predict.R, sienaMargins.r, postestimate.R.
# These tests use synthetic data: no siena07() call required.
# Goal: catch errors early in the computational pipeline before they surface
# at the integration level (sienaAME / predict.sienaFit).

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
    contrib_mat = mat,
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
    contrib_mat = mat,
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
  mat   <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3)
  theta <- c(0.5, -1)
  expect_equal(calculateUtility(mat, theta), as.numeric(mat %*% theta))
})

test_that("calculateUtility: single-column matrix", {
  mat   <- matrix(c(1, -1, 0), nrow = 3)
  theta <- c(-2)
  expect_equal(calculateUtility(mat, theta), c(-2, 2, 0))
})

test_that("calculateUtility: requires matrix input", {
  expect_error(calculateUtility(c(1, 2, 3), c(1)))
})

# ── .computeTieProb ───────────────────────────────────────────────────────────

test_that(".computeTieProb: returns NULL for changeProb type", {
  expect_null(.computeTieProb(c(0.3, 0.7, 0.5), c(1L, -1L, 0L), "changeProb"))
})

test_that(".computeTieProb: density==-1 flips, density==0 → NA, density==1 unchanged", {
  prob        <- c(0.3, 0.7, 0.5)
  density_col <- c(1L, -1L, 0L)
  result <- .computeTieProb(prob, density_col, "tieProb")
  expect_equal(result[1], 0.3)           # density==1:  no flip
  expect_equal(result[2], 1 - 0.7)      # density==-1: flipped
  expect_true(is.na(result[3]))          # density==0:  NA
})

test_that(".computeTieProb: all density==1 → no change", {
  prob <- runif(5)
  result <- .computeTieProb(prob, rep(1L, 5), "tieProb")
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
    int_effectNames = "unspInt", mod_effectNames = "transTrip",
    modContribution = modContrib,
    effectNames = effectNames
  )
  # density * (diff*theta[recip] + diff*density*modContrib*theta[unspInt])
  expected <- densityValue * (1 * theta["recip"] + 1 * densityValue * modContrib * theta["unspInt"])
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
  tieProb_in    <- .computeTieProb(changeProb, densityValue, "tieProb")

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

# ── .attachContribColumns ────────────────────────────────────────────────────

test_that(".attachContribColumns: density==-1 flips other columns; density col unchanged", {
  density <- c(1L, -1L, 0L)
  recip   <- c(2L, 3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density", "recip")))
  df <- data.frame(x = 1:3)

  result <- .attachContribColumns(df, c("density", "recip"), contrib, flip = TRUE)

  expect_equal(result$density, density)             # density itself not flipped
  expect_equal(result$recip,   c(2L, -3L, 4L))     # density==-1 row → negated
})

test_that(".attachContribColumns: flip=FALSE leaves all columns unchanged", {
  density <- c(1L, -1L, -1L)
  recip   <- c(2L, 3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density", "recip")))
  df <- data.frame(x = 1:3)

  result <- .attachContribColumns(df, c("density", "recip"), contrib, flip = FALSE)

  expect_equal(result$recip, c(2L, 3L, 4L))        # no flipping
})

# ── .groupColsList ───────────────────────────────────────────────────────────

test_that(".groupColsList: static (no chain) returns correct names and values", {
  pb <- list(
    group = 1L, period = c(1L, 2L, 3L),
    ego = c(5L, 6L, 7L), choice = c(1L, 0L, 1L)
  )
  result <- .groupColsList(pb)
  expect_named(result, c("group", "period", "ego", "choice"))
  expect_equal(result$ego, c(5L, 6L, 7L))
})

test_that(".groupColsList: static with keep subsets correctly", {
  pb <- list(
    group = 1L, period = c(1L, 2L, 3L),
    ego = c(5L, 6L, 7L), choice = c(1L, 0L, 1L)
  )
  result <- .groupColsList(pb, keep = c(TRUE, FALSE, TRUE))
  expect_equal(result$ego,    c(5L, 7L))
  expect_equal(result$period, c(1L, 3L))
})

test_that(".groupColsList: dynamic (with chain) returns correct names", {
  pb <- list(
    chain = c(1L, 1L, 2L), group = c(1L, 1L, 1L),
    period = c(1L, 2L, 1L), ministep = c(10L, 20L, 30L),
    choice = c(0L, 1L, 1L)
  )
  result <- .groupColsList(pb)
  expect_named(result, c("chain", "group", "period", "ministep", "choice"))
  expect_equal(result$ministep, c(10L, 20L, 30L))
})

# ── alignThetaNoRate ─────────────────────────────────────────────────────────

test_that("alignThetaNoRate: named theta — selects by name", {
  theta <- c(density = -2, recip = 1.5, transTrip = 0.5)
  result <- alignThetaNoRate(theta, c("recip", "transTrip"), ans = NULL)
  expect_named(result, c("recip", "transTrip"))
  expect_equal(unname(result["recip"]), 1.5)
})

test_that("alignThetaNoRate: uses requestedEffects (not effects) for name assignment", {
  theta <- c(2.5, 3.0, -2.38, 3.06, -0.07)
  effectNames <- c("density", "recip", "unspInt")
  mock_ans <- list(
    effects          = list(shortName = c("rateX", "rateX", "density", "recip", "outPop", "unspInt")),  # length 6 ≠ 5
    requestedEffects = list(shortName = c("rateX", "rateX", "density", "recip", "unspInt"))             # length 5 == 5
  )
  result <- alignThetaNoRate(theta, effectNames, mock_ans)
  expect_named(result, effectNames)
  expect_equal(unname(result["density"]), -2.38)
  expect_equal(unname(result["recip"]),    3.06)
  expect_equal(unname(result["unspInt"]), -0.07)
})

test_that("alignThetaNoRate: positional fallback when ans=NULL and theta unnamed", {
  theta <- c(-2, 1.5, 0.5)
  result <- alignThetaNoRate(theta, c("density", "recip", "transTrip"), ans = NULL)
  # returns first length(effectNames) elements positionally
  expect_equal(unname(result), c(-2, 1.5, 0.5))
})

# ── plan_batch ───────────────────────────────────────────────────────────────

test_that("plan_batch: result is always >= 1 and <= nsim", {
  dv <- array(0, dim = c(50, 50, 3))
  mock_data <- list(depvars = list(mynet = dv))
  result <- plan_batch(mock_data, "mynet", nsim = 100, dynamic = FALSE)
  expect_gte(result, 1L)
  expect_lte(result, 100L)
})

test_that("plan_batch: nsim=1 always returns 1", {
  dv <- array(0, dim = c(50, 50, 3))
  mock_data <- list(depvars = list(mynet = dv))
  expect_equal(plan_batch(mock_data, "mynet", nsim = 1L), 1L)
})

test_that("plan_batch: multi-worker result is multiple of nbrNodes", {
  dv <- array(0, dim = c(10, 10, 3))
  mock_data <- list(depvars = list(net = dv))
  result <- plan_batch(mock_data, "net", nsim = 50,
                       useCluster = TRUE, nbrNodes = 4, dynamic = FALSE)
  expect_equal(result %% 4L, 0L)
  expect_gte(result, 4L)
})

test_that("plan_batch: memory_scale reduces batch size", {
  dv <- array(0, dim = c(50, 50, 3))
  mock_data <- list(depvars = list(mynet = dv))
  b1 <- plan_batch(mock_data, "mynet", nsim = 200, dynamic = FALSE)
  b2 <- plan_batch(mock_data, "mynet", nsim = 200, dynamic = FALSE, memory_scale = 4L)
  expect_lte(b2, b1)
})

# ── predictFirstDiff ─────────────────────────────────────────────────────────

test_that("predictFirstDiff: density==0 rows excluded; contrib cols present", {
  cc         <- make_cc2()
  theta_use  <- c(density = -2, recip = 1.5)
  n_valid    <- sum(cc$contrib_mat[, "density"] != 0L)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
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
    interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
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
    interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
    details = TRUE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true("tieProb" %in% names(result))
})

test_that("predictFirstDiff: output is data.table when package available", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  cc        <- make_cc2()
  theta_use <- c(density = -2, recip = 1.5)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true(data.table::is.data.table(result))
})

test_that("predictFirstDiff: base R output is plain data.frame", {
  with_mocked_bindings(
    {
      cc        <- make_cc2()
      theta_use <- c(density = -2, recip = 1.5)

      result <- predictFirstDiff(
        changeContributions = cc, theta_use = theta_use, type = "changeProb",
        effectName = "recip", diff = 1, contrast = NULL,
        interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
        details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
      )
      expect_true(is.data.frame(result))
      expect_false("data.table" %in% class(result))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

# ── predictSecondDiff ────────────────────────────────────────────────────────

test_that("predictSecondDiff: density==0 rows excluded; secondDiff col present", {
  cc         <- make_cc3()
  theta_use  <- c(density = -2, recip = 1.5, transTrip = 0.5)
  n_valid    <- sum(cc$contrib_mat[, "density"] != 0L)

  result <- predictSecondDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName1 = "recip",     diff1 = 1, contrast1 = NULL,
    interaction1 = FALSE, int_effectNames1 = NULL, mod_effectNames1 = NULL,
    effectName2 = "transTrip", diff2 = 1, contrast2 = NULL,
    interaction2 = FALSE, int_effectNames2 = NULL, mod_effectNames2 = NULL,
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
    interaction1 = FALSE, int_effectNames1 = NULL, mod_effectNames1 = NULL,
    effectName2 = "transTrip", diff2 = 1, contrast2 = NULL,
    interaction2 = FALSE, int_effectNames2 = NULL, mod_effectNames2 = NULL,
    details = TRUE, calcRiskRatio = FALSE, mainEffect = "riskDifference"
  )
  expect_true(all(c("changeUtil", "changeProb") %in% names(result)))
})
