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

# ── attachContribColumns ────────────────────────────────────────────────────

test_that("attachContribColumns: density==-1 flips other columns; density col unchanged", {
  density <- c(1L, -1L, 0L)
  recip   <- c(2L, 3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density", "recip")))
  df <- data.frame(x = 1:3)

  result <- attachContribColumns(df, c("density", "recip"), contrib, flip = TRUE)

  expect_equal(result$density, density)             # density itself not flipped
  expect_equal(result$recip,   c(2L, -3L, 4L))     # density==-1 row → negated
})

test_that("attachContribColumns: flip=FALSE leaves all columns unchanged", {
  density <- c(1L, -1L, -1L)
  recip   <- c(2L, 3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density", "recip")))
  df <- data.frame(x = 1:3)

  result <- attachContribColumns(df, c("density", "recip"), contrib, flip = FALSE)

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
  result <- alignThetaNoRate(theta, c("recip", "transTrip"), ans = NULL)
  expect_named(result, c("recip", "transTrip"))
  expect_equal(unname(result["recip"]), 1.5)
})

test_that("alignThetaNoRate: uses requestedEffects (not effects) for name assignment", {
  # theta includes rate params — this is what sienaMargins passes
  theta <- c(-2.38, 3.06, -0.07)   # rate1, rate2, density, recip, unspInt
  # effectNames are composite (rate effects filtered out — these are the names to extract)
  effectNames <- c("density_eval", "recip_eval", "unspInt_eval")

  mock_ans <- list(
    requestedEffects = data.frame(
      shortName = c("rateX", "rateX", "density", "recip", "unspInt"),
      type      = c("rate",  "rate",  "eval",    "eval",  "eval"),
      include   = c(TRUE,    TRUE,    TRUE,       TRUE,    TRUE),
      stringsAsFactors = FALSE
    )
  )
  result <- alignThetaNoRate(theta, effectNames, mock_ans)
  expect_named(result, effectNames)
  expect_equal(unname(result[["density_eval"]]), -2.38)
  expect_equal(unname(result[["recip_eval"]]),    3.06)
  expect_equal(unname(result[["unspInt_eval"]]), -0.07)
})

test_that("alignThetaNoRate: positional fallback when ans=NULL and theta unnamed", {
  theta <- c(-2, 1.5, 0.5)
  result <- alignThetaNoRate(theta, c("density", "recip", "transTrip"), ans = NULL)
  # returns first length(effectNames) elements positionally
  expect_equal(unname(result), c(-2, 1.5, 0.5))
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

test_that("planBatch: multi-worker result is multiple of nbrNodes", {
  dv <- array(0, dim = c(10, 10, 3))
  mock_data <- list(depvars = list(net = dv))
  result <- planBatch(mock_data, "net", nsim = 50,
                       useCluster = TRUE, nbrNodes = 4, dynamic = FALSE)
  expect_equal(result %% 4L, 0L)
  expect_gte(result, 4L)
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

test_that("predictFirstDiff: output is data.table when package available", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  cc        <- make_cc2()
  theta_use <- c(density = -2, recip = 1.5)

  result <- predictFirstDiff(
    changeContributions = cc, theta_use = theta_use, type = "changeProb",
    effectName = "recip", diff = 1, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
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
        interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
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

# ── makeGroupKey ─────────────────────────────────────────────────────────────

test_that("makeGroupKey: composite column name indexes correctly", {
  df <- data.frame(period = c(1L, 1L, 2L), `transTrip_eval` = c(0, 0, 1),
                   check.names = FALSE)
  keys <- makeGroupKey(df, c("period", "transTrip_eval"))
  expect_length(keys, 3L)
  expect_equal(keys[1], keys[2])     # same period AND same val → same key
  expect_false(keys[1] == keys[3])   # different period AND different val → different key
})

test_that("makeGroupKey: empty group_vars returns __all__", {
  keys <- makeGroupKey(data.frame(x = 1:3), character(0))
  expect_true(all(keys == "__all__"))
})

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
  result <- with_mocked_bindings(
    agg("changeProb", df, level = "period", condition = "transTrip_eval"),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
  expect_true("transTrip_eval" %in% names(result))
  expect_equal(nrow(result), 3L)   # period + transTrip_eval + changeProb value
})

test_that("agg: creation/endow composite condition aggregates correctly", {
  df <- data.frame(period = c(1L, 1L, 2L, 2L),
                   `recip_creation` = c(0L, 1L, 0L, 1L),
                   changeProb = c(0.2, 0.8, 0.3, 0.7),
                   check.names = FALSE)
  result <- with_mocked_bindings(
    agg("changeProb", df, level = "period", condition = "recip_creation"),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
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
  agg_a <- with_mocked_bindings(
    agg("changeProb", df_a, level = "period", condition = "density_eval"),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
  agg_b <- with_mocked_bindings(
    agg("changeProb", df_b, level = "period", condition = "density_eval"),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
  expect_false(identical(agg_a$changeProb, agg_b$changeProb))
})

# ── updateStream / finalizeStream ─────────────────────────────────────────────

test_that("updateStream + finalizeStream: composite condition column works end-to-end", {
  sim_df <- data.frame(period = c(1L, 1L, 2L, 2L),
                       `transTrip.eval` = c(0, 1, 0, 1),
                       changeProb = c(0.3, 0.7, 0.4, 0.6),
                       check.names = FALSE)
  group_vars   <- c("period", "transTrip.eval")
  stream_state <- new.env(parent = emptyenv(), hash = TRUE)
  group_state  <- new.env(parent = emptyenv(), hash = TRUE)

  # Two draws of the same result (simple smoke test)
  with_mocked_bindings(
    {
      updateStream(stream_state, group_state, sim_df, "changeProb", group_vars)
      updateStream(stream_state, group_state, sim_df, "changeProb", group_vars)
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )

  result <- with_mocked_bindings(
    finalizeStream(stream_state, group_state, group_vars,
                   uncertainty_summary_fun = function(x, na.rm) list(Mean = mean(x))),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )

  expect_true(is.data.frame(result))
  expect_true("transTrip.eval" %in% names(result))
  expect_true("Mean"           %in% names(result))
  # 2 periods × 2 condition values = 4 groups
  expect_equal(nrow(result), 4L)
  # Mean of two identical draws == draw itself
  expect_true(all(result$Mean %in% c(0.3, 0.7, 0.4, 0.6)))
})

test_that("updateStream + finalizeStream: creation/endow condition column works", {
  sim_df <- data.frame(period = c(1L, 1L),
                       `recip.creation` = c(0L, 1L),
                       changeProb = c(0.25, 0.75),
                       check.names = FALSE)
  group_vars   <- c("period", "recip.creation")
  stream_state <- new.env(parent = emptyenv(), hash = TRUE)
  group_state  <- new.env(parent = emptyenv(), hash = TRUE)

  with_mocked_bindings(
    updateStream(stream_state, group_state, sim_df, "changeProb", group_vars),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )

  result <- with_mocked_bindings(
    finalizeStream(stream_state, group_state, group_vars,
                   uncertainty_summary_fun = function(x, na.rm) list(Mean = mean(x))),
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )

  expect_true("recip.creation" %in% names(result))
  expect_equal(nrow(result), 2L)
})

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

# ── attachContribColumns: composite density name ─────────────────────────────

test_that("attachContribColumns: composite 'density_eval' column triggers flip", {
  density <- c(1L, -1L, 0L)
  recip   <- c(2L,  3L, 4L)
  contrib <- matrix(c(density, recip), nrow = 3,
                    dimnames = list(NULL, c("density_eval", "recip_eval")))
  result <- attachContribColumns(data.frame(x = 1:3),
                                 c("density_eval", "recip_eval"),
                                 contrib, flip = TRUE)
  # grepl("density", "density_eval") is TRUE → flip applies
  expect_equal(result[["density_eval"]], density)        # density col unchanged
  expect_equal(result[["recip_eval"]],   c(2L, -3L, 4L))# density==-1 row negated
})

test_that("attachContribColumns: 'recip_creation' not confused with density trigger", {
  density  <- c(1L, -1L, 1L)
  recip_e  <- c(1L,  1L, 1L)
  recip_c  <- c(1L,  0L, 1L)   # 0 on deletion rows by construction
  contrib  <- matrix(c(density, recip_e, recip_c), nrow = 3,
                     dimnames = list(NULL,
                       c("density_eval", "recip_eval", "recip_creation")))
  result <- attachContribColumns(data.frame(x = 1:3),
                                 colnames(contrib), contrib, flip = TRUE)
  expect_equal(result[["recip_eval"]],     c(1L, -1L, 1L))
  expect_equal(result[["recip_creation"]], c(1L,  0L, 1L))  # 0 negated == 0
})
