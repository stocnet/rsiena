# testthat::skip_on_cran()

# -- Load fixtures from cache (built by helper-models.R) --
ans        <- load_fixture("ans")
mydata     <- load_fixture("mydata")
mymodel    <- load_fixture("mymodel")
mycontrols <- load_fixture("mycontrols")
# Full-mode fixtures (NULL in quick mode; those tests are guarded by skip_slow())
ans_int        <- load_fixture("ans_int")
mydata_int     <- load_fixture("mydata_int")
mymodel_int    <- load_fixture("mymodel_int")
mycontrols_int <- load_fixture("mycontrols_int")
ans_int_uncond <- load_fixture("ans_int_uncond")
ans_co         <- load_fixture("ans_co")
mydata_co      <- load_fixture("mydata_co")
mymodel_co     <- load_fixture("mymodel_co")
mycontrols_co  <- load_fixture("mycontrols_co")
ans_cm         <- load_fixture("ans_cm")
mydata_cm      <- load_fixture("mydata_cm")
mymodel_cm     <- load_fixture("mymodel_cm")

test_that("predict.sienaFit static: no uncertainty, egoChoice level", {
  pred_df <- predict(
    object = ans,
    newdata = mydata,
    type = "tieProb",
    nsim = 10,
    condition = "transTrip_eval",
    level = "egoChoice"
  )
  expect_true(is.data.frame(pred_df))
})

test_that("predict.sienaFit static: MCSE columns, no CI", {
  skip_slow()
  pred_df <- predict(
    object = ans,
    newdata = mydata,
    type = "tieProb",
    nsim = 20,
    condition = "transTrip",
    level = "period",
    uncertainty = TRUE,
    uncertaintyMean = TRUE,
    uncertaintyMcse = TRUE,
    uncertaintymcseBatches = 4,
    uncertaintySd = TRUE,
    uncertaintyCi = FALSE,
    verbose = FALSE
  )
  expect_true(is.data.frame(pred_df))
  expect_true("mcse_Mean" %in% names(pred_df))
  expect_true("mcse_SE" %in% names(pred_df))
  expect_false("q_025" %in% names(pred_df))
  expect_false("q_975" %in% names(pred_df))
  expect_false("Median" %in% names(pred_df))
})

test_that("predict.sienaFit static: sienaPrediction class", {
  pred_dt <- predict.sienaFit(
    object = ans,
    newdata = mydata,
    type = "tieProb",
    nsim = 10,
    condition = "transTrip_eval",
    level = "period",
    uncertainty = TRUE
  )
  expect_true(is.data.frame(pred_dt))
  expect_s3_class(pred_dt, "sienaPrediction")
})

test_that("predict.sienaFit static: FORK cluster", {
  skip_slow()
  pred_dt <- predict(
    object = ans,
    newdata = mydata,
    type = "tieProb",
    nsim = 10,
    condition = "transTrip",
    level = "period",
    uncertainty = TRUE,
    useCluster = TRUE,
    clusterType = "FORK",
    nbrNodes = 2
  )
  expect_null(parallel:::getDefaultCluster())
  expect_true(is.data.frame(pred_dt))
})

test_that("predict.sienaFit dynamic: FORK cluster", {
  skip_slow()
  pred_df <- predict(
        object = ans,
        newdata = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        dynamic = TRUE,
        type = "tieProb",
        n3 = 60,
        nsim = 2,
        condition = "transTrip",
        level = "period",
        uncertainty = TRUE,
        useCluster = TRUE,
        clusterType = "FORK",
    nbrNodes = 2,
    silent = FALSE
  )
  expect_null(parallel:::getDefaultCluster())
  expect_true(is.data.frame(pred_df))
})

test_that("predict.sienaFit dynamic: basic", {
  skip_slow()
  pred_df <- predict(
        object = ans,
        newdata = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        type = "tieProb",
        n3 = 60,
        nsim = 4,
    condition = "density_eval"
  )
  expect_true(is.data.frame(pred_df))
})

test_that("predict.sienaFit dynamic: MCSE, no SD", {
  skip_slow()
  pred_df <- predict(
        object = ans,
        dynamic = TRUE,
        newdata = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        type = "tieProb",
        n3 = 60,
        nsim = 20,
        condition = "density_eval",
        uncertainty = TRUE,
        uncertaintyMean = TRUE,
        uncertaintyMcse = TRUE,
        uncertaintymcseBatches = 4,
        uncertaintySd = FALSE,
        uncertaintyCi = TRUE,
    verbose = FALSE
  )
  expect_true(is.data.frame(pred_df))
  expect_true("mcse_Mean" %in% names(pred_df))
  expect_false("SE" %in% names(pred_df))
  expect_false("mcse_SE" %in% names(pred_df))
  expect_true("q_025" %in% names(pred_df))
  expect_true("q_975" %in% names(pred_df))
})

test_that("predict.sienaFit dynamic: returns data.frame with uncertainty", {
  skip_slow()
  pred_dt <- predict(
    object = ans,
    dynamic = TRUE,
    newdata = mydata,
    effects = mymodel,
    algorithm = mycontrols,
    type = "tieProb",
    n3 = 60,
    nsim = 4,
    condition = "density_eval"
  )
  expect_true(is.data.frame(pred_dt))
})

# Interaction without main effect — model from helper-models.R (full mode):
#   ans_int, mydata_int, mymodel_int, mycontrols_int  (depvar = "mynet_int")
test_that("predict.sienaFit: interactions without main effect, no uncertainty", {
  skip_slow()
  pred_dt <- predict.sienaFit(ans_int,
    newdata =  mydata_int,
    effects = mymodel_int,
    level = "period",
    condition = "recip_eval",
    uncertainty = FALSE
  )
  expect_true(is.data.frame(pred_dt))
})

test_that("predict.sienaFit: interactions without main effect, with uncertainty", {
  skip_slow()
  pred_dt <- predict.sienaFit(ans_int,
    newdata =  mydata_int,
    effects = mymodel_int,
    level = "period",
    condition = "recip_eval",
    uncertainty = TRUE,
    nsim = 5
  )
  expect_true(is.data.frame(pred_dt))
})

# cond=FALSE + interaction without main effect: this is the critical case where
# alignThetaNoRate must use requestedEffects (not effects) to avoid picking up
# injected base-effect slots in theta, which would give wrong (uniform) probs.
# Model from helper-models.R (full mode): ans_int_uncond, mydata_int, mymodel_int

test_that("predict.sienaFit: interaction without main effect, cond=FALSE, no NaN", {
  skip_slow()
  pred_dt <- predict.sienaFit(ans_int_uncond, newdata = mydata_int,
    effects = mymodel_int, uncertainty = FALSE, level = "egoChoice")
  expect_true(is.data.frame(pred_dt))
  expect_true(nrow(pred_dt) > 0)
  expect_false(any(is.nan(pred_dt$changeProb)))
})

# Unit test for nameThetaFromEffects + alignThetaNoRate: verify that the
# two-step naming approach handles cond=FALSE, interaction-without-main correctly.
test_that("nameThetaFromEffects + alignThetaNoRate: cond=FALSE interaction scenario", {
  # Reproduce the cond=FALSE + interaction-without-main-effect scenario:
  #   ans$requestedEffects has 5 rows (2 rate, 3 non-rate — no injected outPop).
  #   ans$theta has 5 elements (with rate params for cond=FALSE).
  #   After nameThetaFromEffects, theta carries shortName_type names.
  #   alignThetaNoRate then selects by effectNames from contributions.
  theta_uncond <- c(5.1, 4.9, -2.38, 3.06, -0.07)  # rate1, rate2, density, recip, unspInt
  effectNames  <- c("density_eval", "recip_eval", "unspInt_eval")
  req_effs <- data.frame(
    shortName = c("rateX", "rateX", "density", "recip", "unspInt"),
    type      = c("rate",  "rate",  "eval",    "eval",  "eval"),
    include   = rep(TRUE, 5L),
    stringsAsFactors = FALSE
  )
  theta_named <- expect_warning(
    nameThetaFromEffects(theta_uncond, req_effs),
    "theta has no names"
  )
  result <- alignThetaNoRate(theta_named, effectNames)
  # Must return the three eval params by composite name
  expect_named(result, effectNames)
  expect_equal(unname(result[["density_eval"]]),  -2.38)
  expect_equal(unname(result[["recip_eval"]]),     3.06)
  expect_equal(unname(result[["unspInt_eval"]]),  -0.07)
})

test_that("numberIntShortNames + getNamesFromEffects handle same-identity interactions", {
  e_single <- data.frame(
    name = rep("mynet", 3),
    shortName = rep("unspInt", 3),
    interaction1 = rep("recip", 3),
    interaction2 = rep("inPop", 3),
    type = c("eval", "creation", "endow"),
    include = rep(TRUE, 3),
    stringsAsFactors = FALSE
  )
  # same logical interaction should not be numbered
  names_single <- getNamesFromEffects(e_single)
  expect_true(any(grepl("unspInt_recipinPop_eval$", names_single)))
  expect_true(any(grepl("unspInt_recipinPop_creation$", names_single)))
  expect_true(any(grepl("unspInt_recipinPop_endow$", names_single)))
  expect_false(any(grepl("unspInt1", names_single)))

  e_double <- data.frame(
    name = rep("mynet", 6),
    shortName = rep("unspInt", 6),
    interaction1 = c(rep("recip", 3), rep("outPop", 3)),
    interaction2 = rep("inPop", 6),
    type = rep(c("eval", "creation", "endow"), 2),
    include = rep(TRUE, 6),
    stringsAsFactors = FALSE
  )
  names_double <- getNamesFromEffects(e_double)
  expect_true(any(grepl("unspInt1_recipinPop_eval$", names_double)))
  expect_true(any(grepl("unspInt2_outPopinPop_eval$", names_double)))
  expect_true(any(grepl("unspInt1_recipinPop_creation$", names_double)))
  expect_true(any(grepl("unspInt2_outPopinPop_endow$", names_double)))
})

# ---------------------------------------------------------------------------
# Two-network co-evolution: verify predict.sienaFit gives separate predictions #
# per network, and that alignThetaNoRate correctly matches parameters to
# effects by name (not position) when cond=FALSE + interaction without main effect
# Model from helper-models.R (full mode): ans_co, mydata_co, mymodel_co, mycontrols_co
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Static predict: co-evolution
# ---------------------------------------------------------------------------
test_that("predict.sienaFit static: co-evolution networks get separate predictions", {
  skip_slow()
  pred_a <- predict.sienaFit(ans_co, newdata = mydata_co, depvar = "mynet_a",
                              effects = mymodel_co, uncertainty = FALSE)
  pred_b <- predict.sienaFit(ans_co, newdata = mydata_co, depvar = "mynet_b",
                              effects = mymodel_co, uncertainty = FALSE)
  # Different networks → different predicted probabilities
  expect_false(identical(pred_a$changeProb, pred_b$changeProb))
  # Basic output structure present
  expect_true("period"     %in% names(pred_a))
  expect_true("changeProb" %in% names(pred_a))
})

# ---------------------------------------------------------------------------
# Dynamic predict: co-evolution
# ---------------------------------------------------------------------------
test_that("predict.sienaFit dynamic: co-evolution networks get separate predictions", {
  skip_slow()
  pred_a <- predict.sienaFit(ans_co, newdata = mydata_co, depvar = "mynet_a",
                              dynamic = TRUE, algorithm = mycontrols_co,
                              effects = mymodel_co, n3 = 50, nsim = 5,
                              uncertainty = FALSE)
  pred_b <- predict.sienaFit(ans_co, newdata = mydata_co, depvar = "mynet_b",
                              dynamic = TRUE, algorithm = mycontrols_co,
                              effects = mymodel_co, n3 = 50, nsim = 5,
                              uncertainty = FALSE)
  expect_false(identical(pred_a$changeProb, pred_b$changeProb))
})

# ---------------------------------------------------------------------------
# Creation / endowment effects: non-unique shortNames *within* a single depvar
#
# recip (eval, included by default) + recip (creation) → shortName "recip" twice
# transTrip (eval) + transTrip (endow) → shortName "transTrip" twice
#
# This tests that alignThetaNoRate correctly matches parameters to effects by name
# (not position) when cond=FALSE + interaction without main effect, so that the
# right parameters get applied to the right effects for utility and probability calcs.
# Model from helper-models.R (full mode): ans_cm, mydata_cm, mymodel_cm, mycontrols_cm
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Static predict: creation/maintenance — all effect types present, distinct
# ---------------------------------------------------------------------------
test_that("predict.sienaFit static: creation/endow effects correctly separated", {
  skip_slow()
  # Test via contributions struct (static)
  wide_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm,
    depvar = "mynet_cm", returnWide = TRUE
  )
  # All five effect slots appear as distinct columns (depvar-prefixed)
  expect_true("mynet_cm_density_eval"    %in% wide_cm$effectNames)
  expect_true("mynet_cm_recip_eval"      %in% wide_cm$effectNames)
  expect_true("mynet_cm_recip_creation"  %in% wide_cm$effectNames)
  expect_true("mynet_cm_transTrip_eval"  %in% wide_cm$effectNames)
  expect_true("mynet_cm_transTrip_endow" %in% wide_cm$effectNames)
  # eval and creation columns must differ (different zero-masks)
  recip_eval_col  <- which(wide_cm$effectNames == "mynet_cm_recip_eval")
  recip_creat_col <- which(wide_cm$effectNames == "mynet_cm_recip_creation")
  expect_false(identical(wide_cm$contribMat[, recip_eval_col],
                         wide_cm$contribMat[, recip_creat_col]))
  pred_cm <- predictProbability(wide_cm, theta = c(ans_cm$rate, ans_cm$theta), type = "tieProb")
  
})

# ---------------------------------------------------------------------------
# Dynamic predict: creation/maintenance
# ---------------------------------------------------------------------------
test_that("predict.sienaFit dynamic: creation/endow effect columns present", {
  skip_slow()
  # Test via the contributions struct (which carries the column names directly)
  # rather than predict.sienaFit whose aggregated output does not include them.
  theta_cm <- c(ans_cm$rate, ans_cm$theta)
  contrib <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = theta_cm, effects = mymodel_cm, depvar = "mynet_cm",
    returnWide = TRUE, useChangeContributions = TRUE
  )
  expect_true("mynet_cm_recip_eval"      %in% contrib$effectNames)
  expect_true("mynet_cm_recip_creation"  %in% contrib$effectNames)
  expect_true("mynet_cm_transTrip_eval"  %in% contrib$effectNames)
  expect_true("mynet_cm_transTrip_endow" %in% contrib$effectNames)
  # eval and creation columns must differ (creation is 0 for existing ties)
  recip_eval_col    <- which(contrib$effectNames == "mynet_cm_recip_eval")
  recip_creat_col   <- which(contrib$effectNames == "mynet_cm_recip_creation")
  expect_false(identical(contrib$contribMat[, recip_eval_col],
                         contrib$contribMat[, recip_creat_col]))
})

# ---------------------------------------------------------------------------
# Dynamic: R-recomputed probabilities match C++ chain-stored probabilities
# ---------------------------------------------------------------------------
test_that("Dynamic R-recomputed probs match C++ chain-stored probs", {
  skip_slow()
  # Rerun siena07 with returnChangeContributions=TRUE to get chain-stored
  # utilities and probabilities (the fixture was built with returnDataFrame=TRUE
  # which doesn't preserve chain attributes).
  theta_all <- c(ans$rate, ans$theta)
  contrib <- getDynamicChangeContributions(
    ans = ans, data = mydata, algorithm = mycontrols,
    theta = theta_all, effects = mymodel, depvar = "mynet",
    returnWide = TRUE, useChangeContributions = FALSE,
    n3 = 10, batch = TRUE, silent = TRUE
  )
  # Chain-stored values
  chain_util <- contrib$changeUtility
  chain_prob <- contrib$changeProbability

  # R-recomputed values
  theta_use <- theta_all[contrib$effectNames]
  r_util <- calculateUtility(contrib$contribMat, theta_use, contrib$permitted)
  r_prob <- calculateChangeProb(r_util, contrib$group_id)

  # Compare on permitted rows only
  perm <- contrib$permitted
  expect_equal(length(chain_util), length(r_util))
  expect_equal(length(chain_prob), length(r_prob))

  # Utilities should match closely on permitted rows
  expect_true(all(is.finite(r_util[perm])))
  max_util_diff <- max(abs(chain_util[perm] - r_util[perm]))
  expect_lt(max_util_diff, 1e-8,
    label = paste("max utility diff =", max_util_diff))

  # Probabilities should match closely on permitted rows
  max_prob_diff <- max(abs(chain_prob[perm] - r_prob[perm]))
  expect_lt(max_prob_diff, 1e-8,
    label = paste("max prob diff =", max_prob_diff))
})

# ---------------------------------------------------------------------------
# Integration: theta sign, utility diff, and predicted probability direction
#
# Uses ans / mydata / mymodel from helper-models.R (base, always available).
# No additional model fitting — getStaticChangeContributions just restructures
# data already stored in ans$changeContributions, so this is always fast.
#
# What is tested:
#   (1) Rate parameters are NOT used in utility computation  — theta is aligned
#       by name (composite effectName), so rate params (unnamed / differently
#       named) must not appear in the dot product.
#   (2) Utility difference is exact: ΔU = Δθ_transTrip × c_transTrip per row.
#       This would silently fail if the wrong theta element is used.
#   (3) Sign is correct: higher θ_transTrip → higher changeProb for additions
#       (density==1) where c_transTrip > 0.
# ---------------------------------------------------------------------------

test_that("predict: theta sign, utility diff, and probability direction are correct", {
  # Build the contribution struct once (pure data restructuring, no simulation).
  contrib <- getStaticChangeContributions(
    ans = ans, data = mydata, effects = mymodel,
    depvar = "mynet", returnWide = TRUE
  )

  tt_col      <- grep("transTrip", contrib$effectNames, fixed = TRUE)
  density_col <- grep("density",   contrib$effectNames, fixed = TRUE)[1L]
  stopifnot(length(tt_col) == 1L)   # model has exactly one transTrip effect

  theta_base  <- ans$theta
  names(theta_base) <- contrib$effectNames  # composite names expected by predictProbability

  # Perturb transTrip upward by a large, clearly detectable amount.
  delta_tt    <- 5.0
  theta_high  <- theta_base
  theta_high[tt_col] <- theta_base[tt_col] + delta_tt

  out_base <- predictProbability(contrib, theta_base)
  out_high <- predictProbability(contrib, theta_high)

  # (1) utility difference is exactly delta_tt * c_transTrip for every row
  util_diff_expected <- delta_tt * contrib$contribMat[, tt_col]
  util_diff_actual   <- out_high$changeUtil - out_base$changeUtil
  expect_equal(util_diff_actual, util_diff_expected, tolerance = 1e-12,
    info = "ΔU must equal Δθ * contribution — fails if wrong theta element is used")

  # (2) Softmax monotonicity: E[c_transTrip] (changeProb-weighted) must be
  #     >= per ego group under theta_high, and strictly > in aggregate.
  #     This is the fundamental property dE[c]/dtheta = Var_theta(c) >= 0.
  #     Per-ego equality happens when all choices have identical transTrip
  #     contributions (Var = 0); the aggregate must strictly increase because
  #     at least some egos have non-constant transTrip contributions.
  #     Individual row changeProb can go either way due to softmax competition.
  tt_contrib <- contrib$contribMat[, tt_col]
  for (gid in unique(contrib$group_id)) {
    rows   <- contrib$group_id == gid
    E_base <- sum(out_base$changeProb[rows] * tt_contrib[rows])
    E_high <- sum(out_high$changeProb[rows] * tt_contrib[rows])
    expect_gte(E_high, E_base - 1e-12,
      label = paste0("E[c_transTrip] must not decrease under theta_high (ego group ", gid, ")"))
  }
  # Aggregate must be strictly higher (at least one ego has Var(c_transTrip) > 0)
  expect_gt(
    sum(out_high$changeProb * tt_contrib),
    sum(out_base$changeProb * tt_contrib),
    label = "aggregate E[c_transTrip] must strictly increase under theta_high"
  )

  # (3) tieProb direction: mean tieProb over rows with positive transTrip
  #     contribution must be higher under theta_high.
  #     changeProb for individual rows can be lower (softmax competition), but
  #     tieProb aggregated over pairs where transitivity helps must increase.
  out_base_tp <- predictProbability(contrib, theta_base, type = "tieProb")
  out_high_tp <- predictProbability(contrib, theta_high, type = "tieProb")
  pos_tt_rows <- tt_contrib > 0 & !is.na(out_base_tp$tieProb)
  expect_gt(length(which(pos_tt_rows)), 0L,
    label = "must have non-NA tieProb rows with positive transTrip contribution")
  expect_gt(
    mean(out_high_tp$tieProb[pos_tt_rows]),
    mean(out_base_tp$tieProb[pos_tt_rows]),
    label = "mean tieProb for positive-transTrip pairs must increase under theta_high"
  )

  # (4) probabilities still sum to 1 per ego under both theta values
  for (gid in unique(contrib$group_id)) {
    rows <- contrib$group_id == gid
    expect_equal(sum(out_base$changeProb[rows]), 1, tolerance = 1e-10)
    expect_equal(sum(out_high$changeProb[rows]), 1, tolerance = 1e-10)
  }
})