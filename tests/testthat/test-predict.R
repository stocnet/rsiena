# testthat::skip_on_cran()

# library(testthat)
# library(RSiena)

# Minimal RSiena setup for reproducible test
mynet <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata <- sienaDataCreate(mynet)
mymodel <- getEffects(mydata)
mymodel <- includeEffects(mymodel, transTrip, name = "mynet")
mycontrols <- sienaAlgorithmCreate(projname = NULL, n3 = 60, cond = FALSE)
ans <- siena07(
  mycontrols,
  data = mydata,
  effects = mymodel,
  returnDeps = TRUE,
  silent = TRUE)

test_that("predict.sienaFit (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans,
        newdata = mydata,
        type = "tieProb",
        nsim = 10,
        condition = "transTrip_eval",
        level = "egoChoice"
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predict.sienaFit (base R fallback, optional MCSE + no CI)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans,
        newdata = mydata,
        type = "tieProb",
        nsim = 20,
        condition = "transTrip_eval",
        level = "period",
        uncertainty = TRUE,
        uncertaintyMcse = TRUE,
        uncertaintymcseBatches = 4,
        uncertaintySd = TRUE,
        uncertaintyCi = FALSE,
        verbose = FALSE
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
      expect_true("mcse_Mean" %in% names(pred_df))
      expect_true("mcse_SE" %in% names(pred_df))
      expect_false("q_025" %in% names(pred_df))
      expect_false("q_975" %in% names(pred_df))
      expect_false("Median" %in% names(pred_df))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predict.sienaFit (base R fallback, streaming uncertainty)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans,
        newdata = mydata,
        type = "tieProb",
        nsim = 10,
        condition = "transTrip",
        level = "period",
        uncertainty = TRUE,
        uncertaintyMode = "stream",
        verbose = FALSE
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
      expect_true("Mean" %in% names(pred_df))
      expect_true("SE" %in% names(pred_df))
      expect_true("q_025" %in% names(pred_df))
      expect_true("q_975" %in% names(pred_df))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predict.sienaFit (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- predict.sienaFit(
    object = ans,
    newdata = mydata,
    type = "tieProb",
    nsim = 10,
    condition = "transTrip_eval",
    level = "period",
    uncertainty = TRUE
  )
  expect_true("data.table" %in% class(pred_dt) && is.data.frame(pred_dt))
})

# does not work here for unknown reasons
# test_that("Test predict.sienaFit with PSOCK clustertype (data.table)", {
#   skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
#   library(data.table)
#   pred_dt <- predict(
#     object = ans,
#     newdata = mydata,
#     type = "tieProb",
#     nsim = 1000,
#     condition = "transTrip.eval",
#     level = "period",
#     uncertainty = TRUE,
#     useCluster = TRUE,
#     clusterType = "PSOCK",
#     nbrNodes = 2
#   )
#   expect_null(parallel:::getDefaultCluster())
#   expect_true("data.table" %in% class(pred_dt) && is.data.frame(pred_dt))
# })

# does not work here because mock bindings are not exported to cluster workers
# test_that("predictDynamic with PSOCK clustertype (base R fallback)", {
#   with_mocked_bindings(
#     {
#       pred_df <- predict.sienaFit(
#         object = ans,
#         newdata = mydata,
#         type = "tieProb",
#         nsim = 10,
#         condition = "transTrip.eval",
#         level = "period",
#         uncertainty = TRUE,
#         useCluster = TRUE,
#         clusterType = "PSOCK",
#         nbrNodes = 2
#       )
#       expect_null(parallel:::getDefaultCluster())
#       expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
#     },
#     requireNamespace = function(pkg, ...) FALSE,
#     .package="base"
#   )
# })

test_that("Test predict.sienaFit with FORK clustertype (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- predict(
    object = ans,
    newdata = mydata,
    type = "tieProb",
    nsim = 10,
    condition = "transTrip_eval",
    level = "period",
    uncertainty = TRUE,
    useCluster = TRUE,
    clusterType = "FORK",
    nbrNodes = 2
  )
  expect_null(parallel:::getDefaultCluster())
  expect_true("data.table" %in% class(pred_dt) && is.data.frame(pred_dt))
})

test_that("predictDynamic with FORK clustertype (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans,
        newdata = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        dynamic = TRUE,
        type = "tieProb",
        n3 = 60,
        nsim = 2,
        condition = "transTrip_eval",
        level = "period",
        uncertainty = TRUE,
        useCluster = TRUE,
        clusterType = "FORK",
        nbrNodes = 2,
        silent = FALSE
      )
      expect_null(parallel:::getDefaultCluster())
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})


test_that("predictDynamic (base R fallback)", {
  with_mocked_bindings(
    {
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
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predictDynamic (base R fallback, optional MCSE + no SD)", {
  with_mocked_bindings(
    {
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
        uncertaintyMcse = TRUE,
        uncertaintymcseBatches = 4,
        uncertaintySd = FALSE,
        uncertaintyCi = TRUE,
        verbose = FALSE
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
      expect_true("mcse_Mean" %in% names(pred_df))
      expect_false("SE" %in% names(pred_df))
      expect_false("mcse_SE" %in% names(pred_df))
      expect_true("q_025" %in% names(pred_df))
      expect_true("q_975" %in% names(pred_df))
      expect_true("Median" %in% names(pred_df))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predictDynamic (base R fallback, streaming uncertainty)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans,
        dynamic = TRUE,
        newdata = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        type = "tieProb",
        n3 = 60,
        nsim = 6,
        condition = "density_eval",
        uncertainty = TRUE,
        uncertaintyMode = "stream",
        verbose = FALSE
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
      expect_true("Mean" %in% names(pred_df))
      expect_true("SE" %in% names(pred_df))
      expect_true("q_025" %in% names(pred_df))
      expect_true("q_975" %in% names(pred_df))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predictDynamic (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
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
  expect_true("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))
})

# Minimal RSiena setup with manual interactions & not using one of the main effects
  mynet2 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
  mydata2 <- sienaDataCreate(mynet2)
  mymodel2 <- getEffects(mydata2)
  # Intentionally do NOT include outPop effect
  # mymodel2 <- includeEffects(mymodel2, outPop, name = "mynet2")
  mymodel2 <- includeInteraction(mymodel2, recip, outPop, name = "mynet2")
  mycontrols2 <- sienaAlgorithmCreate(projname = NULL, seed = 42, n3 = 60)
  ans2 <- siena07(
    mycontrols2,
    data = mydata2,
    effects = mymodel2,
    returnChangeContributions = TRUE,
    returnDataFrame = TRUE
  )
  
test_that("predict.sienaFit with custom interactions & without main effect (data.table),
  cond=FALSE", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  pred_dt <- predict.sienaFit(ans2,
    newdata =  mydata2,
    effects = mymodel2,
    level = "period",
    condition = "recip_eval",
    uncertainty = FALSE
  )
  expect_true("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))

})

test_that("predict.sienaFit with custom interactions & without main effect (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans2,
        newdata =  mydata2,
        effects = mymodel2,
        level = "period",
        condition = "recip_eval",
        uncertainty = FALSE
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("predict.sienaFit with custom interactions & uncertainty & without main effect (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  pred_dt <- predict.sienaFit(ans2,
    newdata =  mydata2,
    effects = mymodel2,
    level = "period",
    condition = "recip_eval",
    uncertainty = TRUE,
    nsim = 5
  )
  expect_true("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))

})

test_that("predict.sienaFit with custom interactions & uncertainty & without main effect (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- predict(
        object = ans2,
        newdata = mydata2,
        effects = mymodel2,
        uncertainty = TRUE,
        nsim = 5,
        condition = "recip_eval"
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

# cond=FALSE + interaction without main effect: this is the critical case where
# alignThetaNoRate must use requestedEffects (not effects) to avoid picking up
# injected base-effect slots in theta, which would give wrong (uniform) probs.
mycontrols2_uncond <- sienaAlgorithmCreate(projname = NULL, seed = 42, n3 = 60, cond = FALSE)
ans2_uncond <- siena07(mycontrols2_uncond, data = mydata2, effects = mymodel2,
  silent = TRUE)

test_that("predict.sienaFit interaction without main effect, cond=FALSE (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- predict.sienaFit(ans2_uncond, newdata = mydata2,
    effects = mymodel2, uncertainty = FALSE, level = "egoChoice")
  expect_true("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))
  expect_true(nrow(pred_dt) > 0)
  expect_false(any(is.nan(pred_dt$changeProb)))
})

# Unit test for alignThetaNoRate: verify it uses requestedEffects (not effects)
# so that injected base effects for interactions don't cause a length mismatch.
test_that("alignThetaNoRate uses requestedEffects$shortName, not effects$shortName", {
  # Reproduce the cond=FALSE + interaction-without-main-effect scenario:
  #   ans$effects has 6 rows (2 rate, 4 non-rate including injected outPop).
  #   ans$requestedEffects has 5 rows (2 rate, 3 non-rate — no outPop).
  #   ans$theta has 3 elements (one per non-rate in requestedEffects, no rates).
  #   Using effects:          non-rate count = 4 ≠ 3 = length(theta) → names not set (fallback)
  #   Using requestedEffects: non-rate count = 3 = 3 = length(theta)  → names set correctly
  theta_uncond <- c(-2.38, 3.06, -0.07)  # density_eval, recip_eval, unspInt_eval
  effectNames  <- c("density_eval", "recip_eval", "unspInt_eval")
  mock_ans <- list(
    effects = data.frame(
      shortName = c("rateX", "rateX", "density", "recip", "outPop", "unspInt"),
      type      = c("rate",  "rate",  "eval",    "eval",  "eval",   "eval"),
      include   = rep(TRUE, 6L),
      stringsAsFactors = FALSE
    ),
    requestedEffects = data.frame(
      shortName = c("rateX", "rateX", "density", "recip", "unspInt"),
      type      = c("rate",  "rate",  "eval",    "eval",  "eval"),
      include   = rep(TRUE, 5L),
      stringsAsFactors = FALSE
    )
  )
  result <- alignThetaNoRate(theta_uncond, effectNames, mock_ans)
  # Must return the three eval params by composite name
  expect_named(result, effectNames)
  expect_equal(unname(result[["density_eval"]]),  -2.38)
  expect_equal(unname(result[["recip_eval"]]),     3.06)
  expect_equal(unname(result[["unspInt_eval"]]),  -0.07)
})

# ---------------------------------------------------------------------------
# Two-network co-evolution: verify predict.sienaFit gives separate predictions #
# per network, and that alignThetaNoRate correctly matches parameters to 
# effects by name (not position) when cond=FALSE + interaction without main effect
# ---------------------------------------------------------------------------
mynet_a <- sienaDependent(array(c(s501, s502), dim=c(50, 50, 2)))
mynet_b <- sienaDependent(array(c(s502, s503), dim=c(50, 50, 2)))
mydata_co <- sienaDataCreate(mynet_a, mynet_b)
mymodel_co <- getEffects(mydata_co)
mymodel_co <- includeEffects(mymodel_co, transTrip, name = "mynet_a")
mymodel_co <- includeEffects(mymodel_co, transTrip, name = "mynet_b")
mycontrols_co <- sienaAlgorithmCreate(projname = NULL, seed = 1, n3 = 60,
                                      cond = FALSE)
ans_co <- siena07(
  mycontrols_co,
  data    = mydata_co,
  effects = mymodel_co,
  returnChangeContributions = TRUE,
  returnDataFrame = FALSE,
  silent  = TRUE
)

# ---------------------------------------------------------------------------
# Static predict: co-evolution
# ---------------------------------------------------------------------------
test_that("predict.sienaFit static: co-evolution networks get separate predictions", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE))
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
  skip_if_not(requireNamespace("data.table", quietly = TRUE))
  pred_a <- predict.sienaFit(ans_co, newdata = mydata_co, depvar = "mynet_a",
                              dynamic = TRUE, algorithm = mycontrols_co,
                              effects = mymodel_co, n3 = 60, nsim = 5,
                              uncertainty = FALSE)
  pred_b <- predict.sienaFit(ans_co, newdata = mydata_co, depvar = "mynet_b",
                              dynamic = TRUE, algorithm = mycontrols_co,
                              effects = mymodel_co, n3 = 60, nsim = 5,
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
# ---------------------------------------------------------------------------
mynet_cm <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata_cm <- sienaDataCreate(mynet_cm)
mymodel_cm <- getEffects(mydata_cm)
# recip (eval) is already included by default; add the creation variant
mymodel_cm <- includeEffects(mymodel_cm, recip,     type = "creation")
# transTrip in both evaluation and endowment flavours
mymodel_cm <- includeEffects(mymodel_cm, transTrip, type = "eval")
mymodel_cm <- includeEffects(mymodel_cm, transTrip, type = "endow")
mycontrols_cm <- sienaAlgorithmCreate(projname = NULL, seed = 3, n3 = 60, cond = FALSE)
ans_cm <- siena07(
  mycontrols_cm,
  data    = mydata_cm,
  effects = mymodel_cm,
  returnChangeContributions = TRUE,
  returnDataFrame = FALSE,
  silent  = TRUE
)

# ---------------------------------------------------------------------------
# Static predict: creation/maintenance — all effect types present, distinct
# ---------------------------------------------------------------------------
test_that("predict.sienaFit static: creation/endow effects correctly separated", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE))
  # Test via contributions struct (static)
  wide_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm,
    depvar = "mynet_cm", returnWide = TRUE
  )
  # All five effect slots appear as distinct columns
  expect_true("density_eval"    %in% wide_cm$effectNames)
  expect_true("recip_eval"      %in% wide_cm$effectNames)
  expect_true("recip_creation"  %in% wide_cm$effectNames)
  expect_true("transTrip_eval"  %in% wide_cm$effectNames)
  expect_true("transTrip_endow" %in% wide_cm$effectNames)
  # eval and creation columns must differ (different zero-masks)
  recip_eval_col  <- which(wide_cm$effectNames == "recip_eval")
  recip_creat_col <- which(wide_cm$effectNames == "recip_creation")
  expect_false(identical(wide_cm$contribMat[, recip_eval_col],
                         wide_cm$contribMat[, recip_creat_col]))
  pred_cm <- predictProbability(wide_cm, theta = c(ans_cm$rate, ans_cm$theta), type = "tieProb")
  
})

# ---------------------------------------------------------------------------
# Dynamic predict: creation/maintenance
# ---------------------------------------------------------------------------
test_that("predict.sienaFit dynamic: creation/endow effect columns present", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE))
  # Test via the contributions struct (which carries the column names directly)
  # rather than predict.sienaFit whose aggregated output does not include them.
  theta_cm <- c(ans_cm$rate, ans_cm$theta)
  contrib <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = theta_cm, effects = mymodel_cm, depvar = "mynet_cm",
    returnWide = TRUE, useChangeContributions = TRUE
  )
  expect_true("recip_eval"      %in% contrib$effectNames)
  expect_true("recip_creation"  %in% contrib$effectNames)
  expect_true("transTrip_eval"  %in% contrib$effectNames)
  expect_true("transTrip_endow" %in% contrib$effectNames)
  # eval and creation columns must differ (creation is 0 for existing ties)
  recip_eval_col    <- which(contrib$effectNames == "recip_eval")
  recip_creat_col   <- which(contrib$effectNames == "recip_creation")
  expect_false(identical(contrib$contribMat[, recip_eval_col],
                         contrib$contribMat[, recip_creat_col]))
})