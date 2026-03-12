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
        condition = "transTrip",
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
        condition = "transTrip",
        level = "period",
        uncertainty = TRUE,
        uncertainty_mcse = TRUE,
        uncertainty_mcse_batches = 4,
        uncertainty_sd = TRUE,
        uncertainty_ci = FALSE,
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
        uncertainty_mode = "stream",
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
    condition = "transTrip",
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
#     condition = "transTrip",
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
#         condition = "transTrip",
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
    condition = "transTrip",
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
        condition = "transTrip",
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
        condition = "density"
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
        condition = "density",
        uncertainty = TRUE,
        uncertainty_mcse = TRUE,
        uncertainty_mcse_batches = 4,
        uncertainty_sd = FALSE,
        uncertainty_ci = TRUE,
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
        condition = "density",
        uncertainty = TRUE,
        uncertainty_mode = "stream",
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
    condition = "density"
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
    condition = "recip",
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
        condition = "recip",
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
    condition = "recip",
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
        nsim = 5
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
  #   ans$effects has an extra "outPop" row injected for the interaction,
  #   making effects length 6 while theta length is only 5.
  #   The fix: use requestedEffects (length 5) for name assignment.
  theta_uncond <- c(2.5, 3.0, -2.38, 3.06, -0.07)  # rate1, rate2, density, recip, unspInt
  effectNames <- c("density", "recip", "unspInt")
  mock_ans <- list(
    effects           = list(shortName = c("rateX", "rateX", "density", "recip", "outPop", "unspInt")),  # 6 != 5
    requestedEffects  = list(shortName = c("rateX", "rateX", "density", "recip", "unspInt"))             # 5 == 5
  )
  result <- alignThetaNoRate(theta_uncond, effectNames, mock_ans)
  # Must return the three eval params by name — not rate params by position
  expect_named(result, effectNames)
  expect_equal(unname(result[["density"]]),  -2.38)
  expect_equal(unname(result[["recip"]]),     3.06)
  expect_equal(unname(result[["unspInt"]]),  -0.07)
})
