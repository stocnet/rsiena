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
  returnChangeContributions = TRUE,
  returnDataFrame = TRUE,
  silent = TRUE
)

test_that("sienaPredict (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- sienaPredict(
        ans = ans,
        data = mydata,
        useTieProb = TRUE,
        condition = "transTrip",
        level = "period"
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("sienaPredict (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- sienaPredict(
    ans = ans,
    data = mydata,
    useTieProb = TRUE,
    nsim = 10,
    condition = "transTrip",
    level = "period",
    uncertainty = TRUE
  )
  expect_true("data.table" %in% class(pred_dt) && is.data.frame(pred_dt))
})

test_that("Test sienaPredict with PSOCK clustertype (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- sienaPredict(
    ans = ans,
    data = mydata,
    useTieProb = TRUE,
    nsim = 1000,
    condition = "transTrip",
    level = "period",
    uncertainty = TRUE,
    useCluster = TRUE,
    clusterType = "PSOCK",
    nbrNodes = 6
  )
  expect_null(parallel:::getDefaultCluster())
  expect_true("data.table" %in% class(pred_dt) && is.data.frame(pred_dt))
})

# does not work here because mock bindings are not exported to cluster workers
# test_that("sienaPredictDynamic with PSOCK clustertype (base R fallback)", {
#   with_mocked_bindings(
#     {
#       pred_df <- sienaPredict(
#         ans = ans,
#         data = mydata,
#         useTieProb = TRUE,
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

test_that("Test sienaPredict with FORK clustertype (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- sienaPredict(
    ans = ans,
    data = mydata,
    useTieProb = TRUE,
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

test_that("sienaPredictDynamic with FORK clustertype (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- sienaPredictDynamic(
        ans = ans,
        data = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        useTieProb = TRUE,
        n3 = 60,
        nsim = 24,
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


test_that("sienaPredictDynamic (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- sienaPredictDynamic(
        ans = ans,
        data = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        useTieProb = TRUE,
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

test_that("sienaPredictDynamic (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  pred_dt <- sienaPredictDynamic(
    ans = ans,
    data = mydata,
    effects = mymodel,
    algorithm = mycontrols,
    useTieProb = TRUE,
    n3 = 60,
    nsim = 10,
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
  mycontrols2 <- sienaAlgorithmCreate(projname = "marginal_example", seed = 42, n3 = 1000)
  ans2 <- siena07(
    mycontrols2,
    data = mydata2,
    effects = mymodel2,
    returnChangeContributions = TRUE,
    returnDataFrame = TRUE
  )
  
test_that("sienaPredict with custom interactions & without main effect (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  pred_dt <- sienaPredict(ans2,
    data =  mydata2,
    effects = mymodel2,
    uncertainty = FALSE
  )
  expect_true("data.table" %in% class(pred_dt) || is.data.frame(pred_dt))

})

test_that("sienaPredict with custom interactions & without main effect (base R fallback)", {
  with_mocked_bindings(
    {
      pred_df <- sienaPredict(
        ans = ans2,
        data = mydata2,
        effects = mymodel2,
        uncertainty = FALSE
      )
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

## add dynamic checks


