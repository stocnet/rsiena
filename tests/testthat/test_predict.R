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
      pred <- sienaPredict(
        ans = ans,
        data = mydata,
        useTieProb = TRUE,
        condition = "transTrip",
        level = "period"
      )
      expect_true(is.data.frame(pred) && !("data.table" %in% class(pred)))
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
    nsim = 10,
    condition = "transTrip",
    level = "period",
    uncertainty = TRUE,
    useCluster = TRUE,
    clusterType = "PSOCK",
    nbrNodes = 2
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
      pred_df <- sienaPredict(
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
      expect_true(is.data.frame(pred_df) && !("data.table" %in% class(pred_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})


test_that("sienaPredictDynamic (base R fallback)", {
  with_mocked_bindings(
    {
      pred <- sienaPredictDynamic(
        ans = ans,
        data = mydata,
        effects = mymodel,
        algorithm = mycontrols,
        useTieProb = TRUE,
        n3 = 60,
        nsim = 10,
        condition = "density"
      )
      expect_true(is.data.frame(pred) && !("data.table" %in% class(pred)))
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

