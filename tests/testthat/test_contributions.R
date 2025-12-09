# testthat::skip_on_cran()

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

# --- Test static contribution extraction (pure list)
test_that("Static contributions (pure list)", {
    stat_ls <- getStaticChangeContributions(ans = ans, data = mydata, 
      effects = mymodel)
    expect_true(is.list(stat_ls))
})


# --- Test static contribution extraction (base R fallback)
test_that("Static contributions (base R fallback)", {
  with_mocked_bindings(
    {
      stat_df <- getStaticChangeContributions(ans = ans, data = mydata, 
        effects = mymodel, returnDataFrame = TRUE)
      expect_true(is.data.frame(stat_df))
      expect_false("data.table" %in% class(stat_df))
    },
    requireNamespace = function(package, quietly=TRUE) FALSE,
    .package="base"
  )
})

# --- Test static contribution extraction (data.table)
test_that("Static contributions (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
        stat_dt <- getStaticChangeContributions(ans = ans, 
        data = mydata,
          returnDataFrame = TRUE)
  expect_true("data.table" %in% class(stat_dt))
})

# # --- Test dynamic contribution extraction (pure list)
test_that("Dynamic contributions (pure list)", {
      dyn_ls <- getDynamicChangeContributions(
        ans = ans,
        data = mydata,
        algorithm = mycontrols,
        theta = c(ans$rate, ans$theta),
        effects = mymodel,
        depvar = "mynet",
        returnDataFrame = FALSE
      )
      expect_true(is.list(dyn_ls))
})

# --- Test dynamic contribution extraction (base R fallback))
test_that("Dynamic contributions (base R fallback)", {
    with_mocked_bindings(
    {
      dyn_df <- getDynamicChangeContributions(
        ans = ans,
        data = mydata,
        theta = c(ans$rate, ans$theta),
        algorithm = mycontrols,
        effects = mymodel,
        depvar = "mynet",
        returnDataFrame = TRUE
      )
      expect_true(is.data.frame(dyn_df) && !("data.table" %in% class(dyn_df)))
    },
    requireNamespace = function(package, quietly=TRUE) FALSE,
    .package="base"
  )
})

# --- Test dynamic contribution extraction (data.table)
test_that("Dynamic contributions (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn_dt <- getDynamicChangeContributions(
    ans = ans,
    data = mydata,
    theta = c(ans$rate, ans$theta),
    algorithm = mycontrols,
    effects = mymodel,
    depvar = "mynet",
    returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_dt))
})

test_that("Dynamic contributions (base R fallback, useChangeContributions = FALSE)", {
  with_mocked_bindings(
    {
      dyn_df <- getDynamicChangeContributions(
        ans = ans,
        data = mydata,
        theta = c(ans$rate, ans$theta),
        algorithm = mycontrols,
        effects = mymodel,
        depvar = "mynet",
        useChangeContributions = FALSE,
        n3 = 60,
        returnDataFrame = TRUE
      )
      expect_true(is.matrix(dyn_df) || is.data.frame(dyn_df))
      expect_false("data.table" %in% class(dyn_df))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("Dynamic contributions (data.table, useChangeContributions = FALSE)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn_dt <- getDynamicChangeContributions(
    ans = ans,
    data = mydata,
    theta = c(ans$rate, ans$theta),
    algorithm = mycontrols,
    effects = mymodel,
    depvar = "mynet",
    useChangeContributions = FALSE,
    n3 = 60,
    returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_dt))
})

test_that("Dynamic contributions (base R fallback, new theta values)", {
  with_mocked_bindings(
    {
      new_theta <- MASS::mvrnorm(n = 1, mu = ans$theta, Sigma = ans$covtheta)
      dyn_df <- getDynamicChangeContributions(
        ans = NULL,
        data = mydata,
        theta = c(ans$rate, new_theta),
        algorithm = mycontrols,
        effects = mymodel,
        depvar = "mynet",
        useChangeContributions = FALSE,
        n3 = 60,
        returnDataFrame = TRUE
      )
      expect_true(is.matrix(dyn_df) || is.data.frame(dyn_df))
      expect_false("data.table" %in% class(dyn_df))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("Dynamic contributions (data.table, new theta values)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  new_theta <- MASS::mvrnorm(n = 1, mu = ans$theta, Sigma = ans$covtheta)
  dyn_dt <- getDynamicChangeContributions(
    ans = NULL,
    data = mydata,
    theta = c(ans$rate, new_theta),
    algorithm = mycontrols,
    effects = mymodel,
    depvar = "mynet",
    useChangeContributions = FALSE,
    n3 = 60,
    returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_dt))
})
