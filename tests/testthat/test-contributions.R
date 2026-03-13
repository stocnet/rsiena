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


# Minimal RSiena setup with manual interactions & not using one of the main effects
mynet2 <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
mydata2 <- sienaDataCreate(mynet2)
mymodel2 <- getEffects(mydata2)
## outdegree recip model
# Intentionally do NOT include outPop effect
# mymodel2 <- includeEffects(mymodel2, outPop, name = "mynet2")
mymodel2 <- includeInteraction(mymodel2, recip, outPop, name = "mynet2")

mycontrols2 <- sienaAlgorithmCreate(projname = NULL, 
  seed = 42, 
  n3 = 60)

ans2 <- siena07(
  mycontrols2,
  data = mydata2,
  effects = mymodel2,
  returnChangeContributions = TRUE,
  returnDataFrame = FALSE
)

test_that("getStaticChangeContributions with custom interactions & without main effect (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  stat_dt <- getStaticChangeContributions(
      ans = ans2,
      data = mydata2,
      effects = mymodel2,
      depvar = "mynet2",
      returnDataFrame = TRUE
  )
  expect_true(!any(stat_dt$effectname == "outPop"))
  expect_true(any(stat_dt$effectname == "unspInt"))
  expect_true("data.table" %in% class(stat_dt) || is.data.frame(stat_dt))

})

test_that("getStaticChangeContributions with custom interactions & without main effect (base R fallback)", {
  with_mocked_bindings(
    {
      stat_df <- getStaticChangeContributions(
          ans = ans2,
          data = mydata2,
          effects = mymodel2,
          depvar = "mynet2",
          returnDataFrame = TRUE
      )
      expect_true(!any(stat_df$effectname == "outPop"))
      expect_true(any(stat_df$effectname == "unspInt"))
      expect_true(is.data.frame(stat_df) && !("data.table" %in% class(stat_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

test_that("getDynamicChangeContributions with custom interactions & without main effect (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  stat_dt <- getDynamicChangeContributions(
      ans = ans2,
      data = mydata2,
      algorithm = mycontrols2,
      effects = mymodel2,
      depvar = "mynet2",
      useChangeContributions = FALSE,
      n3 = 60,
      returnDataFrame = TRUE
  )
  expect_true(!any(stat_dt$effectname == "outPop"))
  expect_true(any(stat_dt$effectname == "unspInt"))
  expect_true("data.table" %in% class(stat_dt) || is.data.frame(stat_dt))

})

test_that("getDynamicChangeContributions with custom interactions & without main effect (base R fallback)", {
  with_mocked_bindings(
    {
      stat_df <- getDynamicChangeContributions(
          ans = ans2,
          data = mydata2,
          algorithm = mycontrols2,
          effects = mymodel2,
          depvar = "mynet2",
          useChangeContributions = FALSE,
          n3 = 60,
          returnDataFrame = TRUE
      )
      expect_true(!any(stat_df$effectname == "outPop"))
      expect_true(any(stat_df$effectname == "unspInt"))
      expect_true(is.data.frame(stat_df) && !("data.table" %in% class(stat_df)))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package="base"
  )
})

# ── thetaValues / simOnly batch approach ─────────────────────────────────────
# sienaAlgorithmCreate(simOnly=TRUE, thetaValues=...) runs a single siena07
# call with one initializeFRAN + terminateFRAN for all theta draws. Phase3
# assigns thetaValues[nit,] to zsmall$theta at each iteration and stores the
# chain in z$changeContributions[[nit]]. Verify:
#   (a) all N3 * NSIM chains are returned and none are NULL,
#   (b) thetaUsed records the correct theta per block, and
#   (c) chains from different theta blocks have different values.

test_that("thetaValues batch returns all chains non-NULL", {
  N3   <- 5L
  NSIM <- 2L
  set.seed(7)
  theta_draws <- MASS::mvrnorm(NSIM, mu = ans$theta, Sigma = ans$covtheta)

  tv <- do.call(rbind, lapply(seq_len(NSIM), function(i) {
    matrix(rep(theta_draws[i, ], each = N3), nrow = N3)
  }))

  alg_batch <- sienaAlgorithmCreate(projname = NULL, nsub = 0,
                                     n3 = nrow(tv), cond = FALSE,
                                     simOnly = TRUE, seed = 11)
  eff_batch <- mymodel
  eff_batch$initialValue[eff_batch$include] <- theta_draws[1, ]

  ans_batch <- siena07(alg_batch, data = mydata, effects = eff_batch,
                        batch = TRUE, silent = TRUE, thetaValues = tv,
                        returnChangeContributions = TRUE)

  cc <- ans_batch$changeContributions
  expect_length(cc, N3 * NSIM)
  expect_true(all(!sapply(cc, is.null)))
})

test_that("thetaValues batch records correct theta per block", {
  N3   <- 5L
  NSIM <- 2L
  set.seed(7)
  theta_draws <- MASS::mvrnorm(NSIM, mu = ans$theta, Sigma = ans$covtheta)

  tv <- do.call(rbind, lapply(seq_len(NSIM), function(i) {
    matrix(rep(theta_draws[i, ], each = N3), nrow = N3)
  }))

  alg_batch <- sienaAlgorithmCreate(projname = NULL, nsub = 0,
                                     n3 = nrow(tv), cond = FALSE,
                                     simOnly = TRUE, seed = 11)
  eff_batch <- mymodel
  eff_batch$initialValue[eff_batch$include] <- theta_draws[1, ]

  ans_batch <- siena07(alg_batch, data = mydata, effects = eff_batch,
                        batch = TRUE, silent = TRUE, thetaValues = tv,
                        returnChangeContributions = TRUE)

  for (i in seq_len(NSIM)) {
    rows <- ((i - 1) * N3 + 1):(i * N3)
    block <- ans_batch$thetaUsed[rows, , drop = FALSE]
    expected <- matrix(rep(theta_draws[i, ], each = N3), nrow = N3)
    expect_equal(block, expected, ignore_attr = TRUE,
                 label = paste("thetaUsed block", i))
  }
})

test_that("thetaValues batch: chains from different theta blocks differ", {
  N3   <- 5L
  NSIM <- 2L
  set.seed(7)
  # Use theta draws that are far apart so chains are almost certainly different
  theta_draws <- rbind(
    ans$theta + 2,
    ans$theta - 2
  )

  tv <- do.call(rbind, lapply(seq_len(NSIM), function(i) {
    matrix(rep(theta_draws[i, ], each = N3), nrow = N3)
  }))

  alg_batch <- sienaAlgorithmCreate(projname = NULL, nsub = 0,
                                     n3 = nrow(tv), cond = FALSE,
                                     simOnly = TRUE, seed = 11)
  eff_batch <- mymodel
  eff_batch$initialValue[eff_batch$include] <- theta_draws[1, ]

  ans_batch <- siena07(alg_batch, data = mydata, effects = eff_batch,
                        batch = TRUE, silent = TRUE, thetaValues = tv,
                        returnChangeContributions = TRUE)

  cc <- ans_batch$changeContributions
  # Extract first ministep matrix from block 1 and block 2
  get_mat <- function(chain) chain[[1]][[1]][[1]]
  mat1 <- get_mat(cc[[1]])
  mat2 <- get_mat(cc[[N3 + 1]])

  expect_identical(dim(mat1), dim(mat2))          # same structure
  expect_false(identical(mat1, mat2))             # different values
})
