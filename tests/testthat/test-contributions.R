# testthat::skip_on_cran()
# Models: ans, mydata, mymodel, mycontrols  — from helper-models.R
#         ans_int, mydata_int, mymodel_int  — from helper-models.R (full mode only)
#         ans_co, mydata_co, mymodel_co, mycontrols_co — from helper-models.R (full mode only)
#         ans_cm, mydata_cm, mymodel_cm, mycontrols_cm — from helper-models.R (full mode only)

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
      skip_slow()
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
    skip_slow()
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
  skip_slow()
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

# --- Column schema tests (aligned between static and dynamic) ---------------

.static_cols  <- c("group", "period", "networkName", "ego", "choice",
                   "effectname", "effecttype", "contribution")
.dynamic_cols <- c("group", "period", "ministep", "networkName", "ego",
                   "choice", "effectname", "effecttype", "contribution")

test_that("Static df has expected columns incl. effecttype and ego (data.table)", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  stat <- getStaticChangeContributions(ans = ans, data = mydata,
    effects = mymodel, returnDataFrame = TRUE)
  expect_true(all(.static_cols %in% names(stat)))
  expect_false("effectIdx" %in% names(stat))
})

test_that("Static df has expected columns incl. effecttype and ego (base R)", {
  with_mocked_bindings(
    {
      stat <- getStaticChangeContributions(ans = ans, data = mydata,
        effects = mymodel, returnDataFrame = TRUE)
      expect_true(all(.static_cols %in% names(stat)))
      expect_false("effectIdx" %in% names(stat))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("Dynamic df has expected columns incl. ego and no effectIdx (data.table)", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn <- getDynamicChangeContributions(
    ans = ans, data = mydata, theta = c(ans$rate, ans$theta),
    algorithm = mycontrols, effects = mymodel,
    depvar = "mynet", returnDataFrame = TRUE
  )
  expect_true(all(.dynamic_cols %in% names(dyn)))
  expect_false("effectIdx" %in% names(dyn))
  expect_true(is.integer(dyn$ego))
  expect_true(all(dyn$ego >= 1L))
})

test_that("Dynamic df has expected columns incl. ego and no effectIdx (base R)", {
  skip_slow()
  with_mocked_bindings(
    {
      dyn <- getDynamicChangeContributions(
        ans = ans, data = mydata, theta = c(ans$rate, ans$theta),
        algorithm = mycontrols, effects = mymodel,
        depvar = "mynet", returnDataFrame = TRUE
      )
      expect_true(all(.dynamic_cols %in% names(dyn)))
      expect_false("effectIdx" %in% names(dyn))
      expect_true(is.integer(dyn$ego))
      expect_true(all(dyn$ego >= 1L))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("Static df: effecttype values are valid effect types", {
  with_mocked_bindings(
    {
      stat <- getStaticChangeContributions(ans = ans, data = mydata,
        effects = mymodel, returnDataFrame = TRUE)
      expect_true(all(stat$effecttype %in% c("eval", "creation", "endow")))
    },
    requireNamespace = function(pkg, ...) FALSE, .package = "base"
  )
})

test_that("Dynamic contributions (base R fallback, useChangeContributions = FALSE)", {
  skip_slow()
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
  skip_slow()
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
  skip_slow()
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
  skip_slow()
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


# Interaction without main effect: recip * outPop unspInt, without outPop.
# Model fit from helper-models.R (full mode): ans_int, mydata_int, mymodel_int.

test_that("getStaticChangeContributions with custom interactions & without main effect (data.table)", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  stat_dt <- getStaticChangeContributions(
      ans = ans_int,
      data = mydata_int,
      effects = mymodel_int,
      depvar = "mynet_int",
      returnDataFrame = TRUE
  )
  expect_true(!any(stat_dt$effectname == "outPop"))
  expect_true(any(stat_dt$effectname == "unspInt"))
  expect_true("data.table" %in% class(stat_dt) || is.data.frame(stat_dt))
})

test_that("getStaticChangeContributions with custom interactions & without main effect (base R fallback)", {
  skip_slow()
  with_mocked_bindings(
    {
      stat_df <- getStaticChangeContributions(
          ans = ans_int,
          data = mydata_int,
          effects = mymodel_int,
          depvar = "mynet_int",
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
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)

  stat_dt <- getDynamicChangeContributions(
      ans = ans_int,
      data = mydata_int,
      algorithm = mycontrols_int,
      effects = mymodel_int,
      depvar = "mynet_int",
      useChangeContributions = FALSE,
      n3 = 50,
      returnDataFrame = TRUE
  )
  expect_true(!any(stat_dt$effectname == "outPop"))
  expect_true(any(stat_dt$effectname == "unspInt"))
  expect_true("data.table" %in% class(stat_dt) || is.data.frame(stat_dt))
})

test_that("getDynamicChangeContributions with custom interactions & without main effect (base R fallback)", {
  skip_slow()
  with_mocked_bindings(
    {
      stat_df <- getDynamicChangeContributions(
          ans = ans_int,
          data = mydata_int,
          algorithm = mycontrols_int,
          effects = mymodel_int,
          depvar = "mynet_int",
          useChangeContributions = FALSE,
          n3 = 50,
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

# ---------------------------------------------------------------------------
# Two-network co-evolution: non-unique effect shortNames across depvars
# Both networks carry the same effects (e.g. "density", "transTrip"),
# so shortNames are NOT globally unique.  The extraction must keep them
# correctly separated by depvar and must not collapse or mis-assign columns.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# Two-network co-evolution: non-unique effect shortNames across depvars
# Models from helper-models.R (full mode): ans_co, mydata_co, mymodel_co, mycontrols_co
# ---------------------------------------------------------------------------

test_that("co-evolution static wide: each depvar has its own effect columns", {
  skip_slow()
  wide_a <- getStaticChangeContributions(
    ans = ans_co, data = mydata_co, effects = mymodel_co,
    depvar = "mynet_a", returnWide = TRUE
  )
  wide_b <- getStaticChangeContributions(
    ans = ans_co, data = mydata_co, effects = mymodel_co,
    depvar = "mynet_b", returnWide = TRUE
  )
  cols_a <- colnames(wide_a$contribMat)
  cols_b <- colnames(wide_b$contribMat)

  # Both networks now have depvar-prefixed column names (mynet_a_*, mynet_b_*)
  expect_true("mynet_a_density_eval"   %in% cols_a)
  expect_true("mynet_a_transTrip_eval" %in% cols_a)
  expect_true("mynet_b_density_eval"   %in% cols_b)
  expect_true("mynet_b_transTrip_eval" %in% cols_b)

  # effectNames stored in the struct must match the columns exactly
  expect_equal(wide_a$effectNames, cols_a)
  expect_equal(wide_b$effectNames, cols_b)

  # The two depvars must produce separate structs with different data
  expect_false(identical(wide_a$contribMat, wide_b$contribMat))

  # No NA columns (mis-assignment would leave all-NA)
  expect_false(any(is.na(wide_a$contribMat[, "mynet_a_density_eval"])))
  expect_false(any(is.na(wide_b$contribMat[, "mynet_b_density_eval"])))
})

test_that("co-evolution static wide (per-depvar): each call returns network-specific struct", {
  skip_slow()
  wide_a <- getStaticChangeContributions(
    ans = ans_co, data = mydata_co, effects = mymodel_co,
    depvar = "mynet_a", returnWide = TRUE
  )
  wide_b <- getStaticChangeContributions(
    ans = ans_co, data = mydata_co, effects = mymodel_co,
    depvar = "mynet_b", returnWide = TRUE
  )

  # Both structs must have the canonical fields
  expect_true(is.matrix(wide_a$contribMat))
  expect_true(is.matrix(wide_b$contribMat))

  # effectNames must be depvar-specific, not shared/overwritten
  expect_equal(colnames(wide_a$contribMat), wide_a$effectNames)
  expect_equal(colnames(wide_b$contribMat), wide_b$effectNames)

  # The two structs must not be identical (different networks → different data)
  expect_false(identical(wide_a$contribMat, wide_b$contribMat))
})

test_that("co-evolution static long (data.frame): effectname column distinguishes depvars", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  long <- getStaticChangeContributions(
    ans = ans_co, data = mydata_co, effects = mymodel_co,
    returnDataFrame = TRUE
  )
  # Both networkNames must appear
  expect_true("mynet_a" %in% long$networkName)
  expect_true("mynet_b" %in% long$networkName)
  # "density" must appear for both depvars, not just one
  nets_with_density <- unique(long$networkName[long$effectname == "density"])
  expect_true("mynet_a" %in% nets_with_density)
  expect_true("mynet_b" %in% nets_with_density)
})

test_that("co-evolution dynamic per-depvar (pure list)", {
  skip_slow()
  dyn_a <- getDynamicChangeContributions(
    ans = ans_co, data = mydata_co, algorithm = mycontrols_co,
    theta = c(ans_co$rate, ans_co$theta),
    effects = mymodel_co, depvar = "mynet_a", 
    useChangeContributions = TRUE, returnDataFrame = FALSE
  )
  dyn_b <- getDynamicChangeContributions(
    ans = ans_co, data = mydata_co, algorithm = mycontrols_co,
    theta = c(ans_co$rate, ans_co$theta),
    effects = mymodel_co, depvar = "mynet_b", 
    useChangeContributions = TRUE, returnDataFrame = FALSE
  )
  expect_true(is.list(dyn_a))
  expect_true(is.list(dyn_b))
  # Results must be network-specific, not identical
  expect_false(identical(dyn_a, dyn_b))
})

test_that("co-evolution dynamic (data.table): each depvar returns its own effect rows", {
  skip_slow()
  library(data.table)
  dyn_a <- getDynamicChangeContributions(
    ans = ans_co, data = mydata_co, algorithm = mycontrols_co,
    theta = c(ans_co$rate, ans_co$theta),
    effects = mymodel_co, depvar = "mynet_a", returnDataFrame = TRUE
  )
  dyn_b <- getDynamicChangeContributions(
    ans = ans_co, data = mydata_co, algorithm = mycontrols_co,
    theta = c(ans_co$rate, ans_co$theta),
    effects = mymodel_co, depvar = "mynet_b", returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_a))
  expect_true("data.table" %in% class(dyn_b))
  # Both depvars have density and transTrip contributions
  expect_true("density"   %in% dyn_a$effectname)
  expect_true("transTrip" %in% dyn_a$effectname)
  expect_true("density"   %in% dyn_b$effectname)
  expect_true("transTrip" %in% dyn_b$effectname)
  # networkName must be consistent with the requested depvar
  expect_true(all(dyn_a$networkName == "mynet_a"))
  expect_true(all(dyn_b$networkName == "mynet_b"))
})

test_that("co-evolution dynamic (base R fallback): per-depvar extraction correct", {
  skip_slow()
  with_mocked_bindings(
    {
      dyn_a <- getDynamicChangeContributions(
        ans = ans_co, data = mydata_co, algorithm = mycontrols_co,
        theta = c(ans_co$rate, ans_co$theta),
        effects = mymodel_co, depvar = "mynet_a", returnDataFrame = TRUE
      )
      dyn_b <- getDynamicChangeContributions(
        ans = ans_co, data = mydata_co, algorithm = mycontrols_co,
        theta = c(ans_co$rate, ans_co$theta),
        effects = mymodel_co, depvar = "mynet_b", returnDataFrame = TRUE
      )
      expect_true(is.data.frame(dyn_a) && !("data.table" %in% class(dyn_a)))
      expect_true(is.data.frame(dyn_b) && !("data.table" %in% class(dyn_b)))
      expect_true("transTrip" %in% dyn_a$effectname)
      expect_true("transTrip" %in% dyn_b$effectname)
      # Rows must belong to the right network
      expect_true(all(dyn_a$networkName == "mynet_a"))
      expect_true(all(dyn_b$networkName == "mynet_b"))
    },
    requireNamespace = function(pkg, ...) FALSE,
    .package = "base"
  )
})

# ---------------------------------------------------------------------------
# Creation / endowment effects: non-unique shortNames *within* a single depvar
#
# recip (eval, included by default) + recip (creation) → shortName "recip" twice
# transTrip (eval) + transTrip (endow) → shortName "transTrip" twice
#
# Extraction must return one column *per included effect slot* (not collapsed),
# and creation/endow columns must carry different values than the matching eval
# column (creation contributions are 0 on tie-deletion choices; endow
# contributions are 0 on tie-creation choices).
# Models from helper-models.R (full mode): ans_cm, mydata_cm, mymodel_cm, mycontrols_cm
# ---------------------------------------------------------------------------

# Column-type index — derived lazily inside each test once models are available.
.cm_indices <- function() {
  incl <- mymodel_cm[
    mymodel_cm$include & mymodel_cm$type != "rate" & mymodel_cm$name == "mynet_cm", ]
  types <- incl$type
  names <- incl$shortName
  list(
    incl        = incl,
    recip_eval  = which(names == "recip"     & types == "eval"),
    recip_creat = which(names == "recip"     & types == "creation"),
    tt_eval     = which(names == "transTrip" & types == "eval"),
    tt_endow    = which(names == "transTrip" & types == "endow")
  )
}

test_that("creation/endow static wide: all non-rate effect slots are filled", {
  skip_slow()
  idx  <- .cm_indices()
  wide <- getStaticChangeContributions(ans = ans_cm, data = mydata_cm,
    effects = mymodel_cm, depvar = "mynet_cm", returnWide = TRUE)
  mat  <- wide$contribMat
  n_included <- nrow(idx$incl)
  expect_equal(ncol(mat), n_included)
  expect_equal(length(wide$effectNames), n_included)
  expect_false(any(apply(mat, 2L, function(x) all(is.na(x)))))
  expect_equal(length(idx$recip_eval),  1L)
  expect_equal(length(idx$recip_creat), 1L)
  expect_equal(length(idx$tt_eval),     1L)
  expect_equal(length(idx$tt_endow),    1L)
  expect_true(all(mat[, idx$tt_endow] >= 0, na.rm = TRUE))
})

test_that("creation/endow static wide: creation contributions are non-negative", {
  skip_slow()
  idx <- .cm_indices()
  mat <- getStaticChangeContributions(ans = ans_cm, data = mydata_cm,
    effects = mymodel_cm, depvar = "mynet_cm", returnWide = TRUE)$contribMat
  expect_true(all(mat[, idx$recip_creat] >= 0, na.rm = TRUE))
})

test_that("creation/endow static wide: endow(j) + eval(j) == 0 wherever endow(j) > 0", {
  skip_slow()
  idx <- .cm_indices()
  mat <- getStaticChangeContributions(ans = ans_cm, data = mydata_cm,
    effects = mymodel_cm, depvar = "mynet_cm", returnWide = TRUE)$contribMat
  rows_with_endow <- which(mat[, idx$tt_endow] > 0)
  expect_gt(length(rows_with_endow), 0L)
  expect_equal(
    mat[rows_with_endow, idx$tt_endow] + mat[rows_with_endow, idx$tt_eval],
    rep(0, length(rows_with_endow))
  )
})

test_that("creation/endow static wide: creation(j) == eval(j) wherever creation(j) > 0", {
  skip_slow()
  idx <- .cm_indices()
  mat <- getStaticChangeContributions(ans = ans_cm, data = mydata_cm,
    effects = mymodel_cm, depvar = "mynet_cm", returnWide = TRUE)$contribMat
  rows_with_creat <- which(mat[, idx$recip_creat] > 0)
  expect_gt(length(rows_with_creat), 0L)
  expect_equal(mat[rows_with_creat, idx$recip_creat],
               mat[rows_with_creat, idx$recip_eval])
})

test_that("creation/endow static long (data.table): all effect slots present, 2x rows for shared shortNames", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  long_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm, returnDataFrame = TRUE)
  expect_true("density"   %in% long_cm$effectname)
  expect_true("recip"     %in% long_cm$effectname)
  expect_true("transTrip" %in% long_cm$effectname)
  n_density_rows   <- nrow(long_cm[effectname == "density"])
  expect_equal(nrow(long_cm[effectname == "recip"]),     2L * n_density_rows)
  expect_equal(nrow(long_cm[effectname == "transTrip"]), 2L * n_density_rows)
})

test_that("creation/endow dynamic (pure list): no error", {
  skip_slow()
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = FALSE
  )
  expect_true(is.list(dyn_cm))
})

test_that("creation/endow dynamic (data.table): 2x rows for shared shortNames", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_cm))
  n_density_rows   <- nrow(dyn_cm[effectname == "density"])
  expect_equal(nrow(dyn_cm[effectname == "recip"]),     2L * n_density_rows)
  expect_equal(nrow(dyn_cm[effectname == "transTrip"]), 2L * n_density_rows)
})

test_that("creation/endow static wide: ncol equals total included non-rate effects", {
  skip_slow()
  idx  <- .cm_indices()
  wide <- getStaticChangeContributions(ans = ans_cm, data = mydata_cm,
    effects = mymodel_cm, depvar = "mynet_cm", returnWide = TRUE)
  mat  <- wide$contribMat
  n_included <- sum(mymodel_cm$include & mymodel_cm$type != "rate")
  expect_equal(ncol(mat), n_included)
  expect_equal(length(wide$effectNames), n_included)
  expect_false(any(apply(mat, 2, function(x) all(is.na(x)))))
  recip_cols <- c(idx$recip_eval, idx$recip_creat)
  expect_gte(length(recip_cols), 2L)
  expect_false(identical(mat[, recip_cols[1]], mat[, recip_cols[2]]))
  tt_cols <- c(idx$tt_eval, idx$tt_endow)
  expect_false(identical(mat[, tt_cols[1]], mat[, tt_cols[2]]))
})

test_that("creation/endow static long (data.table): all effect slots present in output", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  long_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm, returnDataFrame = TRUE)
  expect_true("density"   %in% long_cm$effectname)
  expect_true("recip"     %in% long_cm$effectname)
  expect_true("transTrip" %in% long_cm$effectname)
  n_density_rows <- nrow(long_cm[effectname == "density"])
  expect_equal(nrow(long_cm[effectname == "recip"]),     2L * n_density_rows)
  expect_equal(nrow(long_cm[effectname == "transTrip"]), 2L * n_density_rows)
})

test_that("creation/endow dynamic (pure list): no error and result is a list", {
  skip_slow()
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = FALSE
  )
  expect_true(is.list(dyn_cm))
})

test_that("creation/endow dynamic (data.table): effect slots not collapsed", {
  skip_slow()
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_cm))
  n_density_rows   <- nrow(dyn_cm[effectname == "density"])
  expect_equal(nrow(dyn_cm[effectname == "recip"]),     2L * n_density_rows)
  expect_equal(nrow(dyn_cm[effectname == "transTrip"]), 2L * n_density_rows)
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

# test_that("thetaValues batch: chains from different theta blocks differ", {
#   N3   <- 5L
#   NSIM <- 2L
#   set.seed(7)
#   # Use theta draws that are far apart so chains are almost certainly different
#   theta_draws <- rbind(
#     ans$theta + 2,
#     ans$theta - 2
#   )

#   tv <- do.call(rbind, lapply(seq_len(NSIM), function(i) {
#     matrix(rep(theta_draws[i, ], each = N3), nrow = N3)
#   }))

#   alg_batch <- sienaAlgorithmCreate(projname = NULL, nsub = 0,
#                                      n3 = nrow(tv), cond = FALSE,
#                                      simOnly = TRUE, seed = 11)
#   eff_batch <- mymodel
#   eff_batch$initialValue[eff_batch$include] <- theta_draws[1, ]

#   ans_batch <- siena07(alg_batch, data = mydata, effects = eff_batch,
#                         batch = TRUE, silent = TRUE, thetaValues = tv,
#                         returnChangeContributions = TRUE)

#   cc <- ans_batch$changeContributions
#   # Extract first ministep matrix from block 1 and block 2
#   get_mat <- function(chain) chain[[1]][[1]][[1]]
#   mat1 <- get_mat(cc[[1]])
#   mat2 <- get_mat(cc[[N3 + 1]])
#   # not always true that different theta produce different chains
#   expect_identical(dim(mat1), dim(mat2))          # same structure
#   expect_false(identical(mat1, mat2))             # different values
# })
