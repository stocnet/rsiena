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

# ---------------------------------------------------------------------------
# Two-network co-evolution: non-unique effect shortNames across depvars
# Both networks carry the same effects (e.g. "density", "transTrip"),
# so shortNames are NOT globally unique.  The extraction must keep them
# correctly separated by depvar and must not collapse or mis-assign columns.
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

test_that("co-evolution static wide: each depvar has its own effect columns", {
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

  # Both networks have the same shortNames, i.e. non-unique globally
  expect_true("density_eval"   %in% cols_a)
  expect_true("transTrip_eval" %in% cols_a)
  expect_true("density_eval"   %in% cols_b)
  expect_true("transTrip_eval" %in% cols_b)

  # effectNames stored in the struct must match the columns exactly
  expect_equal(wide_a$effectNames, cols_a)
  expect_equal(wide_b$effectNames, cols_b)

  # The two depvars must produce separate structs with different data
  expect_false(identical(wide_a$contribMat, wide_b$contribMat))

  # No NA columns (mis-assignment would leave all-NA)
  expect_false(any(is.na(wide_a$contribMat[, "density_eval"])))
  expect_false(any(is.na(wide_b$contribMat[, "density_eval"])))
})

test_that("co-evolution static wide (per-depvar): each call returns network-specific struct", {
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
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
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
  returnDataFrame = TRUE,
  silent  = TRUE
)

# Derive column-type mapping from the R effects table (same order as contribMat columns)
.cm_incl <- mymodel_cm[
  mymodel_cm$include & mymodel_cm$type != "rate" & mymodel_cm$name == "mynet_cm", ]
.cm_col_types  <- .cm_incl$type
.cm_col_names  <- .cm_incl$shortName
.recip_eval  <- which(.cm_col_names == "recip"     & .cm_col_types == "eval")
.recip_creat <- which(.cm_col_names == "recip"     & .cm_col_types == "creation")
.tt_eval     <- which(.cm_col_names == "transTrip" & .cm_col_types == "eval")
.tt_endow    <- which(.cm_col_names == "transTrip" & .cm_col_types == "endow")

wide_cm <- getStaticChangeContributions(
  ans = ans_cm, data = mydata_cm, effects = mymodel_cm,
  depvar = "mynet_cm", returnWide = TRUE
)
mat <- wide_cm$contribMat

test_that("creation/endow static wide: all non-rate effect slots are filled", {

  n_included <- nrow(.cm_incl)   # density + recip(eval) + recip(creation) + transTrip(eval) + transTrip(endow)

  # Structure
  expect_equal(ncol(mat), n_included)
  expect_equal(length(wide_cm$effectNames), n_included)
  # No column is entirely NA (would mean a slot was never filled)
  expect_false(any(apply(mat, 2L, function(x) all(is.na(x)))))
  # Each type has its own column
  expect_equal(length(.recip_eval),  1L)
  expect_equal(length(.recip_creat), 1L)
  expect_equal(length(.tt_eval),     1L)
  expect_equal(length(.tt_endow),    1L)
  # Endowment: +contrib(j) for existing ties, 0 for absent ties
  expect_true(all(mat[, .tt_endow] >= 0, na.rm = TRUE))
})

test_that("creation/endow static wide: creation contributions are non-negative", {
  # Creation: +contrib(j) for absent ties, 0 for existing ties
  expect_true(all(mat[, .recip_creat] >= 0, na.rm = TRUE))
})

test_that("creation/endow static wide: endow(j) + eval(j) == 0 wherever endow(j) > 0", {
  # Both endow and eval see the same existing-tie withdrawal:
  #   eval = -calculateContribution(j),  endow = +calculateContribution(j)
  #   → their sum must be exactly zero.
  rows_with_endow <- which(mat[, .tt_endow] > 0)
  expect_gt(length(rows_with_endow), 0L)  # must exist in s501 (real network)
  expect_equal(
    mat[rows_with_endow, .tt_endow] + mat[rows_with_endow, .tt_eval],
    rep(0, length(rows_with_endow))
  )
})

test_that("creation/endow static wide: creation(j) == eval(j) wherever creation(j) > 0", {
  # Both creation and eval see the same absent-tie creation:
  #   eval = +calculateContribution(j),  creation = +calculateContribution(j)
  #   → they must be equal.
  rows_with_creat <- which(mat[, .recip_creat] > 0)
  expect_gt(length(rows_with_creat), 0L)
  expect_equal(
    mat[rows_with_creat, .recip_creat],
    mat[rows_with_creat, .recip_eval]
  )
})

test_that("creation/endow static long (data.table): all effect slots present, 2x rows for shared shortNames", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  long_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm,
    returnDataFrame = TRUE
  )
  expect_true("density"   %in% long_cm$effectname)
  expect_true("recip"     %in% long_cm$effectname)
  expect_true("transTrip" %in% long_cm$effectname)
  n_density_rows   <- nrow(long_cm[effectname == "density"])
  n_recip_rows     <- nrow(long_cm[effectname == "recip"])
  n_transtrip_rows <- nrow(long_cm[effectname == "transTrip"])
  # recip appears twice (eval + creation), transTrip appears twice (eval + endow)
  expect_equal(n_recip_rows,     2L * n_density_rows)
  expect_equal(n_transtrip_rows, 2L * n_density_rows)
})

test_that("creation/endow dynamic (pure list): no error", {
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = FALSE
  )
  expect_true(is.list(dyn_cm))
})

test_that("creation/endow dynamic (data.table): 2x rows for shared shortNames", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_cm))
  n_density_rows   <- nrow(dyn_cm[effectname == "density"])
  n_recip_rows     <- nrow(dyn_cm[effectname == "recip"])
  n_transtrip_rows <- nrow(dyn_cm[effectname == "transTrip"])
  expect_equal(n_recip_rows,     2L * n_density_rows)
  expect_equal(n_transtrip_rows, 2L * n_density_rows)
})

test_that("creation/endow static wide: ncol equals total included non-rate effects", {
  wide_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm,
    depvar = "mynet_cm", returnWide = TRUE
  )
  n_included <- sum(mymodel_cm$include & mymodel_cm$type != "rate")
  mat <- wide_cm$contribMat
  expect_equal(ncol(mat), n_included)
  expect_equal(length(wide_cm$effectNames), n_included)
  # No column may be entirely NA (would indicate a mis-matched index)
  expect_false(any(apply(mat, 2, function(x) all(is.na(x)))))
  # The two "recip" slots must carry different values:
  # recip(creation) is 0 when a tie already exists (choice = existing tie), while
  # recip(eval) is non-zero in both directions.
  recip_cols <- c(.recip_eval, .recip_creat)
  expect_gte(length(recip_cols), 2L)
  expect_false(identical(mat[, recip_cols[1]], mat[, recip_cols[2]]))
  # Same check for transTrip
  tt_cols <- c(.tt_eval, .tt_endow)
  expect_gte(length(tt_cols), 2L)
  expect_false(identical(mat[, tt_cols[1]], mat[, tt_cols[2]]))
})

test_that("creation/endow static long (data.table): all effect slots present in output", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  long_cm <- getStaticChangeContributions(
    ans = ans_cm, data = mydata_cm, effects = mymodel_cm,
    returnDataFrame = TRUE
  )
  expect_true("density"   %in% long_cm$effectname)
  expect_true("recip"     %in% long_cm$effectname)
  expect_true("transTrip" %in% long_cm$effectname)
  # Row count: 5 effect slots × nPeriods × nEgos × nChoices
  # Verify we have more rows than a model with 3 unique shortNames would produce
  # (because the duplicate slots are NOT collapsed)
  n_slots <- sum(mymodel_cm$include & mymodel_cm$type != "rate")
  n_recip_rows    <- nrow(long_cm[effectname == "recip"])
  n_transtrip_rows <- nrow(long_cm[effectname == "transTrip"])
  # Each shortName that has 2 slots should appear with double the rows
  n_density_rows <- nrow(long_cm[effectname == "density"])
  expect_equal(n_recip_rows,    2L * n_density_rows)
  expect_equal(n_transtrip_rows, 2L * n_density_rows)
})

test_that("creation/endow dynamic (pure list): no error and result is a list", {
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = FALSE
  )
  expect_true(is.list(dyn_cm))
})

test_that("creation/endow dynamic (data.table): effect slots not collapsed", {
  skip_if_not(requireNamespace("data.table", quietly = TRUE), "data.table not available")
  library(data.table)
  dyn_cm <- getDynamicChangeContributions(
    ans = ans_cm, data = mydata_cm, algorithm = mycontrols_cm,
    theta = c(ans_cm$rate, ans_cm$theta),
    effects = mymodel_cm, depvar = "mynet_cm", returnDataFrame = TRUE
  )
  expect_true("data.table" %in% class(dyn_cm))
  n_density_rows  <- nrow(dyn_cm[effectname == "density"])
  n_recip_rows    <- nrow(dyn_cm[effectname == "recip"])
  n_transtrip_rows <- nrow(dyn_cm[effectname == "transTrip"])
  expect_equal(n_recip_rows,    2L * n_density_rows)
  expect_equal(n_transtrip_rows, 2L * n_density_rows)
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
