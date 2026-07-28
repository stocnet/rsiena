# Tests for the postestimation configuration-object constructors
# (R/postestConfig.R): set_postest_uncertainty_saom(), set_postest_control_saom(),
# and their print methods.
#
# These constructors are purely additive: they bundle the flat arguments of
# marginalEffects.sienaFit into validated, printable objects. They are not
# (yet) wired into marginalEffects.sienaFit itself, so these tests exercise
# the constructors and print methods directly and compare their defaults
# against the flat-argument defaults for exact equivalence.

flatFormals <- formals(RSiena:::marginalEffects.sienaFit)

# ── A. Class and version attribute ─────────────────────────────────────────

test_that("set_postest_uncertainty_saom returns a sienaPostestUncertainty object with version", {
  u <- set_postest_uncertainty_saom()
  expect_s3_class(u, "sienaPostestUncertainty")
  expect_false(is.null(attr(u, "version")))
  expect_type(attr(u, "version"), "character")
})

test_that("set_postest_control_saom returns a sienaPostestControl object with version", {
  co <- set_postest_control_saom()
  expect_s3_class(co, "sienaPostestControl")
  expect_false(is.null(attr(co, "version")))
  expect_type(attr(co, "version"), "character")
})

# ── B. Defaults match the flat API exactly ─────────────────────────────────

test_that("set_postest_uncertainty_saom defaults match marginalEffects.sienaFit flat defaults", {
  u <- set_postest_uncertainty_saom()

  checks <- list(
    list(field = "enabled",    formal = "uncertainty"),
    list(field = "mode",       formal = "uncertaintyMode"),
    list(field = "nsim",       formal = "nsim"),
    list(field = "sd",         formal = "uncertaintySd"),
    list(field = "ci",         formal = "uncertaintyCi"),
    list(field = "ciInterval", formal = "ciInterval"),
    list(field = "simMean",    formal = "uncertaintyMean"),
    list(field = "simMedian",  formal = "uncertaintyMedian")
  )

  for (chk in checks) {
    formalVal <- eval(flatFormals[[chk$formal]])
    # match.arg-style defaults are given as a character vector; the flat
    # function calls match.arg() internally, which picks the first element.
    if (chk$formal == "uncertaintyMode" && length(formalVal) > 1L) {
      formalVal <- formalVal[1]
    }
    objVal <- u[[chk$field]]
    if (chk$field == "nsim") {
      # nsim is coerced to integer by the constructor; compare numerically.
      expect_equal(as.numeric(objVal), as.numeric(formalVal),
                   info = sprintf("field '%s' vs formal '%s'", chk$field, chk$formal))
    } else {
      expect_equal(objVal, formalVal,
                   info = sprintf("field '%s' vs formal '%s'", chk$field, chk$formal))
    }
  }
})

test_that("set_postest_control_saom defaults match marginalEffects.sienaFit flat defaults", {
  co <- set_postest_control_saom()

  checks <- list(
    list(field = "dynamic",                formal = "dynamic"),
    list(field = "algorithm",               formal = "algorithm"),
    list(field = "n3",                      formal = "n3"),
    list(field = "n3PointEst",              formal = "n3PointEst"),
    list(field = "n3BatchSize",             formal = "n3BatchSize"),
    list(field = "chainStoreMode",          formal = "chainStoreMode"),
    list(field = "useChangeContributions",  formal = "useChangeContributions"),
    list(field = "chainStorePath",          formal = "chainStorePath"),
    list(field = "useCluster",              formal = "useCluster"),
    list(field = "nbrNodes",                formal = "nbrNodes"),
    list(field = "clusterType",             formal = "clusterType"),
    list(field = "cl",                      formal = "cl"),
    list(field = "batchDir",                formal = "batchDir"),
    list(field = "prefix",                  formal = "prefix"),
    list(field = "combineBatch",            formal = "combineBatch"),
    list(field = "batchSize",               formal = "batchSize"),
    list(field = "keepBatch",               formal = "keepBatch"),
    list(field = "verbose",                 formal = "verbose"),
    list(field = "details",                 formal = "details"),
    list(field = "memoryScale",             formal = "memoryScale"),
    list(field = "batchUnitBudget",         formal = "batchUnitBudget"),
    list(field = "dynamicMinistepFactor",   formal = "dynamicMinistepFactor"),
    list(field = "saveDir",                 formal = "saveDir"),
    list(field = "gcEachBatch",             formal = "gcEachBatch"),
    list(field = "gcEachSim",               formal = "gcEachSim")
  )

  intFields <- c("n3", "nbrNodes", "n3BatchSize")

  for (chk in checks) {
    formalVal <- eval(flatFormals[[chk$formal]])
    if (chk$formal %in% c("chainStoreMode", "clusterType") && length(formalVal) > 1L) {
      formalVal <- formalVal[1]
    }
    objVal <- co[[chk$field]]
    if (chk$field %in% intFields) {
      expect_equal(as.numeric(objVal), as.numeric(formalVal),
                   info = sprintf("field '%s' vs formal '%s'", chk$field, chk$formal))
    } else {
      expect_equal(objVal, formalVal,
                   info = sprintf("field '%s' vs formal '%s'", chk$field, chk$formal))
    }
  }

  # batchUnitBudget must not be rounded away from 2.5e8, and n3BatchSize
  # must remain an integer (100L), exactly as specified.
  expect_equal(co$batchUnitBudget, 2.5e8)
  expect_true(is.integer(co$n3BatchSize))
  expect_equal(co$n3BatchSize, 100L)
})

# ── C. Validation rules: set_postest_uncertainty_saom ──────────────────────

test_that("set_postest_uncertainty_saom: mode is matched via match.arg", {
  # match.arg()'s error text is localized (e.g. German under a non-English
  # locale), so match on the (untranslated) choice list rather than the
  # surrounding English sentence.
  expect_error(set_postest_uncertainty_saom(mode = "bogus"),
               regexp = "bootstrap")
})

test_that("set_postest_uncertainty_saom: enabled must be single non-NA logical", {
  expect_error(set_postest_uncertainty_saom(enabled = NA),
               regexp = "'enabled'.*logical")
  expect_error(set_postest_uncertainty_saom(enabled = c(TRUE, FALSE)),
               regexp = "'enabled'.*logical")
  expect_error(set_postest_uncertainty_saom(enabled = "yes"),
               regexp = "'enabled'.*logical")
})

test_that("set_postest_uncertainty_saom: sd/ci/simMean/simMedian must be single non-NA logical", {
  expect_error(set_postest_uncertainty_saom(sd = NA), regexp = "'sd'.*logical")
  expect_error(set_postest_uncertainty_saom(ci = NA), regexp = "'ci'.*logical")
  expect_error(set_postest_uncertainty_saom(mode = "bootstrap", simMean = NA),
               regexp = "'simMean'.*logical")
  expect_error(set_postest_uncertainty_saom(mode = "bootstrap", simMedian = NA),
               regexp = "'simMedian'.*logical")
})

test_that("set_postest_uncertainty_saom: nsim must be single finite positive numeric", {
  expect_error(set_postest_uncertainty_saom(nsim = 0), regexp = "'nsim'.*> 0")
  expect_error(set_postest_uncertainty_saom(nsim = -5), regexp = "'nsim'.*> 0")
  expect_error(set_postest_uncertainty_saom(nsim = Inf), regexp = "'nsim'.*finite")
  expect_error(set_postest_uncertainty_saom(nsim = c(1, 2)), regexp = "'nsim'.*finite")
})

test_that("set_postest_uncertainty_saom: nsim is coerced to integer", {
  u <- set_postest_uncertainty_saom(nsim = 250.9)
  expect_true(is.integer(u$nsim))
  expect_equal(u$nsim, 250L)
})

test_that("set_postest_uncertainty_saom: ciInterval must be length-2 finite numeric", {
  expect_error(set_postest_uncertainty_saom(ciInterval = 0.5),
               regexp = "'ciInterval'.*length 2")
  expect_error(set_postest_uncertainty_saom(ciInterval = c(0.025, NA)),
               regexp = "'ciInterval'.*finite")
  expect_error(set_postest_uncertainty_saom(ciInterval = c(0.025, Inf)),
               regexp = "'ciInterval'.*finite")
})

test_that("set_postest_uncertainty_saom: ciInterval must be strictly within (0, 1)", {
  expect_error(set_postest_uncertainty_saom(ciInterval = c(0, 0.975)),
               regexp = "'ciInterval'.*between 0 and 1")
  expect_error(set_postest_uncertainty_saom(ciInterval = c(0.025, 1)),
               regexp = "'ciInterval'.*between 0 and 1")
})

test_that("set_postest_uncertainty_saom: ciInterval must be strictly increasing", {
  expect_error(set_postest_uncertainty_saom(ciInterval = c(0.975, 0.025)),
               regexp = "'ciInterval'.*increasing")
  expect_error(set_postest_uncertainty_saom(ciInterval = c(0.5, 0.5)),
               regexp = "'ciInterval'.*increasing")
})

# ── D. simMean/simMedian vs. mode interaction ───────────────────────────────

test_that("simMean/simMedian with mode='delta' error", {
  expect_error(set_postest_uncertainty_saom(mode = "delta", simMean = TRUE),
               regexp = "simMean = TRUE requires mode = 'bootstrap'")
  expect_error(set_postest_uncertainty_saom(mode = "delta", simMedian = TRUE),
               regexp = "simMedian = TRUE requires mode = 'bootstrap'")
})

test_that("simMean/simMedian with mode='deltaFull' error", {
  expect_error(set_postest_uncertainty_saom(mode = "deltaFull", simMean = TRUE),
               regexp = "simMean = TRUE requires mode = 'bootstrap'")
  expect_error(set_postest_uncertainty_saom(mode = "deltaFull", simMedian = TRUE),
               regexp = "simMedian = TRUE requires mode = 'bootstrap'")
})

test_that("simMean/simMedian with mode='bootstrap' succeed", {
  expect_silent(u1 <- set_postest_uncertainty_saom(mode = "bootstrap", simMean = TRUE))
  expect_true(u1$simMean)
  expect_silent(u2 <- set_postest_uncertainty_saom(mode = "bootstrap", simMedian = TRUE))
  expect_true(u2$simMedian)
  expect_silent(u3 <- set_postest_uncertainty_saom(mode = "bootstrap",
                                                     simMean = TRUE, simMedian = TRUE))
  expect_true(u3$simMean && u3$simMedian)
})

# ── C/F/G. Validation rules: set_postest_control_saom ──────────────────────

test_that("set_postest_control_saom: chainStoreMode/clusterType are matched via match.arg", {
  # match.arg()'s error text is localized; match on the (untranslated)
  # choice list rather than the surrounding English sentence.
  expect_error(set_postest_control_saom(chainStoreMode = "bogus"),
               regexp = "auto")
  expect_error(set_postest_control_saom(clusterType = "bogus"),
               regexp = "PSOCK")
})

test_that("set_postest_control_saom: dynamic = TRUE without algorithm errors", {
  expect_error(set_postest_control_saom(dynamic = TRUE),
               regexp = "'algorithm' must be provided when dynamic = TRUE")
})

test_that("set_postest_control_saom: dynamic = TRUE with a dummy algorithm succeeds", {
  dummyAlgo <- list(dummy = TRUE)
  expect_silent(co <- set_postest_control_saom(dynamic = TRUE, algorithm = dummyAlgo))
  expect_true(co$dynamic)
  expect_identical(co$algorithm, dummyAlgo)
})

test_that("set_postest_control_saom: logical flags must be single non-NA logical", {
  expect_error(set_postest_control_saom(dynamic = NA), regexp = "'dynamic'.*logical")
  expect_error(set_postest_control_saom(useChangeContributions = NA),
               regexp = "'useChangeContributions'.*logical")
  expect_error(set_postest_control_saom(useCluster = NA), regexp = "'useCluster'.*logical")
  expect_error(set_postest_control_saom(combineBatch = NA), regexp = "'combineBatch'.*logical")
  expect_error(set_postest_control_saom(keepBatch = NA), regexp = "'keepBatch'.*logical")
  expect_error(set_postest_control_saom(details = NA), regexp = "'details'.*logical")
  expect_error(set_postest_control_saom(gcEachBatch = NA), regexp = "'gcEachBatch'.*logical")
  expect_error(set_postest_control_saom(gcEachSim = NA), regexp = "'gcEachSim'.*logical")
})

test_that("set_postest_control_saom: n3/nbrNodes/n3BatchSize must be single finite positive numeric", {
  expect_error(set_postest_control_saom(n3 = 0), regexp = "'n3'.*> 0")
  expect_error(set_postest_control_saom(n3 = -1), regexp = "'n3'.*> 0")
  expect_error(set_postest_control_saom(nbrNodes = Inf), regexp = "'nbrNodes'.*finite")
  expect_error(set_postest_control_saom(n3BatchSize = 0), regexp = "'n3BatchSize'.*> 0")
})

test_that("set_postest_control_saom: n3/n3BatchSize are coerced to integer", {
  co <- set_postest_control_saom(n3 = 250.9, n3BatchSize = 50.1)
  expect_true(is.integer(co$n3))
  expect_equal(co$n3, 250L)
  expect_true(is.integer(co$n3BatchSize))
  expect_equal(co$n3BatchSize, 50L)
})

test_that("set_postest_control_saom: n3PointEst/batchSize allow NULL, else same rule", {
  expect_silent(co <- set_postest_control_saom(n3PointEst = NULL, batchSize = NULL))
  expect_null(co$n3PointEst)
  expect_null(co$batchSize)
  expect_error(set_postest_control_saom(n3PointEst = -1), regexp = "'n3PointEst'.*> 0")
  expect_error(set_postest_control_saom(batchSize = 0), regexp = "'batchSize'.*> 0")
  co2 <- set_postest_control_saom(n3PointEst = 300.7, batchSize = 20.2)
  expect_equal(co2$n3PointEst, 300L)
  expect_equal(co2$batchSize, 20L)
})

test_that("set_postest_control_saom: batchUnitBudget/dynamicMinistepFactor must be finite positive, not coerced to integer", {
  expect_error(set_postest_control_saom(batchUnitBudget = 0), regexp = "'batchUnitBudget'.*> 0")
  expect_error(set_postest_control_saom(dynamicMinistepFactor = -1),
               regexp = "'dynamicMinistepFactor'.*> 0")
  co <- set_postest_control_saom(batchUnitBudget = 1.23e7, dynamicMinistepFactor = 7.5)
  expect_false(is.integer(co$batchUnitBudget))
  expect_equal(co$batchUnitBudget, 1.23e7)
  expect_false(is.integer(co$dynamicMinistepFactor))
  expect_equal(co$dynamicMinistepFactor, 7.5)
})

test_that("set_postest_control_saom: memoryScale allows NULL, else single finite positive numeric", {
  expect_silent(co <- set_postest_control_saom(memoryScale = NULL))
  expect_null(co$memoryScale)
  expect_error(set_postest_control_saom(memoryScale = 0), regexp = "'memoryScale'.*> 0")
  co2 <- set_postest_control_saom(memoryScale = 2.5)
  expect_equal(co2$memoryScale, 2.5)
})

test_that("set_postest_control_saom: batchDir/prefix must be single non-NA character", {
  expect_error(set_postest_control_saom(batchDir = NA), regexp = "'batchDir'.*character")
  expect_error(set_postest_control_saom(prefix = c("a", "b")), regexp = "'prefix'.*character")
})

test_that("set_postest_control_saom: chainStorePath/saveDir allow NULL, else single non-NA character", {
  expect_silent(co <- set_postest_control_saom(chainStorePath = NULL, saveDir = NULL))
  expect_null(co$chainStorePath)
  expect_null(co$saveDir)
  expect_error(set_postest_control_saom(chainStorePath = NA),
               regexp = "'chainStorePath'.*character")
  expect_error(set_postest_control_saom(saveDir = NA), regexp = "'saveDir'.*character")
  co2 <- set_postest_control_saom(chainStorePath = "/tmp/chains", saveDir = "/tmp/save")
  expect_equal(co2$chainStorePath, "/tmp/chains")
  expect_equal(co2$saveDir, "/tmp/save")
})

# ── E. Cluster normalisation ────────────────────────────────────────────────

test_that("passing cl with useCluster = FALSE turns useCluster on and sets nbrNodes", {
  dummyCl <- list(1, 2, 3)  # stand-in for a real cluster object; length is all that matters
  co <- set_postest_control_saom(useCluster = FALSE, cl = dummyCl)
  expect_true(co$useCluster)
  expect_equal(co$nbrNodes, 3L)
})

test_that("cluster normalisation does not override an already-TRUE useCluster / explicit nbrNodes", {
  dummyCl <- list(1, 2, 3)
  co <- set_postest_control_saom(useCluster = TRUE, nbrNodes = 8, cl = dummyCl)
  expect_true(co$useCluster)
  expect_equal(co$nbrNodes, 8L)
})

# ── G. verbose is a level, not a flag ───────────────────────────────────────

test_that("verbose = 2 is accepted and stored as 2, not coerced to logical", {
  co <- set_postest_control_saom(verbose = 2)
  expect_equal(co$verbose, 2)
  expect_false(is.logical(co$verbose))
})

test_that("verbose still rejects non-scalar / NA input", {
  expect_error(set_postest_control_saom(verbose = NA), regexp = "'verbose'")
  expect_error(set_postest_control_saom(verbose = c(1, 2)), regexp = "'verbose'")
  expect_error(set_postest_control_saom(verbose = "loud"), regexp = "'verbose'")
})

# ── H. Print methods ─────────────────────────────────────────────────────────

test_that("print.sienaPostestUncertainty prints key settings (enabled, bootstrap)", {
  u <- set_postest_uncertainty_saom(mode = "bootstrap", nsim = 777)
  expect_output(print(u), "RSiena postestimation uncertainty")
  expect_output(print(u), "bootstrap")
  expect_output(print(u), "777")
  expect_output(print(u), "SE")
  expect_output(print(u), "95%")
  expect_output(print(u), "\\[0.025, 0.975\\]")
})

test_that("print.sienaPostestUncertainty reports disabled uncertainty plainly", {
  u <- set_postest_uncertainty_saom(enabled = FALSE)
  expect_output(print(u), "disabled")
  expect_output(print(u), "point estimates only")
})

test_that("print.sienaPostestUncertainty omits nsim and shows 'no draws' for delta modes", {
  u <- set_postest_uncertainty_saom(mode = "delta")
  out <- capture.output(print(u))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "delta")
  expect_match(txt, "no simulation draws")
  expect_false(grepl("nsim", txt))
})

test_that("print.sienaPostestControl prints key settings and shows single-core by default", {
  co <- set_postest_control_saom()
  expect_output(print(co), "RSiena postestimation control")
  expect_output(print(co), "static")
  expect_output(print(co), "single-core")
  expect_output(print(co), "n3 = 200")
})

test_that("print.sienaPostestControl shows cluster details when useCluster = TRUE", {
  co <- set_postest_control_saom(useCluster = TRUE, nbrNodes = 4, clusterType = "PSOCK")
  out <- capture.output(print(co))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "4 nodes")
  expect_match(txt, "PSOCK")
  expect_false(grepl("single-core", txt))
})

test_that("print.sienaPostestControl shows saveDir only when non-NULL", {
  co1 <- set_postest_control_saom()
  expect_false(grepl("saveDir", paste(capture.output(print(co1)), collapse = "\n")))

  co2 <- set_postest_control_saom(saveDir = "/tmp/saveHere")
  expect_output(print(co2), "/tmp/saveHere")
})

test_that("print methods return their argument invisibly", {
  u <- set_postest_uncertainty_saom()
  ret <- withVisible(print(u))
  expect_false(ret$visible)
  expect_identical(ret$value, u)

  co <- set_postest_control_saom()
  ret2 <- withVisible(print(co))
  expect_false(ret2$visible)
  expect_identical(ret2$value, co)
})
