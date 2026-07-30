# Tests for the postestimation configuration-object constructors
# (R/postestConfig.R): set_postest_uncertainty_saom(), set_postest_algo_saom(),
# and their print methods.
#
# These constructors ARE the interface: since step 5c removed the flat
# arguments, marginalEffects.sienaFit takes nothing but the objects, and each
# constructor's defaults are what a call gets when it says nothing.
#
# Section B used to assert those defaults matched the flat formals one for
# one. That referent is gone, so it now pins the values themselves -- the
# regression it was really guarding against is a default changing silently,
# and that guard does not need a second implementation to compare with.

# ── A. Class and version attribute ─────────────────────────────────────────

test_that("set_postest_uncertainty_saom returns a sienaPostestUncertainty object with version", {
  u <- set_postest_uncertainty_saom()
  expect_s3_class(u, "sienaPostestUncertainty")
  expect_false(is.null(attr(u, "version")))
  expect_type(attr(u, "version"), "character")
})

test_that("set_postest_algo_saom returns a sienaPostestControl object with version", {
  co <- set_postest_algo_saom()
  expect_s3_class(co, "sienaPostestControl")
  expect_false(is.null(attr(co, "version")))
  expect_type(attr(co, "version"), "character")
})

# ── B. Defaults are what we think they are ─────────────────────────────────
#
# A default that changes silently changes every call that did not override it.
# Spelled out as literals rather than derived from anything, so that a change
# has to be made here as well as in the constructor -- which is the point.

test_that("set_postest_uncertainty_saom defaults are unchanged", {
  u <- set_postest_uncertainty_saom()
  expect_true(u$enabled)
  expect_equal(u$mode, "bootstrap")
  expect_equal(as.numeric(u$nsim), 1000)
  expect_true(u$sd)
  expect_true(u$ci)
  expect_equal(u$ciInterval, c(0.025, 0.975))
  expect_false(u$simMean)
  expect_false(u$simMedian)
})

test_that("set_postest_algo_saom defaults are unchanged", {
  co <- set_postest_algo_saom()
  expect_null(co$algorithm)
  expect_equal(as.numeric(co$n3), 200)
  expect_null(co$n3PointEst)
  expect_equal(as.numeric(co$n3BatchSize), 100)
  expect_equal(co$chainStoreMode, "auto")
  expect_false(co$useChangeContributions)
  expect_null(co$chainStorePath)
  expect_false(co$useCluster)
  expect_equal(as.numeric(co$nbrNodes), 1)
  expect_equal(co$clusterType, "PSOCK")
  expect_null(co$cl)
  expect_equal(co$batchDir, "temp")
  expect_equal(co$prefix, "simBatch_b")
  expect_true(co$combineBatch)
  expect_null(co$batchSize)
  expect_false(co$keepBatch)
  expect_true(co$verbose)
  expect_null(co$memoryScale)
  ## Must not be rounded away from 2.5e8 by integer coercion.
  expect_equal(co$batchUnitBudget, 2.5e8)
  expect_equal(as.numeric(co$dynamicMinistepFactor), 10)
  expect_null(co$saveDir)
  expect_false(co$gcEachBatch)
  expect_false(co$gcEachSim)
})

test_that("set_postest_output_saom defaults are unchanged", {
  o <- set_postest_output_saom()
  expect_equal(o$format, "long")
  expect_true(o$combineSameLevel)
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

# ── C/F/G. Validation rules: set_postest_algo_saom ──────────────────────

test_that("set_postest_algo_saom: chainStoreMode/clusterType are matched via match.arg", {
  # match.arg()'s error text is localized; match on the (untranslated)
  # choice list rather than the surrounding English sentence.
  expect_error(set_postest_algo_saom(chainStoreMode = "bogus"),
               regexp = "auto")
  expect_error(set_postest_algo_saom(clusterType = "bogus"),
               regexp = "PSOCK")
})

## `dynamic` moved to make_postest_targets(): it selects the estimand (static vs
## model-implied), not the route to it.  Its consistency rules (accumulated
## requires dynamic; rateWeight inert under dynamic) are therefore checked at
## target construction -- see test-postest-targets-dynamic.R.  The
## `dynamic -> algorithm` requirement is now a call-time check, since the two
## live in different objects.

test_that("set_postest_algo_saom: logical flags must be single non-NA logical", {
    expect_error(set_postest_algo_saom(useChangeContributions = NA),
               regexp = "'useChangeContributions'.*logical")
  expect_error(set_postest_algo_saom(useCluster = NA), regexp = "'useCluster'.*logical")
  expect_error(set_postest_algo_saom(combineBatch = NA), regexp = "'combineBatch'.*logical")
  expect_error(set_postest_algo_saom(keepBatch = NA), regexp = "'keepBatch'.*logical")
  expect_error(set_postest_algo_saom(gcEachBatch = NA), regexp = "'gcEachBatch'.*logical")
  expect_error(set_postest_algo_saom(gcEachSim = NA), regexp = "'gcEachSim'.*logical")
})

test_that("set_postest_algo_saom: n3/nbrNodes/n3BatchSize must be single finite positive numeric", {
  expect_error(set_postest_algo_saom(n3 = 0), regexp = "'n3'.*> 0")
  expect_error(set_postest_algo_saom(n3 = -1), regexp = "'n3'.*> 0")
  expect_error(set_postest_algo_saom(nbrNodes = Inf), regexp = "'nbrNodes'.*finite")
  expect_error(set_postest_algo_saom(n3BatchSize = 0), regexp = "'n3BatchSize'.*> 0")
})

test_that("set_postest_algo_saom: n3/n3BatchSize are coerced to integer", {
  co <- set_postest_algo_saom(n3 = 250.9, n3BatchSize = 50.1)
  expect_true(is.integer(co$n3))
  expect_equal(co$n3, 250L)
  expect_true(is.integer(co$n3BatchSize))
  expect_equal(co$n3BatchSize, 50L)
})

test_that("set_postest_algo_saom: n3PointEst/batchSize allow NULL, else same rule", {
  expect_silent(co <- set_postest_algo_saom(n3PointEst = NULL, batchSize = NULL))
  expect_null(co$n3PointEst)
  expect_null(co$batchSize)
  expect_error(set_postest_algo_saom(n3PointEst = -1), regexp = "'n3PointEst'.*> 0")
  expect_error(set_postest_algo_saom(batchSize = 0), regexp = "'batchSize'.*> 0")
  co2 <- set_postest_algo_saom(n3PointEst = 300.7, batchSize = 20.2)
  expect_equal(co2$n3PointEst, 300L)
  expect_equal(co2$batchSize, 20L)
})

test_that("set_postest_algo_saom: batchUnitBudget/dynamicMinistepFactor must be finite positive, not coerced to integer", {
  expect_error(set_postest_algo_saom(batchUnitBudget = 0), regexp = "'batchUnitBudget'.*> 0")
  expect_error(set_postest_algo_saom(dynamicMinistepFactor = -1),
               regexp = "'dynamicMinistepFactor'.*> 0")
  co <- set_postest_algo_saom(batchUnitBudget = 1.23e7, dynamicMinistepFactor = 7.5)
  expect_false(is.integer(co$batchUnitBudget))
  expect_equal(co$batchUnitBudget, 1.23e7)
  expect_false(is.integer(co$dynamicMinistepFactor))
  expect_equal(co$dynamicMinistepFactor, 7.5)
})

test_that("set_postest_algo_saom: memoryScale allows NULL, else single finite positive numeric", {
  expect_silent(co <- set_postest_algo_saom(memoryScale = NULL))
  expect_null(co$memoryScale)
  expect_error(set_postest_algo_saom(memoryScale = 0), regexp = "'memoryScale'.*> 0")
  co2 <- set_postest_algo_saom(memoryScale = 2.5)
  expect_equal(co2$memoryScale, 2.5)
})

test_that("set_postest_algo_saom: batchDir/prefix must be single non-NA character", {
  expect_error(set_postest_algo_saom(batchDir = NA), regexp = "'batchDir'.*character")
  expect_error(set_postest_algo_saom(prefix = c("a", "b")), regexp = "'prefix'.*character")
})

test_that("set_postest_algo_saom: chainStorePath/saveDir allow NULL, else single non-NA character", {
  expect_silent(co <- set_postest_algo_saom(chainStorePath = NULL, saveDir = NULL))
  expect_null(co$chainStorePath)
  expect_null(co$saveDir)
  expect_error(set_postest_algo_saom(chainStorePath = NA),
               regexp = "'chainStorePath'.*character")
  expect_error(set_postest_algo_saom(saveDir = NA), regexp = "'saveDir'.*character")
  co2 <- set_postest_algo_saom(chainStorePath = "/tmp/chains", saveDir = "/tmp/save")
  expect_equal(co2$chainStorePath, "/tmp/chains")
  expect_equal(co2$saveDir, "/tmp/save")
})

# ── E. Cluster normalisation ────────────────────────────────────────────────

test_that("passing cl with useCluster = FALSE turns useCluster on and sets nbrNodes", {
  dummyCl <- list(1, 2, 3)  # stand-in for a real cluster object; length is all that matters
  co <- set_postest_algo_saom(useCluster = FALSE, cl = dummyCl)
  expect_true(co$useCluster)
  expect_equal(co$nbrNodes, 3L)
})

test_that("cluster normalisation does not override an already-TRUE useCluster / explicit nbrNodes", {
  dummyCl <- list(1, 2, 3)
  co <- set_postest_algo_saom(useCluster = TRUE, nbrNodes = 8, cl = dummyCl)
  expect_true(co$useCluster)
  expect_equal(co$nbrNodes, 8L)
})

# ── G. verbose is a level, not a flag ───────────────────────────────────────

test_that("verbose = 2 is accepted and stored as 2, not coerced to logical", {
  co <- set_postest_algo_saom(verbose = 2)
  expect_equal(co$verbose, 2)
  expect_false(is.logical(co$verbose))
})

test_that("verbose still rejects non-scalar / NA input", {
  expect_error(set_postest_algo_saom(verbose = NA), regexp = "'verbose'")
  expect_error(set_postest_algo_saom(verbose = c(1, 2)), regexp = "'verbose'")
  expect_error(set_postest_algo_saom(verbose = "loud"), regexp = "'verbose'")
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
  co <- set_postest_algo_saom()
  expect_output(print(co), "RSiena postestimation control")
  ## No longer prints static/dynamic: `dynamic` selects the estimand and lives
  ## on the targets object, not here.
  expect_output(print(co), "single-core")
  expect_output(print(co), "n3 = 200")
})

test_that("print.sienaPostestControl shows cluster details when useCluster = TRUE", {
  co <- set_postest_algo_saom(useCluster = TRUE, nbrNodes = 4, clusterType = "PSOCK")
  out <- capture.output(print(co))
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "4 nodes")
  expect_match(txt, "PSOCK")
  expect_false(grepl("single-core", txt))
})

test_that("print.sienaPostestControl shows saveDir only when non-NULL", {
  co1 <- set_postest_algo_saom()
  expect_false(grepl("saveDir", paste(capture.output(print(co1)), collapse = "\n")))

  co2 <- set_postest_algo_saom(saveDir = "/tmp/saveHere")
  expect_output(print(co2), "/tmp/saveHere")
})

test_that("print methods return their argument invisibly", {
  u <- set_postest_uncertainty_saom()
  ret <- withVisible(print(u))
  expect_false(ret$visible)
  expect_identical(ret$value, u)

  co <- set_postest_algo_saom()
  ret2 <- withVisible(print(co))
  expect_false(ret2$visible)
  expect_identical(ret2$value, co)
})

# ── I. Printing a targets object ─────────────────────────────────────────────
#
# The printed targets object is the only place a user sees, in one view, what
# quantity they have asked for.  Three things it has to get right, each of
# which it previously got wrong:
#
#   * the model-level settings are the ones every target inherits, not a
#     record of what the software would have chosen -- so they have to be
#     labelled as the estimand and its reporting, not as "defaults"
#   * "egoX" does not identify a target in a model carrying egoX on several
#     covariates, so the covariate has to appear
#   * a second difference showed only its first component's perturbation,
#     which named neither the statistic being stepped nor the fact that the
#     partner may be a contrast rather than a step

targets_print <- function(tg) paste(capture.output(print(tg)), collapse = "\n")

test_that("printed targets state the estimand rather than listing 'defaults'", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             type = "tieProb", condition = "recip")
  txt <- targets_print(tg)
  expect_match(txt, "risk difference in tie probability")
  expect_match(txt, "static")
  expect_match(txt, "conditional on recip")
  ## "defaults" invited the reading "what would have been used anyway".
  expect_false(grepl("defaults", txt))

  dyn <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                              dynamic = TRUE, mainEffect = "riskRatio")
  dtxt <- targets_print(dyn)
  expect_match(dtxt, "risk ratio in change probability")
  expect_match(dtxt, "dynamic")
  expect_false(grepl("static", dtxt))
})

test_that("printed targets name the covariate an effect refers to", {
  skip_on_cran()
  ans_ego <- load_fixture("ans_ego"); mymodel_ego <- load_fixture("mymodel_ego")
  skip_if(is.null(ans_ego) || is.null(mymodel_ego), "ego fixtures unavailable")

  tg <- make_postest_targets(ans_ego, effects = mymodel_ego)
  txt <- targets_print(tg)
  ## The short name stays -- it is what set_target() is called with -- but on
  ## its own it does not say which covariate.
  expect_match(txt, "egoX \\(mybeh_ego\\)")
  ## Effects without a covariate must not grow an empty parenthesis.
  expect_match(txt, "transTrip \\+1")
  expect_false(grepl("\\(\\)", txt))
})

test_that("printed second differences show both perturbations, not just the first", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_second_diff(tg,
                                         list(transTrip = list(diff = 1),
                                              recip     = list(contrast = c(0, 1))),
                                         name = "interaction_sd"))
  txt <- targets_print(tg)
  expect_match(txt, "second difference")
  ## Both components, each attributed to its effect, and the contrast shown
  ## as two levels rather than as a step.
  expect_match(txt, "transTrip \\+1")
  expect_match(txt, "recip 0 -> 1")
})

test_that("printed first differences distinguish a step from a contrast", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2")
  txt <- targets_print(tg)
  ## density defaults to a contrast, everything else to a unit step.
  expect_match(txt, "density -1 -> 1")
  expect_match(txt, "transTrip \\+1")
  expect_match(txt, "first difference")
})

test_that("printed targets show override VALUES, and only claim overriding when it happens", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  plain <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2")
  ## Nothing overrides anything, so the sentence describing overrides would
  ## be describing a distinction this object does not make.
  expect_false(grepl("bracketed", targets_print(plain)))

  ov <- suppressMessages(set_target(plain, transTrip, diff = 1,
                                    level = "actor"))
  txt <- targets_print(ov)
  expect_match(txt, "bracketed")
  ## The field name alone is not actionable; the value is the point.
  expect_match(txt, "\\[level = actor\\]")
})

test_that("printed target order is the order results come back in", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  ## Built up in an order that differs from the model's effect order, which is
  ## the order the table itself is in.
  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
  tg <- suppressMessages(set_target(tg, density, contrast = c(-1, 1)))

  printed <- grep("first difference", capture.output(print(tg)), value = TRUE)
  printed <- sub("^\\s*(\\S+).*$", "\\1", printed)
  expect_identical(printed, names(RSiena:::.targetsToEffectList(tg)$effectList))
})

# ── J. The per-effect perturbation list ──────────────────────────────────────
#
# A second difference perturbs two effects.  The superseded numbered form
# (diff1 = 1, contrast2 = c(0, 1)) made the reader carry the mapping from
# suffix to position; the list form states each perturbation next to the
# effect it applies to, using the same vocabulary set_target() uses:
#
#   list(transTrip = list(diff = 1), recip = list(contrast = c(0, 1)))
#
# What these check is that each entry's settings reach the component it
# names -- the mapping that used to be the caller's job -- and that what is
# left of the numbered form fails loudly rather than being ignored.

test_that("each list entry lowers onto the component it names", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  lst <- suppressMessages(set_second_diff(tg,
      list(transTrip = list(diff = 1), recip = list(contrast = c(0, 1))),
      name = "sd"))
  sp <- RSiena:::.targetsToEffectList(lst)$effectList[["sd"]]

  ## List order is component order, and each entry's settings land on ITS
  ## component -- the mapping the numbered form made the reader do by hand.
  expect_true(isTRUE(sp$second))
  expect_equal(sp$effectName1, "transTrip")
  expect_equal(sp$diff1, 1)
  expect_null(sp$contrast1)
  expect_equal(sp$effectName2, "recip")
  expect_equal(sp$contrast2, c(0, 1))
  expect_null(sp$diff2)
})

test_that("the perturbation list carries the compound-effect settings too", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  lst <- suppressMessages(set_second_diff(tg, list(
      recip     = list(contrast = c(0, 1), interaction = TRUE,
                       intEffectNames = "transRecTrip",
                       modEffectNames = "transTrip"),
      transTrip = list(diff = 1, interaction = TRUE,
                       intEffectNames = "transRecTrip",
                       modEffectNames = "recip")), name = "sd"))
  sp <- RSiena:::.targetsToEffectList(lst)$effectList[["sd"]]

  ## The moderator pairing is the error the numbered form invited: getting
  ## modEffectNames1/2 the wrong way round was silent and plausible.
  expect_true(isTRUE(sp$interaction1)); expect_true(isTRUE(sp$interaction2))
  expect_equal(sp$modEffectNames1, "transTrip")
  expect_equal(sp$modEffectNames2, "recip")
  expect_equal(sp$intEffectNames1, "transRecTrip")
  expect_equal(sp$intEffectNames2, "transRecTrip")
})

test_that("a NULL entry gives that effect its default perturbation", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  ## density defaults to a contrast, transTrip to a unit step -- the same
  ## defaults make_postest_targets() would have applied.
  lst <- suppressMessages(set_second_diff(tg, list(density = NULL,
                                                   transTrip = list())))
  sp <- RSiena:::.targetsToEffectList(lst)$effectList[[1L]]
  expect_equal(sp$contrast1, c(-1, 1))
  expect_null(sp$diff1)
  expect_equal(sp$diff2, 1)
  expect_null(sp$contrast2)
})

test_that("an effect may be crossed with itself in the list form", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  ## Repeated names are meaningful: the same base effect can enter through
  ## two different compound interactions.  Rejecting duplicates would make
  ## the list form less expressive than the numbered one it replaces.
  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  lst <- suppressMessages(set_second_diff(tg, list(
      recip = list(contrast = c(0, 1), interaction = TRUE,
                   intEffectNames = "transRecTrip",
                   modEffectNames = "transTrip"),
      recip = list(contrast = c(0, 1), interaction = TRUE,
                   intEffectNames = "transRecTrip",
                   modEffectNames = "density"))))
  sp <- RSiena:::.targetsToEffectList(lst)$effectList[[1L]]
  expect_equal(sp$effectName1, "recip")
  expect_equal(sp$effectName2, "recip")
  expect_equal(sp$modEffectNames1, "transTrip")
  expect_equal(sp$modEffectNames2, "density")
})

test_that("malformed perturbation lists are rejected with the effect named", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  sd <- function(...) suppressMessages(set_second_diff(tg, ...))

  ## A typo in a setting name must not be silently ignored -- that would
  ## compute a different quantity than the one written down.
  expect_error(sd(list(recip = list(diff = 1),
                       transTrip = list(constrast = c(0, 1)))),
               "Unknown perturbation setting")
  expect_error(sd(list(recip = list(diff = 1, contrast = c(0, 1)),
                       transTrip = NULL)),
               "either 'diff' or 'contrast' for 'recip'")
  expect_error(sd(list(recip = 1, transTrip = NULL)),
               "must be a list")
  expect_error(sd(list(list(diff = 1), list(diff = 1))),
               "must be named with an effect short name")
})

test_that("a second difference crosses exactly two effects", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  expect_error(suppressMessages(set_second_diff(tg, list(density = NULL))),
               "exactly two effects")
  ## Three entries are the natural way to ask for a third-order difference,
  ## which the machinery does not do yet; say so rather than silently using
  ## the first two.
  expect_error(suppressMessages(set_second_diff(tg,
                   list(density = NULL, recip = NULL, transTrip = NULL))),
               "not supported yet")
})

test_that("the removed numbered form is refused with a pointer, not ignored", {
  skip_on_cran()
  ans2 <- load_fixture("ans2"); mymodel2 <- load_fixture("mymodel2")
  skip_if(is.null(ans2) || is.null(mymodel2), "ans2 fixtures unavailable")

  tg <- make_postest_targets(ans2, effects = mymodel2, depvar = "mynet2",
                             includeDefaults = FALSE)
  ## Two bare effect names: these are effect names, not objects, so the
  ## default failure is "object 'transTrip' not found" -- true and useless.
  err <- tryCatch(suppressMessages(
             set_second_diff(tg, c(transTrip, recip))),
           error = conditionMessage)
  expect_match(err, "named list")
  expect_match(err, "numbered form")

  ## Numbered arguments now land in `...`, where being ignored in silence
  ## would build the target with default perturbations instead.
  expect_error(suppressMessages(set_second_diff(tg,
                   list(recip = NULL, transTrip = NULL), diff1 = 2)),
               "removed numbered form")
  expect_error(suppressMessages(set_second_diff(tg,
                   list(recip = NULL, transTrip = NULL),
                   contrast2 = c(0, 1), modEffectNames1 = "x")),
               "contrast2, modEffectNames1")
})

# ── K. Identifying WHICH target a short name means ───────────────────────────
#
# A short name is the natural way to name a target and is what set_effect()
# takes, but it does not always identify one: a model can carry egoX on
# several covariates, or two unspInt terms.  Picking the first match would
# perturb an effect the caller did not ask for and report it under the name
# of the one they did -- wrong numbers, not an error.  So an ambiguous short
# name is refused, and is resolved the way set_effect() resolves it: by
# naming the covariate.
#
# Internally the qualified ("long") name is what the engine is given, since
# that is what its own resolver can pin down; the readable short name stays
# for naming and printing.

test_that("an ambiguous short name is refused rather than resolved to the first match", {
  skip_on_cran()
  ans_2int <- load_fixture("ans_2int"); mymodel_2int <- load_fixture("mymodel_2int")
  skip_if(is.null(ans_2int) || is.null(mymodel_2int),
          "2int fixtures unavailable (RSENA_FULL_TESTS not set)")

  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int")
  ## Two unspInt terms: "unspInt" names both.
  expect_error(suppressMessages(set_target(tg, unspInt, diff = 2)),
               "does not identify one")
  ## and the message has to say what to type instead.
  err <- tryCatch(suppressMessages(set_target(tg, unspInt, diff = 2)),
                  error = conditionMessage)
  expect_match(err, "unspInt1")
  expect_match(err, "unspInt2")

  ## The qualified name does identify one.
  expect_no_error(suppressMessages(set_target(tg, unspInt1, diff = 2)))
  expect_no_error(suppressMessages(set_target(tg, "unspInt2", diff = 2)))
})

test_that("colliding short names still get distinct target names", {
  skip_on_cran()
  ans_2int <- load_fixture("ans_2int"); mymodel_2int <- load_fixture("mymodel_2int")
  skip_if(is.null(ans_2int) || is.null(mymodel_2int),
          "2int fixtures unavailable (RSENA_FULL_TESTS not set)")

  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int")
  ## Two targets called "unspInt_fd" would collide in the result list and in
  ## set_target's duplicate-name check.
  expect_false(anyDuplicated(tg$name) > 0L)
  expect_true(all(c("unspInt1_fd", "unspInt2_fd") %in% tg$name))
})

test_that("the engine is given the qualified name, the user sees the short one", {
  skip_on_cran()
  ans_ego <- load_fixture("ans_ego"); mymodel_ego <- load_fixture("mymodel_ego")
  skip_if(is.null(ans_ego) || is.null(mymodel_ego),
          "ego fixtures unavailable (RSENA_FULL_TESTS not set)")

  tg <- make_postest_targets(ans_ego, effects = mymodel_ego)
  lowered <- RSiena:::.targetsToEffectList(tg)$effectList
  nm <- vapply(lowered, function(e) e$effectName1, character(1L))
  ## Internal: the covariate-qualified name, which is what pins down which
  ## egoX is meant among the effects the model was fitted with.
  expect_true("egoX_mybeh_ego" %in% nm)
  ## User-facing: the deliberate short name, in the target name and on screen.
  expect_true("egoX_fd" %in% tg$name)
  expect_match(paste(capture.output(print(tg)), collapse = "\n"),
               "egoX \\(mybeh_ego\\)")
})

test_that("a covariate identifies which target a short name means", {
  skip_on_cran()
  ans_ego <- load_fixture("ans_ego"); mymodel_ego <- load_fixture("mymodel_ego")
  skip_if(is.null(ans_ego) || is.null(mymodel_ego),
          "ego fixtures unavailable (RSENA_FULL_TESTS not set)")

  tg <- make_postest_targets(ans_ego, effects = mymodel_ego)
  ## The set_effect() pattern: short name plus the covariate it is defined on.
  expect_no_error(suppressMessages(
      set_target(tg, egoX, covar1 = "mybeh_ego", diff = 2)))
  ## A covariate that is not there must not fall back to "the egoX I found".
  expect_error(suppressMessages(set_target(tg, egoX, covar1 = "nosuch", diff = 2)),
               "No target for effect 'egoX'")
  ## It identifies ONE target, so it cannot be handed a list of names.
  expect_error(suppressMessages(
      set_target(tg, c(egoX, recip), covar1 = "mybeh_ego")),
      "single short name")
  ## and it works inside a second difference's perturbation list.
  expect_no_error(suppressMessages(set_second_diff(tg,
      list(egoX  = list(covar1 = "mybeh_ego", diff = 1),
           recip = list(diff = 1)))))
})

test_that("a second difference records which target each component was", {
  skip_on_cran()
  ans_2int <- load_fixture("ans_2int"); mymodel_2int <- load_fixture("mymodel_2int")
  skip_if(is.null(ans_2int) || is.null(mymodel_2int),
          "2int fixtures unavailable (RSENA_FULL_TESTS not set)")

  tg <- make_postest_targets(ans_2int, effects = mymodel_2int,
                             depvar = "mynet_2int", includeDefaults = FALSE)
  ## An ambiguous component is caught when the target is added, not when it
  ## is computed several steps later.
  expect_error(suppressMessages(set_second_diff(tg,
                   list(unspInt = list(diff = 1), inPop = list(diff = 1)))),
               "does not identify one")

  sd <- suppressMessages(set_second_diff(tg,
            list(unspInt2 = list(diff = 1), inPop = list(diff = 1))))
  sp <- RSiena:::.targetsToEffectList(sd)$effectList[[1L]]
  expect_equal(sp$effectName1, "unspInt2")
  expect_equal(sp$effectName2, "inPop")
})
