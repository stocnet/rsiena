# test-postest-behaviour-baseline.R
#
# Characterization tests for the behaviours step 5d intends to CHANGE.
#
# Written before step 5c removes the flat interface, because that removal takes
# the equivalence harness's reference with it: after 5c there is nothing left
# to check a behaviour change against. These tests pin what the code does
# TODAY, on the object interface, so that when 5d changes it the diff to this
# file is the record of the change -- deliberate, attributable, and separate
# from the removal commit.
#
# So a failure here after 5d is not necessarily a bug: it is the change
# announcing itself. Each test says what it expects to become. What would be a
# bug is any of these changing during 5c, which touches none of them.
#
# Provenance: the step-5c coverage sweep. The corpus had 19 entries and could
# observe none of these conditions; three entries were added there, and this
# file pins what those entries would otherwise only compare.

me_base <- function(ans, dat, tg)
  suppressMessages(marginalEffects(ans, dat, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(verbose = FALSE)))


# ── 5d item 2: the override asymmetry ────────────────────────────────────────
#
#   eff_level       <- if ("level" %in% names(spec)) spec$level else level
#   eff_accumulated <- isTRUE(spec$accumulated) || accumulated
#
# level/condition are override-if-present, so a target can set them freely.
# accumulated/rateWeight are OR'd with the model-level value, so a target can
# switch them ON but never OFF. Almost certainly accidental; 5d makes all five
# override-if-present.

test_that("BASELINE: a per-target rateWeight = FALSE is ignored when the model sets TRUE", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  mk <- function(...) {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                               type = "tieProb", includeDefaults = FALSE, ...)
    suppressMessages(set_target(tg, transTrip, diff = 1))
  }
  on_model_only <- me_base(ans, mydata, mk(rateWeight = TRUE))
  with_override  <- me_base(ans, mydata, suppressMessages(
      set_target(mk(rateWeight = TRUE), transTrip, diff = 1,
                 rateWeight = FALSE)))
  off_entirely   <- me_base(ans, mydata, mk(rateWeight = FALSE))

  ## The setting demonstrably matters -- otherwise the test below would pass
  ## for the trivial reason that nothing depends on it.
  expect_false(isTRUE(all.equal(on_model_only, off_entirely)),
    info = "rateWeight must change the result for this test to mean anything")

  ## AFTER 5d this should be all.equal(with_override, off_entirely): the
  ## target asked for FALSE and should get it.
  expect_equal(with_override, on_model_only)
})

test_that("BASELINE: a per-target level or condition DOES override the model", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  mk <- function(...) {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                               type = "tieProb", includeDefaults = FALSE, ...)
    suppressMessages(set_target(tg, transTrip, diff = 1))
  }
  conditioned   <- me_base(ans, mydata, mk(condition = "recip"))
  turned_off    <- me_base(ans, mydata, suppressMessages(
      set_target(mk(condition = "recip"), transTrip, diff = 1,
                 condition = NULL)))
  never_set     <- me_base(ans, mydata, mk())

  ## An explicit NULL switches conditioning OFF for this target even though
  ## the model sets one. This is the distinction the .overrides list-column
  ## exists to preserve -- "not set" vs "set to NULL" -- and 5d keeps it.
  ## It is also the half of the asymmetry that already behaves correctly.
  expect_false(isTRUE(all.equal(conditioned, turned_off)))
  expect_equal(turned_off, never_set)
})


# ── 5d item 3: massContrasts resolution ──────────────────────────────────────
#
#   eff_massC <- if (!is.null(spec$massContrasts)) spec$massContrasts
#                else (pt1 == "ego") || (eff_second && pt2 == "ego")
#
# Resolved per spec with an auto-detect fallback. On the object interface the
# model-level value is written into every spec when lowering, so it applies;
# on the flat effectList path it is dropped (see the `diverges` entry in the
# corpus). Two consequences to pin.

test_that("BASELINE: model-level massContrasts reaches every target", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", massContrasts = TRUE,
                             includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
  out <- me_base(ans, mydata, tg)
  ## transTrip is not an ego perturbation, so auto-detection would say FALSE:
  ## these columns are here because the model-level setting was honoured.
  expect_true(all(c("massCreation", "massDissolution") %in% names(out)))
})

test_that("BASELINE: an explicit massContrasts = FALSE suppresses auto-detection", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  mk <- function(mc) {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                               type = "tieProb", includeDefaults = FALSE)
    suppressMessages(set_target(tg, transTrip, diff = 1,
                                perturbType = "ego", massContrasts = mc))
  }
  ## An ego perturbation auto-detects TRUE ...
  expect_true("massCreation" %in% names(me_base(ans, mydata, mk(NULL))))
  ## ... and an explicit FALSE overrides that. Defensible, but currently
  ## inherited rather than chosen; 5d decides it deliberately.
  expect_false("massCreation" %in% names(me_base(ans, mydata, mk(FALSE))))
})


# ── 5d item 4: the uncertainty default ───────────────────────────────────────

test_that("BASELINE: uncertainty is computed unless switched off", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))

  ## The redesign leans towards enabled = FALSE, on the grounds that the
  ## expensive thing should be asked for. Pinned so the flip is visible.
  expect_true(formals(set_postest_uncertainty_saom)$enabled)
  out <- suppressMessages(marginalEffects(ans, mydata, targets = tg,
             control_algo = set_postest_algo_saom(verbose = FALSE)))
  expect_true(any(grepl("^se|Se$|SE", names(out))),
    info = "default should currently produce standard errors")
})


# ── The `details` argument: accepted, cannot be honoured ──────────────────────

test_that("BASELINE: details = TRUE still errors, and is not reachable from control_out", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  ## Pre-existing, verified against HEAD before any of the refactor: this is
  ## not a regression. It went unnoticed because no test covered it.
  ## Pinned so that "fix it" and "remove it" are both deliberate acts, and so
  ## the argument cannot quietly go on accepting a value it cannot honour.
  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))
  expect_error(suppressMessages(marginalEffects(ans, mydata, targets = tg,
      details = TRUE,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(verbose = FALSE))))

  ## Withdrawn from the output constructor while it is broken, so that the
  ## object interface does not advertise it.
  expect_false("details" %in% names(formals(set_postest_output_saom)))
})

test_that("BASELINE: a per-target accumulated = FALSE is ignored when the model sets TRUE", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RSENA_FULL_TESTS"), "1"),
              "dynamic simulation; RSENA_FULL_TESTS not set")
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel"); mycontrols <- load_fixture("mycontrols")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel) ||
          is.null(mycontrols), "dynamic fixtures unavailable")

  ## accumulated and rateWeight are OR'd on two SEPARATE lines. Testing only
  ## rateWeight would let a 5d fix correct one and miss the other, so this
  ## pays for a dynamic run to cover the line users actually set.
  mk <- function(...) {
    tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                               type = "tieProb", dynamic = TRUE,
                               includeDefaults = FALSE, ...)
    suppressMessages(set_target(tg, transTrip, diff = 1))
  }
  run <- function(tg) suppressMessages(marginalEffects(ans, mydata, targets = tg,
      control_uncertainty = set_postest_uncertainty_saom(enabled = FALSE),
      control_algo = set_postest_algo_saom(algorithm = mycontrols, n3 = 20,
                                           verbose = FALSE)))

  set.seed(99L); on_model_only <- run(mk(accumulated = TRUE))
  set.seed(99L); with_override <- run(suppressMessages(
      set_target(mk(accumulated = TRUE), transTrip, diff = 1,
                 accumulated = FALSE)))
  set.seed(99L); off_entirely  <- run(mk(accumulated = FALSE))

  expect_false(isTRUE(all.equal(on_model_only, off_entirely)),
    info = "accumulated must change the result for this test to mean anything")
  ## AFTER 5d this should equal off_entirely.
  expect_equal(with_override, on_model_only)
})
