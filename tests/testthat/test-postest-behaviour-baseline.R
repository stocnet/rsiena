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
# CHANGED IN 5d. level/condition were override-if-present; accumulated and
# rateWeight were OR'd with the model-level value, so a target could switch
# them ON but never OFF. All four are override-if-present now, and the tests
# below assert the new behaviour -- this diff IS the record of the change.

test_that("a per-target rateWeight = FALSE overrides a model that sets TRUE", {
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

  ## The target asked for FALSE and gets it (5d).  Before 5d this equalled
  ## on_model_only instead: the request was silently dropped.
  expect_equal(with_override, off_entirely)
})

test_that("a per-target level or condition overrides the model", {
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

test_that("uncertainty is computed by default, analytically", {
  skip_on_cran()
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel),
          "base fixtures unavailable")

  tg <- make_postest_targets(ans, effects = mymodel, depvar = "mynet",
                             type = "tieProb", includeDefaults = FALSE)
  tg <- suppressMessages(set_target(tg, transTrip, diff = 1))

  ## CHANGED IN 5d, but not the way the redesign proposed. It leaned towards
  ## enabled = FALSE, on the grounds that the expensive thing should be asked
  ## for. Measurement moved the answer: the expense was the bootstrap DEFAULT
  ## (~114 s vs 0.57 s for delta on a toy model), not uncertainty itself. So
  ## uncertainty stays on -- a marginal effect without a standard error is
  ## not reportable -- and the default mode became the analytic one.
  expect_true(formals(set_postest_uncertainty_saom)$enabled)
  expect_equal(set_postest_uncertainty_saom()$mode, "delta")
  out <- suppressMessages(marginalEffects(ans, mydata, targets = tg,
             control_algo = set_postest_algo_saom(verbose = FALSE)))
  ## The SE column is called "SE" whichever method produced it (step 5d).
  ## It used to be named for the method, so the default output changed shape
  ## when the default mode changed. Since the modes are mutually exclusive,
  ## the name carried nothing the metadata does not.
  expect_true(any(grepl("(^|_)([Ss][Ee])$", names(out))),
    info = paste("default must produce standard errors; got:",
                 paste(names(out), collapse = ", ")))
  expect_true("SE" %in% names(out))
  ## ...and how it was derived is recorded, since the column name no longer
  ## says. Without this the number is not interpretable: the delta modes
  ## differ in what they hold fixed.
  expect_equal(attr(out, "uncertaintyMethod"), "delta")
})


# ── The `details` capability: parked, not abandoned ───────────────────────────
#
# WHAT IT IS. `details = TRUE` returns the intermediate quantities behind a
# marginal effect rather than just the effect: `utilDiff` (the utility shift
# the perturbation produces), `oldChangeProb` / `newChangeProb` (the choice
# probability before and after it), and `newTieProb` for tieProb targets.
# That is the arithmetic of the difference laid open -- what you want when a
# number looks wrong and you need to see which step produced it, and what a
# teaching or diagnostic use would ask for.
#
# This is NOT the same as `returnDecisionDetails`, which is per-target, works,
# and is covered elsewhere. Decision-level information is not at risk here.
#
# WHY IT IS PARKED. It errors inside `encodeGroupKeys()` -- a defect that
# predates the interface refactor and was never covered by a test, which is
# how it survived. Step 5c removed the flat argument list, so `details` is now
# unreachable from every entry point: the constructors do not offer it and the
# call errors with "unused argument". The internal plumbing is still in place
# and pinned FALSE.
#
# WHEN IT GETS DECIDED. At step 6b (aggregation unification, staged cache),
# which rewrites `encodeGroupKeys` and the aggregation path the bug lives in.
# Deciding earlier means either fixing code about to be rewritten or deleting
# a capability whose replacement is being designed. Fix it or remove it there,
# deliberately -- do not let it drift on as plumbing nobody can reach.

test_that("BASELINE: details = TRUE is unreachable (parked, see note above)", {
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

test_that("a per-target accumulated = FALSE overrides a model that sets TRUE", {
  skip_on_cran()
  skip_if_not(identical(Sys.getenv("RSENA_FULL_TESTS"), "1"),
              "dynamic simulation; RSENA_FULL_TESTS not set")
  ans <- load_fixture("ans"); mydata <- load_fixture("mydata")
  mymodel <- load_fixture("mymodel"); mycontrols <- load_fixture("mycontrols")
  skip_if(is.null(ans) || is.null(mydata) || is.null(mymodel) ||
          is.null(mycontrols), "dynamic fixtures unavailable")

  ## accumulated and rateWeight were OR'd on two SEPARATE lines. Testing only
  ## rateWeight would have let the 5d fix correct one and miss the other, so
  ## this pays for a dynamic run to cover the line users actually set.
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
  expect_equal(with_override, off_entirely)
})
