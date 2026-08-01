# helper-equivalence.R
#
# Infrastructure for the API-equivalence harness (postestimate_api_redesign.md
# Sec. 6, step 2).
#
# Purpose: when the set_postest_*_saom() config objects land, every call
# expressible in the old flat form must produce a bit-identical result in the
# new object form.  Both forms are run in the SAME session and compared
# directly — no stored golden values, so the harness is immune to fixture drift
# and needs no regeneration when a model is refitted.
#
# This file provides:
#   compare_me_output()  — rigorous comparison of two marginalEffects results
#   me_corpus()          — the corpus of representative calls to compare
#
# The comparison function is deliberately strict.  A harness that reports
# "equivalent" while only checking column names and row counts is worse than no
# harness: it converts an unverified refactor into a falsely verified one.
# test-equivalence-harness.R self-tests it against deliberately corrupted input
# to prove it is not vacuous.


# --------------------------------------------------------------------------
# compare_me_output(a, b, tolerance)
#
# Compares two marginalEffects / predict outputs for practical identity.
# Returns TRUE when equivalent, otherwise a character vector of differences
# (so a failing test can report every problem, not just the first).
#
# Checks, in order:
#   - both are data.frames
#   - identical column names, in the same order
#   - identical row count
#   - per column: same type, and same values (numeric within `tolerance`,
#     everything else exactly), with NA patterns required to match
#
# Column order matters: a reordering is a user-visible change to output layout
# even when the values are unchanged.  Row order matters for the same reason.
# --------------------------------------------------------------------------
compare_me_output <- function(a, b, tolerance = 1e-10,
                              label_a = "flat", label_b = "object") {
  problems <- character(0)
  add <- function(...) problems <<- c(problems, paste0(...))

  ## Multi-effect calls return a NAMED LIST of data.frames (one per effect)
  ## rather than a single frame.  Recurse element-wise so those are comparable
  ## too — without this, multi-effect calls could not be verified at all.
  if (!is.data.frame(a) && is.list(a) && !is.data.frame(b) && is.list(b)) {
    if (!identical(sort(names(a)), sort(names(b)))) {
      add("result list names differ: ", label_a, "=[",
          paste(names(a), collapse = ","), "] ", label_b, "=[",
          paste(names(b), collapse = ","), "]")
      return(problems)
    }
    for (nm in names(a)) {
      sub <- compare_me_output(a[[nm]], b[[nm]], tolerance = tolerance,
                               label_a = label_a, label_b = label_b)
      if (!isTRUE(sub)) add("[", nm, "] ", sub)
    }
    return(if (length(problems) == 0L) TRUE else problems)
  }

  if (!is.data.frame(a) || !is.data.frame(b)) {
    add("not both data.frames: ", label_a, "=", class(a)[1],
        ", ", label_b, "=", class(b)[1])
    return(problems)
  }

  na <- names(a); nb <- names(b)
  if (!identical(na, nb)) {
    only_a <- setdiff(na, nb); only_b <- setdiff(nb, na)
    if (length(only_a)) add("columns only in ", label_a, ": ",
                            paste(only_a, collapse = ", "))
    if (length(only_b)) add("columns only in ", label_b, ": ",
                            paste(only_b, collapse = ", "))
    if (!length(only_a) && !length(only_b))
      add("same columns but different order: ", label_a, "=[",
          paste(na, collapse = ","), "] ", label_b, "=[",
          paste(nb, collapse = ","), "]")
    return(problems)   # column mismatch makes per-column comparison meaningless
  }

  if (nrow(a) != nrow(b)) {
    add("row count differs: ", label_a, "=", nrow(a), ", ",
        label_b, "=", nrow(b))
    return(problems)
  }

  for (cn in na) {
    va <- a[[cn]]; vb <- b[[cn]]

    if (!identical(is.na(va), is.na(vb))) {
      add("column '", cn, "': NA pattern differs (",
          sum(is.na(va)), " vs ", sum(is.na(vb)), " NAs)")
      next
    }

    if (is.numeric(va) && is.numeric(vb)) {
      ka <- va[!is.na(va)]; kb <- vb[!is.na(vb)]
      if (length(ka)) {
        d <- max(abs(ka - kb))
        if (!is.finite(d) || d > tolerance)
          add("column '", cn, "': max abs difference ",
              format(d, digits = 3), " exceeds tolerance ",
              format(tolerance, digits = 3))
      }
    } else if (!identical(va, vb)) {
      nd <- sum(va != vb, na.rm = TRUE)
      add("column '", cn, "': ", nd, " of ", length(va),
          " values differ (non-numeric)")
    }
  }

  if (length(problems) == 0L) TRUE else problems
}


# --------------------------------------------------------------------------
# me_corpus()
#
# Representative marginalEffects calls covering the parameter surface that
# moves into the config objects.  Each entry:
#
#   name  — identifier used in test output
#   args  — named list of arguments to do.call(marginalEffects, .)
#   slow  — TRUE if it needs RSENA_FULL_TESTS=1
#   needs — fixture names required, so a missing fixture skips rather than errors
#
# When the constructors land, each entry gains an `obj_args` field expressing
# the same call in config-object form; the harness then asserts
# compare_me_output(do.call(marginalEffects, args),
#                   do.call(marginalEffects, obj_args)).
#
# Coverage intent: every argument here is one that the refactor relocates.
# A gap in this corpus is an argument that sails through the refactor untested.
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# me_corpus(fixtures) / me_call(entry, fixtures)
#
# A corpus of configurations, expressed as DATA rather than as calls, with
# me_call() turning an entry into the arguments marginalEffects() takes.
#
# It was built for the flat-vs-object equivalence check (Part C), which step 5c
# retired along with the flat form.  It is kept, converted to the object
# interface, for two reasons: Part B still checks that every one of these
# configurations runs, and step 6 needs a set of configurations to compare
# before and after the analytic Jacobian replaces the finite-difference
# fallback -- which is also why compare_me_output() is kept.
#
# Entry fields:
#   name, needs, slow  -- identity, required fixtures, whether to skip when fast
#   fit/data/effects   -- fixture NAMES, resolved against `fixtures`
#   depvar             -- dependent variable
#   model              -- extra make_postest_targets() arguments
#   targets            -- named list: effect short name -> set_target() args
#   sd                 -- optional list(name=, perturb=) for set_second_diff()
#   unc / out          -- set_postest_uncertainty_saom() / _output_saom() args
#
# Coverage intent: every setting the refactor relocated appears in at least one
# entry. A gap here is a setting that sails through untested -- and mutation
# testing has twice shown that entry COUNT is not coverage, so several entries
# below carry a note saying what makes their setting observable at all.
# --------------------------------------------------------------------------
me_corpus_fixtures <- function() {
  nms <- c("ans", "mydata", "mymodel", "ans2", "mydata2", "mymodel2",
           "mycontrols", "ans_2int", "mydata_2int", "mymodel_2int",
           "ans_cm", "mydata_cm", "mymodel_cm")
  setNames(lapply(nms, load_fixture), nms)
}


me_call <- function(entry, fixtures) {
  f   <- fixtures
  fit <- f[[entry$fit]]
  tg  <- do.call(make_postest_targets,
                 c(list(x = fit, effects = f[[entry$effects]],
                        depvar = entry$depvar, includeDefaults = FALSE),
                   entry$model %||% list()))
  for (nm in names(entry$targets %||% list()))
    tg <- suppressMessages(do.call(set_target,
            c(list(x = tg, shortNames = nm), entry$targets[[nm]])))
  if (!is.null(entry$sd))
    tg <- suppressMessages(do.call(set_second_diff,
            c(list(x = tg, perturbations = entry$sd$perturb),
              entry$sd[setdiff(names(entry$sd), "perturb")])))

  list(object = fit, data = f[[entry$data]], targets = tg,
       control_uncertainty = do.call(set_postest_uncertainty_saom,
                                     entry$unc %||% list()),
       control_algo = do.call(set_postest_algo_saom,
                              c(list(verbose = FALSE), entry$algo %||% list())),
       control_out = do.call(set_postest_output_saom, entry$out %||% list()))
}

`%||%` <- function(a, b) if (is.null(a)) b else a

me_corpus <- function(fixtures) {
  base1 <- list(fit = "ans",  data = "mydata",  effects = "mymodel",
                depvar = "mynet",  needs = c("ans", "mydata", "mymodel"))
  base2 <- list(fit = "ans2", data = "mydata2", effects = "mymodel2",
                depvar = "mynet2", needs = c("ans2", "mydata2", "mymodel2"))
  no_unc <- list(enabled = FALSE)

  entries <- list(
    c(base1, list(name = "static_firstDiff_period", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1)), unc = no_unc)),

    c(base1, list(name = "static_firstDiff_egoChoice_conditioned", slow = FALSE,
        model = list(type = "tieProb", level = "egoChoice",
                     condition = c("recip", "density")),
        targets = list(transTrip = list(diff = 1)), unc = no_unc)),

    c(base1, list(name = "static_firstDiff_contrast", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(recip = list(contrast = c(0, 1))), unc = no_unc)),

    c(base1, list(name = "static_density_contrast", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(density = list(contrast = c(-1, 1))), unc = no_unc)),

    c(base1, list(name = "static_changeProb", slow = FALSE,
        model = list(type = "changeProb", level = "period"),
        targets = list(transTrip = list(diff = 1)), unc = no_unc)),

    c(base1, list(name = "static_riskRatio", slow = FALSE,
        model = list(type = "tieProb", level = "period",
                     mainEffect = "riskRatio"),
        targets = list(transTrip = list(diff = 1)), unc = no_unc)),

    c(base1, list(name = "static_egoNormalize_FALSE", slow = FALSE,
        model = list(type = "tieProb", level = "period",
                     egoNormalize = FALSE),
        targets = list(transTrip = list(diff = 1)), unc = no_unc)),

    c(base1, list(name = "static_secondDiff", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        sd = list(perturb = list(transTrip = list(diff = 1),
                                 recip = list(contrast = c(0, 1)))),
        unc = no_unc)),

    c(base2, list(name = "static_interaction", slow = FALSE,
        model = list(type = "tieProb", level = "period",
                     condition = "recip"),
        targets = list(transTrip = list(diff = 1, interaction = TRUE,
                                        intEffectNames = "transRecTrip",
                                        modEffectNames = "recip")),
        unc = no_unc)),

    ## `format` is consumed ONLY inside combinePostestResults(), reached only
    ## when combineSameLevel = TRUE AND there are several targets.  Miss either
    ## condition and the setting is inert -- verified by mutation testing,
    ## which the earlier single-target version of this entry could not detect.
    c(base1, list(name = "static_format_long_combined", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1, name = "tt"),
                       recip = list(contrast = c(0, 1), name = "rp")),
        out = list(combineSameLevel = TRUE, format = "long"),
        unc = list(mode = "delta", sd = TRUE, ci = TRUE, nsim = 1L))),

    c(base1, list(name = "static_multi_targets", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1, name = "tt"),
                       recip = list(contrast = c(0, 1), name = "rp")),
        unc = no_unc)),

    c(base1, list(name = "static_delta_uncertainty", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1)),
        unc = list(mode = "delta", sd = TRUE, ci = TRUE, nsim = 1L))),

    ## The TWO-TIER default mechanism: the model sets level = "period", one
    ## target overrides it to "ego".  Without this, dropping the per-target
    ## override entirely goes undetected -- verified by mutation testing.
    c(base1, list(name = "static_per_target_level_override", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1, name = "tt"),
                       recip = list(contrast = c(0, 1), name = "rp",
                                    level = "ego")),
        unc = no_unc)),

    ## An explicit NULL must switch conditioning OFF for one target while the
    ## model sets one -- the distinction .overrides exists to keep.
    c(base1, list(name = "static_per_target_condition_override", slow = FALSE,
        model = list(type = "tieProb", level = "period",
                     condition = "recip"),
        targets = list(transTrip = list(diff = 1, name = "tt",
                                        condition = NULL),
                       recip = list(contrast = c(0, 1), name = "rp",
                                    condition = "density")),
        unc = no_unc)),

    ## The override ASYMMETRY: accumulated/rateWeight are OR'd with the model
    ## value, so a target can switch them on but never off.  Step 5d makes all
    ## five override-if-present; until then this pins that it does not.
    c(base1, list(name = "static_per_target_rateWeight_override", slow = FALSE,
        model = list(type = "tieProb", level = "period", rateWeight = TRUE),
        targets = list(transTrip = list(diff = 1, name = "tt",
                                        rateWeight = FALSE),
                       recip = list(contrast = c(0, 1), name = "rp")),
        unc = no_unc)),

    ## massContrasts adds massCreation/massDissolution; no other entry sets it,
    ## so dropping it from the defaults was invisible.
    c(base1, list(name = "static_massContrasts", slow = FALSE,
        model = list(type = "tieProb", level = "period",
                     massContrasts = TRUE),
        targets = list(recip = list(contrast = c(0, 1))), unc = no_unc)),

    ## Model-level massContrasts must reach a target that does not set it.
    ## This is where the flat path diverged from the object path before 5c:
    ## it dropped the call-level value and fell back to auto-detection.  The
    ## flat path is gone; the correct behaviour is pinned in
    ## test-postest-behaviour-baseline.R.
    c(base1, list(name = "static_massContrasts_multi", slow = FALSE,
        model = list(type = "tieProb", level = "period",
                     massContrasts = TRUE),
        targets = list(transTrip = list(diff = 1, name = "tt"),
                       recip = list(contrast = c(0, 1), name = "rp")),
        unc = no_unc)),

    ## Bootstrap is the ONLY mode where nsim affects output (delta modes are
    ## analytic and draw nothing).  Without this a wiring bug that drops nsim
    ## passes unnoticed -- verified by mutation testing.
    c(base1, list(name = "static_bootstrap_uncertainty", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1)),
        unc = list(mode = "bootstrap", nsim = 12L, sd = TRUE, ci = TRUE))),

    ## perturbType forces the ego-wide mlogit update instead of the default.
    c(base1, list(name = "static_perturbType_ego", slow = FALSE,
        model = list(type = "tieProb", level = "period"),
        targets = list(transTrip = list(diff = 1, perturbType = "ego")),
        unc = no_unc)),

    ## Compound effects on BOTH sides of a second difference, with different
    ## components -- the most demanding perturbation the interface can express.
    list(name = "static_two_unspInt_secondDiff",
         fit = "ans_2int", data = "mydata_2int", effects = "mymodel_2int",
         depvar = "mynet_2int",
         needs = c("ans_2int", "mydata_2int", "mymodel_2int"), slow = TRUE,
         model = list(type = "tieProb", level = "period",
                      condition = c("inPop", "outPop", "density")),
         sd = list(name = "rr", perturb = list(
             recip = list(contrast = c(0, 1), interaction = TRUE,
                          intEffectNames = "unspInt1",
                          modEffectNames = "inPop"),
             recip = list(contrast = c(0, 1), interaction = TRUE,
                          intEffectNames = "unspInt2",
                          modEffectNames = "outPop"))),
         unc = no_unc),

    ## Creation / endowment.  A target collapses eval+creation (or eval+endow)
    ## into ONE row, but the perturbation has to reach BOTH columns -- and the
    ## creation column is active only on ministeps that create a tie, the
    ## endowment column only on ones that drop it.  So the perturbation is
    ## row-dependent and spans several columns, which is why the analytic
    ## Jacobian declines these and finite differences carry them.
    ##
    ## The corpus had no creation/endow entry at all, so this case was not
    ## merely uncovered -- it was invisible to any measurement taken from it.
    list(name = "static_creation_firstDiff",
         fit = "ans_cm", data = "mydata_cm", effects = "mymodel_cm",
         depvar = "mynet_cm",
         needs = c("ans_cm", "mydata_cm", "mymodel_cm"), slow = FALSE,
         model = list(type = "tieProb", level = "period"),
         targets = list(recip = list(contrast = c(0, 1))),
         unc = no_unc),

    list(name = "static_endow_firstDiff",
         fit = "ans_cm", data = "mydata_cm", effects = "mymodel_cm",
         depvar = "mynet_cm",
         needs = c("ans_cm", "mydata_cm", "mymodel_cm"), slow = FALSE,
         model = list(type = "tieProb", level = "period"),
         targets = list(transTrip = list(diff = 1)),
         unc = no_unc),

    ## Both at once, and crossed: the second difference has to keep each
    ## component's eval/creation and eval/endow pairing straight.
    list(name = "static_creation_endow_secondDiff",
         fit = "ans_cm", data = "mydata_cm", effects = "mymodel_cm",
         depvar = "mynet_cm",
         needs = c("ans_cm", "mydata_cm", "mymodel_cm"), slow = FALSE,
         model = list(type = "tieProb", level = "period"),
         sd = list(perturb = list(recip     = list(contrast = c(0, 1)),
                                  transTrip = list(diff = 1))),
         unc = no_unc),

    c(base1, list(name = "dynamic_firstDiff", slow = TRUE,
        needs = c("ans", "mydata", "mymodel", "mycontrols"),
        model = list(type = "tieProb", level = "period", dynamic = TRUE),
        targets = list(transTrip = list(diff = 1)),
        algo = list(algorithm = fixtures$mycontrols, n3 = 20),
        unc = no_unc)),

    c(base1, list(name = "dynamic_accumulated_deltaFull", slow = TRUE,
        needs = c("ans", "mydata", "mymodel", "mycontrols"),
        model = list(type = "tieProb", level = "period", dynamic = TRUE,
                     accumulated = TRUE),
        targets = list(transTrip = list(diff = 1)),
        algo = list(algorithm = fixtures$mycontrols, n3 = 20),
        unc = list(mode = "deltaFull")))
  )

  names(entries) <- vapply(entries, `[[`, character(1), "name")
  entries
}


# --------------------------------------------------------------------------
# split_flat_args(), as_object_args() and as_targets() were REMOVED in step 5c.
#
# They translated a flat corpus entry into the equivalent config-object call so
# Part C could compare the two.  With the flat form gone there is nothing to
# translate from, and they had no other caller.
#
# compare_me_output() above is kept deliberately: steps 6 and 6b need
# before/after comparison of internals when the analytic Jacobian replaces the
# finite-difference fallback.
# --------------------------------------------------------------------------
