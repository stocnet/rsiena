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
me_corpus <- function(fixtures) {
  f <- fixtures
  base1 <- list(object = f$ans,  data = f$mydata,  depvar = "mynet")
  base2 <- list(object = f$ans2, data = f$mydata2, depvar = "mynet2")

  entries <- list(
    list(name = "static_firstDiff_period",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "tieProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_firstDiff_egoChoice_conditioned",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "tieProb", level = "egoChoice",
                              condition = c("recip", "density"),
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_firstDiff_contrast",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "recip", contrast1 = c(0, 1),
                              type = "tieProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_density_contrast",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "density", contrast1 = c(-1, 1),
                              type = "tieProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_changeProb",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "changeProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_riskRatio",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "tieProb", level = "period",
                              mainEffect = "riskRatio",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_egoNormalize_FALSE",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "tieProb", level = "period",
                              egoNormalize = FALSE,
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_secondDiff",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              second = TRUE,
                              effectName2 = "recip", contrast2 = c(0, 1),
                              type = "tieProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_interaction",
         needs = c("ans2", "mydata2"), slow = FALSE,
         args = c(base2, list(effectName1 = "transTrip", diff1 = 1,
                              interaction1 = TRUE,
                              intEffectNames1 = "transRecTrip",
                              modEffectNames1 = "recip",
                              type = "tieProb", level = "period",
                              condition = "recip",
                              uncertainty = FALSE, verbose = FALSE))),

    ## `format` is consumed ONLY inside combinePostestResults(), which is
    ## reached only when combineSameLevel = TRUE AND there are several effects
    ## (`.single_effect` is FALSE).  Miss either condition and the setting is
    ## inert -- verified by mutation testing, which the earlier single-effect,
    ## no-uncertainty version of this entry could not detect.
    list(name = "static_format_long_combined",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(
           effectList = list(
             tt = list(effectName1 = "transTrip", diff1 = 1),
             rp = list(effectName1 = "recip", contrast1 = c(0, 1))),
           type = "tieProb", level = "period",
           combineSameLevel = TRUE, format = "long",
           uncertainty = TRUE, uncertaintyMode = "delta",
           uncertaintySd = TRUE, uncertaintyCi = TRUE,
           nsim = 1L, verbose = FALSE))),

    list(name = "static_multi_effectList",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(
           effectList = list(
             tt = list(effectName1 = "transTrip", diff1 = 1),
             rp = list(effectName1 = "recip", contrast1 = c(0, 1))),
           type = "tieProb", level = "period",
           uncertainty = FALSE, verbose = FALSE))),

    list(name = "static_delta_uncertainty",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "tieProb", level = "period",
                              uncertainty = TRUE, uncertaintyMode = "delta",
                              uncertaintySd = TRUE, uncertaintyCi = TRUE,
                              nsim = 1L, verbose = FALSE))),

    ## Exercises the TWO-TIER default mechanism: the call sets level="period",
    ## one target overrides it to "ego".  Without this, dropping the per-target
    ## override entirely goes undetected -- verified by mutation testing.
    list(name = "static_per_target_level_override",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(
           effectList = list(
             tt = list(effectName1 = "transTrip", diff1 = 1),
             rp = list(effectName1 = "recip", contrast1 = c(0, 1),
                       level = "ego")),
           type = "tieProb", level = "period",
           uncertainty = FALSE, verbose = FALSE))),

    ## massContrasts adds massCreation/massDissolution columns; no other entry
    ## sets it, so dropping it from the defaults was invisible.
    list(name = "static_massContrasts",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "recip", contrast1 = c(0, 1),
                              type = "tieProb", level = "period",
                              massContrasts = TRUE,
                              uncertainty = FALSE, verbose = FALSE))),

    ## Bootstrap is the ONLY mode where nsim affects output (delta modes are
    ## analytic and draw nothing).  Without this entry a wiring bug that drops
    ## nsim passes the harness unnoticed — verified by mutation testing.
    list(name = "static_bootstrap_uncertainty",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              type = "tieProb", level = "period",
                              uncertainty = TRUE,
                              uncertaintyMode = "bootstrap",
                              nsim = 12L, uncertaintySd = TRUE,
                              uncertaintyCi = TRUE, verbose = FALSE))),

    ## Compound effects on BOTH sides of a second difference, with different
    ## components -- the most demanding thing the flat effectList could express.
    list(name = "static_two_unspInt_secondDiff",
         needs = c("ans_2int", "mydata_2int", "mymodel_2int"), slow = TRUE,
         args = list(object = f$ans_2int, data = f$mydata_2int,
                     effects = f$mymodel_2int, depvar = "mynet_2int",
                     effectList = list(rr = list(
                       effectName1 = "recip", contrast1 = c(0, 1),
                       interaction1 = TRUE, intEffectNames1 = "unspInt1",
                       modEffectNames1 = "inPop",
                       second = TRUE,
                       effectName2 = "recip", contrast2 = c(0, 1),
                       interaction2 = TRUE, intEffectNames2 = "unspInt2",
                       modEffectNames2 = "outPop")),
                     type = "tieProb", level = "period",
                     condition = c("inPop", "outPop", "density"),
                     uncertainty = FALSE, verbose = FALSE)),

    ## perturbType forces the ego-wide mlogit update instead of the default.
    list(name = "static_perturbType_ego",
         needs = c("ans", "mydata"), slow = FALSE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              perturbType1 = "ego",
                              type = "tieProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "dynamic_firstDiff",
         needs = c("ans", "mydata", "mymodel", "mycontrols"), slow = TRUE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              effects = f$mymodel, algorithm = f$mycontrols,
                              dynamic = TRUE, n3 = 20,
                              type = "tieProb", level = "period",
                              uncertainty = FALSE, verbose = FALSE))),

    list(name = "dynamic_accumulated_deltaFull",
         needs = c("ans", "mydata", "mymodel", "mycontrols"), slow = TRUE,
         args = c(base1, list(effectName1 = "transTrip", diff1 = 1,
                              effects = f$mymodel, algorithm = f$mycontrols,
                              dynamic = TRUE, n3 = 20,
                              type = "tieProb", level = "period",
                              accumulated = TRUE,
                              uncertainty = TRUE,
                              uncertaintyMode = "deltaFull",
                              verbose = FALSE)))
  )

  names(entries) <- vapply(entries, `[[`, character(1), "name")
  entries
}


# Load every fixture the corpus can reference; missing ones stay NULL and the
# corresponding entries skip rather than error.
me_corpus_fixtures <- function() {
  nms <- c("ans", "mydata", "ans2", "mydata2", "mymodel", "mycontrols",
           "ans_2int", "mydata_2int", "mymodel_2int")
  setNames(lapply(nms, load_fixture), nms)
}


# --------------------------------------------------------------------------
# split_flat_args(args)
#
# Splits a flat marginalEffects() argument list into the three domains the
# config-object refactor separates:
#
#   $uncertainty — arguments that move into set_postest_uncertainty_saom()
#   $control     — arguments that move into set_postest_algo_saom()
#   $rest        — everything else (effect spec + output shape), untouched
#                  by step 3 and still passed flat
#
# Used by the equivalence harness to mechanically derive the object-form call
# from the flat one, so the corpus does not have to spell out both.
# --------------------------------------------------------------------------
.unc_flat_map <- c(
  enabled    = "uncertainty",
  mode       = "uncertaintyMode",
  nsim       = "nsim",
  sd         = "uncertaintySd",
  ci         = "uncertaintyCi",
  ciInterval = "ciInterval",
  simMean    = "uncertaintyMean",
  simMedian  = "uncertaintyMedian"
)

.ctl_flat_names <- c(
  "algorithm", "n3", "n3PointEst", "n3BatchSize",
  "chainStoreMode", "useChangeContributions", "chainStorePath",
  "useCluster", "nbrNodes", "clusterType", "cl",
  "batchDir", "prefix", "combineBatch", "batchSize", "keepBatch",
  "verbose", "memoryScale", "batchUnitBudget",
  "dynamicMinistepFactor", "saveDir", "gcEachBatch", "gcEachSim"
)

split_flat_args <- function(args) {
  unc_flat <- intersect(unname(.unc_flat_map), names(args))
  ctl_flat <- intersect(.ctl_flat_names, names(args))
  list(
    uncertainty = args[unc_flat],
    control     = args[ctl_flat],
    rest        = args[setdiff(names(args), c(unc_flat, ctl_flat))]
  )
}


# --------------------------------------------------------------------------
# as_object_args(args)
#
# Rewrites a flat argument list into the step-3 object form: uncertainty and
# control arguments are folded into their constructors, everything else is
# passed through unchanged.
#
# Only supplies constructor arguments the caller actually set, so unspecified
# ones fall back to the constructor defaults — which is precisely what the
# equivalence test needs to verify (constructor defaults must reproduce the
# flat defaults exactly).
#
# Returns NULL when the constructors are not available yet, so the harness
# skips rather than errors before step 3 lands.
# --------------------------------------------------------------------------
as_object_args <- function(args) {
  if (!exists("set_postest_uncertainty_saom", mode = "function") ||
      !exists("set_postest_algo_saom",     mode = "function"))
    return(NULL)

  sp <- split_flat_args(args)

  unc_call <- list()
  for (obj_nm in names(.unc_flat_map)) {
    flat_nm <- .unc_flat_map[[obj_nm]]
    if (flat_nm %in% names(sp$uncertainty))
      unc_call[[obj_nm]] <- sp$uncertainty[[flat_nm]]
  }

  out <- sp$rest
  out$control_uncertainty <- do.call(set_postest_uncertainty_saom, unc_call)
  if (length(sp$control) > 0L)
    out$control_algo <- do.call(set_postest_algo_saom, sp$control)

  ## Output domain.
  out_flat <- intersect(.out_flat_names, names(args))
  if (length(out_flat) > 0L) {
    out <- out[setdiff(names(out), out_flat)]
    out$control_out <- do.call(set_postest_output_saom, args[out_flat])
  }

  ## Model domain: fold the effect specification onto a targets object and
  ## drop the flat equivalents.
  tg <- as_targets(args)
  if (!is.null(tg)) {
    out <- out[setdiff(names(out),
                       c(.model_defaults_flat, .model_spec_flat))]
    out$targets <- tg
  }
  out
}


## --------------------------------------------------------------------------
## Model-domain arguments: those that move onto the targets object.
## --------------------------------------------------------------------------
.model_defaults_flat <- c("dynamic", "type", "mainEffect", "level", "condition",
                          "egoNormalize", "accumulated", "rateWeight",
                          "massContrasts")

## Shape of the returned R object -> set_postest_output_saom().
## `details` is not exposed on the object interface -- broken feature, see
## set_postest_output_saom().
.out_flat_names <- c("format", "combineSameLevel")

.model_spec_flat <- c("effectName1", "diff1", "contrast1", "interaction1",
                      "intEffectNames1", "modEffectNames1", "second",
                      "effectName2", "diff2", "contrast2", "interaction2",
                      "intEffectNames2", "modEffectNames2", "perturbType1",
                      "perturbType2", "effectList", "effects", "depvar")


## --------------------------------------------------------------------------
## as_targets(args)
##
## Builds the targets object equivalent to a flat call's effect specification.
## Returns NULL when the call cannot be expressed yet, so the harness skips
## rather than silently comparing something else.
##
## Naming note: a single-effect flat call (effectName1 = ...) is internally
## named "single" and unwrapped from the result list.  The equivalent target is
## therefore named "single" too, so the comparison covers the RESULT WRAPPING
## as well as the values.  Whether a one-target call should unwrap on its own
## merits is a step-5 decision; the harness does not prejudge it.
## --------------------------------------------------------------------------
as_targets <- function(args) {
    if (!exists("make_postest_targets", mode = "function")) return(NULL)
    obj <- args$object
    eff <- if (!is.null(args$effects)) args$effects else obj$requestedEffects
    if (is.null(eff)) return(NULL)

    mk <- list(x = obj, effects = eff, includeDefaults = FALSE)
    if (!is.null(args$depvar)) mk$depvar <- args$depvar
    for (nm in .model_defaults_flat)
        if (nm %in% names(args)) mk[[nm]] <- args[[nm]]

    tg <- tryCatch(do.call(make_postest_targets, mk), error = function(e) NULL)
    if (is.null(tg)) return(NULL)

    ## Normalise both entry forms to one list of spec entries.
    specs <- if (!is.null(args$effectList)) args$effectList
             else if (!is.null(args$effectName1))
                 list(single = args[intersect(.model_spec_flat, names(args))])
             else return(NULL)

    for (snm in names(specs)) {
        sp <- specs[[snm]]
        if (isTRUE(sp$second)) {
            ## Built in the per-effect list form, which is what the corpus
            ## is here to exercise: the flat spec IS the numbered form, so
            ## translating it here means every corpus entry checks the list
            ## form against the flat interface it must reproduce.
            plist <- list(
                list(diff = sp$diff1, contrast = sp$contrast1,
                     perturbType = sp$perturbType1,
                     interaction = sp$interaction1,
                     intEffectNames = sp$intEffectNames1,
                     modEffectNames = sp$modEffectNames1),
                list(diff = sp$diff2, contrast = sp$contrast2,
                     perturbType = sp$perturbType2,
                     interaction = sp$interaction2,
                     intEffectNames = sp$intEffectNames2,
                     modEffectNames = sp$modEffectNames2))
            names(plist) <- c(sp$effectName1, sp$effectName2)
            tg <- tryCatch(suppressMessages(
                     set_second_diff(tg, plist, name = snm)),
                   error = function(e) NULL)
            if (is.null(tg)) return(NULL)
        } else {
            st <- list(x = tg, shortNames = sp$effectName1,
                       diff = sp$diff1, contrast = sp$contrast1,
                       perturbType = sp$perturbType1,
                       interaction = sp$interaction1,
                       intEffectNames = sp$intEffectNames1,
                       modEffectNames = sp$modEffectNames1)
            ## Per-spec overrides must be forwarded, and forwarded ONLY when
            ## present: `st[[f]] <- NULL` would drop an explicit NULL, so use
            ## single-bracket list assignment.
            for (f in c("level", "condition", "accumulated", "rateWeight",
                        "massContrasts"))
                if (f %in% names(sp)) st[f] <- list(sp[[f]])
            tg <- tryCatch(suppressMessages(do.call(set_target, st)),
                   error = function(e) NULL)
            if (is.null(tg)) return(NULL)
            ## Row identity must match the flat call's spec name.
            tg$name[tg$effectName1 == sp$effectName1 & !tg$second] <- snm
        }
    }
    tg
}
