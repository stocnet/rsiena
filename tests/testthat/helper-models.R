# Helper: shared model fits — sourced once per test session.
#
# Quick mode (default): only the three base models (ans, ans2, ans3) are fitted.
# Full mode (RSENA_FULL_TESTS=1): additionally fits predict-specific models
# (interaction-without-main, co-evolution, creation/maintenance).
#
# Each test file loads its fixtures via load_fixture("name") at the top.
# Set RSENA_REBUILD_MODELS=1 to force recomputation and cache refresh.

cache_dir <- tryCatch({
  if (requireNamespace("testthat", quietly = TRUE)) {
    testthat::test_path("cache")
  } else {
    file.path("tests", "testthat", "cache")
  }
}, error = function(e) {
  file.path("tests", "testthat", "cache")
})

load_fixture <- function(name) {
  path <- file.path(cache_dir, paste0(name, ".rds"))
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  # sienaAlgorithm objects store a projname pointing to a temp file from the
  # session that built the cache — replace with a fresh temp path.
  if (inherits(obj, "sienaAlgorithm")) {
    obj$projname <- tempfile("Siena")
  }
  obj
}

# Guard: build fixtures only from a test runner, not standalone load_all().
.in_test_run <- function() {
  any(vapply(sys.calls(), function(x) {
    fn <- tryCatch(deparse(x[[1]], nlines = 1), error = function(e) "")
    grepl("(^|::)(test|test_file|test_dir|test_check)$", fn)
  }, logical(1))) || identical(Sys.getenv("RSENA_LOAD_MODELS"), "1")
}

if (.in_test_run()) {
  message("helper-models.R is being sourced")
  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  rebuild <- identical(Sys.getenv("RSENA_REBUILD_MODELS"), "1")
  FULL    <- identical(Sys.getenv("RSENA_FULL_TESTS"), "1")

  minimal_names <- c(
    "mynet", "mydata", "mymodel", "mycontrols", "ans",
    "mynet2", "mydata2", "mymodel2", "mycontrols2", "ans2",
    "mynet3", "mydata3", "mymodel3", "mycontrols3", "ans3"
  )
  full_names <- c(
    "mydata_ego", "mymodel_ego", "mycontrols_ego", "ans_ego",
    "mydata_int", "mymodel_int", "mycontrols_int", "ans_int", "ans_int_uncond",
    "mydata_co", "mymodel_co", "mycontrols_co", "ans_co",
    "mydata_cm", "mymodel_cm", "mycontrols_cm", "ans_cm",
    "mydata_2int", "mymodel_2int", "mycontrols_2int", "ans_2int"
  )
  snap_names_static <- c(
    "snap_me_static_transTrip",
    "snap_me_static_secondDiff",
    "snap_me_delta_transTrip",
    "snap_me_delta_interaction",
    "snap_me_delta_density",
    "snap_predict_changeProb",
    "snap_predict_tieProb"
  )
  needed <- if (FULL) c(minimal_names, full_names, snap_names_static) else c(minimal_names, snap_names_static)
  needed_paths <- file.path(cache_dir, paste0(needed, ".rds"))

  if (!rebuild && all(file.exists(needed_paths))) {
    message("Cache up to date — skipping fixture rebuild")
  } else {
    message("Rebuilding helper model fixtures")

    # ── Base model: density + recip + transTrip ──────────────────────────────────
    mynet    <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
    mydata   <- sienaDataCreate(mynet)
    mymodel  <- getEffects(mydata)
    mymodel  <- includeEffects(mymodel, transTrip, name = "mynet")
    mycontrols <- sienaAlgorithmCreate(projname = NULL, n3 = 200, cond = FALSE, seed = 42)
    ans <- siena07(
      mycontrols,
      data    = mydata,
      effects = mymodel,
      returnDeps               = FALSE,
      returnChangeContributions = TRUE,
      returnDataFrame          = TRUE
    )

    # ── Model 2: add transRecTrip (interaction/moderator tests in sienaMargins) ──
    mynet2    <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
    mydata2   <- sienaDataCreate(mynet2)
    mymodel2  <- getEffects(mydata2)
    mymodel2  <- includeEffects(mymodel2, transTrip, transRecTrip, name = "mynet2")
    mycontrols2 <- sienaAlgorithmCreate(projname = NULL, n3 = 50, cond = FALSE, seed = 7)
    ans2 <- siena07(
      mycontrols2,
      data    = mydata2,
      effects = mymodel2,
      returnChangeContributions = TRUE,
      returnDataFrame           = TRUE,
      silent                    = TRUE
    )

    # ── Model 3: inPop * recip via unspInt (custom interaction in sienaMargins) ──
    mynet3    <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
    mydata3   <- sienaDataCreate(mynet3)
    mymodel3  <- getEffects(mydata3)
    mymodel3  <- includeEffects(mymodel3, inPop, name = "mynet3")
    mymodel3  <- includeInteraction(mymodel3, recip, inPop, name = "mynet3")
    mycontrols3 <- sienaAlgorithmCreate(projname = NULL, n3 = 50, cond = FALSE, seed = 13)
    ans3 <- siena07(
      mycontrols3,
      data    = mydata3,
      effects = mymodel3,
      returnChangeContributions = TRUE,
      returnDataFrame           = TRUE,
      silent                    = TRUE
    )

    for (nm in minimal_names) {
      saveRDS(get(nm), file.path(cache_dir, paste0(nm, ".rds")))
    }

    # ── Numerical snapshots (built from ans / mydata above) ──────────────────
    # These are golden-value baselines for the postestimation refactor.
    # Rebuilt whenever minimal_names are rebuilt (same RSENA_REBUILD_MODELS=1 flag).
    # All calls use uncertainty=FALSE or deterministic uncertainty (delta) so
    # the cached numbers are exactly reproducible — no set.seed needed.
    snap_names <- c(
      "snap_me_static_transTrip",
      "snap_me_static_secondDiff",
      "snap_me_delta_transTrip",
      "snap_me_delta_interaction",
      "snap_me_delta_density",
      "snap_predict_changeProb",
      "snap_predict_tieProb"
    )
    snap_paths <- file.path(cache_dir, paste0(snap_names, ".rds"))

    if (rebuild || !all(file.exists(snap_paths))) {
      message("Building numerical snapshot fixtures…")

      # The builders below and test-snapshots.R must ask for exactly the same
      # quantity, or the snapshot pins nothing.  Both go through a helper of
      # this shape; keep them in step.
      snap_tg <- function(fit, eff, depvar, condition = NULL,
                          dependency = NULL, ...) {
        tg <- make_postest_targets(fit, effects = eff, depvar = depvar,
                                   type = "tieProb", level = "period",
                                   condition = condition,
                                   includeDefaults = FALSE)
        perturb <- list(...)
        for (nm in names(perturb))
          tg <- suppressMessages(do.call(set_target,
                  c(list(x = tg, shortNames = nm), perturb[[nm]])))
        if (!is.null(dependency))
          tg <- suppressMessages(set_dependency(tg, dependency))
        tg
      }
      no_unc  <- set_postest_uncertainty_saom(enabled = FALSE)
      quiet   <- set_postest_algo_saom(verbose = FALSE)

      # (a) firstDiff static — transTrip, period level, no uncertainty
      snap_me_static_transTrip <- marginalEffects(ans, mydata,
        targets = snap_tg(ans, mymodel, "mynet", condition = "density",
                          transTrip = list(diff = 1)),
        control_uncertainty = no_unc, control_algo = quiet)

      # (b) secondDiff static — transTrip × recip, period level, no uncertainty
      snap_me_static_secondDiff <- marginalEffects(ans, mydata,
        targets = suppressMessages(set_second_diff(
          snap_tg(ans, mymodel, "mynet", condition = "density"),
          list(transTrip = list(diff = 1),
               recip     = list(contrast = c(0, 1))))),
        control_uncertainty = no_unc, control_algo = quiet)

      # (c) delta SE — transTrip, period level, uncertaintyMode="delta"
      # Deterministic: no MCMC draws, pure finite-difference Jacobian.
      snap_me_delta_transTrip <- marginalEffects(ans, mydata,
        targets = snap_tg(ans, mymodel, "mynet", condition = "density",
                          transTrip = list(diff = 1)),
        ## delta needs nsim >= 1 to enable the uncertainty path
        control_uncertainty = set_postest_uncertainty_saom(
            mode = "delta", sd = TRUE, ci = FALSE, simMean = FALSE,
            nsim = 1L),
        control_algo = quiet)

      # (c2) delta SE — INTERACTION spec.  Tripwire for the planned analytic
      # Jacobian generalisation (postestimate_api_redesign.md Sec. 2.2):
      # Recorded while a declared relation still lowered to the old
      # interaction arguments, so this snapshot also pins the two encodings
      # together: the numbers must not move now that only the relation runs.
      snap_me_delta_interaction <- marginalEffects(ans2, mydata2,
        targets = snap_tg(ans2, mymodel2, "mynet2", condition = "recip",
                          dependency = transRecTrip ~ transTrip:recip,
                          transTrip = list(diff = 1)),
        control_uncertainty = set_postest_uncertainty_saom(
            mode = "delta", sd = TRUE, ci = FALSE, nsim = 1L),
        control_algo = quiet)

      # (c3) delta SE — DENSITY spec.  Same tripwire: density also returns NULL
      # from calculateUtilityDiffJacobian (delta_u = -2 * changeUtil makes A
      # dense over all K columns), so this too is on the FD fallback today.
      snap_me_delta_density <- marginalEffects(ans, mydata,
        targets = snap_tg(ans, mymodel, "mynet",
                          density = list(contrast = c(-1, 1))),
        control_uncertainty = set_postest_uncertainty_saom(
            mode = "delta", sd = TRUE, ci = FALSE, nsim = 1L),
        control_algo = quiet)

      # (d) predict changeProb — period level, no uncertainty
      snap_predict_changeProb <- predict(
        object  = ans,
        newdata = mydata,
        type    = "changeProb",
        level   = "period", condition = "density",
        uncertainty = FALSE, verbose = FALSE
      )

      # (e) predict tieProb — period level, no uncertainty
      snap_predict_tieProb <- predict(
        object  = ans,
        newdata = mydata,
        type    = "tieProb",
        level   = "period", condition = "density",
        uncertainty = FALSE, verbose = FALSE
      )

      for (nm in snap_names) {
        saveRDS(get(nm), file.path(cache_dir, paste0(nm, ".rds")))
      }
      message("Saved numerical snapshot fixtures to ", cache_dir)
    }

    if (FULL) {
      # ── ego-effect model: density + recip + egoX(mybeh) + egoX*transTrip ──
      mynet_ego   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
      mybeh_ego   <- coCovar(s50a[, 1])
      mydata_ego  <- sienaDataCreate(mynet_ego, mybeh_ego)
      mymodel_ego <- getEffects(mydata_ego)
      mymodel_ego <- includeEffects(mymodel_ego, egoX, interaction1 = "mybeh_ego",
                                    name = "mynet_ego")
      mymodel_ego <- includeEffects(mymodel_ego, transTrip, name = "mynet_ego")
      mymodel_ego <- includeInteraction(mymodel_ego, egoX, transTrip,
                                        interaction1 = c("mybeh_ego", ""),
                                        name = "mynet_ego")
      mycontrols_ego <- sienaAlgorithmCreate(projname = NULL, n3 = 50,
                                            cond = FALSE, seed = 99)
      ans_ego <- siena07(
        mycontrols_ego,
        data    = mydata_ego,
        effects = mymodel_ego,
        returnDeps               = TRUE,
        returnChangeContributions = TRUE,
        returnDataFrame           = TRUE,
        silent                    = TRUE
      )

      # ── interaction-without-main: recip * outPop unspInt WITHOUT outPop ──
      mynet_int   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
      mydata_int  <- sienaDataCreate(mynet_int)
      mymodel_int <- getEffects(mydata_int)
      mymodel_int <- includeInteraction(mymodel_int, recip, outPop, name = "mynet_int")
      mycontrols_int <- sienaAlgorithmCreate(projname = NULL, seed = 42, n3 = 50)
      ans_int <- siena07(
        mycontrols_int,
        data    = mydata_int,
        effects = mymodel_int,
        returnChangeContributions = TRUE,
        returnDataFrame           = TRUE,
        silent                    = TRUE
      )

      # cond=FALSE version of the same model
      mycontrols_int_uncond <- sienaAlgorithmCreate(
        projname = NULL, seed = 42, n3 = 50, cond = FALSE
      )
      ans_int_uncond <- siena07(
        mycontrols_int_uncond,
        data    = mydata_int,
        effects = mymodel_int,
        silent  = TRUE
      )

      # ── two-network co-evolution (50 actors, 2 periods per network) ──
      mynet_a    <- sienaDependent(array(c(s501, s502), dim = c(50, 50, 2)))
      mynet_b    <- sienaDependent(array(c(s502, s503), dim = c(50, 50, 2)))
      mydata_co  <- sienaDataCreate(mynet_a, mynet_b)
      mymodel_co <- getEffects(mydata_co)
      mymodel_co <- includeEffects(mymodel_co, transTrip, name = "mynet_a")
      mymodel_co <- includeEffects(mymodel_co, transTrip, name = "mynet_b")
      mycontrols_co <- sienaAlgorithmCreate(projname = NULL, seed = 1, n3 = 50, cond = FALSE)
      ans_co <- siena07(
        mycontrols_co,
        data    = mydata_co,
        effects = mymodel_co,
        returnChangeContributions = TRUE,
        returnDataFrame           = FALSE,
        silent                    = TRUE
      )

      # ── creation / maintenance effects ──
      mynet_cm   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
      mydata_cm  <- sienaDataCreate(mynet_cm)
      mymodel_cm <- getEffects(mydata_cm)
      mymodel_cm <- includeEffects(mymodel_cm, recip,     type = "creation")
      mymodel_cm <- includeEffects(mymodel_cm, transTrip, type = "eval")
      mymodel_cm <- includeEffects(mymodel_cm, transTrip, type = "endow")
      mycontrols_cm <- sienaAlgorithmCreate(projname = NULL, seed = 3, n3 = 50, cond = FALSE)
      ans_cm <- siena07(
        mycontrols_cm,
        data    = mydata_cm,
        effects = mymodel_cm,
        returnChangeContributions = TRUE,
        returnDataFrame           = FALSE,
        silent                    = TRUE
      )

      # ── double interaction (unspInt1: recip x inPop, unspInt2: recip x outPop) ──
      mynet_2int   <- sienaDependent(array(c(s501, s502, s503), dim = c(50, 50, 3)))
      mydata_2int  <- sienaDataCreate(mynet_2int)
      mymodel_2int <- getEffects(mydata_2int)
      mymodel_2int <- includeEffects(mymodel_2int, inPop, outPop, name = "mynet_2int")
      mymodel_2int <- includeInteraction(mymodel_2int, recip, inPop,  name = "mynet_2int")
      mymodel_2int <- includeInteraction(mymodel_2int, recip, outPop, name = "mynet_2int")
      mycontrols_2int <- sienaAlgorithmCreate(projname = NULL, n3 = 50,
                                              cond = FALSE, seed = 77)
      ans_2int <- siena07(
        mycontrols_2int,
        data    = mydata_2int,
        effects = mymodel_2int,
        returnChangeContributions = TRUE,
        returnDataFrame           = TRUE,
        silent                    = TRUE
      )

      for (nm in full_names) {
        if (exists(nm, inherits = FALSE)) {
          saveRDS(get(nm), file.path(cache_dir, paste0(nm, ".rds")))
        }
      }
    }
    message("Saved helper model fixtures to ", cache_dir)
  }
}