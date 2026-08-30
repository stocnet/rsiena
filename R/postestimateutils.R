# --------------------------------------------------------------------------
# makeSpec — construct and validate a spec entry.
#
# Required: predictFun, outcomeName.
# All other fields have documented defaults.
# --------------------------------------------------------------------------
makeSpec <- function(predictFun,
                     outcomeName,
                     predictArgs           = list(),
                     level                 = "period",
                     condition             = NULL,
                     accumulated           = FALSE,
                     emitLevel             = "ministepChoice",
                     na.rm                 = TRUE,
                     egoNormalize          = FALSE,
                     rateWeight            = FALSE,
                     dynamic               = FALSE,
                     massContrasts         = FALSE,
                     returnDecisionDetails = FALSE,
                     jacobianFun           = NULL,
                     metadata              = list(),
                     ...) {
  dots <- list(...)
  if (length(dots) > 0L)
    warning("makeSpec: unknown fields ignored: ",
            paste(names(dots), collapse = ", "))
  list(
    predictFun            = predictFun,
    predictArgs           = predictArgs,
    outcomeName           = outcomeName,
    level                 = level,
    condition             = condition,
    accumulated           = accumulated,
    ## emitLevel: the level this spec's predictFun EMITS rows at, named with the
    ## same vocabulary as `level`.  Every row-wise prediction emits at
    ## "ministepChoice" (one row per alternative per ministep), which is why it
    ## has never needed saying; the period-level quantities emit at
    ## "chainEgoChoice".  It governs two things: whether .prepBatchContext()'s
    ## ministep-keyed caches and change-statistic condition resolution apply,
    ## and which `level` values can be aggregated to (.canAggregateTo).
    ##
    ## This is NOT the same axis as `accumulated`, which still emits ministep
    ## rows and only asks for a different aggregation over them.
    emitLevel             = emitLevel,
    na.rm                 = na.rm,
    egoNormalize          = egoNormalize,
    rateWeight            = rateWeight,
    dynamic               = dynamic,
    massContrasts         = massContrasts,
    returnDecisionDetails = returnDecisionDetails,
    jacobianFun           = jacobianFun,
    metadata              = metadata
  )
}

# --------------------------------------------------------------------------
# makeContribFun
#
# Backward compatible entry point used across postestimation call paths.
# Supports both legacy mode-based calls and the newer chainStore interface.
# --------------------------------------------------------------------------
makeContribFun <- function(mode = NULL, store = NULL, effects = NULL,
                           depvar = NULL, keepContribMat = NULL,
                           getContribFun = NULL, dynArgs = NULL,
                           preloadedChains = NULL, data = NULL) {
  # Backward compat: old string-mode interface -> wrap in chainStore
  if (is.character(mode)) {
    if (mode == "static") {
      store <- NULL
    } else if (mode == "preloaded" && !is.null(preloadedChains)) {
      nTotal <- length(preloadedChains)
      store <- chainStore_memory(preloadedChains, nTotal)
    } else if (mode == "per_batch" && !is.null(dynArgs)) {
      n3 <- if (!is.null(dynArgs$n3)) dynArgs$n3 else 100L
      store <- chainStore_simulate(dynArgs, n3, n3)
    }
  }

  is_static <- is.null(store)
  if (is.null(keepContribMat))
    keepContribMat <- is_static

# Force evaluation for use in closure below.
  force(store)
  force(effects)
  force(depvar)
  force(keepContribMat)
  force(getContribFun)
  force(data)

  function(theta, batchIdx, nBatches, useChangeContributions = FALSE) {
    cc <- if (is.null(store)) {
      getContribFun(theta)
    } else if (store$mode == "simulate") {
      store$getBatch(batchIdx, theta = theta)
    } else {
      raw <- store$getBatch(batchIdx)
      result <- flattenAndEnrichWide(raw, effects, depvar, data = data)
      rm(raw)
      result
    }

    if (is.null(cc$changeStats))
      cc$changeStats <- contribToChangeStats(cc$contribMat, cc$effectNames,
                                             direction = cc$direction)
    if (!keepContribMat) cc$contribMat <- NULL
    cc
  }
}

## --------------------------------------------------------------------------
## .postestTheta — the parameter vector and covariance to evaluate at
##
## Which parameters go into theta, and which uncertainty may perturb, depends
## on how the model stored its rates, and that in turn differs between the
## static and dynamic paths.  Four cases, none of them interchangeable, which
## is why this is worth naming rather than reading inline:
##
##   dynamic + conditional    rates fixed, spliced into effects$initialValue
##                            so siena07 has rate rows; theta stays non-rate
##   dynamic + unconditional  rates included, so the chain simulation and the
##                            uncertainty draws see the full vector
##   static  + rateWeight     rates included on unconditional models, so
##                            lambda uncertainty is drawn jointly
##   static  otherwise        rates dropped; nothing downstream needs them
##
## Returns the five things the rest of marginalEffects() needs; `effects` is
## among them because the conditional dynamic case adds rows to it.
## --------------------------------------------------------------------------
## `anyRateWeight` is resolved by the CALLER, because the two entry points
## carry it differently: marginalEffects() has a per-target flag and has to
## reduce over its specs, predict() has one scalar.  Everything downstream of
## that reduction is identical, and it is the part worth having in one place --
## the width of theta/covTheta, the condEffects splice, the rate bookkeeping
## and the non-constant-rate guard.
.postestTheta <- function(object, effects, dynamic, anyRateWeight) {
    # ---- Theta / covariance with proper rate handling ----
    # Dynamic path: getDynamicChangeContributions re-runs siena07, which
    # needs rate parameters in effects$initialValue.  Unconditional models
    # store rates inside theta; conditional models store them separately
    # in object$rate with effects rows in attr(object$f, "condEffects").

    if (dynamic && isTRUE(object$cconditional)) {
        # Conditional estimation: rates are fixed (not estimated jointly).
        # Splice condEffects back into effects so siena07 has rate rows.
        condEff <- attr(object$f, "condEffects")
        if (!is.null(condEff) && nrow(condEff) > 0L) {
            condEff$initialValue <- object$rate
            effects <- rbind(condEff, effects)
        }
        # theta: non-rate only (rates are in effects$initialValue above).
        thetaHat <- coef(object)
        covTheta <- vcov(object)
        # covTheta stays non-rate-only: uncertainty draws do not perturb
        # the fixed rate parameters.
    } else if (dynamic) {
        # Unconditional: include rate parameters in theta/covTheta so that
        # getDynamicChangeContributions and uncertainty draws get the full
        # parameter vector.  Avoids the fragile backfill fallback.
        thetaHat <- coef(object, dropRates = FALSE)
        covTheta <- vcov(object, dropRates = FALSE)
    } else {
        # Static path: include basic rate parameters in the draw distribution
        # when rateWeight is active on an unconditional model, so that lambda
        # uncertainty is captured jointly with behavioral param uncertainty.
        # For conditional models rates are fixed by design: no change needed.
        # predictProbability / predictFirstDiff use name-based theta lookup so
        # extra rate entries in theta are silently ignored for all other specs.
        if (anyRateWeight && !isTRUE(object$cconditional)) {
            thetaHat <- coef(object, dropRates = FALSE)
            covTheta <- vcov(object, dropRates = FALSE)
        } else {
            thetaHat <- coef(object)
            covTheta <- vcov(object)
        }
    }

    # ---- Rate parameters for rateWeight (static path) ----
    # rateParams: point-estimate rates (for hat draw and conditional fallback).
    # rateIdx:    integer positions of basic rate params within the
    #             dropRates=FALSE theta vector; NULL for conditional models.
    #             Used in one_batch to extract per-draw rates from theta_sim.
    rateParams <- NULL
    rateIdx    <- NULL
    if (anyRateWeight && !dynamic) {
        if (isTRUE(object$cconditional)) {
            # Conditional: rates fixed at point estimate; rateIdx stays NULL.
            rateParams <- object$rate
        } else {
            eff_df  <- as.data.frame(object$requestedEffects)
            eff_inc <- eff_df[eff_df$include, ]
            # Guard: non-constant rates (structural / covariate) make the
            # scalar lambda_p multiplication incorrect.  Those require a
            # per-actor rate evaluation that is not yet implemented.
            hasNonConstantRates <- any(
                !eff_inc$basicRate & eff_inc$type == "rate")
            if (hasNonConstantRates)
                stop(
                    "rateWeight = TRUE is not supported when the model ",
                    "includes non-constant rate effects (structural or ",
                    "covariate-dependent). rateWeight multiplies by the ",
                    "basic period rate only and would ignore actor-level ",
                    "rate heterogeneity.")
            rate_idx   <- which(eff_inc$basicRate)
            theta_full <- coef(object, dropRates = FALSE)
            rateParams <- theta_full[rate_idx]
            rateIdx    <- rate_idx
        }
        if (length(rateParams) == 0L)
            stop("rateWeight = TRUE but no basic rate parameters found.")
    }

    list(thetaHat = thetaHat, covTheta = covTheta, effects = effects,
         rateParams = rateParams, rateIdx = rateIdx)
}

# --------------------------------------------------------------------------
# .postestSeeds — resolve the randomness for one postestimation call
#
# One seed, taken from control_algo$seed, else algorithm$randomSeed, else none
# (unreproducible, as in siena07 — set.seed() does not govern either).  It
# feeds the simulations via algorithm$randomSeed, and the theta draws and
# per-batch chains via seeds derived from it.
# --------------------------------------------------------------------------

# Evaluate expr with the RNG seeded, then put the caller's stream back exactly
# as it was found -- including the case where there was no stream at all.
.withSeed <- function(seed, expr) {
    have <- exists(".Random.seed", envir = globalenv())
    old  <- if (have) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
        if (!is.null(old))
            assign(".Random.seed", old, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv()))
            rm(".Random.seed", envir = globalenv())
    }, add = TRUE)
    set.seed(seed)
    expr
}

.postestSeeds <- function(algorithm, algoSeed = NULL) {
    simSeed <- if (!is.null(algoSeed)) {
        as.integer(algoSeed)
    } else if (!is.null(algorithm) && !is.null(algorithm$randomSeed)) {
        rs <- algorithm$randomSeed
        if (is.numeric(rs)) {
            as.integer(rs)
        } else {
            warning("algorithm$randomSeed is not numeric (class '",
                    paste(class(rs), collapse = "/"),
                    "') and cannot be used. This run is not reproducible; ",
                    "pass a number via set_postest_algo_saom(seed = ).",
                    call. = FALSE)
            NULL
        }
    } else {
        NULL
    }

    if (!is.null(simSeed) && !is.null(algorithm))
        algorithm$randomSeed <- simSeed

    # Derived by drawing, not by offsetting: adjacent seeds carry no
    # independence guarantee, and simSeed + 1 collided with a batch seed.
    derived <- if (is.null(simSeed)) NULL
               else .withSeed(simSeed, sample.int(.Machine$integer.max, 2L))

    list(algorithm = algorithm,
         simSeed   = simSeed,
         drawSeed  = if (is.null(derived)) NULL else derived[[1L]],
         chainSeed = if (is.null(derived)) NULL else derived[[2L]])
}

## The level a spec emits rows at, defaulting for specs built before emitLevel
## existed (or by hand) -- everything row-wise emits at ministepChoice.
.specEmitLevel <- function(spec) {
    e <- spec[["emitLevel"]]
    if (is.null(e)) "ministepChoice" else e
}
