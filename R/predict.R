##@predict.sienaFit S3 predict
##
## Six arguments, the same shape as marginalEffects():
##
##   object, data          the fit and the data it was fitted to
##   targets               WHAT to predict   -- make_predict_targets() +
##                         set_condition()
##   control_uncertainty   HOW SURE          -- set_postest_uncertainty_saom()
##   control_algo          HOW to simulate   -- set_postest_algo_saom()
##   control_out           HOW to report     -- set_postest_output_saom()
##
## The estimand (type, dynamic, level, condition, egoNormalize, accumulated,
## rateWeight) lives on the targets object; only the simulation settings are a
## question of how hard to work.  This replaced 44 flat formals, of which the
## controls already existed for marginalEffects() and are reused unchanged.
predict.sienaFit <- function(
    object,
    data,
    targets,
    control_uncertainty = set_postest_uncertainty_saom(),
    control_algo        = set_postest_algo_saom(),
    control_out         = set_postest_output_saom(),
    newParams = NULL,
    ...)
{
    if (!inherits(targets, "sienaPredictTargets"))
        stop("'targets' must come from make_predict_targets().", call. = FALSE)
    if (!inherits(control_uncertainty, "sienaPostestUncertainty"))
        stop("'control_uncertainty' must come from ",
             "set_postest_uncertainty_saom().", call. = FALSE)
    if (!inherits(control_algo, "sienaPostestControl"))
        stop("'control_algo' must come from set_postest_algo_saom().",
             call. = FALSE)
    if (!inherits(control_out, "sienaPostestOutput"))
        stop("'control_out' must come from set_postest_output_saom().",
             call. = FALSE)

    ## One prediction per requested row; a single row returns the bare frame,
    ## matching what marginalEffects() does with one target.
    rows <- order(targets$.seq)
    out  <- lapply(rows, function(i)
        .predictOne(object, data, targets, i,
                    control_uncertainty, control_algo, control_out, newParams))
    names(out) <- targets$name[rows]
    if (length(out) == 1L) out[[1L]] else out
}

## The former flat entry point, now internal: one row of the targets object,
## with every setting already resolved.
.predictOne <- function(object, newdata, targets, i,
                        unc, algo, outctl, newParams)
{
    st     <- .predictRowSettings(targets, i)
    effects <- attr(targets, "effects")
    depvar  <- attr(targets, "depvar")

    type         <- st$type
    dynamic      <- st$dynamic
    level        <- st$level
    condition    <- st$condition
    egoNormalize <- st$egoNormalize
    accumulated  <- st$accumulated
    rateWeight   <- st$rateWeight
    na.rm        <- st$na.rm

    uncertainty       <- isTRUE(unc$enabled)
    uncertaintyMode   <- unc$mode
    nsim              <- unc$nsim
    uncertaintySd     <- unc$sd
    uncertaintyCi     <- unc$ci
    uncertaintyMean   <- unc$simMean
    uncertaintyMedian <- unc$simMedian
    ciInterval        <- unc$ciInterval

    algorithm              <- algo$algorithm
    n3                     <- algo$n3
    n3PointEst             <- algo$n3PointEst
    n3BatchSize            <- algo$n3BatchSize
    useChangeContributions <- algo$useChangeContributions
    useCluster             <- algo$useCluster
    nbrNodes               <- algo$nbrNodes
    clusterType            <- algo$clusterType
    cl                     <- algo$cl
    batchDir               <- algo$batchDir
    prefix                 <- algo$prefix
    combineBatch           <- algo$combineBatch
    batchSize              <- algo$batchSize
    keepBatch              <- algo$keepBatch
    verbose                <- algo$verbose
    memoryScale            <- algo$memoryScale
    batchUnitBudget        <- algo$batchUnitBudget
    dynamicMinistepFactor  <- algo$dynamicMinistepFactor
    batch                  <- TRUE
    silent                 <- NULL

    returnDecisionDetails <- isTRUE(outctl$returnDecisionDetails)
    returnComponents      <- isTRUE(outctl$returnComponents)

  if (inherits(newdata, "sienaGroup"))
    stop("predict.sienaFit does not support multi-group data (sienaGroup).")

  if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]
  if (dynamic && is.null(algorithm)) stop("'algorithm' must be provided when dynamic = TRUE")
  if (dynamic && is.null(silent)) silent <- batch

  # ---- accumulated / rateWeight guards ----
  if (accumulated && !dynamic)
    stop("'accumulated = TRUE' requires 'dynamic = TRUE' ",
         "(accumulated prediction sums over simulated ministep chains).")
  if (rateWeight && dynamic)
    message("Note: for dynamic estimation, rate-weighting is absorbed ",
            "by the simulation. 'rateWeight' has no additional effect ",
            "and is ignored.")
  # add option to never batch? but then should do streaming...?
  if (is.null(batchSize)) {
      batchSize <- planBatch(
        data = newdata,
        depvar = depvar,
        nsim = nsim,
        nbrNodes = nbrNodes,
        useCluster = useCluster,
        dynamic = dynamic,
        n3 = n3,
        unitBudget = batchUnitBudget,
        dynamicMinistepFactor = dynamicMinistepFactor,
        memoryScale = memoryScale
      )
  }

  # ---- Theta / covariance with rate handling for rateWeight ----
  anyRateWeight <- rateWeight
  if (anyRateWeight && !dynamic) {
    if (isTRUE(object$cconditional)) {
      thetaHat <- coef(object)
      covTheta <- vcov(object)
    } else {
      eff_df  <- as.data.frame(object$requestedEffects)
      eff_inc <- eff_df[eff_df$include, ]
      hasNonConstantRates <- any(!eff_inc$basicRate & eff_inc$type == "rate")
      if (hasNonConstantRates)
        stop("rateWeight = TRUE is not supported when the model includes ",
             "non-constant rate effects (structural or covariate-dependent).")
      thetaHat <- coef(object, dropRates = FALSE)
      covTheta <- vcov(object, dropRates = FALSE)
    }
  } else {
    thetaHat <- coef(object)
    covTheta <- vcov(object)
  }

  # ---- Rate parameters for rateWeight (static path) ----
  rateParams <- NULL
  rateIdx    <- NULL
  if (anyRateWeight && !dynamic) {
    if (isTRUE(object$cconditional)) {
      rateParams <- object$rate
    } else {
      eff_df   <- as.data.frame(object$requestedEffects)
      eff_inc  <- eff_df[eff_df$include, ]
      rate_idx <- which(eff_inc$basicRate)
      theta_full <- coef(object, dropRates = FALSE)
      rateParams <- theta_full[rate_idx]
      rateIdx    <- rate_idx
    }
    if (length(rateParams) == 0L)
      stop("rateWeight = TRUE but no basic rate parameters found.")
  }

  # ---- Resolve condition ----
  if (!is.null(condition)) condition <- resolveCondition(condition)
  attachContribs <- !is.null(condition)

  # ---- Build contribFun ----
  if (dynamic) {
    if (is.null(effects)) effects <- object$requestedEffects
    n3Hat <- if (!is.null(n3PointEst)) n3PointEst else n3
    dynArgs <- list(
        ans                    = object,
        data                   = newdata,
        algorithm              = algorithm,
        effects                = effects,
        depvar                 = depvar,
        n3                     = n3Hat,
        batch                  = batch,
        silent                 = silent,
        returnWide             = TRUE
    )
    memCheck <- .checkDynMemory(
        data         = newdata,
        depvar       = depvar,
        effects      = effects,
        n3_per_batch = n3Hat,
        n3_uncert    = if (uncertainty) n3 else 0L,
        useCluster   = useCluster,
        nbrNodes     = nbrNodes,
        clusterType  = clusterType,
        uncertainty  = uncertainty,
        verbose      = verbose
    )
    nbrNodes <- memCheck$nbrNodes
    # contribFun <- makeContribFun("per_batch", dynArgs = dynArgs)
    .chains <- if (useChangeContributions) object$changeContributions else NULL
    if (!is.null(.chains)) {
        n3Hat <- if (!is.null(n3PointEst)) n3PointEst else n3
        n3Batch <- if (!is.null(n3BatchSize)) min(n3BatchSize, n3Hat) else n3Hat
        built <- .buildDynChainStore(.chains, dynArgs, length(.chains), n3Batch,
                                      depvar, effects, newdata, verbose,
                                      chainStoreMode = "auto")
        contribFun    <- built$contribFun
        nChainBatches <- built$nChainBatches
        rm(.chains); gc(verbose = FALSE)
    } else {
        contribFun    <- makeContribFun("per_batch", dynArgs = dynArgs)
        nChainBatches <- 1L
    }
  } else {
    staticContributions <- getStaticChangeContributions(
        ans              = object,
        data             = newdata,
        effects          = effects,
        depvar           = depvar,
        algorithm        = algorithm,
        includePermitted = TRUE,
        returnWide       = TRUE
    )
    contribFun <- makeContribFun("static",
        getContribFun = function(theta) staticContributions)
  }

  # ---- Build specs (N=1) ----
  metadata <- list(method = "predict", type = type,
                   depvar = depvar, dynamic = dynamic)
  specs <- setNames(list(makeSpec(
    predictFun   = predictProbability,
    predictArgs  = list(type = type, attachContribs = attachContribs),
    outcomeName  = type,
    level        = level,
    condition    = condition,
    na.rm        = na.rm,
    egoNormalize = egoNormalize,
    dynamic      = dynamic,
    accumulated  = accumulated,
    rateWeight   = rateWeight,
    jacobianFun  = if (!accumulated) predictProbabilityJac else NULL,
    metadata     = metadata
  )), type)

  # returnDecisionDetails: point-estimate call with attachContribs=TRUE
  decisionDetails <- NULL
  if (returnDecisionDetails) {
    cc_hat <- contribFun(thetaHat, 1L, 1L)
    decisionDetails <- predictProbability(cc_hat, thetaHat, type,
                                          attachContribs = TRUE,
                                          returnComponents = returnComponents)
    rm(cc_hat)
  }

  results <- sienaPostestimate(
    contribFun    = contribFun,
    nChainBatches = 1L,
    type          = type,
    specs        = specs,
    thetaHat     = thetaHat,
    covTheta     = covTheta,
    rateParams   = rateParams,
    rateIdx      = rateIdx,
    dynamic       = dynamic,
    dynArgs       = if (dynamic) dynArgs else NULL,
    n3            = n3,
    n3BatchSize   = n3BatchSize,
    useChangeContributions = if (dynamic) useChangeContributions else FALSE,
    uncertainty  = uncertainty,
    ## The analytic Jacobian below was wired in but unreachable: this call
    ## never passed a mode, so sienaPostestimate() took its "bootstrap"
    ## default and predictProbabilityJac() was supplied on every call and
    ## never consumed.  Exposing the mode is what makes the delta path
    ## reachable at all.
    uncertaintyMode = uncertaintyMode,
    nsim         = nsim,
    batchSize    = batchSize,
    useCluster   = useCluster,
    nbrNodes     = nbrNodes,
    clusterType  = clusterType,
    cl           = cl,
    batchDir     = batchDir,
    prefix       = prefix,
    keepBatch    = keepBatch,
    verbose      = verbose,
    na.rm        = na.rm,
    egoNormalize = egoNormalize,
    uncertaintySd         = uncertaintySd,
    uncertaintyCi         = uncertaintyCi,
    uncertaintyMean       = uncertaintyMean,
    uncertaintyMedian     = uncertaintyMedian,
    ciInterval            = ciInterval,
    decisionDetails = if (!is.null(decisionDetails))
                        setNames(list(decisionDetails), type) else NULL
  )

  # Stamp S3 class on each result; sienaPostestimate returns plain data.frames.
  results <- lapply(results, function(r) {
    class(r) <- c("sienaPrediction", class(r))
    r
  })

  result <- results[[type]]
  if (!is.null(decisionDetails))
    attr(result, "decisionDetails") <- decisionDetails
  result
}

print.sienaPrediction <- function(x, ...) {
  cat("SAOM Prediction\n")
  cat("  Type:      ", if (!is.null(attr(x, "type"))) attr(x, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(x, "depvar"))) attr(x, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(x, "level"))) attr(x, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(x, "dynamic"))) attr(x, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(x, "nsim"))) attr(x, "nsim") else NA, "\n")
  cat("\n")
  print.data.frame(x, ...)
  invisible(x)
}

##@summary.sienaPrediction S3 summary
summary.sienaPrediction <- function(object, ...) {
  cat("SAOM Prediction Summary\n")
  cat("  Type:      ", if (!is.null(attr(object, "type"))) attr(object, "type") else "unknown", "\n")
  cat("  Dep. var.: ", if (!is.null(attr(object, "depvar"))) attr(object, "depvar") else "unknown", "\n")
  cat("  Level:     ", if (!is.null(attr(object, "level"))) attr(object, "level") else "unknown", "\n")
  cat("  Dynamic:   ", if (!is.null(attr(object, "dynamic"))) attr(object, "dynamic") else FALSE, "\n")
  cat("  nsim:      ", if (!is.null(attr(object, "nsim"))) attr(object, "nsim") else NA, "\n")
  if (!is.null(attr(object, "condition")))
    cat("  Condition: ", paste(attr(object, "condition"), collapse = ", "), "\n")
  cat("\n")
  cat("  Rows:      ", nrow(object), "\n")
  cat("  Columns:   ", paste(names(object), collapse = ", "), "\n\n")
  summary.data.frame(object, ...)
}

# Estimate memory (GB) for a single getDynamicChangeContributions call.
#
# Returns a list with:
#   estGB_contrib  - estimated size of the contribMat (largest single object)
#   estGB_peak     - estimated peak memory during effect computation
#   estRows        - estimated total rows
#   nActor, nPer   - network dimensions
estimateDynMemory <- function(data, depvar, effects, n3,
                              algorithm = NULL) {
  dv <- data$depvars[[depvar]]
  dvDim <- dim(dv)
  nActor <- if (!is.null(dvDim) && length(dvDim) >= 1L) dvDim[1] else length(dv)
  nChoice <- if (!is.null(dvDim) && length(dvDim) >= 2L) dvDim[2] else nActor
  nPer <- if (!is.null(dvDim) && length(dvDim) >= 3L) max(1L, dvDim[3] - 1L) else 1L

  inc <- effects[effects$include, ]
  rateEff <- inc[inc$basicRate & inc$name == depvar, ]
  evalEff <- inc[inc$type != "rate" & inc$name == depvar, ]
  nEff <- nrow(evalEff)

  # Mean rate parameter from effects; fallback keeps estimate conservative.
  if (nrow(rateEff) > 0L) {
    meanRate <- mean(rateEff$initialValue, na.rm = TRUE)
  } else {
    meanRate <- nActor
  }

  estMinisteps <- as.numeric(nActor) * meanRate * as.numeric(nPer)
  estRows <- estMinisteps * as.numeric(nChoice) * as.numeric(n3)

  # contribMat stores nEff doubles plus metadata fields per row.
  bytesPerRow <- as.numeric(nEff) * 8 + 6 * 4
  estGB_contrib <- estRows * bytesPerRow / 1024^3

  # Peak includes contribMat + csMat + baseline and temporary vectors.
  estGB_peak <- estGB_contrib * 2 + estRows * 4 * 8 / 1024^3

  list(
    estGB_contrib = estGB_contrib,
    estGB_peak = estGB_peak,
    estRows = estRows,
    nActor = as.integer(nActor),
    nPer = as.integer(nPer),
    nEff = nEff,
    meanRate = meanRate,
    n3 = n3
  )
}

# Unified per-theta prediction from a compact contributions struct.
# Works for both static (ego/period/choice coords) and dynamic
# (chain/ministep/period/choice coords) -- groupColsList handles both.

# --------------------------------------------------------------------------
# predictProbabilityJac — analytical Jacobian paired with predictProbability.
#
# d(changeProb_i) / d(theta_k) = J_softmax[i, k]
# d(tieProb_i)   / d(theta_k) = density_i * J_softmax[i, k]
#
# Calling convention matches evalBatchJacobian:
#   cc, theta, changeProb, density, pa (spec$predictArgs), cs, ...
# Returns: n × K_eff matrix.
# --------------------------------------------------------------------------
predictProbabilityJac <- function(cc, theta, changeProb, density,
                                   pa, cs, ...) {
  Jp <- softmax_jac_rcpp(changeProb, cc$contribMat, cc$group_id)
  if (pa$type == "tieProb") density * Jp else Jp
}

predictProbability <- function(changeContributions, theta, type = "changeProb",
                               attachContribs = FALSE,
                               returnComponents = FALSE,
                               baseline = NULL,
                               outcomesOnly = FALSE,
                               ...) {
    contributions <- changeContributions
    theta_use <- theta[contributions$effectNames]
    names(theta_use) <- contributions$effectNames

    # changeStats is eagerly populated by source functions; fall back to
    # on-demand computation for direct callers with raw contribMat.
    if (is.null(contributions$changeStats))
      contributions$changeStats <- contribToChangeStats(
        contributions$contribMat, contributions$effectNames)
    cs <- contributions$changeStats

    # outcomesOnly: reuse pre-computed baseline vectors, skip full recompute.
    if (outcomesOnly && !is.null(baseline)) {
      if (type == "tieProb") {
        tp <- if (!is.null(baseline$tieProb)) baseline$tieProb
              else calculateTieProb(baseline$changeProb, cs$density)
        return(list(tieProb = tp))
      }
      return(list(changeProb = baseline$changeProb))
    }

    thetaEff <- buildThetaEff(theta_use, cs$changeStatsMap)

    # Use C++ values when available (dynamic), else compute in R (static).
    if (!is.null(contributions$changeUtility) &&
        !all(is.na(contributions$changeUtility))) {
      utility    <- contributions$changeUtility
      changeProb <- contributions$changeProbability
    } else {
      utility    <- calculateUtility(cs$csMat, thetaEff,
                                     contributions$permitted, cs$densityIdx)
      changeProb <- calculateChangeProb(utility, contributions$group_id)
    }

    # Return raw intermediates for downstream consumers (entropy, margins).
    if (returnComponents) {
      tieProb <- if (type == "tieProb")
        calculateTieProb(changeProb, cs$density) else NULL
      return(list(
        theta_use  = theta_use,
        thetaEff   = thetaEff,
        utility    = utility,
        changeProb = changeProb,
        tieProb    = tieProb
      ))
    }

    out <- groupColsList(contributions)
    out[["changeUtil"]] <- utility
    out[["changeProb"]] <- changeProb
    if (type == "tieProb") {
      out[["tieProb"]] <- calculateTieProb(changeProb, cs$density)
    }
    if (attachContribs) {
      out <- attachContributions(out, cs$csNames, cs$csMat, flip = FALSE)
    }
    attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
    class(out) <- "data.frame"
    out
}

# ---- Dynamic -------------------------------------------------------
# ---- Shared helpers -------------------------------------------------

# Calculate utility from change statistics and theta.
#
# Supports two calling conventions:
#   (a) Legacy: calculateUtility(mat, theta, permitted)
#       — simple mat %*% theta (eval-only models or backward-compat callers).
#   (b) changeStats: calculateUtility(mat, thetaEff, permitted, densityIdx)
#       — direction-dependent theta; density is column densityIdx in mat.
#
# thetaEff: either a named numeric vector (legacy) or a nEffects x 2 matrix
#   with columns "creation" and "dissolution" (changeStats).
# densityIdx: integer scalar — column index of density in mat (changeStats path).
calculateUtility <- function(mat, theta, permitted = NULL, densityIdx = NULL) {
  stopifnot(is.matrix(mat))

  if (is.null(densityIdx) || isTRUE(is.na(densityIdx))) {
    # No density column in changeStats; creation and dissolution theta are equal
    # (pure eval effects).  Use creation column when theta is a matrix.
    th <- if (is.matrix(theta)) theta[, "creation"] else theta
    util <- as.numeric(mat %*% th)
    # what about no change?
  } else if (!is.matrix(theta)) {
    # Legacy path: simple matrix-vector multiply.
    util <- as.numeric(mat %*% theta)
  } else {
    # changeStats path: u = d × (θ_density + Δs_rest %*% θ_rest)
    # Density column carries ±1 direction; non-density columns are pure Δs.
    thetaEff <- theta
    d <- as.integer(mat[, densityIdx])
    rest <- seq_len(ncol(mat))[-densityIdx]  # non-density columns
    n <- nrow(mat)
    util <- numeric(n)
    cre <- d ==  1L
    dis <- d == -1L
    # d=0 rows stay at 0 (no-change → no utility contribution)
    if (any(cre)) {
      util[cre] <- d[cre] * (thetaEff[densityIdx, "creation"] +
        as.numeric(mat[cre, rest, drop = FALSE] %*% thetaEff[rest, "creation"]))
    }
    if (any(dis)) {
      util[dis] <- d[dis] * (thetaEff[densityIdx, "dissolution"] +
        as.numeric(mat[dis, rest, drop = FALSE] %*% thetaEff[rest, "dissolution"]))
    }
  }

  if (!is.null(permitted)) {
    stopifnot(length(permitted) == length(util))
    util[!permitted] <- -Inf
  }
  util
}

calculateChangeProb <- function(utility, group_id) {
  # is as.numeric not already part of the rcpp softmax?
  as.numeric(softmax_rcpp_by_group(utility, group_id))
}

# wrapper not necessary anymore - clean up later
calculateTieProb <- function(prob, density) {
  calculate_tie_prob_cpp(prob, as.numeric(density))
}
