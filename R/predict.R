##@predict.sienaFit
predict.sienaFit <- function(
    object,
    newdata,
    type = c("changeProb", "tieProb"),
    newParams = NULL,
    effects = NULL,
    depvar = NULL,
    dynamic = FALSE,
    algorithm = NULL,
    n3 = 1000,
    useChangeContributions = FALSE,
    level = "period",
    condition = NULL,
    sum_fun = mean,
    na.rm = TRUE,
    uncertainty = TRUE,
    useCluster = FALSE,
    nbrNodes = 1,
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyMean = FALSE,
    uncertaintyMedian = FALSE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE,
    uncertaintymcseBatches = NULL,
    clusterType = c("PSOCK", "FORK"),
    cl = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    combineBatch = TRUE,
    batch = TRUE,
    silent = NULL,
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    memoryScale = NULL,
    batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10,
    egoNormalize = TRUE,
    returnDecisionDetails = FALSE,
    ...
) {
  if (inherits(newdata, "sienaGroup"))
    stop("predict.sienaFit does not support multi-group data (sienaGroup).")
  type             <- match.arg(type)

  if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]
  if (dynamic && is.null(algorithm)) stop("'algorithm' must be provided when dynamic = TRUE")
  if (dynamic && is.null(silent)) silent <- batch
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

  thetaHat <- object[["theta"]]
  covTheta <- object[["covtheta"]]

  # ---- Resolve condition ----
  if (!is.null(condition)) condition <- resolveCondition(condition)
  attachContribs <- !is.null(condition)

  # ---- Build contribFun ----
  if (dynamic) {
    dynArgs <- list(
        ans                    = object,
        data                   = newdata,
        algorithm              = algorithm,
        effects                = effects,
        depvar                 = depvar,
        n3                     = n3,
        batch                  = batch,
        silent                 = silent,
        returnWide             = TRUE
    )
    contribFun <- makeContribFun("per_batch", dynArgs = dynArgs)
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
    predictFun   = predictProbFun,
    predictArgs  = list(type = type, attachContribs = attachContribs),
    outcomeName  = type,
    level        = level,
    condition    = condition,
    na.rm        = na.rm,
    egoNormalize = egoNormalize,
    dynamic      = dynamic,
    metadata     = metadata
  )), type)

  # returnDecisionDetails: point-estimate call with attachContribs=TRUE
  decisionDetails <- NULL
  if (returnDecisionDetails) {
    cc_hat <- contribFun(thetaHat, 1L, 1L)
    decisionDetails <- predictProbability(cc_hat, thetaHat, type,
                                          attachContribs = TRUE)
    rm(cc_hat)
  }

  results <- sienaPostestimate(
    contribFun    = contribFun,
    nChainBatches = 1L,
    type          = type,
    specs        = specs,
    thetaHat     = thetaHat,
    covTheta     = covTheta,
    useChangeContributions = if (dynamic) useChangeContributions else FALSE,
    uncertainty  = uncertainty,
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
    uncertaintyProbs      = uncertaintyProbs,
    uncertaintyMcse       = uncertaintyMcse,
    uncertaintymcseBatches = uncertaintymcseBatches,
    decisionDetails = if (!is.null(decisionDetails))
                        setNames(list(decisionDetails), type) else NULL
  )

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

# ---- predictFun for predict.sienaFit specs ------------------------------
# Wraps predictProbability in the (cc, theta, baseline, ...) interface
# expected by makeEstimatorFun.  baseline (returnComponents=TRUE list) is
# not used: predictProbability needs data.frame output, so it re-calls
# with returnComponents=FALSE.  For dynamic contributions the re-call is
# cheap (reads precomputed C++ fields); for static it is deterministic.
#
# outcomesOnly=TRUE: return just the outcome vector(s) as a named list,
# using the pre-computed baseline to avoid redundant work.
predictProbFun <- function(changeContributions, theta, baseline,
                           type, attachContribs, outcomesOnly = FALSE) {
  if (outcomesOnly) {
    # baseline already has changeProb/tieProb — no recomputation needed
    if (type == "tieProb") {
      cs <- changeContributions$changeStats
      if (is.null(cs))
        cs <- contribToChangeStats(changeContributions$contribMat,
                                   changeContributions$effectNames)
      tp <- if (!is.null(baseline$tieProb)) baseline$tieProb
            else calculateTieProb(baseline$changeProb, cs$density)
      return(list(tieProb = tp))
    }
    return(list(changeProb = baseline$changeProb))
  }
  predictProbability(changeContributions, theta, type,
                     attachContribs = attachContribs)
}

predictProbability <- function(contributions, theta, type = "changeProb",
                               attachContribs = FALSE,
                               returnComponents = FALSE) {
    theta_use <- theta[contributions$effectNames]
    names(theta_use) <- contributions$effectNames

    # changeStats is eagerly populated by source functions; fall back to
    # on-demand computation for direct callers with raw contribMat.
    if (is.null(contributions$changeStats))
      contributions$changeStats <- contribToChangeStats(
        contributions$contribMat, contributions$effectNames)
    cs <- contributions$changeStats

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
  as.numeric(softmax_arma_by_group(utility, group_id))
}

# wrapper not necessary anymore - clean up later
calculateTieProb <- function(prob, density) {
  calculate_tie_prob_cpp(prob, as.numeric(density))
}
