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
  if (dynamic) {
    predictFun  <- predictProbabilityDynamic
    predictArgs <- list(
        ans                    = object,
        data                   = newdata,
        algorithm              = algorithm,
        effects                = effects,
        type                   = type,
        depvar                 = depvar,
        n3                     = n3,
        useChangeContributions = useChangeContributions,
        batch                  = batch,
        silent                 = silent,
        attachContribs         = FALSE
    )

  } else {
    # Build the contribution matrix ONCE — reused across all theta draws
    staticContributions <- getStaticChangeContributions(
        ans            = object,
        data           = newdata,
        effects        = effects,
        depvar         = depvar,
        algorithm      = algorithm,
        includePermitted = TRUE,
        returnWide = TRUE
    )
    predictFun  <- predictProbabilityStatic
    predictArgs <- list(
        staticContributions = staticContributions,
        type                = type,
        attachContribs      = FALSE
    )
  }

  sienaPostestimate(
  predictFun               = predictFun,
  predictArgs              = predictArgs,
  outcomeName              = type,
  level                    = level,
  condition                = condition,
  sum_fun                  = sum_fun,
  na.rm                    = na.rm,
  thetaHat                = object[["theta"]], # adjust coef to be used
  covTheta                = object[["covtheta"]], # adjust voctheta to be used
  uncertainty              = uncertainty,
  nsim                     = nsim,
  uncertaintySd           = uncertaintySd,
  uncertaintyCi           = uncertaintyCi,
  uncertaintyMean         = uncertaintyMean,
  uncertaintyMedian       = uncertaintyMedian,
  uncertaintyProbs        = uncertaintyProbs,
  uncertaintyMcse         = uncertaintyMcse,
  uncertaintymcseBatches = uncertaintymcseBatches,
  useCluster               = useCluster,
  nbrNodes                 = nbrNodes,
  clusterType              = clusterType,
  cl                       = cl,
  batchDir                = batchDir,
  prefix                   = prefix,
  combineBatch            = combineBatch,
  batchSize               = batchSize,
  keepBatch               = keepBatch,
  verbose                  = verbose,
  useChangeContributions   = if (dynamic) useChangeContributions else NULL,
  egoNormalize             = egoNormalize,
  returnDecisionDetails    = returnDecisionDetails,
  metadata = list(
    method      = "predict",
    type        = type,
    level       = level,
    condition   = condition,
    depvar      = depvar,
    dynamic     = dynamic,
    nsim        = nsim,
    outcomeName = type
  )
  )
}

# Estimate memory (GB) for a single getDynamicChangeContributions call.
#
# Returns a list with:
#   estGB_contrib  – estimated size of the contribMat (the largest single object)
#   estGB_peak     – estimated peak memory during effect computation
#   estRows        – estimated total rows (n3 * sum_periods(ministeps * nChoice))
#   nActor, nPer   – network dimensions
#
# The estimate uses: rows ≈ n3 * nPer * nActor * meanRate * nChoice,
# where meanRate is the average rate parameter (ministeps/actor/period).
# Each row stores nEff doubles (contribMat) plus integer metadata fields.
# Peak includes csMat + baseline + one effect's temporaries (~3× contribMat).
estimateDynMemory <- function(data, depvar, effects, n3,
                              algorithm = NULL) {
    dv     <- data$depvars[[depvar]]
    dvDim  <- dim(dv)
    nActor <- if (!is.null(dvDim) && length(dvDim) >= 1L) dvDim[1] else length(dv)
    nChoice <- if (!is.null(dvDim) && length(dvDim) >= 2L) dvDim[2] else nActor
    nPer   <- if (!is.null(dvDim) && length(dvDim) >= 3L) max(1L, dvDim[3] - 1L) else 1L

    inc     <- effects[effects$include, ]
    rateEff <- inc[inc$basicRate & inc$name == depvar, ]
    evalEff <- inc[inc$type != "rate" & inc$name == depvar, ]
    nEff    <- nrow(evalEff)

    # Mean rate parameter — use estimated values, or n3 from algorithm
    if (nrow(rateEff) > 0L) {
        meanRate <- mean(rateEff$initialValue, na.rm = TRUE)
    } else {
        meanRate <- nActor  # rough fallback
    }

    # Estimated ministeps per chain = nActor * meanRate * nPer
    # (rate ~ expected ministeps/actor/period for unconditional;
    #  for conditional, actual ministeps per period ≈ observed distance)
    estMinisteps <- as.numeric(nActor) * meanRate * as.numeric(nPer)
    # Each ministep has nChoice candidate alterations
    estRows <- estMinisteps * as.numeric(nChoice) * as.numeric(n3)
    # contribMat: nEff doubles per row + 6 integer metadata fields
    bytesPerRow <- as.numeric(nEff) * 8 + 6 * 4
    estGB_contrib <- estRows * bytesPerRow / 1024^3
    # Peak: contribMat + csMat (≈ same size) + baseline vectors (1× rows × 8)
    # + one effect's temporaries (~3 vectors × rows × 8)
    estGB_peak <- estGB_contrib * 2 + estRows * 4 * 8 / 1024^3

    list(estGB_contrib = estGB_contrib, estGB_peak = estGB_peak,
         estRows = estRows, nActor = as.integer(nActor),
         nPer = as.integer(nPer), nEff = nEff,
         meanRate = meanRate, n3 = n3)
}

planBatch <- function(
  data, 
  depvar, 
  nsim,
  nbrNodes = 1L,
  useCluster = FALSE,
  clusterType = c("PSOCK", "FORK"),
  dynamic = FALSE,
  n3 = NULL,
  unitBudget = 2.5e8,
  dynamicMinistepFactor = 10,
  memoryScale = NULL,
  useChangeContributions = FALSE
) {
  clusterType <- match.arg(clusterType)
  dv     <- data$depvars[[depvar]]
  dvDim <- dim(dv)
  nEgo  <- if (!is.null(dvDim) && length(dvDim) >= 1L) as.integer(dvDim[1]) else as.integer(length(dv))
  nChoice   <- if (!is.null(dvDim) && length(dvDim) >= 2L) as.integer(dvDim[2]) else 1L
  nPer  <- if (!is.null(dvDim) && length(dvDim) >= 3L) max(1L, as.integer(dvDim[3] - 1L)) else 1L
  units  <- as.numeric(nEgo) * as.numeric(nChoice) * as.numeric(nPer)

  effective_workers <- if (isTRUE(useCluster)) max(1L, as.integer(nbrNodes)) else 1L

  dynamic_costly <- dynamic

  if (dynamic_costly) {
    n3Val <- if (is.null(n3)) 1L else max(1L, as.integer(n3))
    unitsPerCall <- units * as.numeric(dynamicMinistepFactor) * as.numeric(n3Val)
  } else {
    n3Val <- 1L
    unitsPerCall <- units
  }

  unitsPerAgg <- max(1.0, units * as.numeric(n3Val) *
                        if (dynamic_costly) as.numeric(dynamicMinistepFactor) else 1.0)

  # For FORK workers, children run in separate processes and their working
  # memory does not add to the parent's budget (CoW means pages are shared
  # until written).  Only the aggregation buffer in the parent matters.
  # For PSOCK, results are serialised back to the parent, so worker cost
  # does count against the budget.
  is_fork <- isTRUE(useCluster) && clusterType == "FORK"
  budgetForAgg <- if (is_fork) {
    as.numeric(unitBudget)
  } else {
    as.numeric(unitBudget) - effective_workers * unitsPerCall
  }

  if (budgetForAgg <= 0) {
    if (effective_workers > 1L) {
      warning(sprintf(
        "Memory budget (%.0f units) may be insufficient for %d parallel worker(s)\n        at %.0f units each. Consider reducing nbrNodes or increasing batchUnitBudget.",
        unitBudget, effective_workers, unitsPerCall
      ))
    }

    maxBatch <- max(1L, as.integer(floor(as.numeric(unitBudget) / unitsPerAgg)))
  } else {
    maxBatch <- max(1L, as.integer(floor(budgetForAgg / unitsPerAgg)))
  }

  batch_raw <- min(as.integer(nsim), maxBatch)

  if (!is.null(memoryScale) && as.integer(memoryScale) > 1L)
    batch_raw <- max(1L, as.integer(floor(batch_raw / as.integer(memoryScale))))

  k <- if (isTRUE(useCluster) && as.integer(nbrNodes) > 1L) as.integer(nbrNodes) else 1L
  if (k > 1L && as.integer(nsim) >= k) {
    b2 <- (batch_raw %/% k) * k
    batch_raw <- if (b2 < k) k else b2
  }
  min(max(1L, batch_raw), as.integer(nsim))
}

# Unified per-theta prediction from a compact contributions struct.
# Works for both static (ego/period/choice coords) and dynamic
# (chain/ministep/period/choice coords) -- groupColsList handles both.
predictProbability <- function(contributions, theta, type = "changeProb",
                               attachContribs = FALSE,
                               returnComponents = FALSE) {
    theta_use <- theta[contributions$effectNames]
    names(theta_use) <- contributions$effectNames

    # Ensure changeStats are cached on the contribution struct.
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

# Per-theta static prediction: build contributions once, reused per theta draw.
predictProbabilityStatic <- function(staticContributions, theta,
                                     type = "changeProb",
                                     attachContribs = FALSE) {
    predictProbability(staticContributions, theta, type,
                       attachContribs = attachContribs)
}

# ---- Dynamic -------------------------------------------------------

# Per-theta dynamic prediction: simulate chains, flatten to compact struct.
predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
                                      type = "changeProb", depvar, n3 = NULL,
                                      useChangeContributions = FALSE,
                                      batch = TRUE, silent = TRUE,
                                      attachContribs = FALSE) {
    contributions <- getDynamicChangeContributions(
      ans = ans, theta = theta, data = data, algorithm = algorithm,
      effects = effects, depvar = depvar, n3 = n3,
      useChangeContributions = useChangeContributions,
      returnWide = TRUE, batch = batch, silent = silent
    )
    predictProbability(contributions, theta, type,
                       attachContribs = attachContribs)
}

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

  if (is.null(densityIdx) || !is.matrix(theta)) {
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

# ===========================================================================
# === V1 BACKUP: predict.sienaFit + static/dynamic per-theta functions    ===
# === (data.frame-internal pipeline, superseded by V2 above)              ===
# ===========================================================================

# ##@predict.sienaFit
# predict.sienaFit <- function(
#     object,
#     newdata,
#     type = c("changeProb", "tieProb"), #behavior not implemented *here*
#     newParams = NULL,
#     effects = NULL,
#     depvar = NULL,
#     dynamic = FALSE,
#     algorithm = NULL,
#     n3 = 1000,
#     useChangeContributions = FALSE,
#     level = "period",
#     condition = NULL,
#     sum_fun = mean,
#     na.rm = TRUE,
#     uncertainty = TRUE,
#     uncertaintyMode = c("batch", "stream"),
#     useCluster = FALSE,
#     nbrNodes = 1,
#     nsim = 1000,
#     uncertaintySd = TRUE,
#     uncertaintyCi = TRUE,
#     uncertaintyProbs = c(0.025, 0.5, 0.975),
#     uncertaintyMcse = FALSE,
#     uncertaintymcseBatches = NULL,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batchDir = "temp",
#     prefix = "simBatch_b",
#     combineBatch = TRUE,
#     batch = TRUE,
#     silent = NULL,
#     batchSize = NULL,
#     keepBatch = FALSE,
#     verbose = TRUE,
#     memoryScale = NULL,
#     batchUnitBudget = 5e6,
#     dynamicMinistepFactor = 10,
#     ...
#   ){
#       uncertaintyMode <- match.arg(uncertaintyMode)
#       type <- match.arg(type)
#       if (is.null(depvar)) depvar <- names(newdata[["depvars"]])[1]
#       if (dynamic && is.null(algorithm)) {
#         stop("'algorithm' must be provided when dynamic = TRUE")
#         # could default now instead
#       }
#       if (dynamic && is.null(silent)) silent <- batch
#       if (is.null(batchSize)) {
#         batchSize <- planBatch(
#           data = newdata, 
#           depvar = depvar, 
#           nsim = nsim,
#           nbrNodes = nbrNodes, 
#           useCluster = useCluster,
#           dynamic = dynamic, 
#           n3 = n3,
#           unitBudget = batchUnitBudget,
#           dynamicMinistepFactor = dynamicMinistepFactor,
#           memoryScale = memoryScale
#         )
#       }
#       if (dynamic) {
#         predictFun <- predictProbabilityDynamic
#         predictArgs <- list(
#           ans = object,
#           data = newdata,
#           algorithm = algorithm,
#           effects = effects,
#           type = type,
#           depvar = depvar,
#           n3 = n3,
#           useChangeContributions = useChangeContributions,
#           batch = batch,
#           silent = silent
#         )
#       } else {
#           staticContributions <- getStaticChangeContributions(
#             ans = object,
#             data = newdata,
#             effects = effects,
#             depvar = depvar,
#             returnDataFrame = TRUE
#           )
#         predictFun <- predictProbabilityStatic
#         predictArgs <- list(
#           ans = object,
#           staticContributions = staticContributions,
#           type = type
#         )
#       }
# 
#       sienaPostestimate(
#         predictFun = predictFun,
#         predictArgs = predictArgs,
#         outcomeName = type,
#         level = level,
#         condition = condition,
#         sum_fun = sum_fun,
#         na.rm = na.rm,
#         thetaHat = object[["theta"]], # change into coef
#         covTheta = object[["covtheta"]], # change into cov
#         uncertainty = uncertainty,
#         uncertaintyMode = uncertaintyMode,
#         nsim = nsim,
#         uncertaintySd = uncertaintySd,
#         uncertaintyCi = uncertaintyCi,
#         uncertaintyProbs = uncertaintyProbs,
#         uncertaintyMcse = uncertaintyMcse,
#         uncertaintymcseBatches = uncertaintymcseBatches,
#         useCluster = useCluster,
#         nbrNodes = nbrNodes,
#         clusterType = clusterType,
#         cluster = cluster,
#         batchDir = batchDir,
#         prefix = prefix,
#         combineBatch = combineBatch,
#         batchSize = batchSize,
#         keepBatch = keepBatch,
#         verbose = verbose,
#         useChangeContributions = if (dynamic) useChangeContributions else NULL
#       )
# }
# 
# 
# 
# predictProbabilityStatic <- function(
#   ans, 
#   staticContributions, 
#   theta, 
#   type = "changeProb") {
#     effectNames <- getEffectNamesNoRate(ans[["requestedEffects"]])
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
#     staticContributions <- widenStaticContribution(staticContributions)
#     df <- predictProbability(staticContributions, effectNames, thetaNoRate,
#         group_vars = c("period", "ego"), type = type)
#     df <- conditionalReplace(df, df[["density"]] == -1,
#         setdiff(effectNames, "density"), function(x) x * -1)
#     df
# }
# 
# predictDynamic <- function(
#     ans,
#     newdata,
#     effects,
#     algorithm,
#     type = c("changeProb", "tieProb"),
#     depvar = NULL,
#     level = "period",
#     condition = NULL,
#     sum_fun = mean,
#     na.rm = TRUE,
#     n3 = 1000,
#     useChangeContributions = FALSE,
#     uncertainty = TRUE,
#     uncertaintyMode = c("batch", "stream"),
#     useCluster = FALSE,
#     nbrNodes = 1,
#     nsim = 100,
#     uncertaintySd = TRUE,
#     uncertaintyCi = TRUE,
#     uncertaintyProbs = c(0.025, 0.5, 0.975),
#     uncertaintyMcse = FALSE,
#     uncertaintymcseBatches = NULL,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batchDir = "temp",
#     prefix = "simBatch_b",
#     combineBatch = TRUE,
#     batchSize = NULL,
#     keepBatch = FALSE,
#     verbose = TRUE,
#     batch = TRUE,
#     silent = NULL,
#     memoryScale = NULL,
#     batchUnitBudget = 5e6,
#     dynamicMinistepFactor = 10
# ){
#   predict.sienaFit(
#     object = ans,
#     newdata = newdata,
#     type = type,
#     effects = effects,
#     depvar = depvar,
#     dynamic = TRUE,
#     algorithm = algorithm,
#     n3 = n3,
#     useChangeContributions = useChangeContributions,
#     level = level,
#     condition = condition,
#     sum_fun = sum_fun,
#     na.rm = na.rm,
#     uncertainty = uncertainty,
#     uncertaintyMode = uncertaintyMode,
#     useCluster = useCluster,
#     nbrNodes = nbrNodes,
#     nsim = nsim,
#     uncertaintySd = uncertaintySd,
#     uncertaintyCi = uncertaintyCi,
#     uncertaintyProbs = uncertaintyProbs,
#     uncertaintyMcse = uncertaintyMcse,
#     uncertaintymcseBatches = uncertaintymcseBatches,
#     clusterType = clusterType,
#     cluster = cluster,
#     batchDir = batchDir,
#     prefix = prefix,
#     combineBatch = combineBatch,
#     batchSize = batchSize,
#     keepBatch = keepBatch,
#     verbose = verbose,
#     batch = batch,
#     silent = silent,
#     memoryScale = memoryScale,
#     batchUnitBudget = batchUnitBudget,
#     dynamicMinistepFactor = dynamicMinistepFactor
#   )
# }
# 
# predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
#     type = "changeProb", depvar, n3 = NULL, useChangeContributions = FALSE,   
#     batch = TRUE, silent = TRUE) {
#     effectNames <- getEffectNamesNoRate(effects)
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
# 
#     df <- getDynamicChangeContributions(
#         ans = ans, 
#         theta = theta, 
#         data = data, 
#         algorithm = algorithm,
#         effects = effects, 
#         depvar = depvar, 
#         n3 = n3, 
#         useChangeContributions = useChangeContributions, 
#         returnDataFrame = TRUE,
#         batch = batch,
#         silent = silent
#     )
#     df <- widenDynamicContribution(df)
#     df <- predictProbability(df, effectNames, thetaNoRate,
#         group_vars = c("chain", "period", "ministep"), type = type)
#     conditionalReplace(df, df[["density"]] == -1,
#         setdiff(effectNames, "density"), function(x) x * -1)
# }

# ===========================================================================
# === V1 BACKUP: predictProbability, addUtilityColumn, addProbabilityColumn ===
# === (V1-only per-row helpers; calculateUtility is now a live function)    ===
# ===========================================================================

# 
# #' Shared core: compute utility and probability columns
# #' Used by predictProbabilityStatic, predictProbabilityDynamic,
# #' predictFirstDiff, and predictSecondDiff.
# predictProbability <- function(df, effectNames, thetaNoRate, group_vars, type) {
#     df <- addUtilityColumn(df, effectNames, thetaNoRate)
#     df <- addProbabilityColumn(df, group_vars = group_vars, type = type)
#     df
# }
# 
# addUtilityColumn <- function(df, effectNames, theta) {
#   availableEffects <- effectNames[effectNames %in% names(df)]
#   if (!length(availableEffects)) {
#     stop("No requested effect columns found in data for utility calculation")
#   }
# 
#   # checked already?
#   if (!is.null(names(theta))) {
#     theta_use <- theta[availableEffects]
#   } else {
#     theta_use <- theta[seq_along(availableEffects)]
#   }
# 
#   if (requireNamespace("data.table", quietly = TRUE) && 
#       data.table::is.data.table(df)) {
#     df[, ("changeUtil") := calculateUtility(
#       as.matrix(.SD), 
#       theta_use), 
#       .SDcols = availableEffects]
#   } else {
#     df[["changeUtil"]] <- calculateUtility(
#       as.matrix(df[, availableEffects, drop = FALSE]), 
#       theta_use)
#     return(df)
#   }
# }
# 
# # for density == 0: can we filter density 0 out?
# addProbabilityColumn <- function(
#     df,
#     group_vars,
#     type = "changeProb"
# ) {
#   stopifnot(is.data.frame(df))
#   # Zero-copy list of grouping column pointers -> single C++ call
#   group_cols <- lapply(group_vars, function(v) as.integer(df[[v]]))
#   df[["changeProb"]] <- as.vector(softmax_rcpp_grouped_lst(df[["changeUtil"]], group_cols))
#   if (type == "tieProb") {
#     df[["tieProb"]] <- df[["changeProb"]]
#     if ("density" %in% names(df)) {
#       idx_neg1 <- which(df[["density"]] == -1)
#       df[["tieProb"]][idx_neg1] <- 1 - df[["tieProb"]][idx_neg1]
#       df[["tieProb"]][df[["density"]] == 0] <- NA
#     }
#   }
#   df
# }
