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
    cluster = NULL,
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
        silent                 = silent
    )
  } else {
    # Build the contribution matrix ONCE — reused across all theta draws
    staticContributions <- getStaticChangeContributions(
        ans     = object,
        data    = newdata,
        effects = effects,
        depvar  = depvar,
        returnWide = TRUE
    )
    predictFun  <- predictProbabilityStatic
    predictArgs <- list(
        staticContributions = staticContributions,
        type     = type
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
  cluster                  = cluster,
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

planBatch <- function(
  data, 
  depvar, 
  nsim,
  nbrNodes = 1L,
  useCluster = FALSE,
  dynamic = FALSE,
  n3 = NULL,
  unitBudget = 2.5e8,
  dynamicMinistepFactor = 10,
  memoryScale = NULL
) {
  dv     <- data$depvars[[depvar]]
  dvDim <- dim(dv)
  nEgo  <- if (!is.null(dvDim) && length(dvDim) >= 1L) as.integer(dvDim[1]) else as.integer(length(dv))
  nChoice   <- if (!is.null(dvDim) && length(dvDim) >= 2L) as.integer(dvDim[2]) else 1L
  nPer  <- if (!is.null(dvDim) && length(dvDim) >= 3L) max(1L, as.integer(dvDim[3] - 1L)) else 1L
  units  <- as.numeric(nEgo) * as.numeric(nChoice) * as.numeric(nPer)

  effective_workers <- if (isTRUE(useCluster)) max(1L, as.integer(nbrNodes)) else 1L

  if (dynamic) {
    n3Val <- if (is.null(n3)) 1L else max(1L, as.integer(n3))
    unitsPerCall <- units * as.numeric(dynamicMinistepFactor) * as.numeric(n3Val)
  } else {
    n3Val <- 1L
    unitsPerCall <- units
  }


  unitsPerAgg <- max(1.0, units * as.numeric(n3Val) *
                        if (dynamic) as.numeric(dynamicMinistepFactor) else 1.0)

  budgetForAgg <- as.numeric(unitBudget) - effective_workers * unitsPerCall

  if (budgetForAgg <= 0) {
    if (effective_workers > 1L) {
      warning(sprintf(
        "Memory budget (%.0f units) may be insufficient for %d parallel worker(s)
        at %.0f units each. Consider reducing nbrNodes or increasing batchUnitBudget.",
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
predictProbability <- function(contributions, theta, type = "changeProb") {
  # should not be necessary?
    theta_use <- theta[contributions$effectNames]
    names(theta_use) <- contributions$effectNames
    utility <- calculateUtility(contributions$contribMat, theta_use)
    changeProb <- calculateChangeProb(utility, contributions$group_id)
    out <- data.frame(groupColsList(contributions),
                      changeUtil = utility, changeProb = changeProb,
                      stringsAsFactors = FALSE)
    if (type == "tieProb") {
      out[["tieProb"]] <- calculateTieProb(
        changeProb,
        contributions$contribMat[, grep("density", contributions$effectNames, fixed = TRUE)[1L]]
      )
    }
    out <- attachContribColumns(out, contributions$effectNames,
                                contributions$contribMat, flip = TRUE)
    # if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out) # data.table removed
    out
}

# Per-theta static prediction: build contributions once, reused per theta draw.
predictProbabilityStatic <- function(staticContributions, theta,
                                     type = "changeProb") {
    predictProbability(staticContributions, theta, type)
}

# ---- Dynamic -------------------------------------------------------

# Per-theta dynamic prediction: simulate chains, flatten to compact struct.
predictProbabilityDynamic <- function(ans, data, theta, algorithm, effects,
                                      type = "changeProb", depvar, n3 = NULL,
                                      useChangeContributions = FALSE,
                                      batch = TRUE, silent = TRUE) {
    contributions <- getDynamicChangeContributions(
      ans = ans, theta = theta, data = data, algorithm = algorithm,
      effects = effects, depvar = depvar, n3 = n3,
      useChangeContributions = useChangeContributions,
      returnWide = TRUE, batch = batch, silent = silent
    )
    predictProbability(contributions, theta, type)
}

# ---- Shared helpers -------------------------------------------------

calculateUtility <- function(mat, theta) {
  stopifnot(is.matrix(mat))
  as.numeric(mat %*% theta)
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
