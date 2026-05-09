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