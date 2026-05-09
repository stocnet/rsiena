##@entropy.sienaFit Entropy with optional uncertainty via drawSimBatch
entropy.sienaFit <- function(object, data, effects = NULL, depvar = NULL,
    dynamic = FALSE, algorithm = NULL, n3 = 500,
    useChangeContributions = FALSE,
    level = "ego", sum_fun = mean, na.rm = TRUE,
    uncertainty = TRUE,
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyMean = FALSE,
    uncertaintyMedian = FALSE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    batchDir = "temp",
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    memoryScale = NULL,
    batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10,
    ...) {

    if (is.null(effects)) effects <- object$requestedEffects
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    if (dynamic && is.null(algorithm))
        stop("'algorithm' must be provided when dynamic = TRUE")

    if (is.null(batchSize)) {
        batchSize <- planBatch(
            data = data, depvar = depvar, nsim = nsim,
            nbrNodes = nbrNodes, useCluster = useCluster,
            dynamic = dynamic, n3 = n3,
            unitBudget = batchUnitBudget,
            dynamicMinistepFactor = dynamicMinistepFactor,
            memoryScale = memoryScale
        )
    }

    thetaHat <- nameThetaFromEffects(object[["theta"]], object$requestedEffects)
    covTheta  <- object[["covtheta"]]

    # ---- Build contribFun ----
    if (dynamic) {
        dynArgs <- list(
            ans                    = object,
            data                   = data,
            algorithm              = algorithm,
            effects                = effects,
            depvar                 = depvar,
            n3                     = n3,
            batch                  = TRUE,
            silent                 = TRUE,
            returnWide             = TRUE
        )
        contribFun <- makeContribFun("per_batch", dynArgs = dynArgs)
    } else {
        staticContribs <- getStaticChangeContributions(
            ans = object, data = data, effects = effects,
            depvar = depvar, returnWide = TRUE
        )
        contribFun <- makeContribFun("static",
            getContribFun = function(theta) staticContribs)
    }

    # ---- Build specs (N=1) ----
    metadata <- list(method = "entropy", depvar = depvar, dynamic = dynamic)
    specs <- list(R_entropy = makeSpec(
        predictFun  = predictEntropyFun,
        outcomeName = "R_entropy",
        level       = level,
        na.rm       = na.rm,
        dynamic     = dynamic,
        metadata    = metadata
    ))

    result <- sienaPostestimate(
        contribFun    = contribFun,
        nChainBatches = 1L,
        type          = "changeProb",
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
        batchDir     = batchDir,
        keepBatch    = keepBatch,
        verbose      = verbose,
        na.rm        = na.rm,
        egoNormalize = FALSE,
        uncertaintySd     = uncertaintySd,
        uncertaintyCi     = uncertaintyCi,
        uncertaintyMean   = uncertaintyMean,
        uncertaintyMedian = uncertaintyMedian,
        uncertaintyProbs  = uncertaintyProbs
    )

    result[["R_entropy"]]
}

# ---- predictFun for entropy specs -------------------------------------------
# Spec-interface wrapper: (cc, theta, baseline, ...) → data.frame with one
# R_entropy value per actor-period (or per ego, depending on grouping).
# Uses baseline$changeProb (pre-computed by makeEstimatorFun) to avoid
# redundant softmax.
predictEntropyFun <- function(changeContributions, theta, baseline, ...) {
    fullProb <- baseline$changeProb
    group_id <- changeContributions$group_id
    uGroups  <- unique(group_id)
    firstIdx <- match(uGroups, group_id)
    coordCols <- groupColsList(changeContributions, firstIdx)
    coordCols[["choice"]] <- NULL
    rh <- numeric(length(uGroups))
    for (i in seq_along(uGroups)) {
        rh[i] <- calculateEntropy(fullProb[group_id == uGroups[i]])
    }
    data.frame(coordCols, R_entropy = rh, stringsAsFactors = FALSE)
}

##@calculateEntropy Core scalar helper: degree of certainty from a probability vector
# Computes R_H = 1 - H/log(J) where H = -sum(p*log(p)) and
# J is the number of valid (non-NA, positive) probabilities.
# Returns NA if fewer than 2 valid choices.
calculateEntropy <- function(p) {
    p <- p[!is.na(p) & p > 0]
    if (length(p) >= 2L) {
        1 - (-sum(p * log(p)) / log(length(p)))
    } else {
        NA_real_
    }
}

## make this into a c++ function was well?



##@sienaEntropy Compute normalised entropy from predicted probabilities
# Compute normalised entropy (degree of certainty) from predicted probabilities.
#
# Can operate on the output of predict.sienaFit (a data.frame containing
# a changeProb column and grouping columns like period, ego),
# or directly on raw probability vectors grouped by group_id.
#
# Returns R_H: the Brillouin-style "degree of certainty"
# R_H = 1 - H / log(J) where H = -sum(p_j log p_j)
# and J is the number of non-NA choices.
#
# x:        A data.frame with a changeProb column and grouping columns
#           (typically output of predict.sienaFit), or a numeric vector of
#           probabilities.
# group_id: Integer vector of group IDs (required if x is numeric).
# level:    Aggregation level for the output (default "ego").
# sum_fun:  Aggregation function (default mean).
# Returns a data.frame with the grouping columns and an R_entropy column.
sienaEntropy <- function(x, group_id = NULL, level = "ego",
                          sum_fun = mean, na.rm = TRUE) {
    if (is.data.frame(x)) {
        if (!"changeProb" %in% names(x))
            stop("data.frame 'x' must contain a 'changeProb' column")
        # Determine group: ego-period for static, chain-ministep for dynamic
        if ("ego" %in% names(x)) {
            grp <- interaction(x$period, x$ego, drop = TRUE)
        } else if ("ministep" %in% names(x)) {
            grp <- interaction(x$chain, x$period, x$ministep, drop = TRUE)
        } else {
            stop("Cannot determine grouping from data.frame columns")
        }
        probs <- x$changeProb
    } else if (is.numeric(x)) {
        if (is.null(group_id))
            stop("'group_id' required when 'x' is a numeric vector")
        probs <- x
        grp   <- group_id
    } else {
        stop("'x' must be a data.frame (from predict.sienaFit) or a numeric vector")
    }
    # Compute per-group entropy
    uGroups    <- unique(grp)
    rh         <- setNames(numeric(length(uGroups)), uGroups)
    for (g in uGroups) {
        rh[as.character(g)] <- calculateEntropy(probs[grp == g])
    }
    data.frame(group = names(rh), R_entropy = unname(rh),
               stringsAsFactors = FALSE)
}

# Per-theta entropy from a static wide struct.
predictEntropyStatic <- function(staticContributions, theta) {
    baseline  <- predictProbability(staticContributions, theta,
                                    returnComponents = TRUE)
    fullProb  <- baseline$changeProb

    group_id <- staticContributions$group_id
    uGroups  <- unique(group_id)
    firstIdx <- match(uGroups, group_id)
    coordCols <- groupColsList(staticContributions, firstIdx)
    coordCols[["choice"]] <- NULL

    rh <- numeric(length(uGroups))
    for (i in seq_along(uGroups)) {
        rh[i] <- calculateEntropy(fullProb[group_id == uGroups[i]])
    }

    out <- data.frame(coordCols, R_entropy = rh, stringsAsFactors = FALSE)
    out
}

## the following function could be changed into rcpp for speed as well

##@entropy Generic for computing entropy-based degree of certainty
entropy <- function(object, ...) UseMethod("entropy", object)

##@entropy.default Default method: degree of certainty from a distribution vector
entropy.default <- function(object = NULL, ...) {
    calculateEntropy(object)
}