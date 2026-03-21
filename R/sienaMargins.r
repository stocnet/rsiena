##@marginalEffects Generic
marginalEffects <- function(object, ...) UseMethod("marginalEffects", object)

##@marginalEffects.sienaFit Method
marginalEffects.sienaFit <- function(
    object,
    data,
    effects = NULL,
    depvar = NULL,
    effectName1,
    diff1 = NULL,
    contrast1 = NULL,
    interaction1 = FALSE,
    intEffectNames1 = NULL,
    modEffectNames1 = NULL,
    effectName2 = NULL,
    diff2 = NULL,
    contrast2 = NULL,
    interaction2 = FALSE,
    intEffectNames2 = NULL,
    modEffectNames2 = NULL,
    type = c("changeProb", "tieProb"),
    second = FALSE,
    level = "period",
    condition = NULL,
    sum_fun = mean,
    na.rm = TRUE,
    dynamic = FALSE,
    algorithm = NULL,
    n3 = 500,
    useChangeContributions = FALSE,
    uncertainty = TRUE,
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
    uncertaintyMean = FALSE,
    uncertaintyMedian = FALSE,
    uncertaintyProbs = c(0.025, 0.5, 0.975),
    uncertaintyMcse = FALSE,
    uncertaintymcseBatches = NULL,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batchDir = "temp",
    prefix = "simBatch_b",
    combineBatch = TRUE,
    batchSize = NULL,
    keepBatch = FALSE,
    verbose = TRUE,
    mainEffect = "riskDifference", # allow utility later as well
    details = FALSE,
    perturbType1 = NULL,
    perturbType2 = NULL,
    massContrasts = NULL,
    memoryScale = NULL,
    batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10,
    egoNormalize = TRUE,
    returnDecisionDetails = FALSE,
    ...
) {
    if (inherits(data, "sienaGroup"))
      stop("marginalEffects does not support multi-group data (sienaGroup).")
    type <- match.arg(type)
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]

    # ---- behaviour DV guard ----
    dvType <- attr(data[["depvars"]][[depvar]], "type")
    if (!is.null(dvType) && dvType == "behavior")
      stop("marginalEffects is currently only implemented for network ",
           "dependent variables. Behaviour variable '", depvar,
           "' was selected.")

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

    # ---- outcome name ----
    diffName <- if (!second) {
    ifelse(mainEffect == "riskDifference", "firstDiff", "firstRiskRatio")
    } else {
    ifelse(mainEffect == "riskDifference", "secondDiff", "secondRiskRatio")
    }

    # still necessary?
    req      <- object$requestedEffects
    effects  <- if (is.null(effects)) req else effects
    # ---- resolve all user-supplied effect names once ----
    if (dynamic) {
      knownEffectNames <- getEffectNamesNoRate(effects, depvar)
    } else {
      staticContributions <- getStaticChangeContributions(
          ans = object, data = data, effects = effects, depvar = depvar,
          returnWide = TRUE
      )
      # Pre-flip to change-statistic space once (like getChangeStatistics).
      # Reused across all theta draws — avoids per-draw matrix copy.
      staticContributions$csMat <- contribToChangeStats(
          staticContributions$contribMat, staticContributions$effectNames)
      knownEffectNames <- staticContributions$effectNames
    }

    resolvedNames <- resolveAMEEffectNames(
      knownEffectNames,
      effectName1, intEffectNames1, modEffectNames1,
      effectName2, intEffectNames2, modEffectNames2,
      second
    )
    effectName1     <- resolvedNames$effectName1
    intEffectNames1 <- resolvedNames$intEffectNames1
    modEffectNames1 <- resolvedNames$modEffectNames1
    effectName2     <- resolvedNames$effectName2
    intEffectNames2 <- resolvedNames$intEffectNames2
    modEffectNames2 <- resolvedNames$modEffectNames2

    # ---- auto-center contrasts for centered covariates ----
    # Users supply contrast on the raw (uncentered) scale, e.g. c(0, 1).
    # Internally the change statistics use centered covariate values, so we
    # shift the contrast by the centering mean.
    if (!is.null(contrast1)) {
      cm1 <- getCovCenteringMean(effectName1, effects, data, depvar)
      if (cm1 != 0) contrast1 <- contrast1 - cm1
    }
    if (second && !is.null(contrast2)) {
      cm2 <- getCovCenteringMean(effectName2, effects, data, depvar)
      if (cm2 != 0) contrast2 <- contrast2 - cm2
    }

    # ---- resolve perturbation type (alter vs ego) per effect ----
    eiTypes <- if (!dynamic) staticContributions$effectInteractionTypes else NULL
    # For the dynamic path, build effectInteractionTypes from the effects object
    if (dynamic && "interactionType" %in% names(effects)) {
        noRate <- effects$type != "rate"
        inc <- effects[noRate & effects$include & effects$name == depvar, ]
        eiTypes <- setNames(inc$interactionType, knownEffectNames)
    }
    pt1 <- resolvePerturbType(effectName1, eiTypes, perturbType1)
    pt2 <- if (second) resolvePerturbType(effectName2, eiTypes, perturbType2) else "alter"

    if (is.null(massContrasts)) {
        massContrasts <- (pt1 == "ego") || (second && pt2 == "ego")
    }

    if (!is.null(condition) && massContrasts) {
      warning("'condition' is applied to 'firstDiff' only. ",
      "Mass contrasts (massCreation, massDissolution) are ego-level ",
      "quantities and are always averaged unconditionally over level = '",
      level, "'.", call. = FALSE)
    }
    # ---- build predictArgs ----
    if (dynamic) {
      predictArgs <- list(
          ans = object, data = data, effects = effects,
          algorithm = algorithm, type = type, depvar = depvar,
          n3 = n3, useChangeContributions = useChangeContributions,
          mainEffect = mainEffect, details = details,
          attachContribs = TRUE
      )
      if (!second) {
          diffFun <- predictFirstDiffDynamic
          predictArgs <- c(predictArgs, list(
          effectName = effectName1, diff = diff1, contrast = contrast1,
          interaction = interaction1, intEffectNames = intEffectNames1,
          modEffectNames = modEffectNames1,
          perturbType = pt1,
          massContrasts = massContrasts
          ))
      } else {
          diffFun <- predictSecondDiffDynamic
          predictArgs <- c(predictArgs, list(
          effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
          interaction1 = interaction1, intEffectNames1 = intEffectNames1,
          modEffectNames1 = modEffectNames1,
          effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
          interaction2 = interaction2, intEffectNames2 = intEffectNames2,
          modEffectNames2 = modEffectNames2,
          perturbType1 = pt1, perturbType2 = pt2,
          massContrasts = massContrasts
          ))
      }
    } else {
      predictArgs <- list(
          staticContributions = staticContributions,
          type = type, mainEffect = mainEffect, details = details,
          attachContribs = TRUE
      )
      if (!second) {
          diffFun <- predictFirstDiffStatic
          predictArgs <- c(predictArgs, list(
          effectName = effectName1, diff = diff1, contrast = contrast1,
          interaction = interaction1, intEffectNames = intEffectNames1,
          modEffectNames = modEffectNames1,
          perturbType = pt1,
          massContrasts = massContrasts
          ))
      } else {
          diffFun <- predictSecondDiffStatic
          predictArgs <- c(predictArgs, list(
          effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
          interaction1 = interaction1, intEffectNames1 = intEffectNames1,
          modEffectNames1 = modEffectNames1,
          effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
          interaction2 = interaction2, intEffectNames2 = intEffectNames2,
          modEffectNames2 = modEffectNames2,
          perturbType1 = pt1, perturbType2 = pt2,
          massContrasts = massContrasts
          ))
      }
    }

    sienaPostestimate(
    predictFun               = diffFun,
    predictArgs              = predictArgs,
    outcomeName              = diffName,
    level                    = level,
    condition                = condition,
    sum_fun                  = sum_fun,
    na.rm                    = na.rm,
    thetaHat                = object[["theta"]], # change to coef later
    covTheta                = object[["covtheta"]], # change to vcov later
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
      method      = "marginalEffects",
      type        = type,
      effectName1 = effectName1,
      contrast1   = contrast1,
      effectName2 = if (second) effectName2 else NULL,
      contrast2   = if (second) contrast2 else NULL,
      second      = second,
      level       = level,
      condition   = condition,
      depvar      = depvar,
      dynamic     = dynamic,
      nsim        = nsim,
      outcomeName = diffName,
      mainEffect  = mainEffect
    )
    )
}



# ---- Thin wrappers (static / dynamic) ----------------------------------------

predictFirstDiffStatic <- function(theta, staticContributions,
    type = "changeProb",
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference",
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE)
{
  effectNames <- staticContributions$effectNames
  theta_use   <- theta[effectNames]
  predictFirstDiff(changeContributions = staticContributions, theta_use, type,
      effectName, diff, contrast, interaction, intEffectNames,
      modEffectNames, details, calcRiskRatio, mainEffect,
      perturbType = perturbType, massContrasts = massContrasts,
      attachContribs = attachContribs)
}

predictSecondDiffStatic <- function(theta, staticContributions,
    type = "changeProb",
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    mainEffect = "riskDifference", details = FALSE,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE)
{
  effectNames <- staticContributions$effectNames
  theta_use   <- theta[effectNames]
  predictSecondDiff(changeContributions = staticContributions, theta_use, type,
      effectName1, diff1, contrast1, interaction1, intEffectNames1, modEffectNames1,
      effectName2, diff2, contrast2, interaction2, intEffectNames2, modEffectNames2,
      details, FALSE, mainEffect,
      perturbType1 = perturbType1, perturbType2 = perturbType2,
      massContrasts = massContrasts, attachContribs = attachContribs)
}

predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference",
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE)
{
  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynContrib <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions, returnWide = TRUE
  )
  dynContrib$csMat <- contribToChangeStats(dynContrib$contribMat,
                                            dynContrib$effectNames)
  theta_use <- theta[dynContrib$effectNames]
  predictFirstDiff(changeContributions = dynContrib, theta_use, type,
      effectName, diff, contrast, interaction, intEffectNames,
      modEffectNames, details, calcRiskRatio, mainEffect,
      perturbType = perturbType, massContrasts = massContrasts,
      attachContribs = attachContribs)
}

predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    calcRiskRatio = FALSE, mainEffect = "riskDifference", details = FALSE,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE)
{
  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynContrib <- getDynamicChangeContributions(
    theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions, returnWide = TRUE
  )
  dynContrib$csMat <- contribToChangeStats(dynContrib$contribMat,
                                            dynContrib$effectNames)
  theta_use <- theta[dynContrib$effectNames]
  predictSecondDiff(changeContributions = dynContrib, theta_use, type,
      effectName1, diff1, contrast1, interaction1, intEffectNames1, modEffectNames1,
      effectName2, diff2, contrast2, interaction2, intEffectNames2, modEffectNames2,
      details, calcRiskRatio, mainEffect,
      perturbType1 = perturbType1, perturbType2 = perturbType2,
      massContrasts = massContrasts, attachContribs = attachContribs)
}

# Core computation shared by both static and dynamic paths.
# `changeContributions` — unified wide struct (contribMat + coord vectors + effectNames)
# `theta_use`           — named numeric vector aligned to changeContributions$effectNames
predictFirstDiff <- function(changeContributions, theta_use, type,
    effectName, diff, contrast, interaction, intEffectNames,
    modEffectNames, details, calcRiskRatio, mainEffect,
    perturbType = "alter", massContrasts = FALSE, attachContribs = TRUE)
{
  contribMat         <- changeContributions$contribMat
  effectNames        <- changeContributions$effectNames
  densityName <- resolveEffectName("density", effectNames)
  density     <- contribMat[, densityName]
  # Use pre-flipped change-statistic matrix if available (static path),
  # otherwise compute it (backward compat / standalone calls).
  csMat <- if (!is.null(changeContributions$csMat)) changeContributions$csMat
           else contribToChangeStats(contribMat, effectNames)
  utility    <- calculateUtility(contribMat, theta_use)
  changeProb <- calculateChangeProb(utility, changeContributions$group_id)
  tieProb    <- if (type == "tieProb") calculateTieProb(changeProb, density) else NULL

  fd <- calculateFirstDiff(
    densityValue       = density,
    changeProb         = changeProb,
    changeUtil         = utility,
    effectName         = effectName,
    effectContribution = csMat[, effectName],
    diff               = diff,
    contrast           = contrast,
    interaction        = interaction,
    intEffectNames     = intEffectNames,
    modEffectNames     = modEffectNames,
    modContribution    = if (!is.null(modEffectNames)) csMat[, modEffectNames] else NULL,
    effectNames        = effectNames,
    theta              = theta_use,
    type               = type,
    tieProb            = tieProb,
    details            = details,
    calcRiskRatio      = calcRiskRatio,
    mainEffect         = mainEffect,
    perturbType        = perturbType,
    group_id           = changeContributions$group_id
  )

  keep <- density != 0L
  out  <- groupColsList(changeContributions, keep)
  out[names(fd)] <- lapply(fd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- changeProb[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  # Mass contrasts for ego-wide perturbations (Delta P_i^+ and Delta P_i^-)
  if (massContrasts) {
    diffColName <- intersect(c("firstDiff", "firstRiskRatio"), names(fd))[1L]
    if (!is.na(diffColName) && diffColName == "firstDiff") {
      mc <- computeMassContrasts(
        firstDiff = fd[["firstDiff"]][keep],
        density   = density[keep],
        ego       = changeContributions$ego[keep],
        period    = changeContributions$period[keep],
        group     = if (!is.null(changeContributions$group))
                      changeContributions$group[keep] else rep(1L, sum(keep)),
        type      = type
      )
      out[["massCreation"]]    <- mc[["massCreation"]]
      out[["massDissolution"]] <- mc[["massDissolution"]]
    }
  }

  if (attachContribs) {
    out <- attachContributions(out, effectNames,
                               csMat[keep, , drop = FALSE], flip = FALSE)
  }
  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}


predictSecondDiff <- function(changeContributions, theta_use, type,
    effectName1, diff1, contrast1, interaction1,
    intEffectNames1, modEffectNames1,
    effectName2, diff2, contrast2, interaction2,
    intEffectNames2, modEffectNames2,
    details, calcRiskRatio, mainEffect,
    perturbType1 = "alter", perturbType2 = "alter",
    massContrasts = FALSE, attachContribs = TRUE)
{
  contribMat         <- changeContributions$contribMat
  effectNames        <- changeContributions$effectNames
  densityName <- resolveEffectName("density", effectNames)
  density     <- contribMat[, densityName]
  # Use pre-flipped change-statistic matrix if available.
  csMat <- if (!is.null(changeContributions$csMat)) changeContributions$csMat
           else contribToChangeStats(contribMat, effectNames)
  utility    <- calculateUtility(contribMat, theta_use)
  changeProb <- calculateChangeProb(utility, changeContributions$group_id)
  tieProb    <- if (type == "tieProb") calculateTieProb(changeProb, density) else NULL

  sd <- calculateSecondDiff(
    densityValue        = density,
    changeProb          = changeProb,
    changeUtil          = utility,
    effectName1         = effectName1,
    effectContribution1 = csMat[, effectName1],
    diff1               = diff1,
    contrast1           = contrast1,
    interaction1        = interaction1,
    intEffectNames1     = intEffectNames1,
    modEffectNames1     = modEffectNames1,
    modContribution1    = if (!is.null(modEffectNames1)) csMat[, modEffectNames1] else NULL,
    effectName2         = effectName2,
    effectContribution2 = if (!is.null(effectName2)) csMat[, effectName2] else NULL,
    diff2               = diff2,
    contrast2           = contrast2,
    interaction2        = interaction2,
    intEffectNames2     = intEffectNames2,
    modEffectNames2     = modEffectNames2,
    modContribution2    = if (!is.null(modEffectNames2)) csMat[, modEffectNames2] else NULL,
    effectNames         = effectNames,
    theta               = theta_use,
    type                = type,
    tieProb             = tieProb,
    details             = details,
    calcRiskRatio       = calcRiskRatio,
    mainEffect          = mainEffect,
    perturbType1        = perturbType1,
    perturbType2        = perturbType2,
    group_id            = changeContributions$group_id
  )

  keep <- density != 0L
  out  <- groupColsList(changeContributions, keep)
  out[names(sd)] <- lapply(sd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- changeProb[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  # Mass contrasts for ego-wide perturbations
  if (massContrasts) {
    diffColName <- intersect(c("secondDiff", "secondRiskRatio"), names(sd))[1L]
    if (!is.na(diffColName) && diffColName == "secondDiff") {
      mc <- computeMassContrasts(
        firstDiff = sd[["secondDiff"]][keep],
        density   = density[keep],
        ego       = changeContributions$ego[keep],
        period    = changeContributions$period[keep],
        group     = if (!is.null(changeContributions$group))
                      changeContributions$group[keep] else rep(1L, sum(keep)),
        type      = type
      )
      out[["massCreation"]]    <- mc[["massCreation"]]
      out[["massDissolution"]] <- mc[["massDissolution"]]
    }
  }

  if (attachContribs) {
    out <- attachContributions(out, effectNames,
                               csMat[keep, , drop = FALSE], flip = FALSE)
  }
  attr(out, "row.names") <- .set_row_names(length(out[[1L]]))
  class(out) <- "data.frame"
  out
}

calculateFirstDiff <- function(densityValue,
                               changeProb,
                               changeUtil,
                               effectName,
                               effectContribution,
                               diff = NULL,
                               contrast = NULL, 
                               interaction = FALSE,
                               intEffectNames = NULL,
                               modEffectNames = NULL,
                               modContribution = NULL,
                               effectNames,
                               theta, 
                               type = "changeProb",
                               tieProb = NULL,
                               details = FALSE,
                               calcRiskRatio = FALSE,
                               mainEffect = "firstDiff",
                               perturbType = "alter",
                               group_id = NULL){

  if (effectName == "density") {
    if((!is.null(diff))) stop("firstDiff for density must be contrast c(-1,1)")
    if(!is.null(contrast)){
      if(any(setdiff(contrast, c(-1,1)))) {stop("firstDiff for density can only be be calculated for c(-1,1)")}
      if(interaction == TRUE) {stop("Interaction with density is not possible")}
      oldChangeStatistic <- densityValue
      newChangeStatistic <- rep(NA, length(oldChangeStatistic))
      newChangeStatistic[oldChangeStatistic == contrast[1]] <- contrast[2]
      newChangeStatistic[oldChangeStatistic == contrast[2]] <- contrast[1]
      utilDiff <- ifelse(is.na(newChangeStatistic), NA, -2 * changeUtil)
      densityValue <- newChangeStatistic
    }
  } else {
    if(!is.null(contrast)){
      # effectContribution is in CS space (delta) — density-independent.
      # Contrast values match delta directly for both ego and alter effects.
      oldChangeStatistic <- effectContribution
      newChangeStatistic <- rep(NA, length(oldChangeStatistic))
      newChangeStatistic[oldChangeStatistic == contrast[1]] <- contrast[2]
      newChangeStatistic[oldChangeStatistic == contrast[2]] <- contrast[1]
      if (sum(!is.na(newChangeStatistic)) == 0L)
        warning("contrast: no rows matched the supplied values ",
                "(after auto-centering, if applicable). Check that the ",
                "contrast values correspond to observed change statistics.",
                call. = FALSE)
      diff <- newChangeStatistic - oldChangeStatistic # this is a vector
    }
    utilDiff <- calculateUtilityDiff(effectName = effectName, diff = diff, 
                                      theta = theta, densityValue = densityValue,
                                      interaction = interaction,
                                      intEffectNames = intEffectNames,
                                      modEffectNames = modEffectNames,
                                      modContribution = modContribution,
                                      effectNames = effectNames)
  }

  changeProb_cf <- mlogit_update_r(changeProb, utilDiff, group_id, perturbType)
  changeProb_cf[densityValue == 0] <- NA
  if (type == "tieProb") {
    tieProb_cf <- changeProb_cf
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf[idx] <- 1 - changeProb_cf[idx]
    firstDiff <- tieProb_cf - tieProb
  } else {
     firstDiff <- changeProb_cf - changeProb
  }

  if(!is.null(contrast)){
    idx_flip <- which(newChangeStatistic == min(contrast))
    firstDiff[idx_flip] <- -firstDiff[idx_flip]
  }

  if (calcRiskRatio || mainEffect == "riskRatio") {
    if (type == "tieProb") {
      firstRiskRatio <- tieProb_cf / tieProb
    } else {
      firstRiskRatio <- changeProb_cf / changeProb
    }
    if (!is.null(contrast)) {
      idx_flip <- which(newChangeStatistic == min(contrast))
      firstRiskRatio[idx_flip] <- 1 / firstRiskRatio[idx_flip]
    }
  }

  if(details){ # mostly for debugging
    out <- data.frame(
      "firstDiff" = firstDiff,
      "utilDiff" = utilDiff,
      "newChangeProb" = changeProb_cf,
      "oldChangeProb" = changeProb
    )
    if(type == "tieProb"){
      out[["newTieProb"]] <- tieProb_cf
      out[["oldTieProb"]] <- tieProb
    }
    if(calcRiskRatio|| mainEffect == "riskRatio"){
      out[["firstRiskRatio"]] <- firstRiskRatio
    }
    return(out)
  } else if (mainEffect == "riskRatio") {
    return(list(firstRiskRatio = firstRiskRatio))
  } else {
    return(list(firstDiff = firstDiff))
  }
}

calculateSecondDiff <- function(densityValue,
                                changeProb,
                                changeUtil,
                                effectName1,
                                effectContribution1,
                                diff1 = NULL, 
                                contrast1 = NULL,
                                interaction1 = FALSE,
                                intEffectNames1 = NULL, # later it would be nice to detect these automatically
                                modEffectNames1 = NULL,
                                modContribution1 = NULL,
                                effectName2, 
                                effectContribution2,
                                diff2 = NULL,
                                contrast2 = NULL,
                                interaction2 = FALSE,
                                intEffectNames2 = NULL, # later it would be nice to detect these automatically
                                modEffectNames2 = NULL,
                                modContribution2 = NULL,
                                effectNames,
                                theta,
                                type = "changeProb",
                                tieProb = NULL,
                                details = FALSE,
                                mainEffect = "riskDifference",
                                calcRiskRatio = FALSE,
                                perturbType1 = "alter",
                                perturbType2 = "alter",
                                group_id = NULL){
  firstDiff <- calculateFirstDiff(
    densityValue = densityValue,
    changeProb = changeProb,
    changeUtil = changeUtil,
    effectName = effectName1,
    effectContribution = effectContribution1,
    diff = diff1,
    contrast = contrast1, 
    interaction = interaction1,
    intEffectNames = intEffectNames1,
    modEffectNames = modEffectNames1,
    modContribution = modContribution1,
    effectNames = effectNames,
    theta = theta, 
    type = type,
    tieProb = tieProb,
    mainEffect = mainEffect,
    details = details,
    perturbType = perturbType1,
    group_id = group_id
  )

  if(!is.null(contrast2)){
    # effectContribution2 is in CS space (delta) — density-independent.
    oldChangeStatistic2 <- effectContribution2
    newChangeStatistic2 <- rep(NA, length(oldChangeStatistic2))
    newChangeStatistic2[oldChangeStatistic2 == contrast2[1]] <- contrast2[2]
    newChangeStatistic2[oldChangeStatistic2 == contrast2[2]] <- contrast2[1]
    if (sum(!is.na(newChangeStatistic2)) == 0L)
      warning("contrast2: no rows matched the supplied values ",
              "(after auto-centering, if applicable). Check that the ",
              "contrast values correspond to observed change statistics.",
              call. = FALSE)
    diff2 <- newChangeStatistic2 - oldChangeStatistic2
  }
  utilDiff21 <- calculateUtilityDiff(effectName = effectName2,
                                      diff = diff2, theta = theta, 
                                      densityValue = densityValue,
                                      interaction = interaction2,
                                      intEffectNames = intEffectNames2,
                                      modEffectNames = modEffectNames2,
                                      modContribution = modContribution2,
                                      effectNames = effectNames)
  changeProb_cf21 <- mlogit_update_r(changeProb, utilDiff21, group_id, perturbType2)
  if (type == "tieProb") {
    tieProb_cf21 <- changeProb_cf21
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf21[idx] <- 1 - changeProb_cf21[idx]
  }
  changeUtil21 <- changeUtil + utilDiff21

  # Update moderator for effectName1's interaction after the step-2 shift:
  # if the shifted effect (effectName2) IS one of the moderators, that moderator
  # must reflect the post-shift state; otherwise the interaction contribution
  # to firstDiff2 uses stale values and the utility-level SD is missed.
  # Both diff2 and modContribution1 are in CS space — just add directly.
  if (interaction1 && !is.null(modEffectNames1)) {
    mod_shift <- if (is.null(diff2)) 1 else diff2
    mod_shift[is.na(mod_shift)] <- 0
    match_idx <- which(modEffectNames1 == effectName2)
    if (length(match_idx) > 0L) {
      if (is.matrix(modContribution1)) {
        for (mi in match_idx) modContribution1[, mi] <- modContribution1[, mi] + mod_shift
      } else {
        modContribution1 <- modContribution1 + mod_shift
      }
    }
  }

  # Calculate first difference of changing effect1 if effect2 was already changed
  firstDiff2 <- calculateFirstDiff(
    densityValue = densityValue, # careful if one was density!
    changeProb = changeProb_cf21,
    changeUtil = changeUtil21,
    effectName = effectName1, 
    diff = diff1,
    contrast = contrast1,
    effectContribution = effectContribution1,
    theta = theta, 
    type = type,
    tieProb = tieProb_cf21,
    interaction = interaction1,
    intEffectNames = intEffectNames1,
    modEffectNames = modEffectNames1,
    modContribution = modContribution1,
    effectNames = effectNames,
    mainEffect = mainEffect,
    details = details,
    perturbType = perturbType1,
    group_id = group_id
  )

  secondDiff <- firstDiff2[["firstDiff"]] - firstDiff[["firstDiff"]]
  if(!is.null(contrast2)){
    secondDiff[which(newChangeStatistic2 == min(contrast2))] <- -secondDiff[which(newChangeStatistic2 == min(contrast2))]
  }

  if (mainEffect == "riskRatio") {
      secondRiskRatio <- firstDiff2[["firstRiskRatio"]] / firstDiff[["firstRiskRatio"]]
    if (!is.null(contrast2)) {
      idx_flip <- which(newChangeStatistic2 == min(contrast2))
      secondRiskRatio[idx_flip] <- 1 / secondRiskRatio[idx_flip]
    }
  }
  ## really use data.frame here?
  if(details){
    out <- data.frame(
      "changeProb_base" = changeProb,        
      "changeProb_main" = firstDiff[["newChangeProb"]],        
      "changeProb_mod" = changeProb_cf21,       
      "changeProb_both" = firstDiff2[["newChangeProb"]],         
      "firstDiff1" = firstDiff[["firstDiff"]],
      "firstDiff2" = firstDiff2[["firstDiff"]],
      "secondDiff" = secondDiff
    )
    if (type == "tieProb") {
      out$tieProb_base  <- tieProb
      out$tieProb_main  <- firstDiff[["newTieProb"]]
      out$tieProb_mod   <- firstDiff2[["newTieProb"]]
      out$tieProb_both  <- tieProb_cf21
    }
    if (mainEffect == "riskRatio" || calcRiskRatio) {
      out$firstRiskRatio1 <- firstDiff[["firstRiskRatio"]]
      out$firstRiskRatio2 <- firstDiff2[["firstRiskRatio"]]
      out$secondRiskRatio <- secondRiskRatio
    }
    return(out)
  } else if (mainEffect == "riskRatio") {
    return(list(secondRiskRatio = secondRiskRatio))
  } else {
    return(list(secondDiff = secondDiff))
  }
}

# Compute the utility shift from perturbing effectName by diff.
# All inputs (diff, modContribution) are in change-statistic space (delta).
# The utility shift is  d * diff * (theta_e + sum_k mod_k * theta_int_k).
calculateUtilityDiff <- function(effectName, diff, 
                                 theta, densityValue,
                                 interaction = FALSE,
                                 intEffectNames = NULL,
                                 modEffectNames = NULL,
                                 modContribution = NULL,
                                 effectNames = NULL){
  if (is.null(diff)) diff <- 1  # NULL means "+1 unit" perturbation
  effectNum <- which(effectNames == effectName)
  if(interaction == TRUE){
    if (is.null(intEffectNames))
      stop("'intEffectNames' must not be NULL when interaction = TRUE.")
    if (is.null(modEffectNames))
      stop("'modEffectNames' must not be NULL when interaction = TRUE.")
    # d * diff * (theta_e + sum_k mod_k * theta_int_k)
    inner <- theta[effectNum]
    nInt <- length(intEffectNames)
    for (k in seq_len(nInt)) {
      mod_k <- if (is.matrix(modContribution)) modContribution[, k] else modContribution
      int_num <- which(effectNames == intEffectNames[k])
      inner <- inner + mod_k * theta[int_num]
    }
    util_diff <- densityValue * diff * inner
  } else {
    util_diff <- densityValue * diff * theta[effectNum]
  }
  util_diff
 }

# Resolve a user-supplied bare effect shortName to the composite "shortName_type"
# column name used in contribMat. If it already contains "_" it passes through.
# If multiple type-variants exist, "_eval" is preferred.
resolveEffectName <- function(effectName, effectNames) {
  if (is.null(effectName)) return(NULL)
  vapply(effectName, function(nm) {
    if (nm %in% effectNames) return(nm)
    # If nm already contains '_' (composite shortName_type or full name), try
    # matching it as a depvar-prefixed suffix: e.g. "recip_eval" → "mynet_recip_eval"
    if (grepl("_", nm, fixed = TRUE)) {
      m <- grep(paste0("_", nm, "$"), effectNames, perl = TRUE, value = TRUE)
      if (length(m) > 0L) return(m[1L])
    }
    # Try shortName_eval (most common default type)
    m <- grep(paste0("(^|_)", nm, "_eval$"), effectNames, perl = TRUE, value = TRUE)
    if (length(m) > 0L) return(m[1L])
    # Try any type suffix
    m <- grep(paste0("(^|_)", nm, "_"), effectNames, perl = TRUE, value = TRUE)
    if (length(m) > 0L) return(m[1L])
    stop("Effect '", nm, "' not found in contribMat columns: ",
         paste(effectNames, collapse = ", "), call. = FALSE)
  }, character(1L), USE.NAMES = FALSE)
}

# Helper: resolve the effect-name arguments of marginalEffects against a known column
# list.  Returns a named list with the same structure used in predictArgs.
resolveAMEEffectNames <- function(effectNames,
                                  effectName1, intEffectNames1, modEffectNames1,
                                  effectName2, intEffectNames2, modEffectNames2,
                                  second) {
  list(
    effectName1     = resolveEffectName(effectName1,      effectNames),
    intEffectNames1 = resolveEffectName(intEffectNames1,  effectNames),
    modEffectNames1 = resolveEffectName(modEffectNames1,  effectNames),
    effectName2     = if (second) resolveEffectName(effectName2,     effectNames) else NULL,
    intEffectNames2 = if (second) resolveEffectName(intEffectNames2, effectNames) else NULL,
    modEffectNames2 = if (second) resolveEffectName(modEffectNames2, effectNames) else NULL
  )
}

