sienaAME <- function(
    ans,
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
    uncertaintyMode = c("batch", "stream"),
    nsim = 1000,
    uncertaintySd = TRUE,
    uncertaintyCi = TRUE,
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
    mainEffect = "riskDifference", #allow utility later as well
    details = FALSE,
    memoryScale = NULL,
    batchUnitBudget = 2.5e8,
    dynamicMinistepFactor = 10
) {
    if (inherits(data, "sienaGroup"))
      stop("sienaAME does not support multi-group data (sienaGroup).")
    uncertaintyMode <- match.arg(uncertaintyMode)
    type <- match.arg(type)
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

    # ---- outcome name ----
    diffName <- if (!second) {
    ifelse(mainEffect == "riskDifference", "firstDiff", "firstRiskRatio")
    } else {
    ifelse(mainEffect == "riskDifference", "secondDiff", "secondRiskRatio")
    }

    # ---- resolve all user-supplied effect names once ----
    if (dynamic) {
      knownEffectNames <- getEffectNamesNoRate(effects, depvar)
    } else {
      staticContributions <- getStaticChangeContributions(
          ans = ans, data = data, effects = effects, depvar = depvar,
          returnWide = TRUE
      )
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

    # ---- build predictArgs ----
    if (dynamic) {
      predictArgs <- list(
          ans = ans, data = data, effects = effects,
          algorithm = algorithm, type = type, depvar = depvar,
          n3 = n3, useChangeContributions = useChangeContributions,
          mainEffect = mainEffect, details = details
      )
      if (!second) {
          diffFun <- predictFirstDiffDynamic
          predictArgs <- c(predictArgs, list(
          effectName = effectName1, diff = diff1, contrast = contrast1,
          interaction = interaction1, intEffectNames = intEffectNames1,
          modEffectNames = modEffectNames1
          ))
      } else {
          diffFun <- predictSecondDiffDynamic
          predictArgs <- c(predictArgs, list(
          effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
          interaction1 = interaction1, intEffectNames1 = intEffectNames1,
          modEffectNames1 = modEffectNames1,
          effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
          interaction2 = interaction2, intEffectNames2 = intEffectNames2,
          modEffectNames2 = modEffectNames2
          ))
      }
    } else {
      predictArgs <- list(
          ans = ans, staticContributions = staticContributions,
          type = type, mainEffect = mainEffect, details = details
      )
      if (!second) {
          diffFun <- predictFirstDiffStatic
          predictArgs <- c(predictArgs, list(
          effectName = effectName1, diff = diff1, contrast = contrast1,
          interaction = interaction1, intEffectNames = intEffectNames1,
          modEffectNames = modEffectNames1
          ))
      } else {
          diffFun <- predictSecondDiffStatic
          predictArgs <- c(predictArgs, list(
          effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
          interaction1 = interaction1, intEffectNames1 = intEffectNames1,
          modEffectNames1 = modEffectNames1,
          effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
          interaction2 = interaction2, intEffectNames2 = intEffectNames2,
          modEffectNames2 = modEffectNames2
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
    thetaHat                = ans[["theta"]],
    covTheta                = ans[["covtheta"]],
    uncertainty              = uncertainty,
    uncertaintyMode         = uncertaintyMode,
    nsim                     = nsim,
    uncertaintySd           = uncertaintySd,
    uncertaintyCi           = uncertaintyCi,
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
    useChangeContributions   = if (dynamic) useChangeContributions else NULL
    )
}

# ---------------------------------------------------------------------------
# Backward-compatibility wrapper: sienaAMEDynamic()
# ---------------------------------------------------------------------------
sienaAMEDynamic <- function(...)
{
  sienaAME(..., dynamic = TRUE)
}

# Core computation shared by both static and dynamic paths.
# `changeContributions` — unified wide struct (contribMat + coord vectors + effectNames + thetaIdx)
# `theta_use`           — named numeric vector aligned to changeContributions$effectNames
predictFirstDiff <- function(changeContributions, theta_use, type,
    effectName, diff, contrast, interaction, intEffectNames,
    modEffectNames, details, calcRiskRatio, mainEffect)
{
  mat         <- changeContributions$contribMat
  densityName <- resolveEffectName("density", changeContributions$effectNames)
  density     <- mat[, densityName]
  utility    <- calculateUtility(mat, theta_use)
  changeProb <- calculateChangeProb(utility, changeContributions$group_id)
  tieProb    <- if (type == "tieProb") calculateTieProb(changeProb, density) else NULL

  fd <- calculateFirstDiff(
    densityValue       = density,
    changeProb         = changeProb,
    changeUtil         = utility,
    effectName         = effectName,
    effectContribution = mat[, effectName],
    diff               = diff,
    contrast           = contrast,
    interaction        = interaction,
    intEffectNames     = intEffectNames,
    modEffectNames     = modEffectNames,
    modContribution    = if (!is.null(modEffectNames)) mat[, modEffectNames] else NULL,
    effectNames        = changeContributions$effectNames,
    theta              = theta_use,
    type               = type,
    tieProb            = tieProb,
    details            = details,
    calcRiskRatio      = calcRiskRatio,
    mainEffect         = mainEffect
  )

  keep   <- density != 0L
  out    <- data.frame(groupColsList(changeContributions, keep), stringsAsFactors = FALSE)
  out    <- attachContribColumns(out, changeContributions$effectNames,
                                 mat[keep, , drop = FALSE], flip = TRUE)
  out[names(fd)] <- lapply(fd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- changeProb[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
  out
}


predictSecondDiff <- function(changeContributions, theta_use, type,
    effectName1, diff1, contrast1, interaction1,
    intEffectNames1, modEffectNames1,
    effectName2, diff2, contrast2, interaction2,
    intEffectNames2, modEffectNames2,
    details, calcRiskRatio, mainEffect)
{
  mat         <- changeContributions$contribMat
  densityName <- resolveEffectName("density", changeContributions$effectNames)
  density     <- mat[, densityName]
  utility    <- calculateUtility(mat, theta_use)
  changeProb <- calculateChangeProb(utility, changeContributions$group_id)
  tieProb    <- if (type == "tieProb") calculateTieProb(changeProb, density) else NULL

  sd <- calculateSecondDiff(
    densityValue        = density,
    changeProb          = changeProb,
    changeUtil          = utility,
    effectName1         = effectName1,
    effectContribution1 = mat[, effectName1],
    diff1               = diff1,
    contrast1           = contrast1,
    interaction1        = interaction1,
    intEffectNames1     = intEffectNames1,
    modEffectNames1     = modEffectNames1,
    modContribution1    = if (!is.null(modEffectNames1)) mat[, modEffectNames1] else NULL,
    effectName2         = effectName2,
    effectContribution2 = if (!is.null(effectName2)) mat[, effectName2] else NULL,
    diff2               = diff2,
    contrast2           = contrast2,
    interaction2        = interaction2,
    intEffectNames2     = intEffectNames2,
    modEffectNames2     = modEffectNames2,
    modContribution2    = if (!is.null(modEffectNames2)) mat[, modEffectNames2] else NULL,
    effectNames         = changeContributions$effectNames,
    theta               = theta_use,
    type                = type,
    tieProb             = tieProb,
    details             = details,
    calcRiskRatio       = calcRiskRatio,
    mainEffect          = mainEffect
  )

  keep <- density != 0L
  out  <- data.frame(groupColsList(changeContributions, keep), stringsAsFactors = FALSE)
  out  <- attachContribColumns(out, changeContributions$effectNames,
                               mat[keep, , drop = FALSE], flip = TRUE)
  out[names(sd)] <- lapply(sd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- changeProb[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
  out
}


# ---- Thin wrappers (static / dynamic) ----------------------------------------

predictFirstDiffStatic <- function(ans, theta, staticContributions,
    type = "changeProb",
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference")
{
  effectNames <- staticContributions$effectNames
  theta_use   <- alignThetaNoRate(theta, effectNames, ans)[effectNames]
  predictFirstDiff(changeContributions = staticContributions, theta_use, type,
      effectName, diff, contrast, interaction, intEffectNames,
      modEffectNames, details, calcRiskRatio, mainEffect)
}

predictSecondDiffStatic <- function(ans, theta, staticContributions,
    type = "changeProb",
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    mainEffect = "riskDifference", details = FALSE)
{
  effectNames <- staticContributions$effectNames
  theta_use   <- alignThetaNoRate(theta, effectNames, ans)[effectNames]
  predictSecondDiff(changeContributions = staticContributions, theta_use, type,
      effectName1, diff1, contrast1, interaction1, intEffectNames1, modEffectNames1,
      effectName2, diff2, contrast2, interaction2, intEffectNames2, modEffectNames2,
      details, FALSE, mainEffect)
}

predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, intEffectNames = NULL, modEffectNames = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference")
{
  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynContrib <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions, returnWide = TRUE
  )
  effectNames <- dynContrib$effectNames
  theta_use <- alignThetaNoRate(theta, effectNames, ans)[effectNames]
  predictFirstDiff(changeContributions = dynContrib, theta_use, type,
      effectName, diff, contrast, interaction, intEffectNames,
      modEffectNames, details, calcRiskRatio, mainEffect)
}

predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, intEffectNames1 = NULL, modEffectNames1 = NULL,
    effectName2 = NULL, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, intEffectNames2 = NULL, modEffectNames2 = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    calcRiskRatio = FALSE, mainEffect = "riskDifference", details = FALSE)
{
  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynContrib <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions, returnWide = TRUE
  )
  effectNames <- dynContrib$effectNames
  theta_use <- alignThetaNoRate(theta, effectNames, ans)[effectNames]
  predictSecondDiff(changeContributions = dynContrib, theta_use, type,
      effectName1, diff1, contrast1, interaction1, intEffectNames1, modEffectNames1,
      effectName2, diff2, contrast2, interaction2, intEffectNames2, modEffectNames2,
      details, calcRiskRatio, mainEffect)
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
                               mainEffect = "firstDiff"){

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
      oldChangeStatistic <- densityValue * effectContribution
      newChangeStatistic <- rep(NA, length(oldChangeStatistic))
      newChangeStatistic[oldChangeStatistic == contrast[1]] <- contrast[2]
      newChangeStatistic[oldChangeStatistic == contrast[2]] <- contrast[1]
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

  expDiff <- exp(utilDiff)
  changeProb_cf <- as.vector(changeProb * expDiff / 
    (1 - changeProb + changeProb * expDiff))
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
                                calcRiskRatio = FALSE){
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
    details = details
  )

  if(!is.null(contrast2)){
    oldChangeStatistic2 <- densityValue * effectContribution2
    newChangeStatistic2 <- rep(NA, length(oldChangeStatistic2))
    newChangeStatistic2[oldChangeStatistic2 == contrast2[1]] <- contrast2[2]
    newChangeStatistic2[oldChangeStatistic2 == contrast2[2]] <- contrast2[1]
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
  expDiff21 <- exp(utilDiff21)
  changeProb_cf21 <- as.vector(changeProb * expDiff21 / (1-changeProb + changeProb * expDiff21))
  if (type == "tieProb") {
    tieProb_cf21 <- changeProb_cf21
    idx <- which(!is.na(densityValue) & densityValue == -1)
    if (length(idx) > 0) tieProb_cf21[idx] <- 1 - changeProb_cf21[idx]
  }
  changeUtil21 <- changeUtil + utilDiff21
  ## dangerous with interactions because other effect values are not "corrected"
  ## should be correct now as long as user provides them correctly and there
  ## are no higher order interactions *or* multiple interactions with same moderator

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
    details = details
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

calculateUtilityDiff <- function(effectName, diff, 
                                 theta, densityValue,
                                 interaction = FALSE,
                                 intEffectNames = NULL,
                                 modEffectNames = NULL,
                                 modContribution = NULL,
                                 effectNames = NULL){
  effectNum <- which(effectNames == effectName)
  if(interaction == TRUE){
    moderator_values <- densityValue * modContribution
    interaction_effectNums <- which(effectNames == intEffectNames)
    util_diff <- densityValue * (diff * theta[effectNum] + 
                              diff * moderator_values * theta[interaction_effectNums])
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

# Helper: resolve the effect-name arguments of sienaAME against a known column
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


# ===========================================================================
# === V1 BACKUP: sienaAME, predictFirstDiffStatic, predictSecondDiffStatic ===
# === (data.frame-internal, superseded by V2 above)                        ===
# ===========================================================================

# sienaAME <- function(
#     ans,
#     data,
#     effects = NULL,
#     depvar = NULL,
#     effectName1,
#     diff1 = NULL,
#     contrast1 = NULL,
#     interaction1 = FALSE,
#     intEffectNames1 = NULL,
#     modEffectNames1 = NULL,
#     effectName2 = NULL,
#     diff2 = NULL,
#     contrast2 = NULL,
#     interaction2 = FALSE,
#     intEffectNames2 = NULL,
#     modEffectNames2 = NULL,
#     type = c("changeProb", "tieProb"),
#     second = FALSE,
#     level = "period",
#     condition = NULL, 
#     sum_fun = mean,
#     na.rm = TRUE,
#     dynamic = FALSE,
#     algorithm = NULL,
#     n3 = 500,
#     useChangeContributions = FALSE,
#     uncertainty = TRUE,
#     uncertaintyMode = c("batch", "stream"),
#     nsim = 1000,
#     uncertaintySd = TRUE,
#     uncertaintyCi = TRUE,
#     uncertaintyProbs = c(0.025, 0.5, 0.975),
#     uncertaintyMcse = FALSE,
#     uncertaintymcseBatches = NULL,
#     useCluster = FALSE,
#     nbrNodes = 1,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batchDir = "temp",
#     prefix = "simBatch_b",
#     combineBatch = TRUE,
#     batchSize = NULL,
#     keepBatch = FALSE,
#     verbose = TRUE,
#     mainEffect = "riskDifference",
#     details = FALSE,
#     memoryScale = NULL,
#     batchUnitBudget = 5e6,
#     dynamicMinistepFactor = 10
# ){
#     uncertaintyMode <- match.arg(uncertaintyMode)
#     type <- match.arg(type)
#     if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
# 
#     if (dynamic && is.null(algorithm)) {
#         stop("'algorithm' must be provided when dynamic = TRUE")
#     }
# 
#     if (is.null(batchSize)) {
#       batchSize <- planBatch(
#         data = data, 
#         depvar = depvar, 
#         nsim = nsim,
#         nbrNodes = nbrNodes, 
#         useCluster = useCluster,
#         dynamic = dynamic, 
#         n3 = n3,
#         unitBudget = batchUnitBudget,
#         dynamicMinistepFactor = dynamicMinistepFactor,
#         memoryScale = memoryScale
#       )
#     }
#     # ---- outcome name ----
#     diffName <- if (!second) {
#         ifelse(mainEffect == "riskDifference", "firstDiff", "firstRiskRatio")
#     } else {
#         ifelse(mainEffect == "riskDifference", "secondDiff", "secondRiskRatio")
#     }
# 
#     # ---- arm-specific predictFun + predictArgs ----
#     if (dynamic) {
#         predictArgs <- list(
#             ans = ans, data = data, effects = effects,
#             algorithm = algorithm, type = type, depvar = depvar,
#             n3 = n3, useChangeContributions = useChangeContributions,
#             mainEffect = mainEffect, details = details
#         )
#         if (!second) {
#             diffFun <- predictFirstDiffDynamic
#             predictArgs <- c(predictArgs, list(
#                 effectName = effectName1, diff = diff1, contrast = contrast1,
#                 interaction = interaction1, intEffectNames = intEffectNames1,
#                 modEffectNames = modEffectNames1
#             ))
#         } else {
#             diffFun <- predictSecondDiffDynamic
#             predictArgs <- c(predictArgs, list(
#                 effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#                 interaction1 = interaction1, intEffectNames1 = intEffectNames1,
#                 modEffectNames1 = modEffectNames1,
#                 effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#                 interaction2 = interaction2, intEffectNames2 = intEffectNames2,
#                 modEffectNames2 = modEffectNames2
#             ))
#         }
#     } else {
#         staticContributions <- getStaticChangeContributions(
#             ans = ans, data = data, effects = effects,
#             depvar = depvar, returnDataFrame = TRUE
#         )
#         predictArgs <- list(
#             ans = ans, staticContributions = staticContributions,
#             type = type, mainEffect = mainEffect, details = details
#         )
#         if (!second) {
#             diffFun <- predictFirstDiffStatic
#             predictArgs <- c(predictArgs, list(
#                 effectName = effectName1, diff = diff1, contrast = contrast1,
#                 interaction = interaction1, intEffectNames = intEffectNames1,
#                 modEffectNames = modEffectNames1
#             ))
#         } else {
#             diffFun <- predictSecondDiffStatic
#             predictArgs <- c(predictArgs, list(
#                 effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#                 interaction1 = interaction1, intEffectNames1 = intEffectNames1,
#                 modEffectNames1 = modEffectNames1,
#                 effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#                 interaction2 = interaction2, intEffectNames2 = intEffectNames2,
#                 modEffectNames2 = modEffectNames2
#             ))
#         }
#     }
# 
#     sienaPostestimate(
#         predictFun = diffFun,
#         predictArgs = predictArgs,
#         outcomeName = diffName,
#         level = level,
#         condition = condition,
#         sum_fun = sum_fun,
#         na.rm = na.rm,
#         thetaHat = ans$theta,
#         covTheta = ans$covtheta,
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
#         useChangeContributions = if (dynamic)
#             useChangeContributions else NULL
#     )
# }
# 
# predictFirstDiffStatic <- function(ans, theta, staticContributions,
#     type = "changeProb",
#     effectName, diff = NULL, contrast = NULL,
#     interaction = FALSE,
#     intEffectNames = NULL,
#     modEffectNames = NULL,
#     details = FALSE,
#     calcRiskRatio = FALSE,
#     mainEffect = "riskDifference"
# ) {
#     effectNames <- getEffectNamesNoRate(ans[["effects"]])
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
#     df <- widenStaticContribution(staticContributions)
#     predictFirstDiff(df, effectNames, thetaNoRate,
#         group_vars = c("period", "ego"),
#         type = type, effectName = effectName, diff = diff,
#         contrast = contrast, interaction = interaction,
#         intEffectNames = intEffectNames,
#         modEffectNames = modEffectNames,
#         details = details, calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect)
# }
# 
# predictSecondDiffStatic <- function(ans, theta, staticContributions,
#     type = "changeProb",
#     effectName1, diff1 = NULL, contrast1 = NULL,
#     interaction1 = FALSE,
#     intEffectNames1 = NULL,
#     modEffectNames1 = NULL,
#     effectName2, diff2 = NULL, contrast2 = NULL,
#     interaction2 = FALSE,
#     intEffectNames2 = NULL,
#     modEffectNames2 = NULL,
#     mainEffect = "riskDifference",
#     details = FALSE
# ) {
#     effectNames <- getEffectNamesNoRate(ans[["effects"]])
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
#     df <- widenStaticContribution(staticContributions)
#     predictSecondDiff(df, effectNames, thetaNoRate,
#         group_vars = c("period", "ego"),
#         type = type,
#         effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#         interaction1 = interaction1,
#         intEffectNames1 = intEffectNames1,
#         modEffectNames1 = modEffectNames1,
#         effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#         interaction2 = interaction2,
#         intEffectNames2 = intEffectNames2,
#         modEffectNames2 = modEffectNames2,
#         details = details,
#         calcRiskRatio = FALSE,
#         mainEffect = mainEffect)
# }

# ===========================================================================
# === V1 BACKUP: sienaAMEDynamic (full), predictFirst/SecondDiffDynamic,   ===
# === predictFirst/SecondDiff (V1 per-row helpers), removeZeroDensity      ===
# ===========================================================================

# ##@sienaAMEDynamic backward-compat wrapper; use sienaAME(dynamic=TRUE) instead
# sienaAMEDynamic <- function(
#     ans,
#     data,
#     effectName1,
#     diff1 = NULL,
#     contrast1 = NULL,
#     interaction1 = FALSE,
#     intEffectNames1 = NULL,
#     modEffectNames1 = NULL,
#     effectName2 = NULL,
#     diff2 = NULL,
#     contrast2 = NULL,
#     interaction2 = FALSE,
#     intEffectNames2 = NULL,
#     modEffectNames2 = NULL,
#     effects,
#     algorithm,
#     type = c("changeProb", "tieProb"),
#     depvar = NULL,
#     second = FALSE,
#     level = "none",
#     condition = NULL,
#     sum_fun = mean,
#     na.rm = TRUE,
#     n3 = 500,
#     useChangeContributions = FALSE,
#     uncertainty = TRUE,
#     uncertaintyMode = c("batch", "stream"),
#     nsim = 100,
#     uncertaintySd = TRUE,
#     uncertaintyCi = TRUE,
#     uncertaintyProbs = c(0.025, 0.5, 0.975),
#     uncertaintyMcse = FALSE,
#     uncertaintymcseBatches = NULL,
#     useCluster = FALSE,
#     nbrNodes = 1,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batchDir = "temp",
#     prefix = "simBatch_b",
#     batchSize = NULL,
#     combineBatch = TRUE,
#     keepBatch = FALSE,
#     verbose = TRUE,
#     batch = TRUE,
#     silent = NULL,
#     mainEffect = "riskDifference",
#     details = FALSE,
#     memoryScale = NULL,
#     batchUnitBudget = 5e6,
#     dynamicMinistepFactor = 10
# ){
#     sienaAME(
#         ans = ans,
#         data = data,
#         effects = effects,
#         depvar = depvar,
#         effectName1 = effectName1,
#         diff1 = diff1,
#         contrast1 = contrast1,
#         interaction1 = interaction1,
#         intEffectNames1 = intEffectNames1,
#         modEffectNames1 = modEffectNames1,
#         effectName2 = effectName2,
#         diff2 = diff2,
#         contrast2 = contrast2,
#         interaction2 = interaction2,
#         intEffectNames2 = intEffectNames2,
#         modEffectNames2 = modEffectNames2,
#         type = type,
#         second = second,
#         level = level,
#         condition = condition,
#         sum_fun = sum_fun,
#         na.rm = na.rm,
#         dynamic = TRUE,
#         algorithm = algorithm,
#         n3 = n3,
#         useChangeContributions = useChangeContributions,
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
#         mainEffect = mainEffect,
#         details = details,
#         memoryScale = memoryScale,
#         batchUnitBudget = batchUnitBudget,
#         dynamicMinistepFactor = dynamicMinistepFactor
#     )
# }
# 
# predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
#     type = "changeProb", depvar = NULL,
#     effectName, diff = NULL, contrast = NULL,
#     interaction = FALSE,
#     intEffectNames = NULL,
#     modEffectNames = NULL,
#     n3 = NULL,
#     useChangeContributions = FALSE,
#     details = FALSE,
#     calcRiskRatio = FALSE,
#     mainEffect = "riskDifference"
# ){
#     if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
#     effectNames <- getEffectNamesNoRate(effects)
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
# 
#     dynamicContributions <- getDynamicChangeContributions(
#         ans = ans,
#         data = data,
#         theta = theta,
#         algorithm = algorithm,
#         effects = effects,
#         depvar = depvar,
#         n3 = n3,
#         useChangeContributions = useChangeContributions,
#         returnDataFrame = TRUE
#     )
#     dynamicContributions <- widenDynamicContribution(dynamicContributions)
# 
#     predictFirstDiff(dynamicContributions, effectNames, thetaNoRate,
#         group_vars = c("chain", "period", "ministep"),
#         type = type, effectName = effectName, diff = diff,
#         contrast = contrast, interaction = interaction,
#         intEffectNames = intEffectNames,
#         modEffectNames = modEffectNames,
#         details = details, calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect)
# }
# 
# predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
#     type = "changeProb", depvar = NULL,
#     effectName1, 
#     diff1 = NULL, contrast1 = NULL,
#     interaction1 = FALSE,
#     intEffectNames1 = NULL,
#     modEffectNames1 = NULL,
#     effectName2,
#     diff2 = NULL, contrast2 = NULL,
#     interaction2 = FALSE,
#     intEffectNames2 = NULL,
#     modEffectNames2 = NULL,
#     n3 = NULL,
#     useChangeContributions = FALSE,
#     calcRiskRatio = FALSE,
#     mainEffect = "riskDifference",
#     details = FALSE
# ){
#     if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
#     effectNames <- getEffectNamesNoRate(effects)
#     thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)
# 
#     dynamicContributions <- getDynamicChangeContributions(
#         ans = ans,
#         data = data,
#         theta = theta,
#         algorithm = algorithm,
#         effects = effects,
#         depvar = depvar,
#         n3 = n3,
#         useChangeContributions = useChangeContributions,
#         returnDataFrame = TRUE
#     )
#     dynamicContributions <- widenDynamicContribution(dynamicContributions)
# 
#     predictSecondDiff(dynamicContributions, effectNames, thetaNoRate,
#         group_vars = c("chain", "period", "ministep"),
#         type = type,
#         effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#         interaction1 = interaction1,
#         intEffectNames1 = intEffectNames1,
#         modEffectNames1 = modEffectNames1,
#         effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#         interaction2 = interaction2,
#         intEffectNames2 = intEffectNames2,
#         modEffectNames2 = modEffectNames2,
#         details = details,
#         calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect)
# }
# 
# #' Shared core: compute first difference (reuses predictProbability)
# predictFirstDiff <- function(df, effectNames, thetaNoRate, group_vars,
#     type, effectName, diff, contrast, interaction, intEffectNames,
#     modEffectNames, details, calcRiskRatio, mainEffect) {
# 
#     df <- predictProbability(df, effectNames, thetaNoRate, group_vars, type)
#     df <- removeZeroDensity(df)
# 
#     fd <- calculateFirstDiff(
#         densityValue = df[["density"]],
#         changeProb = df[["changeProb"]],
#         changeUtil = df[["changeUtil"]],
#         effectName = effectName,
#         effectContribution = df[[effectName]],
#         diff = diff, contrast = contrast,
#         interaction = interaction,
#         intEffectNames = intEffectNames,
#         modEffectNames = modEffectNames,
#         modContribution = df[[modEffectNames]],
#         effectNames = effectNames,
#         theta = thetaNoRate,
#         type = type,
#         tieProb = df[["tieProb"]],
#         details = details,
#         calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect
#     )
# 
#     if (requireNamespace("data.table", quietly = TRUE) &&
#         data.table::is.data.table(df)) {
#         df[, (names(fd)) := fd]
#     } else {
#         df[names(fd)] <- fd
#     }
# 
#     conditionalReplace(df, df[["density"]] == -1,
#         setdiff(effectNames, "density"), function(x) x * -1)
# }
# 
# #' Shared core: compute second difference (reuses predictProbability)
# predictSecondDiff <- function(df, effectNames, thetaNoRate, group_vars,
#     type,
#     effectName1, diff1, contrast1, interaction1,
#     intEffectNames1, modEffectNames1,
#     effectName2, diff2, contrast2, interaction2,
#     intEffectNames2, modEffectNames2,
#     details, calcRiskRatio, mainEffect) {
# 
#     df <- predictProbability(df, effectNames, thetaNoRate, group_vars, type)
#     df <- removeZeroDensity(df)
# 
#     sd <- calculateSecondDiff(
#         densityValue = df[["density"]],
#         changeProb = df[["changeProb"]],
#         changeUtil = df[["changeUtil"]],
#         effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#         effectContribution1 = df[[effectName1]],
#         interaction1 = interaction1,
#         intEffectNames1 = intEffectNames1,
#         modEffectNames1 = modEffectNames1,
#         modContribution1 = df[[modEffectNames1]],
#         effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#         effectContribution2 = df[[effectName2]],
#         interaction2 = interaction2,
#         intEffectNames2 = intEffectNames2,
#         modEffectNames2 = modEffectNames2,
#         modContribution2 = df[[modEffectNames2]],
#         effectNames = effectNames,
#         theta = thetaNoRate,
#         type = type,
#         tieProb = df[["tieProb"]],
#         details = details,
#         calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect
#     )
# 
#     if (requireNamespace("data.table", quietly = TRUE) &&
#         data.table::is.data.table(df)) {
#         df[, (names(sd)) := sd]
#     } else {
#         df[names(sd)] <- sd
#     }
# 
#     conditionalReplace(df, df[["density"]] == -1,
#         setdiff(effectNames, "density"), function(x) x * -1)
# }
# 
# removeZeroDensity <- function(df) {
#     if (requireNamespace("data.table", quietly = TRUE) && 
#         data.table::is.data.table(df)) {
#         df[df[["density"]] != 0]
#     } else {
#         df[df[["density"]] != 0, , drop = FALSE]
#     }
# }
