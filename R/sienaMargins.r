sienaAME <- function(
    ans,
    data,
    effects = NULL,
    depvar = NULL,
    effectName1,
    diff1 = NULL,
    contrast1 = NULL,
    interaction1 = FALSE,
    int_effectNames1 = NULL,
    mod_effectNames1 = NULL,
    effectName2 = NULL,
    diff2 = NULL,
    contrast2 = NULL,
    interaction2 = FALSE,
    int_effectNames2 = NULL,
    mod_effectNames2 = NULL,
    type = c("changeProb", "tieProb"),
    second = FALSE,
    level = "period",
    condition = NULL,
    sum.fun = mean,
    na.rm = TRUE,
    dynamic = FALSE,
    algorithm = NULL,
    n3 = 500,
    useChangeContributions = FALSE,
    uncertainty = TRUE,
    uncertainty_mode = c("batch", "stream"),
    nsim = 1000,
    uncertainty_sd = TRUE,
    uncertainty_ci = TRUE,
    uncertainty_probs = c(0.025, 0.5, 0.975),
    uncertainty_mcse = FALSE,
    uncertainty_mcse_batches = NULL,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch_size = NULL,
    keep_batch = FALSE,
    verbose = TRUE,
    mainEffect = "riskDifference",
    details = FALSE,
    memory_scale = NULL,
    batch_unit_budget = 2.5e8,
    dynamic_ministep_factor = 10
) {
    if (inherits(data, "sienaGroup"))
      stop("sienaAME does not support multi-group data (sienaGroup).")
    uncertainty_mode <- match.arg(uncertainty_mode)
    type <- match.arg(type)
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]

    if (dynamic && is.null(algorithm))
      stop("'algorithm' must be provided when dynamic = TRUE")

    if (is.null(batch_size)) {
        batch_size <- plan_batch(
        data = data, depvar = depvar, nsim = nsim,
        nbrNodes = nbrNodes, useCluster = useCluster,
        dynamic = dynamic, n3 = n3,
        unit_budget = batch_unit_budget,
        dynamic_ministep_factor = dynamic_ministep_factor,
        memory_scale = memory_scale
        )
    }

    # ---- outcome name ----
    diffName <- if (!second) {
    ifelse(mainEffect == "riskDifference", "firstDiff", "firstRiskRatio")
    } else {
    ifelse(mainEffect == "riskDifference", "secondDiff", "secondRiskRatio")
    }

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
          interaction = interaction1, int_effectNames = int_effectNames1,
          mod_effectNames = mod_effectNames1
          ))
      } else {
          diffFun <- predictSecondDiffDynamic
          predictArgs <- c(predictArgs, list(
          effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
          interaction1 = interaction1, int_effectNames1 = int_effectNames1,
          mod_effectNames1 = mod_effectNames1,
          effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
          interaction2 = interaction2, int_effectNames2 = int_effectNames2,
          mod_effectNames2 = mod_effectNames2
          ))
      }
    } else {
      staticContributions <- getStaticChangeContributions(
          ans = ans, data = data, effects = effects, depvar = depvar,
          returnWide = TRUE
      )
      predictArgs <- list(
          ans = ans, staticContributions = staticContributions,
          type = type, mainEffect = mainEffect, details = details
      )
      if (!second) {
          diffFun <- predictFirstDiffStatic
          predictArgs <- c(predictArgs, list(
          effectName = effectName1, diff = diff1, contrast = contrast1,
          interaction = interaction1, int_effectNames = int_effectNames1,
          mod_effectNames = mod_effectNames1
          ))
      } else {
          diffFun <- predictSecondDiffStatic
          predictArgs <- c(predictArgs, list(
          effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
          interaction1 = interaction1, int_effectNames1 = int_effectNames1,
          mod_effectNames1 = mod_effectNames1,
          effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
          interaction2 = interaction2, int_effectNames2 = int_effectNames2,
          mod_effectNames2 = mod_effectNames2
          ))
      }
    }

    sienaPostestimate(
    predictFun               = diffFun,
    predictArgs              = predictArgs,
    outcome                  = diffName,
    level                    = level,
    condition                = condition,
    sum.fun                  = sum.fun,
    na.rm                    = na.rm,
    theta_hat                = ans[["theta"]],
    cov_theta                = ans[["covtheta"]],
    uncertainty              = uncertainty,
    uncertainty_mode         = uncertainty_mode,
    nsim                     = nsim,
    uncertainty_sd           = uncertainty_sd,
    uncertainty_ci           = uncertainty_ci,
    uncertainty_probs        = uncertainty_probs,
    uncertainty_mcse         = uncertainty_mcse,
    uncertainty_mcse_batches = uncertainty_mcse_batches,
    useCluster               = useCluster,
    nbrNodes                 = nbrNodes,
    clusterType              = clusterType,
    cluster                  = cluster,
    batch_dir                = batch_dir,
    prefix                   = prefix,
    combine_batch            = combine_batch,
    batch_size               = batch_size,
    keep_batch               = keep_batch,
    verbose                  = verbose,
    useChangeContributions   = if (dynamic) useChangeContributions else NULL
    )
}

# ---------------------------------------------------------------------------
# Backward-compatibility wrapper: sienaAMEDynamic()
# ---------------------------------------------------------------------------
sienaAMEDynamic <- function(...) {
  sienaAME(..., dynamic = TRUE)
}


predictFirstDiff <- function(changeContributions, theta_use, type,
    effectName, diff, contrast, interaction, int_effectNames,
    mod_effectNames, details, calcRiskRatio, mainEffect) {

  utility     <- calculateUtility(changeContributions$contrib_mat, theta_use)
  prob        <- as.numeric(softmax_arma_by_group(utility, changeContributions$group_id))
  density_col <- changeContributions$contrib_mat[, "density"]
  tieProb     <- .computeTieProb(prob, density_col, type)

  fd <- calculateFirstDiff(
    densityValue       = density_col,
    changeProb         = prob,
    changeUtil         = utility,
    effectName         = effectName,
    effectContribution = changeContributions$contrib_mat[, effectName],
    diff               = diff,
    contrast           = contrast,
    interaction        = interaction,
    int_effectNames    = int_effectNames,
    mod_effectNames    = mod_effectNames,
    modContribution    = if (!is.null(mod_effectNames)) changeContributions$contrib_mat[, mod_effectNames] else NULL,
    effectNames        = changeContributions$effectNames,
    theta              = theta_use,
    type               = type,
    tieProb            = tieProb,
    details            = details,
    calcRiskRatio      = calcRiskRatio,
    mainEffect         = mainEffect
  )

  # Exclude density==0 (no-change) rows from output — calculateFirstDiff already sets them to NA
  keep   <- density_col != 0L
  contrib <- changeContributions$contrib_mat[keep, , drop = FALSE]
  out <- data.frame(.groupColsList(changeContributions, keep), stringsAsFactors = FALSE)
  out <- .attachContribColumns(out, changeContributions$effectNames, contrib, flip = TRUE)
  out[names(fd)] <- lapply(fd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- prob[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
  out
}


predictSecondDiff <- function(changeContributions, theta_use, type,
    effectName1, diff1, contrast1, interaction1,
    int_effectNames1, mod_effectNames1,
    effectName2, diff2, contrast2, interaction2,
    int_effectNames2, mod_effectNames2,
    details, calcRiskRatio, mainEffect) {

  utility     <- calculateUtility(changeContributions$contrib_mat, theta_use)
  prob        <- as.numeric(softmax_arma_by_group(utility, changeContributions$group_id))
  density_col <- changeContributions$contrib_mat[, "density"]
  tieProb     <- .computeTieProb(prob, density_col, type)

  sd <- calculateSecondDiff(
    densityValue        = density_col,
    changeProb          = prob,
    changeUtil          = utility,
    effectName1         = effectName1,
    effectContribution1 = changeContributions$contrib_mat[, effectName1],
    diff1               = diff1,
    contrast1           = contrast1,
    interaction1        = interaction1,
    int_effectNames1    = int_effectNames1,
    mod_effectNames1    = mod_effectNames1,
    modContribution1    = if (!is.null(mod_effectNames1)) changeContributions$contrib_mat[, mod_effectNames1] else NULL,
    effectName2         = effectName2,
    effectContribution2 = changeContributions$contrib_mat[, effectName2],
    diff2               = diff2,
    contrast2           = contrast2,
    interaction2        = interaction2,
    int_effectNames2    = int_effectNames2,
    mod_effectNames2    = mod_effectNames2,
    modContribution2    = if (!is.null(mod_effectNames2)) changeContributions$contrib_mat[, mod_effectNames2] else NULL,
    effectNames         = changeContributions$effectNames,
    theta               = theta_use,
    type                = type,
    tieProb             = tieProb,
    details             = details,
    calcRiskRatio       = calcRiskRatio,
    mainEffect          = mainEffect
  )

  keep   <- density_col != 0L
  contrib <- changeContributions$contrib_mat[keep, , drop = FALSE]
  out <- data.frame(.groupColsList(changeContributions, keep), stringsAsFactors = FALSE)
  out <- .attachContribColumns(out, changeContributions$effectNames, contrib, flip = TRUE)
  out[names(sd)] <- lapply(sd, `[`, keep)

  if (details) {
    out[["changeUtil"]] <- utility[keep]
    out[["changeProb"]] <- prob[keep]
    if (!is.null(tieProb)) out[["tieProb"]] <- tieProb[keep]
  }

  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDT(out)
  out
}


# ---- Static per-theta wrappers -------------------------------------------

predictFirstDiffStatic <- function(ans, theta, staticContributions,
    type = "changeProb",
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference") {

  effectNames <- staticContributions[[1]]$effectNames
  thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)

  results <- vector("list", length(staticContributions))
  for (i in seq_along(staticContributions)) {
    staticCont <- staticContributions[[i]]
    theta_use <- thetaNoRate[staticCont$effectNames]
    results[[i]] <- predictFirstDiff(staticCont, theta_use, type,
        effectName, diff, contrast, interaction, int_effectNames,
        mod_effectNames, details, calcRiskRatio, mainEffect)
  }
  if (requireNamespace("data.table", quietly = TRUE))
    data.table::rbindlist(results, use.names = TRUE)
  else {
    out <- do.call(rbind, results)
    rownames(out) <- NULL
    out
  }
}

predictSecondDiffStatic <- function(ans, theta, staticContributions,
    type = "changeProb",
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, int_effectNames1 = NULL, mod_effectNames1 = NULL,
    effectName2, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, int_effectNames2 = NULL, mod_effectNames2 = NULL,
    mainEffect = "riskDifference", details = FALSE) {

  effectNames <- staticContributions[[1]]$effectNames
  thetaNoRate <- alignThetaNoRate(theta, effectNames, ans)

  results <- vector("list", length(staticContributions))
  for (i in seq_along(staticContributions)) {
    staticCont <- staticContributions[[i]]
    theta_use <- thetaNoRate[staticCont$effectNames]
    results[[i]] <- predictSecondDiff(staticCont, theta_use, type,
        effectName1, diff1, contrast1, interaction1,
        int_effectNames1, mod_effectNames1,
        effectName2, diff2, contrast2, interaction2,
        int_effectNames2, mod_effectNames2,
        details, FALSE, mainEffect)
  }
  if (requireNamespace("data.table", quietly = TRUE))
    data.table::rbindlist(results, use.names = TRUE)
  else {
    out <- do.call(rbind, results)
    rownames(out) <- NULL
    out
  }
}


# ---- Dynamic per-theta wrappers ------------------------------------------

predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE, int_effectNames = NULL, mod_effectNames = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    details = FALSE, calcRiskRatio = FALSE, mainEffect = "riskDifference") {

  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynamicChangeContributions <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions,
    returnWide = TRUE
  )
  thetaNoRate <- alignThetaNoRate(theta, getEffectNamesNoRate(effects, depvar), ans)
  theta_use   <- thetaNoRate[dynamicChangeContributions$effectNames]

  predictFirstDiff(dynamicChangeContributions, theta_use, type,
      effectName, diff, contrast, interaction, int_effectNames,
      mod_effectNames, details, calcRiskRatio, mainEffect)
}

predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
    type = "changeProb", depvar = NULL,
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE, int_effectNames1 = NULL, mod_effectNames1 = NULL,
    effectName2, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE, int_effectNames2 = NULL, mod_effectNames2 = NULL,
    n3 = NULL, useChangeContributions = FALSE,
    calcRiskRatio = FALSE, mainEffect = "riskDifference", details = FALSE) {

  if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
  dynamicChangeContributions <- getDynamicChangeContributions(
    ans = ans, theta = theta, data = data, algorithm = algorithm,
    effects = effects, depvar = depvar, n3 = n3,
    useChangeContributions = useChangeContributions,
    returnWide = TRUE
  )
  thetaNoRate <- alignThetaNoRate(theta, getEffectNamesNoRate(effects, depvar), ans)
  theta_use   <- thetaNoRate[dynamicChangeContributions$effectNames]

  predictSecondDiff(dynamicChangeContributions, theta_use, type,
      effectName1, diff1, contrast1, interaction1,
      int_effectNames1, mod_effectNames1,
      effectName2, diff2, contrast2, interaction2,
      int_effectNames2, mod_effectNames2,
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
                               int_effectNames = NULL,
                               mod_effectNames = NULL,
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
                                      int_effectNames = int_effectNames,
                                      mod_effectNames = mod_effectNames,
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
                                int_effectNames1 = NULL, # later it would be nice to detect these automatically
                                mod_effectNames1 = NULL,
                                modContribution1 = NULL,
                                effectName2, 
                                effectContribution2,
                                diff2 = NULL,
                                contrast2 = NULL,
                                interaction2 = FALSE,
                                int_effectNames2 = NULL, # later it would be nice to detect these automatically
                                mod_effectNames2 = NULL,
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
    int_effectNames = int_effectNames1,
    mod_effectNames = mod_effectNames1,
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
                                      int_effectNames = int_effectNames2,
                                      mod_effectNames = mod_effectNames2,
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
    int_effectNames = int_effectNames1,
    mod_effectNames = mod_effectNames1,
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
                                 int_effectNames = NULL,
                                 mod_effectNames = NULL,
                                 modContribution = NULL,
                                 effectNames = NULL){
  effectNum <- which(effectNames == effectName)
  if(interaction == TRUE){
    moderator_values <- densityValue * modContribution
    interaction_effectNums <- which(effectNames == int_effectNames)
    util_diff <- densityValue * (diff * theta[effectNum] + 
                              diff * moderator_values * theta[interaction_effectNums])
  } else {
    util_diff <- densityValue * diff * theta[effectNum]
  }
  util_diff
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
#     int_effectNames1 = NULL,
#     mod_effectNames1 = NULL,
#     effectName2 = NULL,
#     diff2 = NULL,
#     contrast2 = NULL,
#     interaction2 = FALSE,
#     int_effectNames2 = NULL,
#     mod_effectNames2 = NULL,
#     type = c("changeProb", "tieProb"),
#     second = FALSE,
#     level = "period",
#     condition = NULL, 
#     sum.fun = mean,
#     na.rm = TRUE,
#     dynamic = FALSE,
#     algorithm = NULL,
#     n3 = 500,
#     useChangeContributions = FALSE,
#     uncertainty = TRUE,
#     uncertainty_mode = c("batch", "stream"),
#     nsim = 1000,
#     uncertainty_sd = TRUE,
#     uncertainty_ci = TRUE,
#     uncertainty_probs = c(0.025, 0.5, 0.975),
#     uncertainty_mcse = FALSE,
#     uncertainty_mcse_batches = NULL,
#     useCluster = FALSE,
#     nbrNodes = 1,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batch_dir = "temp",
#     prefix = "simBatch_b",
#     combine_batch = TRUE,
#     batch_size = NULL,
#     keep_batch = FALSE,
#     verbose = TRUE,
#     mainEffect = "riskDifference",
#     details = FALSE,
#     memory_scale = NULL,
#     batch_unit_budget = 5e6,
#     dynamic_ministep_factor = 10
# ){
#     uncertainty_mode <- match.arg(uncertainty_mode)
#     type <- match.arg(type)
#     if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
# 
#     if (dynamic && is.null(algorithm)) {
#         stop("'algorithm' must be provided when dynamic = TRUE")
#     }
# 
#     if (is.null(batch_size)) {
#       batch_size <- plan_batch(
#         data = data, 
#         depvar = depvar, 
#         nsim = nsim,
#         nbrNodes = nbrNodes, 
#         useCluster = useCluster,
#         dynamic = dynamic, 
#         n3 = n3,
#         unit_budget = batch_unit_budget,
#         dynamic_ministep_factor = dynamic_ministep_factor,
#         memory_scale = memory_scale
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
#                 interaction = interaction1, int_effectNames = int_effectNames1,
#                 mod_effectNames = mod_effectNames1
#             ))
#         } else {
#             diffFun <- predictSecondDiffDynamic
#             predictArgs <- c(predictArgs, list(
#                 effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#                 interaction1 = interaction1, int_effectNames1 = int_effectNames1,
#                 mod_effectNames1 = mod_effectNames1,
#                 effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#                 interaction2 = interaction2, int_effectNames2 = int_effectNames2,
#                 mod_effectNames2 = mod_effectNames2
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
#                 interaction = interaction1, int_effectNames = int_effectNames1,
#                 mod_effectNames = mod_effectNames1
#             ))
#         } else {
#             diffFun <- predictSecondDiffStatic
#             predictArgs <- c(predictArgs, list(
#                 effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
#                 interaction1 = interaction1, int_effectNames1 = int_effectNames1,
#                 mod_effectNames1 = mod_effectNames1,
#                 effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#                 interaction2 = interaction2, int_effectNames2 = int_effectNames2,
#                 mod_effectNames2 = mod_effectNames2
#             ))
#         }
#     }
# 
#     sienaPostestimate(
#         predictFun = diffFun,
#         predictArgs = predictArgs,
#         outcome = diffName,
#         level = level,
#         condition = condition,
#         sum.fun = sum.fun,
#         na.rm = na.rm,
#         theta_hat = ans$theta,
#         cov_theta = ans$covtheta,
#         uncertainty = uncertainty,
#         uncertainty_mode = uncertainty_mode,
#         nsim = nsim,
#         uncertainty_sd = uncertainty_sd,
#         uncertainty_ci = uncertainty_ci,
#         uncertainty_probs = uncertainty_probs,
#         uncertainty_mcse = uncertainty_mcse,
#         uncertainty_mcse_batches = uncertainty_mcse_batches,
#         useCluster = useCluster,
#         nbrNodes = nbrNodes,
#         clusterType = clusterType,
#         cluster = cluster,
#         batch_dir = batch_dir,
#         prefix = prefix,
#         combine_batch = combine_batch,
#         batch_size = batch_size,
#         keep_batch = keep_batch,
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
#     int_effectNames = NULL,
#     mod_effectNames = NULL,
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
#         int_effectNames = int_effectNames,
#         mod_effectNames = mod_effectNames,
#         details = details, calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect)
# }
# 
# predictSecondDiffStatic <- function(ans, theta, staticContributions,
#     type = "changeProb",
#     effectName1, diff1 = NULL, contrast1 = NULL,
#     interaction1 = FALSE,
#     int_effectNames1 = NULL,
#     mod_effectNames1 = NULL,
#     effectName2, diff2 = NULL, contrast2 = NULL,
#     interaction2 = FALSE,
#     int_effectNames2 = NULL,
#     mod_effectNames2 = NULL,
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
#         int_effectNames1 = int_effectNames1,
#         mod_effectNames1 = mod_effectNames1,
#         effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#         interaction2 = interaction2,
#         int_effectNames2 = int_effectNames2,
#         mod_effectNames2 = mod_effectNames2,
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
#     int_effectNames1 = NULL,
#     mod_effectNames1 = NULL,
#     effectName2 = NULL,
#     diff2 = NULL,
#     contrast2 = NULL,
#     interaction2 = FALSE,
#     int_effectNames2 = NULL,
#     mod_effectNames2 = NULL,
#     effects,
#     algorithm,
#     type = c("changeProb", "tieProb"),
#     depvar = NULL,
#     second = FALSE,
#     level = "none",
#     condition = NULL,
#     sum.fun = mean,
#     na.rm = TRUE,
#     n3 = 500,
#     useChangeContributions = FALSE,
#     uncertainty = TRUE,
#     uncertainty_mode = c("batch", "stream"),
#     nsim = 100,
#     uncertainty_sd = TRUE,
#     uncertainty_ci = TRUE,
#     uncertainty_probs = c(0.025, 0.5, 0.975),
#     uncertainty_mcse = FALSE,
#     uncertainty_mcse_batches = NULL,
#     useCluster = FALSE,
#     nbrNodes = 1,
#     clusterType = c("PSOCK", "FORK"),
#     cluster = NULL,
#     batch_dir = "temp",
#     prefix = "simBatch_b",
#     batch_size = NULL,
#     combine_batch = TRUE,
#     keep_batch = FALSE,
#     verbose = TRUE,
#     batch = TRUE,
#     silent = NULL,
#     mainEffect = "riskDifference",
#     details = FALSE,
#     memory_scale = NULL,
#     batch_unit_budget = 5e6,
#     dynamic_ministep_factor = 10
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
#         int_effectNames1 = int_effectNames1,
#         mod_effectNames1 = mod_effectNames1,
#         effectName2 = effectName2,
#         diff2 = diff2,
#         contrast2 = contrast2,
#         interaction2 = interaction2,
#         int_effectNames2 = int_effectNames2,
#         mod_effectNames2 = mod_effectNames2,
#         type = type,
#         second = second,
#         level = level,
#         condition = condition,
#         sum.fun = sum.fun,
#         na.rm = na.rm,
#         dynamic = TRUE,
#         algorithm = algorithm,
#         n3 = n3,
#         useChangeContributions = useChangeContributions,
#         uncertainty = uncertainty,
#         uncertainty_mode = uncertainty_mode,
#         nsim = nsim,
#         uncertainty_sd = uncertainty_sd,
#         uncertainty_ci = uncertainty_ci,
#         uncertainty_probs = uncertainty_probs,
#         uncertainty_mcse = uncertainty_mcse,
#         uncertainty_mcse_batches = uncertainty_mcse_batches,
#         useCluster = useCluster,
#         nbrNodes = nbrNodes,
#         clusterType = clusterType,
#         cluster = cluster,
#         batch_dir = batch_dir,
#         prefix = prefix,
#         combine_batch = combine_batch,
#         batch_size = batch_size,
#         keep_batch = keep_batch,
#         verbose = verbose,
#         mainEffect = mainEffect,
#         details = details,
#         memory_scale = memory_scale,
#         batch_unit_budget = batch_unit_budget,
#         dynamic_ministep_factor = dynamic_ministep_factor
#     )
# }
# 
# predictFirstDiffDynamic <- function(ans, data, theta, effects, algorithm,
#     type = "changeProb", depvar = NULL,
#     effectName, diff = NULL, contrast = NULL,
#     interaction = FALSE,
#     int_effectNames = NULL,
#     mod_effectNames = NULL,
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
#         int_effectNames = int_effectNames,
#         mod_effectNames = mod_effectNames,
#         details = details, calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect)
# }
# 
# predictSecondDiffDynamic <- function(ans, data, theta, effects, algorithm,
#     type = "changeProb", depvar = NULL,
#     effectName1, 
#     diff1 = NULL, contrast1 = NULL,
#     interaction1 = FALSE,
#     int_effectNames1 = NULL,
#     mod_effectNames1 = NULL,
#     effectName2,
#     diff2 = NULL, contrast2 = NULL,
#     interaction2 = FALSE,
#     int_effectNames2 = NULL,
#     mod_effectNames2 = NULL,
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
#         int_effectNames1 = int_effectNames1,
#         mod_effectNames1 = mod_effectNames1,
#         effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#         interaction2 = interaction2,
#         int_effectNames2 = int_effectNames2,
#         mod_effectNames2 = mod_effectNames2,
#         details = details,
#         calcRiskRatio = calcRiskRatio,
#         mainEffect = mainEffect)
# }
# 
# #' Shared core: compute first difference (reuses predictProbability)
# predictFirstDiff <- function(df, effectNames, thetaNoRate, group_vars,
#     type, effectName, diff, contrast, interaction, int_effectNames,
#     mod_effectNames, details, calcRiskRatio, mainEffect) {
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
#         int_effectNames = int_effectNames,
#         mod_effectNames = mod_effectNames,
#         modContribution = df[[mod_effectNames]],
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
#     int_effectNames1, mod_effectNames1,
#     effectName2, diff2, contrast2, interaction2,
#     int_effectNames2, mod_effectNames2,
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
#         int_effectNames1 = int_effectNames1,
#         mod_effectNames1 = mod_effectNames1,
#         modContribution1 = df[[mod_effectNames1]],
#         effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
#         effectContribution2 = df[[effectName2]],
#         interaction2 = interaction2,
#         int_effectNames2 = int_effectNames2,
#         mod_effectNames2 = mod_effectNames2,
#         modContribution2 = df[[mod_effectNames2]],
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
