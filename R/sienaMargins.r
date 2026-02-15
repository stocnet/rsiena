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
    uncertainty = TRUE,
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
    batch_unit_budget = 5e6
){
    type <- match.arg(type)
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]

    if (is.null(memory_scale)) {
        memory_scale <- compute_memory_scale(
            data = data,
            depvar = depvar,
            dynamic = FALSE
        )
    }
    if (is.null(batch_size)) {
      auto_batch <- auto_batch_from_budget(
        data = data,
        depvar = depvar,
        nsim = nsim,
        nbrNodes = nbrNodes,
        useCluster = useCluster,
        dynamic = FALSE,
        unit_budget = batch_unit_budget,
        dynamic_ministep_factor = 1
      )
      batch_size <- auto_batch$batch_size
    }

    staticContributions <- getStaticChangeContributions(ans = ans, 
      data = data, 
      effects = effects,
      depvar = depvar,
      returnDataFrame = TRUE)
    if(!second){
        diffName <- ifelse(mainEffect == "riskDifference", 
          "firstDiff", 
          "firstRiskRatio")
        diffFun <- predictFirstDiff
        predictArgs <- list(
            ans = ans,
            staticContributions = staticContributions,
            effectName = effectName1,
            diff = diff1,
            contrast = contrast1,
            interaction = interaction1,
            int_effectNames = int_effectNames1,
            mod_effectNames = mod_effectNames1,
            type = type,
            mainEffect = mainEffect,
            details = details
        )
    }else{
        diffName <- ifelse(mainEffect == "riskDifference", "secondDiff", "secondRiskRatio")
        diffFun <- predictSecondDiff
        predictArgs <- list(
            ans = ans,
            staticContributions = staticContributions,
            effectName1 = effectName1,
            diff1 = diff1,
            contrast1 = contrast1,
            interaction1 = interaction1,
            int_effectNames1 = int_effectNames1,
            mod_effectNames1 = mod_effectNames1,
            effectName2 = effectName2,
            diff2 = diff2,
            contrast2 = contrast2,
            interaction2 = interaction2,
            int_effectNames2 = int_effectNames2,
            mod_effectNames2 = mod_effectNames2,
            type = type,
            mainEffect = mainEffect,
            details = details
        )
    }
    
    sienaPostestimate(
      predictFun = diffFun,
        predictArgs = predictArgs,
      outcome = diffName,
        level = level,
        condition = condition,
        sum.fun = sum.fun,
        na.rm = na.rm,
        theta_hat = ans$theta,
        cov_theta = ans$covtheta,
        uncertainty = uncertainty,
        nsim = nsim,
        uncertainty_sd = uncertainty_sd,
        uncertainty_ci = uncertainty_ci,
        uncertainty_probs = uncertainty_probs,
        uncertainty_mcse = uncertainty_mcse,
        uncertainty_mcse_batches = uncertainty_mcse_batches,
        useCluster = useCluster,
        nbrNodes = nbrNodes,
        clusterType = clusterType,
        cluster = cluster,
        batch_dir = batch_dir,
        prefix = prefix,
        combine_batch = combine_batch,
        batch_size = batch_size,
        keep_batch = keep_batch,
        verbose = verbose,
        memory_scale = memory_scale
    )
}

predictFirstDiff <- function(ans, theta, staticContributions,
    type = "changeProb",
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE,
    int_effectNames = NULL,
    mod_effectNames = NULL,
    details = FALSE,
    calcRiskRatio = FALSE,
    mainEffect = "riskDifference"
) {
    effects <- ans[["effects"]] # provide effects instead?
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]

  if (is.null(names(theta))) {
      # Get all effect names (including rate effects)
      allEffectNames <- ans$effects$shortName
      names(theta) <- allEffectNames
  }
  thetaNoRate <- theta[effectNames]
  df <- widenStaticContribution(staticContributions)
  df <- addUtilityColumn(df, effectNames, thetaNoRate)
  df <- addProbabilityColumn(df, group_vars = c("period", "ego"), type = type)
  
  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df <- df[df[["density"]] != 0]
  } else {
    df <- df[df[["density"]] != 0, , drop = FALSE]
  }

  fd <- calculateFirstDiff(
      densityValue = df[["density"]],
      changeProb = df[["changeProb"]],
      changeUtil = df[["changeUtil"]],
      effectName = effectName,
      effectContribution = df[[effectName]],
      diff = diff, contrast = contrast,
      interaction = interaction,
      int_effectNames = int_effectNames,
      mod_effectNames = mod_effectNames,
      modContribution = df[[mod_effectNames]],
      effectNames = effectNames,
      theta = thetaNoRate,
      type = type,
      tieProb = df[["tieProb"]],
      details = details,
      calcRiskRatio = calcRiskRatio,
      mainEffect = mainEffect
  )

  if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
    df[, (names(fd)) := fd]
  } else {
    df[names(fd)] <- fd
  }
  rm(fd)

  df <- conditionalReplace(df, df[["density"]] == -1, 
  setdiff(effectNames, "density"), function(x) x * -1)
  df
}

predictSecondDiff <- function(ans, theta, staticContributions,
    type = "changeProb",
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE,
    int_effectNames1 = NULL,
    mod_effectNames1 = NULL,
    effectName2, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE,
    int_effectNames2 = NULL,
    mod_effectNames2 = NULL,
    mainEffect = "riskDifference",
    details = FALSE
) {
    effects <- ans[["effects"]] # provide effects instead?
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]

    if (is.null(names(theta))) {
        # Get all effect names (including rate effects)
        allEffectNames <- ans$effects$shortName
        names(theta) <- allEffectNames
    }
    thetaNoRate <- theta[effectNames]

    df <- widenStaticContribution(staticContributions)
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("period", "ego"), type = type)

    if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
      df <- df[df[["density"]] != 0]
    } else {
      df <- df[df[["density"]] != 0, , drop = FALSE]
    }

    sd <- calculateSecondDiff(
        densityValue = df[["density"]], 
        changeProb = df[["changeProb"]],
        changeUtil = df[["changeUtil"]],
        effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
        effectContribution1 = df[[effectName1]],
        interaction1 = interaction1,
        int_effectNames1 = int_effectNames1,
        mod_effectNames1 = mod_effectNames1,
        modContribution1 = df[[mod_effectNames1]],
        effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
        effectContribution2 = df[[effectName2]],
        interaction2 = interaction2,
        int_effectNames2 = int_effectNames2,
        mod_effectNames2 = mod_effectNames2,
        modContribution2 = df[[mod_effectNames2]],
        effectNames = effectNames,
        theta = thetaNoRate,
        type = type, 
        tieProb = df[["tieProb"]], # if not tieProb NULL?
        details = details,
        calcRiskRatio = calcRiskRatio,
        mainEffect = mainEffect
    )
    if (requireNamespace("data.table", quietly = TRUE) && data.table::is.data.table(df)) {
      df[, (names(sd)) := sd]
    } else {
      df[names(sd)] <- sd
    }
    rm(sd)

    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
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
      "newChangeProb" = changeProb_cf, # make clear if change or tieProb
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
    return(as.data.frame(firstRiskRatio))
  } else {
    return(as.data.frame(firstDiff))
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
    diff2 <- newChangeStatistic2 - oldChangeStatistic2 # this is a vector
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
    # print(paste("firstRiskratio 1 is of length:", length(firstDiff[["firstRiskRatio"]]), " and firstRiskratio 2 is of length:", length(firstDiff2[["firstRiskRatio"]])))
      secondRiskRatio <- firstDiff2[["firstRiskRatio"]] / firstDiff[["firstRiskRatio"]]
    if (!is.null(contrast2)) {
      idx_flip <- which(newChangeStatistic2 == min(contrast2))
      secondRiskRatio[idx_flip] <- 1 / secondRiskRatio[idx_flip]
    }
  }
 # print(mainEffect)

  if(details){
    out <- data.frame(
      "changeProb_base" = changeProb,                # unchanged average probability
      "changeProb_main" = firstDiff[["newChangeProb"]],             # after main effect is changed
      "changeProb_mod" = firstDiff2[["newChangeProb"]],             # after moderator alone is changed
      "changeProb_both" = changeProb_cf21,           # after both have changed
      "firstDiff1" = firstDiff[["firstDiff"]],
      "firstDiff2" = firstDiff2[["firstDiff"]],
      "secondDiff" = secondDiff
    )
    if (type == "tieProb") {
      out$tieprob_base  <- tieProb
      out$tieprob_main  <- firstDiff[["newTieProb"]]
      out$tieprob_mod   <- firstDiff2[["newTieProb"]] # If you have tieProb_cf2, use that instead
      out$tieprob_both  <- tieProb_cf21
    }
    if (mainEffect == "riskRatio" || calcRiskRatio) {
      out$firstRiskRatio1 <- firstDiff[["firstRiskRatio"]]
      out$firstRiskRatio2 <- firstDiff2[["firstRiskRatio"]]
      out$secondRiskRatio <- secondRiskRatio
    }
    return(out)
  } else if (mainEffect == "riskRatio") {
    return(as.data.frame(secondRiskRatio))
  } else {
    return(as.data.frame(secondDiff))
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