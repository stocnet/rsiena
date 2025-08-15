sienaAME <- function(
    ans,
    data,
    effectName1,
    diff1 = NULL,
    contrast1 = NULL,
    interaction1 = FALSE, # clashes with what interaction1 means in siena
    int_effectNames1 = NULL,
    mod_effectNames1 = NULL,
    effectName2 = NULL,
    diff2 = NULL,
    contrast2 = NULL,
    interaction2 = FALSE,
    int_effectNames2 = NULL,
    mod_effectNames2 = NULL,
    useTieProb = TRUE,
    depvar = NULL,
    second = FALSE,
    level = "period",
    condition = NULL, 
    sum.fun = mean,
    na.rm = TRUE,
    uncertainty = TRUE,
    nsim = 1000,
    useCluster = FALSE,
    nbrNodes = 1,
    clusterType = c("PSOCK", "FORK"),
    cluster = NULL,
    batch_dir = "temp",
    prefix = "simBatch_b",
    combine_batch = TRUE,
    batch_size = 50,
    keep_batch = FALSE,
    verbose = TRUE
){
    if (is.null(depvar)) depvar <- names(data[["depvars"]])[1]
    staticContributions <- calculateContribution(ans, data, depvar)

    if(!second){
        diffName <- "firstDiff"
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
            useTieProb = useTieProb
        )
    }else{
        diffName <- "secondDiff"
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
            useTieProb = useTieProb
        )
        ME <- "secondDiff"
    }
    sienaPostestimate(
        predictFun = diffFun,
        predictArgs = predictArgs,
        outcome = diffName,
        level = level,
        condition = condition,
        sum.fun = sum.fun,
        na.rm = na.rm,
        theta_hat = ans[["theta"]],
        cov_theta = ans[["covtheta"]],
        uncertainty = uncertainty,
        nsim = nsim,
        useCluster = useCluster,
        nbrNodes = nbrNodes,
        clusterType = clusterType,
        cluster = cluster,
        batch_dir = batch_dir,
        prefix = prefix,
        combine_batch = combine_batch,
        batch_size = batch_size,
        keep_batch = keep_batch,
        verbose = verbose
    )
}

predictFirstDiff <- function(ans, theta, staticContributions,
    useTieProb = TRUE,
    effectName, diff = NULL, contrast = NULL,
    interaction = FALSE,
    int_effectNames = NULL,
    mod_effectNames = NULL
) {
    effects <- ans[["effects"]] # provide effects instead?
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    if(length(theta) > length(includedEffects)) {
        thetaNoRate <- theta[noRateIncluded]
    } else {
        thetaNoRate <- theta
    }
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]


    df <- staticContributions
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("period", "ego"), useTieProb = useTieProb)
    df[["firstDiff"]] <- calculateFirstDiff(
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
        useTieProb = useTieProb, 
        tieProb = df[["tieProb"]]
    )
    df <- subset(df, df[["density"]] != 0)
    df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
    df
}

predictSecondDiff <- function(ans, theta, staticContributions,
    useTieProb = TRUE,
    effectName1, diff1 = NULL, contrast1 = NULL,
    interaction1 = FALSE,
    int_effectNames1 = NULL,
    mod_effectNames1 = NULL,
    effectName2, diff2 = NULL, contrast2 = NULL,
    interaction2 = FALSE,
    int_effectNames2 = NULL,
    mod_effectNames2 = NULL
) {
    effects <- ans[["effects"]] # provide effects instead?
    include <- effects[["include"]]
    includedEffects <- effects[include, ]
    noRateIncluded <- includedEffects[["type"]] != "rate"
    if(length(theta) > length(includedEffects)) {
        thetaNoRate <- theta[noRateIncluded]
    } else {
        thetaNoRate <- theta
    }
    effectNames  <- includedEffects[["shortName"]][noRateIncluded]


    df <- staticContributions
    df <- addUtilityColumn(df, effectNames, thetaNoRate)
    df <- addProbabilityColumn(df, group_vars = c("period", "ego"), useTieProb = useTieProb)
    df[["secondDiff"]] <- calculateSecondDiff(
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
        useTieProb = useTieProb, 
        tieProb = df[["tieProb"]]# if not tieProb NULL?
    )
    df <- subset(df, df[["density"]] != 0)
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
                               useTieProb = TRUE,
                               tieProb = NULL,
                               details = FALSE){
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
    # is the same for all without interaction except for *density
    utilDiff <- calculateUtilityDiff(effectName = effectName, diff = diff, 
                                      theta = theta, densityValue = densityValue,
                                      interaction = interaction,
                                      int_effectNames = int_effectNames,
                                      mod_effectNames = mod_effectNames,
                                      modContribution = modContribution,
                                      effectNames = effectNames)
  }
  expDiff <- exp(utilDiff)
  changeProb_cf <- as.vector(changeProb * expDiff / (1 - changeProb + changeProb * expDiff))
  if (useTieProb == TRUE) {
    tieProb_cf <- changeProb_cf
    tieProb_cf[densityValue == -1] <- 1 - changeProb_cf[densityValue == -1]
    firstDiff <- tieProb_cf - tieProb
  } else {
     firstDiff <- changeProb_cf - changeProb
  }

  # Ensure correct sign for changes in contrast
  if(!is.null(contrast)){
    firstDiff[which(newChangeStatistic == min(contrast))] <- -firstDiff[which(newChangeStatistic == min(contrast))]
  }

  if(details){ # mostly for debugging
    out <- data.frame(
      "firstDiff" = firstDiff,
      "utilDiff" = utilDiff,
      "newChangeProb" = changeProb_cf, # make clear if change or tieProb
      "oldChangeProb" = changeProb
    )
    if(useTieProb){
      out[["newTieProb"]] <- tieProb_cf
      out[["oldTieProb"]] <- tieProb
    }
    return(out)
  } else {
    return(firstDiff) # returns vector!
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
                                useTieProb = TRUE,
                                tieProb = NULL,
                                details = FALSE){
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
    useTieProb = useTieProb,
    tieProb = tieProb
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
  if (useTieProb == TRUE) {
    tieProb_cf21 <- changeProb_cf21
    tieProb_cf21[densityValue == -1] <- 1 - changeProb_cf21[densityValue == -1]
  }
  changeUtil21 <- changeUtil + utilDiff21
  ## dangerous with interactions because other effect values are not "corrected"

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
    useTieProb = useTieProb,
    tieProb = tieProb_cf21,
    interaction = interaction1,
    int_effectNames = int_effectNames1,
    mod_effectNames = mod_effectNames1,
    modContribution = modContribution1,
    effectNames = effectNames
  )

  secondDiff <- firstDiff2 - firstDiff
  
  if(!is.null(contrast2)){
    secondDiff[which(newChangeStatistic2 == min(contrast2))] <- -secondDiff[which(newChangeStatistic2 == min(contrast2))]
  }

  if(details){
    out <- data.frame(
      "utilDiff21" = utilDiff21,
      "newChangeProb21" = changeProb_cf21,
      "firstDiff1" = firstDiff,
      "firstDiff2" = firstDiff2,
      "secondDiff" = secondDiff)
     if(useTieProb){
      out[["newTieProb21"]] <- tieProb_cf21
      out[["oldTieProb"]] <- tieProb
    }
    return(out)
  } else {
    return(secondDiff)
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