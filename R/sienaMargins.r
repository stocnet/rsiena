#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaMargins.r
# *
# * Description: Calculates predicted edge probabilities and
# * (average) marginal effects
# *****************************************************************************/

##@sienaAME
sienaAME <- function(ans,
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
                            effects,
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
                            keep_batch = FALSE,
                            verbose = TRUE){

  ## add check if neither diff nor cont is set
  if(!second){
    sim_fun <- simFirstDiff
    sim_args <- list(
      ans = ans,
      data = data,
      effectName = effectName1,
      diff = diff1,
      contrast = contrast1,
      interaction = interaction1,
      int_effectNames = int_effectNames1,
      mod_effectNames = mod_effectNames1,
      effects = effects,
      useTieProb = useTieProb,
      depvar = depvar,
      sim_theta = TRUE,
      level = level,
      condition = condition,
      sum.fun = sum.fun,
      na.rm = na.rm
    )
    ME <- "firstDiff"
  }else{
    sim_fun <- simSecondDiff
    sim_args <- list(
      ans = ans,
      data = data,
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
      effects = effects,
      useTieProb = useTieProb,
      depvar = depvar,
      sim_theta = TRUE,
      level = level,
      condition = condition,
      sum.fun = sum.fun,
      na.rm = na.rm
    )
    ME <- "secondDiff"
  }

  # Create marginal function depending on non rate theta
  static_ame_fun <- function(theta = NULL) {
    # Is also done in calculateContribution
    # Make a copy and set only the non-rate effects
    sim_args_theta <- sim_args
    #ans_tmp <- ans
    if(!is.null(theta))
    {
      #ans_tmp[["theta"]][["ans"]][noRate] <- theta
      sim_args_theta[["ans"]][["theta"]][noRate] <- theta
    }
    sim_args_theta[["sim_theta"]] <- FALSE
    
    AME <- do.call(sim_fun, sim_args_theta)
    AME

  }
  effectsIncluded <- ans$effects[ans$effects$include == TRUE]
  noRate <- effectsIncluded$type != "rate"
  effectsIncludedNoRate <- effectsIncluded[noRate, ]
  thetaNoRate <- ans[["theta"]][noRate]

  point_estimate2 <- static_ame_fun(thetaNoRate)

  if(!useCluster){
    # not nice, but works since mclapply reduces to lapply for nbrNodes = 1
    clusterType <- "FORK" 
    nbrNodes <- 1
  }
  
  # Point estimate
  sim_args_theta <- sim_args
  sim_args_theta$sim_theta <- FALSE
  
  point_AME <- do.call(sim_fun, sim_args_theta)

  if(!uncertainty){
    return(point_AME)
  }else{

    uncert_AME <- drawSim(
      sim_fun = sim_fun,
      sim_args = sim_args,
      nbrNodes = nbrNodes,
      nsim = nsim,
      clusterType = clusterType,
      cluster = cluster,
      batch_dir = batch_dir,
      prefix = prefix,
      combine_batch = combine_batch,
      keep_batch = keep_batch,
      verbose = verbose
    )
    
    uncert_AME <- agg(ME, uncert_AME, level = level, condition = condition, sum.fun = summarizeValue)

    # ## delta method starting here
    # library(numDeriv)
    # ## function of static_ame_fun and thetaNoRate (and point_estimate2 which should be unneccesary)
    # k <- length(point_estimate2)
    # p <- length(thetaNoRate)
    # J <- jacobian(static_ame_fun, thetaNoRate)
    # cov_theta <- ans[["covtheta"]][noRate, noRate, drop=FALSE]
    # var_ame <- J %*% cov_theta %*% t(J)
    # se_ame <- sqrt(diag(var_ame))
  
    # hessians <- vector("list", k)
    # for(i in seq_len(k)) {
    #     hessians[[i]] <- hessian(function(th) static_ame_fun(th)[i], thetaNoRate)
    # }

    # Tmat <- matrix(NA_real_, nrow=k, ncol=k)

    # for(i in seq_len(k)) {
    #   for(j in i:k) { # only upper triangle
    #     val <- sum(hessians[[i]] * (cov_theta %*% hessians[[j]] %*% cov_theta))
    #     Tmat[i, j] <- val
    #     if (i != j) Tmat[j, i] <- val
    #   }
    # }

    # var_ame2 <- var_ame + 0.5 * Tmat
    # se_ame2 <- sqrt(diag(var_ame2))

    # delta_df <- data.frame(
    #   ame = point_estimate2,
    #   se = se_ame,
    #   se2 = se_ame2
    # )
  
    if(is.null(condition)){
      AME <- cbind(point_AME, uncert_AME) # merge would be safer, but does not handle results without id vars well
    } else {
      AME <- merge(point_AME, uncert_AME)
    }
    # return(list(margin_sim = AME1, 
    #   margin_delta = delta_df))
    AME
  }
}

simFirstDiff <- function(ans, data, 
                         effectName, diff = NULL, contrast = NULL,
                         interaction = FALSE,
                         int_effectNames = NULL,
                         mod_effectNames = NULL,
                         effects = NULL, # currently unused
                         useTieProb = TRUE, 
                         depvar = NULL, # curently unused
                         sim_theta = TRUE,
                         aggregateValues = TRUE,
                         level = "period",
                         condition = NULL,
                         sum.fun = mean,
                         na.rm = TRUE){
  changeProb <- changeUtil <- chain <- period <- ministep <- NULL # To resolve R CMD checks not understanding data.table syntax

  ## effectNames can just be extracted
  if(sim_theta){
    ## might not work with estimated rate effects
    theta <- MASS::mvrnorm(n=1,
                           mu = ans$theta,
                           Sigma = ans$covtheta)
    ## add option to change algorithm -> number of chains -> 
  }else{
    theta <- ans$theta
  }
  df <- calculateChoiceProbability(ans, data, useTieProb = useTieProb, theta = theta)
  
  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  df <- subset(df, density != 0)
  ## change for data.table

  include <- effects$include
  includedEffects <- effects[include, ]
  noRateIncluded <- includedEffects[["type"]] != "rate"
  thetaNoRate <- theta[noRateIncluded]
  effectNames  <- includedEffects$shortName[noRateIncluded]

  df[["firstDiff"]] <- calculateFirstDiff(densityValue = df[["density"]],
                                      changeProb = df[["changeProb"]],
                                      effectName = effectName, 
                                      effectContribution = df[[effectName]],
                                      diff = diff, contrast = contrast,
                                      effectNames = effectNames,
                                      theta = thetaNoRate,
                                      useTieProb = useTieProb, 
                                      tieProb = df[["tieProb"]],
                                      interaction = interaction,
                                      int_effectNames = int_effectNames,
                                      mod_effectNames = mod_effectNames,
                                      modContribution = df[[mod_effectNames]])
  # Transform contributions to change statistics for aggregation and output
  df  <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)
  
  if(aggregateValues) {
    df <- agg("firstDiff",
                           df,
                           level = level,
                           condition = condition,
                           sum.fun = sum.fun,
                           na.rm = TRUE)
  }
  df
}

simSecondDiff <- function(ans, data, 
                          effectName1, diff1 = NULL, contrast1 = NULL,
                          interaction1 = FALSE,
                          int_effectNames1 = NULL,
                          mod_effectNames1 = NULL,
                          effectName2, diff2 = NULL, contrast2 = NULL,
                          interaction2 = FALSE,
                          int_effectNames2 = NULL,
                          mod_effectNames2 = NULL,
                          effects = NULL, # currently unused
                          useTieProb = TRUE,
                          depvar = NULL, # currently unused
                          sim_theta = TRUE,
                          aggregateValues = TRUE,
                          level = "period",
                          condition = NULL,
                          sum.fun = mean,
                          na.rm = TRUE){
  ## effectNames can just be extracted
  if(sim_theta){
    ## might not work with estimated rate effects
    theta <- MASS::mvrnorm(n=1,
                           mu = ans$theta,
                           Sigma = ans$covtheta)
  }else{
    theta <- ans$theta
  }
  ## add option to provide df OR combine firstDiff and secondDiff in one call
  df <- calculateChoiceProbability(ans, data, useTieProb = useTieProb, theta = theta)
  
  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  df <- subset(df, density != 0)
  
  include <- effects$include
  includedEffects <- effects[include, ]
  noRateIncluded <- includedEffects[["type"]] != "rate"
  thetaNoRate <- theta[noRateIncluded]
  effectNames  <- includedEffects$shortName[noRateIncluded]

  df[["secondDiff"]] <- calculateSecondDiff(densityValue = df[["density"]], 
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
  # Transform contributions to change statistics for aggregation and output
  conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)  
  
  if(aggregateValues) {
    df <- agg("secondDiff",
                            df,
                           level = level,
                           condition = condition,
                           sum.fun = sum.fun,
                           na.rm = TRUE)
  }
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
  # density <- df[[ "density"]]
  
  if (effectName == "density") {
    if((!is.null(diff))) stop("firstDiff for density must be contrast c(-1,1)")
    if(!is.null(contrast)){
      if(any(setdiff(contrast, c(-1,1)))) {stop("firstDiff for density can only be be calculated for c(-1,1)")}
      if(interaction == TRUE) {stop("Interaction with density is not possible")}
      oldChangeStatistic <- densityValue
      newChangeStatistic <- rep(NA, length(oldChangeStatistic))
      newChangeStatistic[oldChangeStatistic == contrast[1]] <- contrast[2]
      newChangeStatistic[oldChangeStatistic == contrast[2]] <- contrast[1]
      # changeUtils <- df[["changeUitil"]]
      utilDiff <- ifelse(is.na(newChangeStatistic), NA, -2 * changeUtils)
      densityValue <- newChangeStatistic
    }
  } else {
    if(!is.null(contrast)){
      # effectContributions <- df[["effectName"]]
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
  # changeProb <- df[["changeProb"]]
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

  # which ... not used anymore?
  # effectNum1 <- which(effectNames == effectName1)
  # effectNum2 <- which(effectNames == effectName2)

  # density <- df[["density"]]
  
  # Calculate Probabilitiy if only effect2 has been changed
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
  # changeProb <- df[["changeProb"]]
  changeProb_cf21 <- as.vector(changeProb * expDiff21 / (1-changeProb + changeProb * expDiff21))
  if (useTieProb == TRUE) {
    tieProb_cf21 <- changeProb_cf21
    tieProb_cf21[densityValue == -1] <- 1 - changeProb_cf21[densityValue == -1]
  }
  # changeUtil <- df[["changeUtil"]]
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