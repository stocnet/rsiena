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
                            sienaData,
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
                            effectNames,
                            effects,
                            tieProb = TRUE,
                            depvar = NULL,
                            second = FALSE,
                            level = "Period",
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
  if(!second){
    sim_fun <- simFirstDiff
    sim_args <- list(
      ans = ans,
      sienaData = sienaData,
      effectName = effectName1,
      diff = diff1,
      contrast = contrast1,
      interaction = interaction1,
      int_effectNames = int_effectNames1,
      mod_effectNames = mod_effectNames1,
      effectNames = effectNames,
      effects = effects,
      tieProb = tieProb,
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
      sienaData = sienaData,
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
      effectNames = effectNames,
      effects = effects,
      tieProb = tieProb,
      depvar = depvar,
      sim_theta = TRUE,
      level = level,
      condition = condition,
      sum.fun = sum.fun,
      na.rm = na.rm
    )
    ME <- "secondDiff"
  }
  
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
    
    if(is.null(condition)){
      return(cbind(point_AME, uncert_AME)) # merge would be safer, but does not handle results without id vars well
    } else {
      return(merge(point_AME, uncert_AME))
    }
  }
}

simFirstDiff <- function(ans, sienaData, 
                         effectName, diff = NULL, contrast = NULL,
                         interaction = FALSE,
                         int_effectNames = NULL,
                         mod_effectNames = NULL,
                         effectNames, 
                         effects = NULL, # currently unused
                         tieProb = TRUE, 
                         depvar = NULL, # curently unused
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
    ## add option to change algorithm -> number of chains -> 
  }else{
    theta <- ans$theta
  }
  prob_sim <- calculateChoiceProbability(ans, sienaData, tieProb = tieProb, theta = theta)
  
  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  prob_sim <- subset(prob_sim, density != 0)

  firstDiff_sim <- calculateFirstDiff(prob_sim, effectName = effectName, diff = diff, contrast = contrast,
                                      theta = theta,
                                      tieProb = tieProb, 
                                      interaction = interaction,
                                      int_effectNames = int_effectNames,
                                      mod_effectNames = mod_effectNames,
                                      effectNames = effectNames)
  # Transform contributions to change statistics for aggregation and output
  prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] <- prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] * -1
  firstDiff_sim <- cbind(prob_sim, firstDiff_sim) ## probably inefficient
  if(aggregateValues) {
    firstDiff_sim <- agg("firstDiff",
                           firstDiff_sim,
                           level = level,
                           condition = condition,
                           sum.fun = sum.fun,
                           na.rm = TRUE)
  }
  firstDiff_sim
}

simSecondDiff <- function(ans, sienaData, 
                          effectName1, diff1 = NULL, contrast1 = NULL,
                          interaction1 = FALSE,
                          int_effectNames1 = NULL,
                          mod_effectNames1 = NULL,
                          effectName2, diff2 = NULL, contrast2 = NULL,
                          interaction2 = FALSE,
                          int_effectNames2 = NULL,
                          mod_effectNames2 = NULL,
                          effectNames, 
                          effects = NULL, # currently unused
                          tieProb = TRUE,
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
    ## add option to change algorithm -> number of chains -> 
  }else{
    theta <- ans$theta
  }
  ## add option to provide prob_sim OR combine firstDiff and secondDiff in one call
  prob_sim <- calculateChoiceProbability(ans, sienaData, tieProb = tieProb, theta = theta)
  
  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  prob_sim <- subset(prob_sim, density != 0)
  
  secondDiff_sim <- calculateSecondDiff(prob_sim, 
                                        effectName1 = effectName1, diff1 = diff1, contrast1 = contrast1,
                                        interaction1 = interaction1,
                                        int_effectNames1 = int_effectNames1,
                                        mod_effectNames1 = mod_effectNames1,
                                        effectName2 = effectName2, diff2 = diff2, contrast2 = contrast2,
                                        interaction2 = interaction2,
                                        int_effectNames2 = int_effectNames2,
                                        mod_effectNames2 = mod_effectNames2,
                                        theta = theta,
                                        tieProb = tieProb, 
                                        effectNames = effectNames)
  
  # Transform contributions to change statistics for aggregation and output
  prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] <- prob_sim[prob_sim[, "density"] == -1, setdiff(effectNames, "density")] * -1
  secondDiff_sim <- cbind(prob_sim, secondDiff_sim)
  
  if(aggregateValues) {
    secondDiff_sim <- agg("secondDiff",
                            secondDiff_sim,
                           level = level,
                           condition = condition,
                           sum.fun = sum.fun,
                           na.rm = TRUE)
  }
  secondDiff_sim
}

calculateFirstDiff <- function(df, effectName, 
                               diff=NULL,
                               contrast=NULL, 
                               theta, 
                               tieProb = TRUE, 
                               effectNames,
                               interaction = FALSE,
                               int_effectNames = NULL, # later it would be nice to detect these automatically
                               mod_effectNames = NULL){
  density <- df[, "density"]
  
  if (effectName == "density") {
    if((!is.null(diff))) stop("firstDiff for density must be contrast c(-1,1)")
    if(!is.null(contrast)){
      if(any(setdiff(contrast, c(-1,1)))) {stop("firstDiff for density can only be be calculated for c(-1,1)")}
      if(interaction == TRUE) {stop("Interaction with density is not possible")}
      old_value <- density
      new_value <- rep(NA, length(old_value))
      new_value[old_value == contrast[1]] <- contrast[2]
      new_value[old_value == contrast[2]] <- contrast[1]
      util_diff <- ifelse(is.na(new_value), NA, -1 * df[,"changeUtil"] - df[,"changeUtil"])
      density <- new_value
    }
  } else {
    if(!is.null(contrast)){
      old_value <- density * df[, effectName]
      new_value <- rep(NA, length(old_value))
      new_value[old_value == contrast[1]] <- contrast[2]
      new_value[old_value == contrast[2]] <- contrast[1]
      diff <- new_value - old_value
    }
    ## is the same for all without interaction except for *density
    util_diff <- calculateUtilityDiff(df = df, effectName = effectName, diff = diff, 
                                      theta = theta, density = density,
                                      interaction = interaction,
                                      int_effectNames = int_effectNames,
                                      mod_effectNames = mod_effectNames,
                                      effectNames = effectNames)
  }
  
  exp_diff <- exp(util_diff)
  prob <- df[,"changeProb"]
  prob_cf <- as.vector(prob * exp_diff / (1-prob + prob * exp_diff))
  if (tieProb == TRUE) {
    prob_cf[density == -1] <- 1 - prob_cf[density == -1]
    prob <- df[,"tieProb"]
  } else {
    prob <- df[,"changeProb"]
  }
  
  firstDiff <- prob_cf - prob
  if(!is.null(contrast)){
    firstDiff[which(new_value == min(contrast))] <- -firstDiff[which(new_value == min(contrast))]
  }
  
  data.frame("utilDiff" = util_diff,
        "newProb" = prob_cf,
        "oldProb" = prob,
        "firstDiff" = firstDiff) # could also return new density values and changeUtil
}

calculateSecondDiff <- function(df, effectName1, 
                                diff1 = NULL, 
                                contrast1 = NULL,
                                interaction1 = FALSE,
                                int_effectNames1 = NULL, # later it would be nice to detect these automatically
                                mod_effectNames1 = NULL,
                                effectName2, 
                                diff2 = NULL,
                                contrast2 = NULL,
                                interaction2 = FALSE,
                                int_effectNames2 = NULL, # later it would be nice to detect these automatically
                                mod_effectNames2 = NULL,
                                effectNames,
                                theta, 
                                tieProb = TRUE){
  # Calculate firstDiff if it is not already present in dataframe
  if(!("firstDiff" %in% names(df))) {
    df[, "firstDiff"] <- calculateFirstDiff(df, effectName = effectName1, 
                                            diff = diff1,
                                            contrast = contrast1, 
                                            theta = theta, 
                                            tieProb = tieProb,
                                            effectNames = effectNames,
                                            interaction = interaction1,
                                            int_effectNames = int_effectNames1,
                                            mod_effectNames = mod_effectNames1)[,"firstDiff"]
  }
  firstDiff <- df[, "firstDiff"]
  # might be more efficient to use tieProb when we already have it
  prob <- df[,"changeProb"]
  # which ...
  effectNum1 <- which(effectNames == effectName1)
  effectNum2 <- which(effectNames == effectName2)

  density <- df[, "density"]
  
  
  if(!is.null(contrast2)){
    old_value <- density * df[, effectName2]
    new_value <- rep(NA, length(old_value))
    new_value[old_value == contrast2[1]] <- contrast2[2]
    new_value[old_value == contrast2[2]] <- contrast2[1]
    diff2 <- new_value - old_value
  }
  
  
  util_diff21 <- calculateUtilityDiff(df = df, effectName = effectName2,
                                      diff = diff2, theta = theta, 
                                      density = density,
                                      interaction = interaction2,
                                      int_effectNames = int_effectNames2,
                                      mod_effectNames = mod_effectNames2,
                                      effectNames = effectNames)
  
  exp_diff21 <- exp(util_diff21)
  prob_cf21 <- as.vector(prob * exp_diff21 / (1-prob + prob * exp_diff21))
  
  util_orig <- df[,"changeUtil"]
  temp <- df
  temp[,"changeUtil"] <- util_orig + util_diff21
  temp[, "changeProb"] <- prob_cf21
  ## dangerous with interactions because other effect values are not "corrected"
  temp <- calculateFirstDiff(temp, effectName = effectName1, 
                                          diff = diff1,
                                          contrast = contrast1, 
                                          theta = theta, 
                                          tieProb = tieProb,
                                          effectNames = effectNames,
                                          interaction = interaction1,
                                          int_effectNames = int_effectNames1,
                                          mod_effectNames = mod_effectNames1)
  firstDiff2 <- temp[,"firstDiff"]
  util_diff22 <- temp[,"utilDiff"]
  prob_cf22 <- temp[,"newProb"]
  prob21 <- temp[,"oldProb"]
  
  ## should not be necessary, since already set in call to calculateFirstDiff
  
  ## Check if correct
  # if (tieProb == TRUE) {
  #   prob_cf21 <- ifelse(df[, "density"] == -1, 1 - prob_cf1, prob_cf21)
  #   prob_cf22 <- ifelse(df[, "density"] == -1, 1 - prob_cf2, prob_cf22)
  #   prob <- ifelse(df[, "density"] == -1, 1 - prob, prob)
  # }
  
  secondDiff <- firstDiff2 - firstDiff
  if(!is.null(contrast2)){
    secondDiff[which(new_value == min(contrast2))] <- -secondDiff[which(new_value == min(contrast2))]
  }
  ## what do we actually want to return and how?
  data.frame("utilDiff21" = util_diff21,
        "newProb21" = prob_cf21,
        "utilDiff22" = util_diff22,
        "newProb22" = prob_cf22,
        "oldProb" = prob21,
        "firstDiff2_test" = firstDiff2,
        "firstDiff_test" = firstDiff,
        "secondDiff" = secondDiff)
}

calculateUtilityDiff <- function(df, effectName, diff, 
                                 theta, density,
                                 interaction = FALSE,
                                 int_effectNames = NULL,
                                 mod_effectNames = NULL,
                                 effectNames = NULL){
  effectNum <- which(effectNames == effectName)
  if(interaction == TRUE){
    moderator_values <- density * df[, mod_effectNames]
    interaction_effectNums <- which(effectNames == int_effectNames)
    util_diff <- density * (diff * theta[effectNum] + 
                              diff * moderator_values * theta[interaction_effectNums])
  } else {
    util_diff <- density * diff * theta[effectNum]
  }
  return(util_diff)
}