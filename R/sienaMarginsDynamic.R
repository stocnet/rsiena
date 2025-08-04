#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: https://www.stats.ox.ac.uk/~snijders/siena/
# *
# * File: sienaMargins.r
# *
# * Description: Calculates predicted edge probabilities and
# * (average) marginal effects for simulated chains
# *****************************************************************************/

##@sienaAMEDynamic
sienaAMEDynamic <- function(ans,
                            data,
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
                            effects,
                            algorithm,
                            useTieProb = TRUE,
                            depvar = NULL,
                            second = FALSE,
                            level = "none",
                            condition = NULL,
                            sum.fun = mean,
                            na.rm = TRUE,
                            n3 = 500,
                            useChangeContributions = TRUE,
                            uncertainty = TRUE,
                            nsim = 100,
                            useCluster = FALSE,
                            nbrNodes = 1,
                            clusterType = c("PSOCK", "FORK"),
                            cluster = NULL,
                            batch_dir = "temp",
                            prefix = "simBatch_b",
                            batch_size = 1,
                            combine_batch = TRUE,
                            keep_batch = FALSE,
                            verbose = TRUE){
  if(!second){
    sim_fun <- simFirstDiffDynamic
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
      algorithm = algorithm,
      useTieProb = useTieProb,
      depvar = depvar,
      sim_theta = TRUE,
      level = level,
      condition = condition,
      sum.fun = sum.fun,
      na.rm = na.rm,
      n3 = n3,
      useChangeContributions = useChangeContributions
    )
    ME <- "firstDiff"
  }else{
    sim_fun <- simSecondDiffDynamic
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
      algorithm = algorithm,
      useTieProb = useTieProb,
      depvar = depvar,
      sim_theta = TRUE,
      level = level,
      condition = condition,
      sum.fun = sum.fun,
      na.rm = na.rm,
      n3 = n3,
      useChangeContributions = useChangeContributions
    )
    ME <- "secondDiff"
  }

# Create marginal function depending on non rate theta
  dynamic_ame_fun <- function(theta = NULL) {
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

  point_estimate2 <- dynamic_ame_fun(thetaNoRate)

  if(!useCluster){
    # not nice, but works since mclapply reduces to lapply for nbrNodes = 1
    clusterType <- "FORK"
    nbrNodes <- 1
  }

  # Point estimate
  sim_args_theta <- sim_args
  sim_args_theta[["sim_theta"]] <- FALSE

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
      batch_size = batch_size,
      combine_batch = combine_batch,
      keep_batch = keep_batch,
      verbose = verbose
    )

    uncert_AME <- agg(ME, uncert_AME, level = level, condition = condition, sum.fun = summarizeValue)

# ## delta method starting here
    # library(numDeriv)
    # ## function of dynamic_ame_fun and thetaNoRate (and point_estimate2 which should be unneccesary)
    # k <- length(point_estimate2)
    # p <- length(thetaNoRate)
    # J <- jacobian(dynamic_ame_fun, thetaNoRate)
    # cov_theta <- ans[["covtheta"]][noRate, noRate, drop=FALSE]
    # var_ame <- J %*% cov_theta %*% t(J)
    # se_ame <- sqrt(diag(var_ame))
  
    # hessians <- vector("list", k)
    # for(i in seq_len(k)) {
    #     hessians[[i]] <- hessian(function(th) dynamic_ame_fun(th)[i], thetaNoRate)
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

simFirstDiffDynamic <- function(ans, data, effectName,
                                diff = NULL, contrast = NULL,
                                interaction = FALSE,
                                int_effectNames = NULL,
                                mod_effectNames = NULL,
                                effects,
                                algorithm,
                                useTieProb = TRUE,
                                depvar = NULL,
                                sim_theta = TRUE,
                                aggregateValues = TRUE,
                                level = "none",
                                condition = NULL,
                                sum.fun = mean,
                                na.rm = TRUE,
                                n3 = NULL,
                                useChangeContributions = FALSE){
  changeProb <- changeUtil <- chain <- period <- ministep <- NULL # To resolve R CMD checks not understanding data.table syntax



  if(sim_theta){
    ## with estimated rate effects, this draws from the joint distribution!
    theta <- MASS::mvrnorm(n=1,
                           mu = ans$theta,
                           Sigma = ans$covtheta)
    useChangeContributions <- FALSE
    ans <- NULL

  }else{
    theta <- ans$theta
  }

  include <- effects$include
  includedEffects <- effects[include, ]
  noRateIncluded <- includedEffects[["type"]] != "rate"
  thetaNoRate <- theta[noRateIncluded]
  effectNames  <- includedEffects$shortName[noRateIncluded]

  if(is.null(depvar)){
    depvar <- names(data$depvars)[1]
  }
  df <- getChangeContributionsDynamic(
    ans = ans,
    data = data,
    theta = theta,
    algorithm = algorithm,
    effects = effects,
    depvar = depvar,
    n3 = n3,
    useChangeContributions = useChangeContributions,
    returnDataFrame = TRUE
  )

  df <- widenContribution(df)
  df <- addUtilityColumn(df, effectNames, thetaNoRate)
  df <- addProbabilityColumn(df, group_vars=c("chain", "period", "ministep"), useTieProb = useTieProb)

  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  df <- subset(df, df[["density"]] != 0)

  df[["firstDiff"]] <- calculateFirstDiff(
    densityValue = df[["density"]],
    changeProb = df[["changeProb"]],
    effectName = effectName, diff = diff, contrast = contrast,
    effectContribution = df[[effectName]],
    theta = thetaNoRate,
    interaction = interaction,
    int_effectNames = int_effectNames,
    mod_effectNames = mod_effectNames,
    modContribution = df[[mod_effectNames]],
    useTieProb = useTieProb,
    tieProb = df[["tieProb"]],
    effectNames = effectNames
  )

  # Transform contributions to change statistics for aggregation and output
  df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)

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

simSecondDiffDynamic <- function(ans, data, effectName1,
                                diff1 = NULL, contrast1 = NULL,
                                interaction1 = FALSE,
                                int_effectNames1 = NULL,
                                mod_effectNames1 = NULL,
                                effectName2,
                                diff2 = NULL, contrast2 = NULL,
                                interaction2 = FALSE,
                                int_effectNames2 = NULL,
                                mod_effectNames2 = NULL,
                                effects,
                                algorithm,
                                useTieProb = TRUE,
                                depvar = NULL,
                                sim_theta = TRUE,
                                aggregateValues = TRUE,
                                level = "none",
                                condition = NULL,
                                sum.fun = mean,
                                na.rm = TRUE,
                                n3 = NULL,
                                useChangeContributions = FALSE){
  if (is.null(depvar)) {
    depvar <- names(data[["depvars"]])[1]
  }

  include <- effects[["include"]]
  includedEffects <- effects[include, ]
  noRateIncluded <- includedEffects[["type"]] != "rate"
  effectNames <- includedEffects[["shortName"]][noRateIncluded]

  if(sim_theta){
    ## might not work with estimated rate effects
    theta <- MASS::mvrnorm(n=1,
                               mu = ans$theta,
                               Sigma = ans$covtheta)
    useChangeContributions <- FALSE
    ans <- NULL
  }else{
    theta <- ans$theta
  }

  thetaNoRate <- theta[noRateIncluded]
  df <- getChangeContributionsDynamic(
    ans = ans,
    data = data,
    theta = theta,
    algorithm = algorithm,
    effects = effects,
    depvar = depvar,
    n3 = n3,
    useChangeContributions = useChangeContributions,
    returnDataFrame = TRUE
  )  
  df <- widenContribution(df)
  df <- addUtilityColumn(df, effectNames, thetaNoRate)
  df <- addProbabilityColumn(df, group_vars = c("chain", "period", "ministep"), useTieProb = useTieProb)

  # density == 0 is only relevant for the normalizing constant, not needed for difference calculation
  df <- subset(df, df[["density"]] != 0)

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

  # Transform contributions to change statistics for aggregation and output
  df <- conditionalReplace(df, df[["density"]] == -1, setdiff(effectNames, "density"), function(x) x * -1)

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
